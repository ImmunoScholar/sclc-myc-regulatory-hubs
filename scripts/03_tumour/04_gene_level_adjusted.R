# -----------------------------------------------------------------------------
# 04_gene_level_adjusted.R — per-gene paralog association, lineage-adjusted.
#
#   Rscript scripts/03_tumour/04_gene_level_adjusted.R
#
# WHY THIS EXISTS. The regulon-level test is null: lineage explains 39-49% of
# regulon score variance against 0.0-1.6% uniquely attributable to the paralog,
# and 0/3 associations survive adjustment (R-01, D-032).
#
# But MOES operates PER GENE, not per regulon. A 500-gene aggregate can be null
# while individual genes carry lineage-independent paralog signal — averaging over
# 370 measurable genes will bury a subset of 20 real ones. This tests that
# directly rather than assuming it either way.
#
# METHOD. Partial Spearman: rank-transform each gene across samples, regress out
# the lineage design (NE score + ASCL1 + NEUROD1 + POU2F3 + YAP1), then correlate
# residuals. df = n - 2 - k = 79 - 2 - 5 = 72.
#
# THE HONEST TEST IS ENRICHMENT, NOT COUNT. With ~33,700 genes tested, some will
# survive FDR by chance. The question is whether REGULON genes survive at a higher
# rate than non-regulon genes. If they do not, the surviving genes are not
# regulon-related and the transcriptional domain has no paralog-specific content
# to contribute.
#
# OUTCOME DETERMINES M7-M8:
#   enrichment present -> lineage-adjusted per-gene evidence feeds MOES
#   enrichment absent   -> transcriptional domain DROPS OUT; MOES runs on 3 domains
#
# Output: data/processed/tumour/gene_level_adjusted.rds
#         data/metadata/m6_gene_level_survivors.csv
#         results/tables/m6_gene_level_adjusted.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml) })
CFG <- yaml::read_yaml("config/params.yml")
set.seed(CFG$project$seed)

tum <- readRDS("data/processed/tumour/gse60052.rds")
lin <- readRDS("data/processed/tumour/lineage_confounding.rds")
reg <- readRDS("data/processed/regions/regulons.rds")$regulons
E   <- tum$expr_filtered[, !tum$is_normal, drop = FALSE]
n   <- ncol(E)
cat("genes: ", format(nrow(E), big.mark = ","), "   tumours: ", n, "\n", sep = "")

TFS <- intersect(c("ASCL1","NEUROD1","POU2F3","YAP1"), rownames(E))
PARA <- intersect(c("MYC","MYCN","MYCL1"), rownames(E))

# --- lineage design -----------------------------------------------------------
tfz <- t(scale(t(E[TFS, , drop = FALSE])))
X <- cbind(1, NE = lin$ne_score, t(tfz))
k <- ncol(X) - 1
df <- n - 2 - k
cat("covariates adjusted: NE score + ", paste(TFS, collapse = " + "),
    "  (k = ", k, ", df = ", df, ")\n", sep = "")
cat("detectable |partial rho| at p<0.05: ", round(sqrt(qf(0.95, 1, df) / (qf(0.95, 1, df) + df)), 3), "\n\n", sep = "")

# --- rank-transform, then residualise ----------------------------------------
# Rank first so the result is a partial SPEARMAN rather than partial Pearson —
# these are log2 values with 51% zeros, so a rank-based measure is the right one.
Rk <- t(apply(E, 1, rank, ties.method = "average"))
qrX <- qr(X)
Res <- t(qr.resid(qrX, t(Rk)))          # gene x sample residuals
rownames(Res) <- rownames(E)

partial_assoc <- function(paralog) {
  yp <- Res[paralog, ]
  sy <- sqrt(sum(yp^2))
  sx <- sqrt(rowSums(Res^2))
  r <- as.numeric(Res %*% yp) / (sx * sy)
  r[!is.finite(r)] <- NA_real_
  tstat <- r * sqrt(df / (1 - r^2))
  p <- 2 * stats::pt(abs(tstat), df, lower.tail = FALSE)
  data.frame(gene = rownames(Res), rho = r, p = p,
             fdr = stats::p.adjust(p, "BH"), stringsAsFactors = FALSE)
}

cat("=========== per-gene partial association ===========\n")
out <- list(); summ <- data.frame()
for (p in PARA) {
  a <- partial_assoc(p)
  a <- a[a$gene != p, ]                       # exclude the paralog itself
  a$in_regulon <- a$gene %in% reg[[p]]
  out[[p]] <- a

  n_sig <- sum(a$fdr < 0.05, na.rm = TRUE)
  n_sig_pos <- sum(a$fdr < 0.05 & a$rho > 0, na.rm = TRUE)
  in_r <- a$in_regulon
  sig  <- a$fdr < 0.05 & a$rho > 0
  # enrichment: do regulon genes survive at a higher rate than non-regulon genes?
  tb <- table(factor(in_r, c(FALSE,TRUE)), factor(sig, c(FALSE,TRUE)))
  ft <- if (all(dim(tb) == c(2,2)) && sum(tb[,2]) > 0)
          stats::fisher.test(tb) else NULL
  rate_in  <- if (sum(in_r) > 0) mean(sig[in_r], na.rm = TRUE) else NA
  rate_out <- if (sum(!in_r) > 0) mean(sig[!in_r], na.rm = TRUE) else NA

  summ <- rbind(summ, data.frame(
    paralog = p, n_tested = nrow(a),
    regulon_genes_tested = sum(in_r),
    n_fdr05 = n_sig, n_fdr05_positive = n_sig_pos,
    regulon_survivors = sum(sig & in_r, na.rm = TRUE),
    rate_in_regulon = round(100 * rate_in, 2),
    rate_outside = round(100 * rate_out, 2),
    enrich_or = if (is.null(ft)) NA else round(as.numeric(ft$estimate), 2),
    enrich_p = if (is.null(ft)) NA else signif(ft$p.value, 3),
    stringsAsFactors = FALSE))

  cat(sprintf("  %-6s %5d regulon genes | %5d genes survive FDR<0.05 positive | %3d are regulon members\n",
              p, sum(in_r), n_sig_pos, sum(sig & in_r, na.rm = TRUE)))
  cat(sprintf("         survival rate: in regulon %.2f%%  outside %.2f%%  OR %s  p %s\n",
              100 * rate_in, 100 * rate_out,
              if (is.null(ft)) "NA" else sprintf("%.2f", ft$estimate),
              if (is.null(ft)) "NA" else format.pval(ft$p.value, digits = 3)))
}

cat("\n=========== enrichment summary ===========\n")
print(summ[, c("paralog","regulon_genes_tested","n_fdr05_positive",
               "regulon_survivors","rate_in_regulon","rate_outside",
               "enrich_or","enrich_p")], row.names = FALSE)

enriched <- !is.na(summ$enrich_p) & summ$enrich_p < 0.05 & summ$enrich_or > 1
cat("\nregulons enriched for lineage-independent survivors: ",
    sum(enriched), "/", nrow(summ), "\n", sep = "")

# --- top survivors, the candidate MOES evidence ------------------------------
cat("\n=========== top lineage-independent regulon genes ===========\n")
top_all <- data.frame()
for (p in PARA) {
  a <- out[[p]]
  s <- a[a$in_regulon & a$fdr < 0.05 & a$rho > 0, ]
  if (!nrow(s)) { cat("  ", p, ": none\n", sep = ""); next }
  s <- s[order(-s$rho), ]
  cat("  ", p, " (", nrow(s), " genes): ",
      paste(utils::head(s$gene, 12), collapse = ", "),
      if (nrow(s) > 12) " ..." else "", "\n", sep = "")
  top_all <- rbind(top_all, data.frame(paralog = p, s[, c("gene","rho","fdr")],
                                       stringsAsFactors = FALSE))
}
write.csv(top_all, "data/metadata/m6_gene_level_survivors.csv", row.names = FALSE)
write.csv(summ, "data/metadata/m6_gene_level_summary.csv", row.names = FALSE)

# --- verdict ------------------------------------------------------------------
cat("\n=========== VERDICT: does MOES get a transcriptional domain? ===========\n")
if (sum(enriched) >= 2) {
  cat("YES. Regulon genes survive lineage adjustment at a higher rate than\n")
  cat("background in ", sum(enriched), "/", nrow(summ), " paralogs. The per-gene\n", sep = "")
  cat("lineage-adjusted association is admissible as MOES transcriptional evidence.\n")
  cat("The regulon AGGREGATE remains null — only the per-gene layer is usable, and\n")
  cat("that distinction must be stated wherever this domain is used.\n")
} else if (sum(enriched) == 1) {
  cat("PARTIAL. Only ", summ$paralog[enriched], " shows enrichment. Admit the\n", sep = "")
  cat("transcriptional domain for that paralog only; the others contribute none.\n")
} else {
  cat("NO. Regulon genes do not survive lineage adjustment at above background\n")
  cat("rate. The transcriptional domain has no paralog-specific content to add and\n")
  cat("DROPS OUT of MOES, which then runs on three domains (cis-regulatory,\n")
  cat("network, functional). This must be declared in the methods, not omitted.\n")
}

saveRDS(list(per_gene = out, summary = summ, survivors = top_all,
             covariates = c("NE_score", TFS), df = df),
        "data/processed/tumour/gene_level_adjusted.rds")

md <- c("# M6 — per-gene paralog association, lineage-adjusted", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste0("GSE60052, n = ", n, " tumours. Adjusted for NE score + ",
               paste(TFS, collapse = ", "), " (df = ", df, ")."), "",
        "The regulon-level test is null (R-01). MOES operates per gene, so this asks",
        "whether individual genes carry lineage-independent paralog signal. The test",
        "is ENRICHMENT relative to non-regulon genes, not a raw survivor count —",
        "with ~33,700 genes some survive FDR by chance.", "",
        "| paralog | regulon genes | survivors (FDR<0.05, +) | in regulon | rate in | rate out | OR | p |",
        "|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(summ)))
  md <- c(md, sprintf("| %s | %d | %d | %d | %s%% | %s%% | %s | %s |",
                      summ$paralog[i], summ$regulon_genes_tested[i],
                      summ$n_fdr05_positive[i], summ$regulon_survivors[i],
                      summ$rate_in_regulon[i], summ$rate_outside[i],
                      summ$enrich_or[i], summ$enrich_p[i]))
writeLines(md, "results/tables/m6_gene_level_adjusted.md")
cat("\nwrote results/tables/m6_gene_level_adjusted.md\n")

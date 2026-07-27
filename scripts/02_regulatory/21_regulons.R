# -----------------------------------------------------------------------------
# 21_regulons.R — paralog regulons and their internal-validity gate.
#
#   Rscript scripts/02_regulatory/21_regulons.R
#
# CONSTRUCTION (config regulons, as clarified in D-028)
# A paralog regulon is the TOP-RANKED genes by aggregate link score to that
# paralog's active regions — not every gene with a link. MYC links to thousands
# of genes and singscore/GSVA need focused sets, so the selection rule is part of
# the definition and is stated, not applied silently.
#
# VALIDITY (D-028). The original spec was not computable and was circular:
# it wanted an amplified-vs-non-amplified AUC (MYC amplification is UNVERIFIED,
# D-026) and expression correlation (no cell-line RNA-seq). And a regulon built
# from a paralog's own lines scores high in those lines BY CONSTRUCTION.
#
# Three checks that can each fail:
#   1. LEAVE-ONE-LINE-OUT. Build the paralog's active regions from its remaining
#      line(s); test whether the HELD-OUT line's own ChIP signal is enriched
#      there against background. The held-out ChIP never informed the regions, so
#      this is not circular.
#   2. CROSS-PARALOG DISTINCTNESS. Regulons more similar than max_regulon_jaccard
#      are not separable, and that must be reported rather than glossed.
#   3. BENCHMARK. The MYC regulon should be enriched for HALLMARK_MYC_TARGETS_V1.
#      A MYC regulon with no excess of known MYC targets is not a MYC regulon.
#
# Output: data/processed/regions/regulons.rds
#         data/metadata/regulon_validity.csv
#         results/tables/m5_regulons.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(yaml); library(msigdbr)
})

CFG <- yaml::read_yaml("config/params.yml")
RG  <- CFG$regulons; RV <- CFG$regulon_validity
FOLD <- CFG$active_regions$primary_fold
set.seed(CFG$project$seed)

nrm  <- readRDS("data/processed/signal/region_signal_normalised.rds")
p2g  <- readRDS("data/processed/regions/peak_to_gene.rds")
M    <- nrm$mean_fob; meta <- nrm$meta; u <- nrm$regions
links <- p2g$retained
cat("links (high+moderate): ", format(nrow(links), big.mark = ","), "\n", sep = "")

PARALOG_LINES <- list(MYC = c("H1048","H211","H524","H847","SHP77"),
                      MYCN = c("H526","H69"), MYCL1 = c("COLO668","H889"))
col_for <- function(a, l) { j <- which(meta$assay == a & meta$line == l); if (length(j)==1) j else NA }

active_in_line <- function(p, l, f = FOLD) {
  jp <- col_for(p, l); jk <- col_for("H3K27ac", l)
  if (is.na(jp) || is.na(jk)) return(NULL)
  (M[, jp] >= f) & (M[, jk] >= f)
}
active_set <- function(p, lines = PARALOG_LINES[[p]], min_rep = 2) {
  mats <- Filter(Negate(is.null), lapply(lines, active_in_line, p = p))
  if (!length(mats)) return(rep(FALSE, nrow(M)))
  rowSums(do.call(cbind, mats)) >= min(min_rep, length(mats))
}

build_regulon <- function(region_mask) {
  idx <- which(region_mask)
  L <- links[links$region %in% idx, ]
  if (!nrow(L)) return(character(0))
  agg <- tapply(L$score, L$gene, sum)
  agg <- sort(agg, decreasing = TRUE)
  n <- min(length(agg), RG$max_size)
  if (n < RG$min_size) return(character(0))
  names(agg)[seq_len(n)]
}

cat("\n=========== regulons ===========\n")
reg <- list(); regions_used <- list()
for (p in names(PARALOG_LINES)) {
  mask <- active_set(p)
  g <- build_regulon(mask)
  reg[[p]] <- g; regions_used[[p]] <- mask
  cat(sprintf("  %-6s active regions %6d -> regulon %4d genes\n", p, sum(mask), length(g)))
}

# ---- 1. leave-one-line-out ---------------------------------------------------
cat("\n=========== validity 1: leave-one-line-out ===========\n")
cat("Regions built WITHOUT the held-out line; the held-out line's own ChIP is\n")
cat("then tested for enrichment there. Its data never informed the regions.\n\n")
loo <- list()
for (p in names(PARALOG_LINES)) {
  lns <- PARALOG_LINES[[p]]
  lns <- lns[!vapply(lns, function(l) is.null(active_in_line(p, l)), logical(1))]
  for (lh in lns) {
    rest <- setdiff(lns, lh)
    if (!length(rest)) next
    mask <- active_set(p, lines = rest, min_rep = 1)
    jh <- col_for(p, lh)
    if (is.na(jh) || !sum(mask)) next
    inr <- M[mask,  jh]; out <- M[!mask, jh]
    fold <- stats::median(inr) / max(stats::median(out), 1e-9)
    wt <- suppressWarnings(stats::wilcox.test(inr, out, alternative = "greater"))
    # AUC from the rank-sum statistic
    auc <- as.numeric(wt$statistic) / (length(inr) * length(out))
    loo[[length(loo)+1L]] <- data.frame(
      paralog = p, held_out = lh, n_train_lines = length(rest),
      n_regions = sum(mask), median_fold = round(fold, 2),
      auc = round(auc, 3), p = signif(wt$p.value, 3),
      pass = auc >= 0.70 && wt$p.value < 0.05, stringsAsFactors = FALSE)
    cat(sprintf("  %-6s hold out %-8s regions %6d  fold %5.2f  AUC %.3f  p %8.2g  %s\n",
                p, lh, sum(mask), fold, auc, wt$p.value,
                if (auc >= 0.70 && wt$p.value < 0.05) "PASS" else "FAIL"))
  }
}
loo <- do.call(rbind, loo)
v1 <- mean(loo$pass) >= 0.70
cat(sprintf("\nleave-one-line-out: %d/%d held-out lines pass -> %s\n",
            sum(loo$pass), nrow(loo), if (v1) "PASS" else "FAIL"))

# ---- 2. cross-paralog distinctness ------------------------------------------
cat("\n=========== validity 2: cross-paralog distinctness ===========\n")
jac <- function(a, b) length(intersect(a,b)) / max(length(union(a,b)), 1)
J <- outer(names(reg), names(reg), Vectorize(function(a,b) round(jac(reg[[a]], reg[[b]]), 3)))
dimnames(J) <- list(names(reg), names(reg))
print(J)
offd <- J[upper.tri(J)]
v2 <- all(offd <= RV$max_regulon_jaccard)
cat(sprintf("\nmax off-diagonal Jaccard %.3f (limit %.2f) -> %s\n",
            max(offd), RV$max_regulon_jaccard, if (v2) "PASS" else "FAIL"))
if (!v2) cat("Regulons overlap too heavily to be treated as separable programmes.\n")

# ---- 3. benchmark enrichment -------------------------------------------------
cat("\n=========== validity 3: HALLMARK MYC targets ===========\n")
hm <- msigdbr(species = "Homo sapiens", collection = "H")
hm <- hm[hm$gs_name == RV$benchmark_enrichment, ]
hs <- unique(hm$gene_symbol)
universe_genes <- unique(links$gene)
cat("benchmark set: ", RV$benchmark_enrichment, " (", length(hs), " genes; ",
    length(intersect(hs, universe_genes)), " in our linkable universe)\n\n", sep = "")
bench <- list()
for (p in names(reg)) {
  g <- reg[[p]]
  a <- length(intersect(g, hs)); b <- length(setdiff(g, hs))
  c_ <- length(setdiff(intersect(hs, universe_genes), g))
  d <- length(setdiff(universe_genes, union(g, hs)))
  ft <- stats::fisher.test(matrix(c(a,b,c_,d), 2, byrow = TRUE), alternative = "greater")
  bench[[length(bench)+1L]] <- data.frame(
    paralog = p, n_regulon = length(g), n_hallmark_hit = a,
    odds_ratio = round(as.numeric(ft$estimate), 2),
    p = signif(ft$p.value, 3), stringsAsFactors = FALSE)
  cat(sprintf("  %-6s %4d genes, %3d Hallmark MYC targets, OR %5.2f, p %8.2g\n",
              p, length(g), a, ft$estimate, ft$p.value))
}
bench <- do.call(rbind, bench)
v3 <- bench$p[bench$paralog == "MYC"] < 0.05 && bench$odds_ratio[bench$paralog == "MYC"] > 1
cat("\nMYC regulon enriched for known MYC targets: ", if (v3) "PASS" else "FAIL", "\n", sep = "")

# ---- verdict -----------------------------------------------------------------
res <- data.frame(
  check = c("leave_one_line_out","cross_paralog_distinctness","hallmark_myc_benchmark"),
  detail = c(sprintf("%d/%d held-out lines pass", sum(loo$pass), nrow(loo)),
             sprintf("max Jaccard %.3f (limit %.2f)", max(offd), RV$max_regulon_jaccard),
             sprintf("MYC OR %.2f, p %.2g", bench$odds_ratio[bench$paralog=="MYC"],
                     bench$p[bench$paralog=="MYC"])),
  pass = c(v1, v2, v3), stringsAsFactors = FALSE)
write.csv(rbind(
  data.frame(type="loo", loo[, c("paralog","held_out","auc","p","pass")]),
  data.frame(type="bench", paralog=bench$paralog, held_out=NA,
             auc=bench$odds_ratio, p=bench$p, pass=bench$p<0.05)),
  "data/metadata/regulon_validity.csv", row.names = FALSE)

cat("\n=========== REGULON VALIDITY GATE ===========\n")
print(res, row.names = FALSE)
saveRDS(list(regulons = reg, regions = regions_used, loo = loo,
             jaccard = J, benchmark = bench, verdict = res),
        "data/processed/regions/regulons.rds")

md <- c("# M5 — paralog regulons", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        paste0("Top-ranked genes by aggregate link score, capped at ", RG$max_size,
               ". Links are high+moderate tier only, from an H3K27ac activity",
               " proxy (D-027) — not expression."), "",
        "| paralog | active regions | regulon genes |", "|---|---|---|")
for (p in names(reg))
  md <- c(md, sprintf("| %s | %s | %s |", p,
                      format(sum(regions_used[[p]]), big.mark = ","), length(reg[[p]])))
md <- c(md, "", "## Validity gate", "", "| check | detail | result |", "|---|---|---|")
for (i in seq_len(nrow(res)))
  md <- c(md, sprintf("| %s | %s | %s |", res$check[i], res$detail[i],
                      ifelse(res$pass[i], "PASS", "**FAIL**")))
writeLines(md, "results/tables/m5_regulons.md")

cat("\nwrote data/processed/regions/regulons.rds\n")
if (!all(res$pass)) {
  cat("\nRESULT: FAIL — regulons must not be scored in tumours until this passes.\n")
  quit(status = 1)
}
cat("\nRESULT: PASS — regulons cleared for M6 tumour scoring.\n")
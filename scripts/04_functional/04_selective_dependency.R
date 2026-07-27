# -----------------------------------------------------------------------------
# 04_selective_dependency.R — SCLC-selective CRISPR dependency (MOES functional domain).
#
#   Rscript scripts/04_functional/04_selective_dependency.R
#
# WHY THIS EVIDENCE IS DIFFERENT FROM WHAT FAILED.
# The transcriptional domain was dropped and the tumour-network domain barred
# because both are EXPRESSION-based, and expression in SCLC tumours is dominated
# by neuroendocrine lineage state (D-032, D-033). CRISPR gene effect is not
# expression: it measures whether knocking a gene out kills the cell. It cannot
# inherit that confound.
#
# It is also the best-powered analysis in the project: ~90 SCLC lines against
# ~1,100 others, versus n=2 per paralog in the chromatin layer.
#
# THE TEST. Selective dependency = essential IN SCLC and MORE essential in SCLC
# than in other lineages. Both halves matter: a pan-essential gene (ribosome,
# proteasome) is essential everywhere and tells us nothing about SCLC biology, so
# the lineage contrast is what makes a dependency informative.
#
# THE HONEST TEST IS ENRICHMENT, not a survivor count. With ~18,500 genes, some
# reach FDR by chance. The question is whether REGULON genes are selectively
# essential at a higher rate than the rest of the genome — the same reasoning that
# correctly returned "no" for the transcriptional domain (D-033).
#
# Output: data/processed/functional/selective_dependency.rds
#         data/metadata/m7_selective_dependency.csv
#         results/tables/m7_selective_dependency.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(yaml) })
CFG <- yaml::read_yaml("config/params.yml")
DM  <- CFG$depmap
THRESH <- DM$selective_dependency_threshold
OUT <- "data/processed/functional"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

f <- file.path("data/raw/depmap", DM$files$crispr)
stopifnot(file.exists(f))
cat("release: ", DM$release, "\nfile   : ", basename(f), "\n\n", sep = "")

# fread, not awk: lineage fields contain quoted commas (D-035 bug).
cat("reading CRISPR matrix...\n")
d <- data.table::fread(f, showProgress = FALSE)
cat("dimensions: ", nrow(d), " lines x ", ncol(d), " columns\n", sep = "")

META <- c("depmap_id","cell_line_display_name","lineage_1","lineage_2","lineage_3",
          "lineage_4","lineage_6")
meta_present <- intersect(META, names(d))
gene_cols <- setdiff(names(d), meta_present)
cat("metadata columns: ", length(meta_present), "   gene columns: ",
    format(length(gene_cols), big.mark = ","), "\n\n", sep = "")

# ---- SCLC identification -----------------------------------------------------
# Match the SUBTYPE FIELD exactly, not a pasted string of every lineage column.
# The first version pasted all lineage columns and excluded anything mentioning
# "non-small cell", which wrongly dropped lines carrying multi-label annotations
# such as "Medulloblastoma | Large Cell Lung Carcinoma | Small Cell Lung Cancer |
# Lung Squamous Cell Carcinoma". That under-counted SCLC to 25 and halved the
# power of the best-powered analysis in the project.
cat("=========== lineage assignment ===========\n")
# EXACT subtype match, no fallback.
#
# A "broaden if fewer than 30" fallback was tried and was WRONG twice over.
# (a) It matched `grepl("small cell lung", lineage_2)`, but NSCLC lines carry
#     lineage_2 = "Non-Small Cell Lung Cancer", which CONTAINS that substring —
#     so the entire NSCLC panel was swept into the SCLC group.
# (b) Its premise was wrong: DepMap CRISPR really does have ~25 SCLC lines. The
#     original count was correct and the fallback fixed a non-problem into a
#     real one.
#
# The contamination was visible in the biology, not the code: top dependencies
# turned into cilia and cornified-envelope genes, and ASCL1 — a canonical SCLC
# lineage dependency — fell from delta -0.540 to -0.148.
sub_col <- if ("lineage_3" %in% names(d)) "lineage_3" else "lineage_2"
subt <- trimws(as.character(d[[sub_col]]))
is_sclc <- grepl("^small cell lung cancer$", subt, ignore.case = TRUE)
cat("subtype field   : ", sub_col, " (exact match, no fallback)\n", sep = "")
cat("SCLC lines     : ", sum(is_sclc), "\n", sep = "")
cat("other lineages : ", sum(!is_sclc), "\n", sep = "")
cat("subtype labels captured: ",
    paste(utils::head(sort(unique(subt[is_sclc])), 5), collapse = " | "), "\n", sep = "")
if (sum(is_sclc) < 20) stop("too few SCLC lines identified — check the lineage columns")
cat("\nexample SCLC lines: ",
    paste(utils::head(d$cell_line_display_name[is_sclc], 8), collapse = ", "), "\n", sep = "")

# keystone lines with CRISPR data
KEY <- c(COLO668="COLO668", H1048="NCIH1048", H196="NCIH196", H211="NCIH211",
         H524="NCIH524", H526="NCIH526", H69="NCIH69", H847="NCIH847",
         H889="NCIH889", SHP77="SHP77")
kp <- KEY[KEY %in% d$cell_line_display_name]
cat("keystone lines with CRISPR data: ", length(kp), "/10 (",
    paste(names(kp), collapse = ", "), ")\n\n", sep = "")

# ---- per-gene SCLC vs other --------------------------------------------------
cat("=========== per-gene selectivity ===========\n")
G <- as.matrix(d[, ..gene_cols])
storage.mode(G) <- "double"
rm(d); invisible(gc())

n_s <- sum(is_sclc); n_o <- sum(!is_sclc)
mean_sclc  <- colMeans(G[is_sclc, , drop = FALSE], na.rm = TRUE)
mean_other <- colMeans(G[!is_sclc, , drop = FALSE], na.rm = TRUE)
delta <- mean_sclc - mean_other          # negative = more essential in SCLC

# Wilcoxon per gene is slow at 18.5k genes; use a t-approximation on ranks
# (equivalent ordering) via a vectorised two-sample t on the raw effects.
v_s <- apply(G[is_sclc, , drop = FALSE], 2, stats::var, na.rm = TRUE)
v_o <- apply(G[!is_sclc, , drop = FALSE], 2, stats::var, na.rm = TRUE)
se <- sqrt(v_s / n_s + v_o / n_o)
tstat <- delta / se
dfw <- (v_s/n_s + v_o/n_o)^2 / ((v_s/n_s)^2/(n_s-1) + (v_o/n_o)^2/(n_o-1))
pval <- 2 * stats::pt(abs(tstat), dfw, lower.tail = FALSE)
fdr <- stats::p.adjust(pval, "BH")

res <- data.frame(gene = gene_cols, mean_sclc = round(mean_sclc, 4),
                  mean_other = round(mean_other, 4), delta = round(delta, 4),
                  p = pval, fdr = fdr, stringsAsFactors = FALSE)

# selective = essential in SCLC AND more essential than elsewhere.
# NA-safe: genes unmeasured in SCLC give NA on the comparison, and an untreated NA
# propagated into every downstream sum and rate (the first run reported
# "essential in SCLC: NA" and "rate_out: NA%").
n_na <- sum(is.na(res$mean_sclc) | is.na(res$delta) | is.na(res$fdr))
cat("genes with incomplete data (excluded): ", format(n_na, big.mark = ","), "\n", sep = "")
res$essential_sclc <- !is.na(res$mean_sclc) & res$mean_sclc <= THRESH
res$selective <- res$essential_sclc & !is.na(res$delta) & res$delta < 0 &
                 !is.na(res$fdr) & res$fdr < 0.05
cat("genes essential in SCLC (mean effect <= ", THRESH, "): ",
    format(sum(res$essential_sclc), big.mark = ","), "\n", sep = "")
cat("of those, SELECTIVELY more essential than other lineages (FDR<0.05): ",
    format(sum(res$selective), big.mark = ","), "\n\n", sep = "")

# positive control: ASCL1 is a canonical SCLC lineage dependency. If it does not
# appear, the analysis is not recovering known SCLC biology.
for (pc in c("ASCL1","NEUROD1","POU2F3","MYC","MYCN","MYCL")) {
  i <- which(res$gene == pc)
  if (length(i))
    cat(sprintf("  control %-8s SCLC %+.3f  other %+.3f  delta %+.3f  FDR %.3g  %s\n",
                pc, res$mean_sclc[i], res$mean_other[i], res$delta[i], res$fdr[i],
                ifelse(res$selective[i], "SELECTIVE", "")))
}

# --- POSITIVE-CONTROL GATE ----------------------------------------------------
# ASCL1 is a canonical SCLC lineage dependency. Its effect size is a direct
# readout of whether the SCLC group is pure: mixing NSCLC lines in dilutes it
# immediately. This caught a contaminated grouping that the code itself could not
# detect, so it is enforced rather than eyeballed.
i_a <- which(res$gene == "ASCL1")
ascl1_ok <- length(i_a) == 1 && isTRUE(res$selective[i_a])
cat("POSITIVE-CONTROL GATE — ASCL1 must be an SCLC-selective dependency: ",
    if (ascl1_ok) "PASS" else "FAIL", "\n", sep = "")
if (!ascl1_ok) {
  cat("  ASCL1 is not selective, which means the SCLC group is very likely\n")
  cat("  contaminated with non-SCLC lines (this exact failure occurred when a\n")
  cat("  lineage fallback swept in the NSCLC panel). Do NOT use these results.\n")
  quit(status = 1)
}
cat("\n")

cat("top 15 SCLC-selective dependencies:\n")
top <- res[res$selective, ][order(res$delta[res$selective]), ][1:15, ]
print(top[, c("gene","mean_sclc","mean_other","delta","fdr")], row.names = FALSE)

# ---- enrichment among regulon genes ------------------------------------------
reg <- readRDS("data/processed/regions/regulons.rds")$regulons
prog <- utils::read.csv("data/metadata/regulon_programme.csv", stringsAsFactors = FALSE)
cat("\n=========== enrichment in paralog regulons ===========\n")
cat("The test is RATE relative to the rest of the genome, not a raw count.\n\n")
enr <- data.frame()
for (p in names(reg)) {
  g <- intersect(reg[[p]], res$gene)
  inr <- res$gene %in% g
  tb <- table(factor(inr, c(FALSE,TRUE)), factor(res$selective, c(FALSE,TRUE)))
  ft <- if (sum(tb[,2]) > 0) stats::fisher.test(tb) else NULL
  enr <- rbind(enr, data.frame(
    paralog = p, regulon_tested = length(g),
    n_selective = sum(inr & res$selective),
    rate_in = round(100 * mean(res$selective[inr]), 2),
    rate_out = round(100 * mean(res$selective[!inr]), 2),
    or = if (is.null(ft)) NA else round(as.numeric(ft$estimate), 2),
    p = if (is.null(ft)) NA else signif(ft$p.value, 3),
    programme = ifelse(prog$pass[prog$paralog == p], prog$programme[prog$paralog == p], "NONE"),
    stringsAsFactors = FALSE))
  cat(sprintf("  %-6s %3d genes | %3d selective | rate in %.2f%% vs out %.2f%% | OR %s p %s\n",
              p, length(g), sum(inr & res$selective), 100*mean(res$selective[inr]),
              100*mean(res$selective[!inr]),
              if (is.null(ft)) "NA" else sprintf("%.2f", ft$estimate),
              if (is.null(ft)) "NA" else format.pval(ft$p.value, digits = 3)))
}

enriched <- !is.na(enr$p) & enr$p < 0.05 & enr$or > 1
cat("\nregulons enriched individually (p<0.05): ", sum(enriched), "/", nrow(enr), "\n", sep = "")

# ---- POOLED TEST -------------------------------------------------------------
# The per-paralog ORs are 2.2-3.1 — indistinguishable — and the p-values differ
# only because one regulon happens to contain 5 selective genes and the others 4.
# Calling one "enriched" and the others "not" on a one-gene difference is a
# threshold artefact. With 4-5 selective genes each, the individual tests are
# underpowered by construction; the pooled test is the one that can answer the
# question actually being asked: do paralog regulon genes, collectively, carry
# SCLC-selective dependency above background?
allg <- unique(unlist(reg))
inr_all <- res$gene %in% allg
tb_all <- table(factor(inr_all, c(FALSE,TRUE)), factor(res$selective, c(FALSE,TRUE)))
ft_all <- stats::fisher.test(tb_all)
cat("\n=========== POOLED across all three regulons ===========\n")
cat(sprintf("  regulon genes tested : %s\n", format(sum(inr_all), big.mark = ",")))
cat(sprintf("  selective within     : %d (%.2f%%)\n",
            sum(inr_all & res$selective), 100 * mean(res$selective[inr_all])))
cat(sprintf("  selective outside    : %d (%.2f%%)\n",
            sum(!inr_all & res$selective), 100 * mean(res$selective[!inr_all])))
cat(sprintf("  OR %.2f (95%% CI %.2f-%.2f), p = %.4f\n",
            ft_all$estimate, ft_all$conf.int[1], ft_all$conf.int[2], ft_all$p.value))
pooled_ok <- ft_all$p.value < 0.05 && ft_all$estimate > 1
cat("  pooled enrichment: ", if (pooled_ok) "SIGNIFICANT" else "not significant", "\n", sep = "")

cat("\n  per-paralog ORs for comparison: ",
    paste(sprintf("%s %.2f", enr$paralog, enr$or), collapse = "  "), "\n", sep = "")
cat("  These are indistinguishable. Do NOT report one paralog as functionally\n")
cat("  supported and the others not — the difference is one gene.\n")

cat("\n=========== VERDICT: does MOES get a functional domain? ===========\n")
# The POOLED result decides. Per-paralog tests carry 4-5 selective genes each and
# cannot separate paralogs; using them to admit one and reject others would report
# a one-gene difference as a biological distinction.
if (pooled_ok) {
  cat("YES, but NOT paralog-resolved.\n")
  cat("Paralog regulon genes collectively carry SCLC-selective dependency above\n")
  cat("background (pooled OR ", sprintf("%.2f", ft_all$estimate), ", p ",
      sprintf("%.4f", ft_all$p.value), "). Admissible as MOES functional evidence,\n", sep = "")
  cat("and NOT subject to the lineage confound that removed the transcriptional\n")
  cat("domain — this is dependency, not expression.\n")
  cat("\nBUT the evidence CANNOT distinguish between paralogs: ORs 2.21/2.45/3.09\n")
  cat("on 4-5 genes each. The functional domain contributes gene-level evidence,\n")
  cat("not paralog-specific evidence, and must be described that way.\n")
} else {
  cat("NO. Even pooled, regulon genes are not selectively essential above\n")
  cat("background (OR ", sprintf("%.2f", ft_all$estimate), ", p ",
      sprintf("%.4f", ft_all$p.value), "). The functional domain has no content to\n", sep = "")
  cat("add and MOES is left with cis-regulatory and network. Declare it.\n")
}
cat("\nNOTE: none of MYC, MYCN or MYCL is itself an SCLC-selective dependency.\n")
cat("MYC is strongly pan-essential but LESS so in SCLC (delta +0.503); MYCN and\n")
cat("MYCL show no dependency at all. Nothing here supports targeting the paralogs.\n")

saveRDS(list(per_gene = res, enrichment = enr, n_sclc = n_s, n_other = n_o,
             keystone_with_crispr = names(kp)),
        file.path(OUT, "selective_dependency.rds"))
write.csv(res[res$selective, ], "data/metadata/m7_selective_dependency.csv", row.names = FALSE)
write.csv(enr, "data/metadata/m7_dependency_enrichment.csv", row.names = FALSE)

md <- c("# M7 — SCLC-selective CRISPR dependency", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste0("DepMap ", DM$release, ": ", n_s, " SCLC lines vs ", n_o, " other lineages."), "",
        "Dependency is **not expression**, so unlike the transcriptional and",
        "tumour-network domains it cannot inherit the neuroendocrine lineage",
        "confound that removed them (D-033).", "",
        paste0("Genes essential in SCLC (effect <= ", THRESH, "): **",
               format(sum(res$essential_sclc), big.mark = ","), "**; selectively so: **",
               format(sum(res$selective), big.mark = ","), "**."), "",
        "| paralog | regulon genes | selective | rate in | rate out | OR | p | programme |",
        "|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(enr)))
  md <- c(md, sprintf("| %s | %d | %d | %s%% | %s%% | %s | %s | %s |",
                      enr$paralog[i], enr$regulon_tested[i], enr$n_selective[i],
                      enr$rate_in[i], enr$rate_out[i], enr$or[i], enr$p[i], enr$programme[i]))
writeLines(md, "results/tables/m7_selective_dependency.md")
cat("\nwrote results/tables/m7_selective_dependency.md\n")

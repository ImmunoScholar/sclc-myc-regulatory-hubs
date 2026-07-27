# -----------------------------------------------------------------------------
# 06_support_rule.R — does a cross-line support rule rescue agreement with MACS2?
#
#   Rscript scripts/02_regulatory/06_support_rule.R
#
# Calibration (04_) gave a best F1 of only 0.32 against H524's independent MACS2
# peaks. The suspicion is that single-line regions — 65.6% of the universe, and
# only 39.9% externally corroborated — are dragging precision down, and that
# requiring reproducibility across keystone lines is itself a background control:
# noise does not recur at the same locus in independent samples.
#
# This script TESTS that rather than assuming it. It re-scores each candidate
# support rule against the same H524 chr1 reference used for calibration.
#
# It does NOT choose the rule. It reports.
#
# Output: data/metadata/support_rule_comparison.csv
#         results/tables/m5_support_rule.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(yaml)
})
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
RP  <- CFG$regions
CAL_CHROM <- "chr1"

universe <- readRDS("data/processed/regions/universe.rds")
cat("universe: ", format(length(universe), big.mark = ","), " regions\n", sep = "")

# --- reference: H524 external MACS2 peaks, chr1 -------------------------------
lf  <- list.files("data/processed/liftover", pattern = "H524.*\\.hg19\\.rds$", full.names = TRUE)
ref <- GenomicRanges::reduce(do.call(c, unname(lapply(lf, readRDS))),
                             min.gapwidth = RP$merge_gap)
ref <- ref[as.character(seqnames(ref)) == CAL_CHROM]
cat("reference (H524 MACS2, chr1): ", length(ref), " peaks, ",
    format(sum(width(ref)), big.mark = ","), " bp\n\n", sep = "")

genome_bp <- sum(load_hg19_sizes()[ANALYSIS_CHROMS_UCSC])
nl  <- mcols(universe)$n_lines_supporting
ext <- mcols(universe)$n_external_support > 0

rules <- list(
  "all (>=1 line)"                = rep(TRUE, length(universe)),
  ">=2 lines"                     = nl >= 2,
  ">=3 lines"                     = nl >= 3,
  ">=2 lines OR ext-corroborated" = (nl >= 2) | ext,
  ">=1 line AND ext-corroborated" = ext
)

# H524 is one of the nine keystone lines, so restrict the comparison to regions
# our own H524 track supports. Scoring regions that H524 never called against an
# H524 reference would penalise the rule for regions it was never claiming.
h524_supported <- mcols(universe)$atac_H524

res <- list()
for (nm in names(rules)) {
  keep <- rules[[nm]]
  g    <- universe[keep]
  n    <- length(g)
  bp   <- sum(as.numeric(width(g)))

  gc1  <- g[as.character(seqnames(g)) == CAL_CHROM &
            mcols(g)$atac_H524]                       # H524-claimed, chr1
  if (length(gc1)) {
    inter <- sum(width(GenomicRanges::intersect(gc1, ref)))
    prec  <- inter / sum(width(gc1))
    rec   <- inter / sum(width(ref))
    f1    <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
    jac   <- inter / sum(width(GenomicRanges::union(gc1, ref)))
  } else { prec <- rec <- f1 <- jac <- NA_real_ }

  res[[length(res) + 1L]] <- data.frame(
    rule = nm, n_regions = n,
    pct_of_all = round(100 * n / length(universe), 1),
    total_mb = round(bp / 1e6, 1),
    pct_genome = round(100 * bp / genome_bp, 3),
    pct_ext_corrob = round(100 * mean(mcols(g)$n_external_support > 0), 1),
    chr1_h524_n = length(gc1),
    precision = round(prec, 4), recall = round(rec, 4),
    f1 = round(f1, 4), jaccard = round(jac, 4),
    stringsAsFactors = FALSE)
}

cmp <- do.call(rbind, res)
write.csv(cmp, "data/metadata/support_rule_comparison.csv", row.names = FALSE)

cat("=========== support rule comparison ===========\n")
print(cmp[, c("rule","n_regions","pct_genome","pct_ext_corrob","precision","recall","f1","jaccard")],
      row.names = FALSE)

cat("\nreference for scale: calibration on the unfiltered set gave F1 = 0.320\n")

# Which rule maximises agreement?
best <- cmp[which.max(cmp$f1), ]
cat("\nbest F1: ", best$rule, "  (F1 ", best$f1, ", precision ", best$precision,
    ", recall ", best$recall, ")\n", sep = "")

md <- c("# M5 — accessible-region support rule", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "Each candidate rule scored against H524's independent MACS2 peak set on",
        "chr1, restricted to regions our own H524 ATAC track supports.", "",
        "| rule | regions | % genome | ext. corroborated | precision | recall | F1 | Jaccard |",
        "|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(cmp)))
  md <- c(md, sprintf("| %s | %s | %s | %s%% | %s | %s | **%s** | %s |",
                      cmp$rule[i], format(cmp$n_regions[i], big.mark = ","),
                      cmp$pct_genome[i], cmp$pct_ext_corrob[i],
                      cmp$precision[i], cmp$recall[i], cmp$f1[i], cmp$jaccard[i]))
md <- c(md, "",
        "Agreement with MACS2 is a calibration target, not ground truth (D-023).",
        "Recall is bounded by the fact that the reference is an independent ATAC",
        "experiment in transduced rather than parental H524; published replicate",
        "Jaccards between independent ATAC experiments run roughly 0.3-0.6.")
writeLines(md, "results/tables/m5_support_rule.md")

cat("\nwrote data/metadata/support_rule_comparison.csv\n")
cat("wrote results/tables/m5_support_rule.md\n")
cat("\nRESULT: comparison complete. The rule is NOT yet decided.\n")

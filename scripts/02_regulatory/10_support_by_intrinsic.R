# -----------------------------------------------------------------------------
# 10_support_by_intrinsic.R — choose the support rule using INTRINSIC quality.
#
#   Rscript scripts/02_regulatory/10_support_by_intrinsic.R
#
# Step 06 asked this question against the wrong yardstick (agreement with another
# laboratory's peak calls). Same question, right metric.
#
# The problem: individual ATAC tracks are good (H524/chr1 TSS enrichment 8.9x),
# but the nine-line UNION gives 6.99% of the genome and only 5.3% promoter-proximal
# regions — implying TSS enrichment near 3.5x. Real ATAC peak sets are 20-30%
# promoter-proximal. With 65% of regions supported by a single line, the union is
# accumulating nine independent tracks' worth of noise.
#
# So: does requiring cross-line reproducibility restore intrinsic quality? Noise
# should not recur at the same locus in independent samples, so a support filter
# ought to raise TSS enrichment back toward the single-line value. If it does not,
# the problem is elsewhere.
#
# Reports TSS enrichment, promoter fraction and external corroboration at every
# support level, for each threshold. Chooses nothing.
#
# Output: data/metadata/support_by_intrinsic.csv
#         results/tables/m5_support_intrinsic.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
})
source("R/genome_utils.R")

TSS_WIN <- 1000L
sizes     <- load_hg19_sizes()
genome_bp <- sum(sizes[ANALYSIS_CHROMS_UCSC])

tx  <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tss <- promoters(tx, upstream = 0, downstream = 1)
tss <- tss[as.character(seqnames(tss)) %in% ANALYSIS_CHROMS_UCSC]
tss_win <- GenomicRanges::reduce(resize(tss, width = 2 * TSS_WIN, fix = "center"))
tss_frac <- sum(as.numeric(width(tss_win))) / genome_bp
cat("TSS +/-", TSS_WIN, " covers ", sprintf("%.2f%%", 100 * tss_frac),
    " of the analysis genome\n", sep = "")
cat("random-expectation promoter overlap rate for small regions: ~",
    sprintf("%.1f%%", 100 * tss_frac), "\n\n", sep = "")

files <- list.files("data/processed/regions", pattern = "^universe_x[0-9.]+\\.rds$",
                    full.names = TRUE)
stopifnot(length(files) > 0)

res <- list()
for (f in files) {
  m <- as.numeric(sub("^universe_x([0-9.]+)\\.rds$", "\\1", basename(f)))
  u <- readRDS(f)
  cat("=== threshold x", m, " (", format(length(u), big.mark = ","), " regions) ===\n", sep = "")
  for (k in 1:9) {
    g <- u[mcols(u)$n_lines_supporting >= k]
    if (!length(g)) next
    prom <- mean(mcols(g)$is_promoter)
    res[[length(res) + 1L]] <- data.frame(
      multiple = m, min_lines = k, n_regions = length(g),
      total_mb = round(sum(as.numeric(width(g))) / 1e6, 1),
      pct_genome = round(100 * sum(as.numeric(width(g))) / genome_bp, 3),
      pct_promoter = round(100 * prom, 1),
      tss_enrichment = round(prom / tss_frac, 2),
      pct_ext_corrob = round(100 * mean(mcols(g)$n_external_support > 0), 1),
      median_width = median(width(g)),
      stringsAsFactors = FALSE)
    cat(sprintf("  >=%d lines: %8s regions  %5.2f%% genome  prom %4.1f%%  TSS %5.2fx  ext %4.1f%%\n",
                k, format(length(g), big.mark = ","),
                100 * sum(as.numeric(width(g))) / genome_bp,
                100 * prom, prom / tss_frac,
                100 * mean(mcols(g)$n_external_support > 0)))
  }
  cat("\n")
}

v <- do.call(rbind, res)
write.csv(v, "data/metadata/support_by_intrinsic.csv", row.names = FALSE)

cat("=============================================================\n")
cat("TARGETS for a credible accessible-region universe:\n")
cat("  genome fraction  1-3%   (accessible chromatin)\n")
cat("  promoter frac    15-30% (typical ATAC peak sets)\n")
cat("  TSS enrichment   >5x, ideally >8x\n")
cat("=============================================================\n\n")

v$ok <- v$pct_genome >= 0.8 & v$pct_genome <= 3.5 &
        v$pct_promoter >= 12 & v$tss_enrichment >= 5
good <- v[v$ok, ]
if (nrow(good)) {
  cat("Combinations meeting all three targets:\n")
  print(good[order(-good$n_regions),
             c("multiple","min_lines","n_regions","pct_genome","pct_promoter",
               "tss_enrichment","pct_ext_corrob")], row.names = FALSE)
  best <- good[which.max(good$n_regions), ]
  cat("\nLargest such universe (maximum coverage subject to quality):\n")
  cat("  x", best$multiple, ", >=", best$min_lines, " lines  ->  ",
      format(best$n_regions, big.mark = ","), " regions, ",
      best$pct_genome, "% genome, ", best$pct_promoter, "% promoter, TSS ",
      best$tss_enrichment, "x\n", sep = "")
} else {
  cat("NO combination meets all three targets. The universe cannot be rescued by\n")
  cat("a support rule and the approach needs reconsidering.\n")
}

md <- c("# M5 — support rule chosen on intrinsic quality", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "Individual ATAC tracks are good (H524/chr1 TSS enrichment 8.9x) but the",
        "nine-line union dilutes to ~3.5x because 65% of regions are single-line.",
        "This tests whether cross-line reproducibility restores intrinsic quality.", "",
        paste0("Random-expectation promoter overlap: ",
               sprintf("%.1f%%", 100 * tss_frac), "."), "",
        "| x mean | >= lines | regions | % genome | promoter % | TSS enrich | ext corrob % |",
        "|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(v)))
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | **%s** | %s |",
                      v$multiple[i], v$min_lines[i],
                      format(v$n_regions[i], big.mark = ","), v$pct_genome[i],
                      v$pct_promoter[i], v$tss_enrichment[i], v$pct_ext_corrob[i]))
writeLines(md, "results/tables/m5_support_intrinsic.md")
cat("\nwrote data/metadata/support_by_intrinsic.csv\n")
cat("wrote results/tables/m5_support_intrinsic.md\n")

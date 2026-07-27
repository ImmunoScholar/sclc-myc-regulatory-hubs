# -----------------------------------------------------------------------------
# 08_intrinsic_validation.R — validate the ATAC universe on its OWN merits.
#
#   Rscript scripts/02_regulatory/08_intrinsic_validation.R [--chrom chr1]
#
# WHY THIS REPLACES THE EXTERNAL CALIBRATION
#
# Steps 04-07 tuned the ATAC threshold to agree with GSE269424's H524 MACS2 peaks.
# That target was wrong. It is an independent experiment — different lab, protocol
# and pipeline, in EGFP-transduced rather than parental H524 — and cross-lab ATAC
# concordance is Jaccard 0.3-0.5 at best. Optimising against it was measuring
# inter-laboratory reproducibility, not whether our grid is fit for purpose.
#
# The universe exists to be a quantification grid for MYC and H3K27ac ChIP signal
# in THESE nine lines. So it is validated on properties intrinsic to these data:
#
#   1. TSS ENRICHMENT — the standard ATAC quality metric. Real accessible regions
#      concentrate at promoters. Computed as observed/expected against the
#      genome-wide expectation for randomly placed regions of the same total width.
#
#   2. MATCHED-CELL H3K27ac SIGNAL — fold-enrichment of H3K27ac pileup inside the
#      regions versus the chromosome background, using the SAME cell line. Matched
#      samples, and directly relevant because active regions require H3K27ac.
#      Uses raw signal, so it does not depend on thresholding H3K27ac too.
#
#   3. INCREMENTAL ENRICHMENT — the stopping rule. As the threshold loosens, are
#      the regions NEWLY admitted still TSS-enriched and H3K27ac-enriched? When the
#      marginal additions stop looking like regulatory elements, loosening further
#      is adding noise. This gives an optimum without any external peak set.
#
# Output: data/metadata/intrinsic_validation.csv
#         results/tables/m5_intrinsic_validation.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(yaml)
})
source("R/genome_utils.R")

args  <- commandArgs(trailingOnly = TRUE)
CHROM <- if (length(args) >= 2 && args[1] == "--chrom") args[2] else "chr1"
CHROM_ENS <- sub("^chr", "", CHROM)

CFG <- yaml::read_yaml("config/params.yml")
RP  <- CFG$regions
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

CAL_LINE  <- "H524"
TSS_WIN   <- 1000L
MULTIPLES <- c(1.5, 2, 2.5, 3, 4, 5, 6, 8)

sizes  <- load_hg19_sizes()
chr_bp <- sizes[[CHROM]]

# --- TSS reference ------------------------------------------------------------
tx  <- TxDb.Hsapiens.UCSC.hg19.knownGene
tss <- promoters(genes(tx), upstream = 0, downstream = 1)
tss <- tss[as.character(seqnames(tss)) == CHROM]
tss_win <- GenomicRanges::reduce(resize(tss, width = 2 * TSS_WIN, fix = "center"))
tss_bp  <- sum(width(tss_win))
cat("TSS on ", CHROM, ": ", length(tss), " genes, +/-", TSS_WIN,
    " covers ", format(tss_bp, big.mark = ","), " bp (",
    sprintf("%.2f%%", 100 * tss_bp / chr_bp), " of the chromosome)\n\n", sep = "")

man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
man$dest <- file.path("data/raw", man$dataset_id, man$file_name)
atac_f <- man$dest[grepl(sprintf("GSM[0-9]+_%s_treat_pileup", CAL_LINE), man$file_name)][1]
k27_f  <- man$dest[grepl(sprintf("GSM[0-9]+_H3K27Ac_%s_", CAL_LINE), man$file_name)][1]
stopifnot(file.exists(atac_f), file.exists(k27_f))
cat("ATAC   : ", basename(atac_f), "\n", sep = "")
cat("H3K27ac: ", basename(k27_f),  "\n\n", sep = "")

hist_mean <- function(line) {
  raw <- readLines(file.path("data/processed/atac_profile", paste0(line, ".hist.tsv")))
  h <- utils::read.delim(text = paste(raw[!grepl("^#", raw)], collapse = "\n"),
                         header = FALSE, col.names = c("bin","bp"))
  h$signal <- h$bin / 2; cov <- h[h$signal > 0, ]
  sum(cov$signal * cov$bp) / sum(cov$bp)
}
gm <- hist_mean(CAL_LINE)
cat("ATAC mean signal over covered bp: ", round(gm, 3), "\n\n", sep = "")

regions_at <- function(thr) {
  awk <- sprintf("awk -F'\\t' '$1==\"%s\" && $4+0 >= %f { print $2\"\\t\"$3 }'", CHROM_ENS, thr)
  con <- pipe(paste("zcat", shQuote(atac_f), "|", awk), "r"); on.exit(close(con))
  d <- try(utils::read.delim(con, header = FALSE, col.names = c("start","end")), silent = TRUE)
  if (inherits(d, "try-error") || !nrow(d)) return(GRanges())
  gr <- GenomicRanges::reduce(GRanges(CHROM, IRanges(d$start + 1L, d$end)),
                              min.gapwidth = RP$merge_gap)
  gr[width(gr) >= RP$min_width & width(gr) <= RP$max_width]
}

# H3K27ac signal, binned once, then averaged over any region set.
BIN <- 1000L
cat("binning H3K27ac signal (one pass)...\n")
awk_bin <- sprintf(
  "awk -F'\\t' '$1==\"%s\" { w=$3-$2; v=$4+0; b=int($2/%d); s[b]+=v*w } END { for (k in s) print k\"\\t\"s[k] }'",
  CHROM_ENS, BIN)
k27 <- utils::read.delim(pipe(paste("zcat", shQuote(k27_f), "|", awk_bin)),
                         header = FALSE, col.names = c("bin","sumsw"))
k27_vec <- numeric(ceiling(chr_bp / BIN) + 1L)
k27_vec[k27$bin + 1L] <- k27$sumsw / BIN
k27_bg <- sum(k27$sumsw) / chr_bp
cat("H3K27ac chromosome-wide mean signal: ", round(k27_bg, 4), "\n\n", sep = "")

k27_fold <- function(gr) {
  if (!length(gr)) return(NA_real_)
  b1 <- floor(start(gr) / BIN) + 1L
  b2 <- floor(end(gr)   / BIN) + 1L
  v  <- vapply(seq_along(gr), function(i) mean(k27_vec[b1[i]:b2[i]]), numeric(1))
  mean(v, na.rm = TRUE) / k27_bg
}

tss_enrich <- function(gr) {
  if (!length(gr)) return(NA_real_)
  hit <- mean(IRanges::overlapsAny(gr, tss_win))
  exp <- tss_bp / chr_bp          # expected hit rate for randomly placed regions
  hit / exp
}

res  <- list()
prev <- NULL
for (m in rev(MULTIPLES)) {            # strict -> loose, so "new" regions are well defined
  gr <- regions_at(m * gm)
  if (!length(gr)) next
  newg <- if (is.null(prev)) gr else IRanges::subsetByOverlaps(gr, prev, invert = TRUE)
  res[[length(res) + 1L]] <- data.frame(
    multiple = m, threshold = round(m * gm, 3),
    n_regions = length(gr), mb = round(sum(as.numeric(width(gr))) / 1e6, 2),
    pct_chrom = round(100 * sum(as.numeric(width(gr))) / chr_bp, 3),
    tss_enrich = round(tss_enrich(gr), 2),
    k27_fold   = round(k27_fold(gr), 2),
    n_new = length(newg),
    tss_enrich_new = round(tss_enrich(newg), 2),
    k27_fold_new   = round(k27_fold(newg), 2),
    stringsAsFactors = FALSE)
  prev <- gr
  cat(sprintf("  x%-4g n=%6d  TSS=%5.2f  K27=%5.2f   | new n=%6d TSS=%5.2f K27=%5.2f\n",
              m, length(gr), tss_enrich(gr), k27_fold(gr),
              length(newg), tss_enrich(newg), k27_fold(newg)))
}

v <- do.call(rbind, res)
v <- v[order(v$multiple), ]
write.csv(v, "data/metadata/intrinsic_validation.csv", row.names = FALSE)

cat("\n=========== intrinsic validation (", CAL_LINE, ", ", CHROM, ") ===========\n", sep = "")
print(v, row.names = FALSE)

cat("\nHow to read this:\n")
cat("  tss_enrich     : fold over random expectation. Real ATAC gives >5, good data >10.\n")
cat("  k27_fold       : H3K27ac signal in regions vs chromosome background, SAME cells.\n")
cat("  *_new          : the SAME metrics for regions admitted by loosening one step.\n")
cat("                   The stopping rule: loosen while new regions still look\n")
cat("                   regulatory; stop when they stop.\n")

ok <- v[!is.na(v$tss_enrich_new) & v$tss_enrich_new >= 2 & v$k27_fold_new >= 1.2, ]
if (nrow(ok)) {
  pick <- ok[which.min(ok$multiple), ]
  cat("\nLoosest threshold whose NEW regions remain enriched (TSS>=2, K27>=1.2):\n")
  cat("  x", pick$multiple, "  ->  ", format(pick$n_regions, big.mark = ","),
      " regions, ", pick$pct_chrom, "% of ", CHROM,
      ", TSS ", pick$tss_enrich, "x, K27 ", pick$k27_fold, "x\n", sep = "")
} else {
  cat("\nNo threshold has enriched marginal regions — the signal may be too shallow.\n")
}

md <- c("# M5 — intrinsic validation of the ATAC universe", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        paste0("Line **", CAL_LINE, "**, chromosome **", CHROM, "**."), "",
        "Validated on properties intrinsic to these data rather than on agreement",
        "with another laboratory's peak calls. Prior calibration against GSE269424",
        "was measuring inter-laboratory ATAC reproducibility, not fitness for purpose.", "",
        "| x mean | regions | % chrom | TSS enrich | H3K27ac fold | new regions | new TSS | new K27 |",
        "|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(v)))
  md <- c(md, sprintf("| %s | %s | %s | **%s** | **%s** | %s | %s | %s |",
                      v$multiple[i], format(v$n_regions[i], big.mark = ","), v$pct_chrom[i],
                      v$tss_enrich[i], v$k27_fold[i], format(v$n_new[i], big.mark = ","),
                      v$tss_enrich_new[i], v$k27_fold_new[i]))
writeLines(md, "results/tables/m5_intrinsic_validation.md")
cat("\nwrote data/metadata/intrinsic_validation.csv\n")
cat("wrote results/tables/m5_intrinsic_validation.md\n")

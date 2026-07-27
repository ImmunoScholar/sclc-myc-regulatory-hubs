# -----------------------------------------------------------------------------
# 04_calibrate_threshold.R — fit the ATAC signal threshold against real peak calls.
#
#   Rscript scripts/02_regulatory/04_calibrate_threshold.R
#
# THE POINT OF THIS SCRIPT
#
# The keystone ATAC ships as bedGraph signal with no control track, so regions
# must be called by thresholding. A hand-picked threshold would be exactly the
# kind of free parameter that made the original M5 gate circular (D-020).
#
# H524 is the one cell line with BOTH keystone ATAC signal AND an independent
# MACS2 peak set (GSE269424, lifted to hg19). So the threshold is chosen by
# fitting it to reproduce a real peak caller IN THE SAME CELLS, and the fit
# quality is reported rather than assumed.
#
# Calibration runs on chr1 only — ~8% of the genome and several thousand peaks,
# which is ample to locate an optimum and keeps the sweep to minutes. The chosen
# threshold is then verified genome-wide in 05_.
#
# Honest limits: agreement with MACS2 is a calibration target, not ground truth.
# MACS2 with input has information we do not have. A high Jaccard means our
# regions behave like real peak calls; it does not make them peak calls.
#
# Output: data/metadata/threshold_calibration.csv
#         results/tables/m5_threshold_calibration.md
#         config note: chosen value written to data/metadata/chosen_threshold.txt
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(yaml)
})
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
RP  <- CFG$regions
CAL_CHROM_ENS <- "1"        # Ensembl naming in the keystone files
CAL_CHROM_UCSC <- "chr1"

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# --- reference: external MACS2 peaks for H524, already lifted to hg19 ---------
lift_files <- list.files("data/processed/liftover", pattern = "H524.*\\.hg19\\.rds$",
                         full.names = TRUE)
stopifnot(length(lift_files) >= 1)
ref <- GenomicRanges::reduce(do.call(c, unname(lapply(lift_files, readRDS))),
                             min.gapwidth = RP$merge_gap)
ref <- ref[as.character(seqnames(ref)) == CAL_CHROM_UCSC]
cat("reference MACS2 peaks (H524, ", CAL_CHROM_UCSC, "): ", length(ref), "\n", sep = "")
cat("reference bp: ", format(sum(width(ref)), big.mark = ","), "\n\n", sep = "")

# --- candidate thresholds from the H524 signal histogram ----------------------
hist_file <- "data/processed/atac_profile/H524.hist.tsv"
stopifnot(file.exists(hist_file))
raw <- readLines(hist_file)
h <- utils::read.delim(text = paste(raw[!grepl("^#", raw)], collapse = "\n"),
                       header = FALSE, col.names = c("bin", "bp"))
h$signal <- h$bin / 2
h <- h[order(h$signal), ]
cov <- h[h$signal > 0, ]
mean_cov <- sum(cov$signal * cov$bp) / sum(cov$bp)
cat("H524 mean signal over covered bp: ", round(mean_cov, 3), "\n", sep = "")

# Thresholds as multiples of mean covered signal — scale-free, so the same
# multiples transfer sensibly to the other eight lines despite differing depth.
MULTIPLES <- c(2, 3, 4, 5, 6, 8, 10, 12, 15, 20)
cands <- round(mean_cov * MULTIPLES, 3)
cat("candidate thresholds: ", paste(cands, collapse = ", "), "\n\n", sep = "")

# --- extract chr1 intervals above a threshold, streaming ----------------------
atac_file <- {
  man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
  p <- file.path("data/raw", man$dataset_id, man$file_name)
  p[grepl("GSM7230515_H524_treat_pileup", man$file_name)][1]
}
stopifnot(!is.na(atac_file), file.exists(atac_file))

regions_at <- function(thr) {
  awk <- sprintf(
    "awk -F'\\t' '$1==\"%s\" && $4+0 >= %f { print $2\"\\t\"$3 }'",
    CAL_CHROM_ENS, thr)
  con <- pipe(paste("zcat", shQuote(atac_file), "|", awk), "r")
  on.exit(close(con))
  d <- try(utils::read.delim(con, header = FALSE, col.names = c("start", "end")),
           silent = TRUE)
  if (inherits(d, "try-error") || !nrow(d)) return(GRanges())
  gr <- GRanges(CAL_CHROM_UCSC, IRanges(d$start + 1L, d$end))
  gr <- GenomicRanges::reduce(gr, min.gapwidth = RP$merge_gap)
  gr[width(gr) >= RP$min_width & width(gr) <= RP$max_width]
}

jaccard <- function(a, b) {
  if (!length(a) || !length(b)) return(0)
  i <- sum(width(GenomicRanges::intersect(a, b)))
  u <- sum(width(GenomicRanges::union(a, b)))
  if (u == 0) 0 else i / u
}

res <- list()
for (k in seq_along(cands)) {
  thr <- cands[k]
  t0  <- Sys.time()
  gr  <- regions_at(thr)
  if (!length(gr)) { cat(sprintf("  x%-3g thr=%-8.2f no regions\n", MULTIPLES[k], thr)); next }

  # bp-level precision/recall against the MACS2 reference
  inter <- sum(width(GenomicRanges::intersect(gr, ref)))
  prec  <- inter / sum(width(gr))
  rec   <- inter / sum(width(ref))
  f1    <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  jac   <- jaccard(gr, ref)

  res[[length(res) + 1L]] <- data.frame(
    multiple = MULTIPLES[k], threshold = thr,
    n_regions = length(gr), median_width = median(width(gr)),
    total_mb = round(sum(as.numeric(width(gr))) / 1e6, 2),
    precision = round(prec, 4), recall = round(rec, 4),
    f1 = round(f1, 4), jaccard = round(jac, 4),
    stringsAsFactors = FALSE)

  cat(sprintf("  x%-3g thr=%-8.2f n=%6d  prec=%.3f rec=%.3f F1=%.3f J=%.3f  (%.1fs)\n",
              MULTIPLES[k], thr, length(gr), prec, rec, f1, jac,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

cal <- do.call(rbind, res)
write.csv(cal, "data/metadata/threshold_calibration.csv", row.names = FALSE)

best <- cal[which.max(cal$f1), ]
cat("\n=========== calibration result ===========\n")
print(cal, row.names = FALSE)
cat("\nbest by F1: multiple x", best$multiple, "  threshold ", best$threshold,
    "  (F1 ", best$f1, ", Jaccard ", best$jaccard, ")\n", sep = "")

writeLines(as.character(best$multiple), "data/metadata/chosen_threshold.txt")

md <- c("# M5 — ATAC signal threshold calibration", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        paste0("Calibrated on **H524 / ", CAL_CHROM_UCSC,
               "** against its independent MACS2 peak set (GSE269424, lifted to hg19)."),
        paste0("Reference peaks: ", length(ref), " covering ",
               format(sum(width(ref)), big.mark = ","), " bp."), "",
        "Thresholds are expressed as multiples of the line's mean signal over",
        "covered bases, so the same multiple transfers across lines of differing depth.", "",
        "| x mean | threshold | regions | median width | Mb | precision | recall | F1 | Jaccard |",
        "|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(cal)))
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | **%s** | %s |",
                      cal$multiple[i], cal$threshold[i],
                      format(cal$n_regions[i], big.mark = ","), cal$median_width[i],
                      cal$total_mb[i], cal$precision[i], cal$recall[i],
                      cal$f1[i], cal$jaccard[i]))
md <- c(md, "",
        paste0("**Chosen: x", best$multiple, " mean covered signal** (F1 = ", best$f1, ")."), "",
        "Agreement with MACS2 is a calibration target, not ground truth. MACS2 with",
        "an input control uses information these deposits do not contain. A high F1",
        "means our regions behave like real peak calls in matched cells; it does not",
        "make them peak calls, and the residual weakness is recorded in D-023.")
writeLines(md, "results/tables/m5_threshold_calibration.md")

cat("\nwrote data/metadata/threshold_calibration.csv\n")
cat("wrote results/tables/m5_threshold_calibration.md\n")
cat("wrote data/metadata/chosen_threshold.txt\n")
cat("\nNext: Rscript scripts/02_regulatory/05_derive_regions.R\n")

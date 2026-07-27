# -----------------------------------------------------------------------------
# 01_liftover_intervals.R — lift hg38 interval files to hg19, with loss accounting.
#
#   Rscript scripts/02_regulatory/01_liftover_intervals.R
#
# Policy D-014: intervals lift, continuous signal never does. Only interval files
# are handled here. The hg38 bigWigs (GSE281523, GSE281524) are deliberately NOT
# lifted — regions are called in their native build first, later.
#
# Every lift reports its loss rate. A lift is not a free coordinate change: the
# builds differ by real sequence, and losses are non-uniform across the genome.
# A lifted region set is never presented as equivalent to a natively-hg19 one.
#
# Inputs  : data/raw/GSE269424/*_EGFP_*_peaks.narrowPeak.gz   (hg38)
#           data/raw/GSE256345/*_Peaks_merge500.bed.gz        (hg38)
# Outputs : data/processed/liftover/<dataset>__<file>.hg19.rds
#           data/metadata/liftover_report.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(yaml)
})
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
OUT <- "data/processed/liftover"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CHAIN_GZ <- "data/raw/ucsc_chain_hg38_to_hg19/hg38ToHg19.over.chain.gz"
stopifnot(file.exists(CHAIN_GZ))

# rtracklayer::import.chain needs an uncompressed file.
chain_file <- file.path(tempdir(), "hg38ToHg19.over.chain")
if (!file.exists(chain_file)) {
  system2("gunzip", c("-c", shQuote(CHAIN_GZ)), stdout = chain_file)
}
chain <- import.chain(chain_file)
cat("chain loaded:", length(chain), "chromosome chains\n\n")

# --- what to lift -------------------------------------------------------------
man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
man$dest <- file.path("data/raw", man$dataset_id, man$file_name)

to_lift <- man[
  file.exists(man$dest) &
  as.logical(man$liftover_required) %in% TRUE &
  grepl("\\.(narrowPeak|bed)\\.gz$", man$file_name, ignore.case = TRUE), ]

cat("interval files requiring lift:", nrow(to_lift), "\n")
if (!nrow(to_lift)) stop("nothing to lift — check the manifest")

read_intervals <- function(path) {
  d <- utils::read.delim(gzfile(path), header = FALSE, stringsAsFactors = FALSE,
                         comment.char = "#", quote = "")
  # BED/narrowPeak: chrom, start(0-based), end
  gr <- GRanges(
    seqnames = to_ucsc_seqnames(d[[1]]),          # chain files are UCSC-named
    ranges   = IRanges(start = as.integer(d[[2]]) + 1L,   # BED is 0-based
                       end   = as.integer(d[[3]]))
  )
  if (ncol(d) >= 7) mcols(gr)$signal <- suppressWarnings(as.numeric(d[[7]]))
  gr
}

report <- list()

for (i in seq_len(nrow(to_lift))) {
  r  <- to_lift[i, ]
  gr <- tryCatch(read_intervals(r$dest), error = function(e) NULL)
  if (is.null(gr) || !length(gr)) {
    cat("  SKIP (unreadable):", r$file_name, "\n"); next
  }

  # Restrict to analysis chromosomes before lifting; scaffolds are out of scope
  # and would inflate the apparent loss rate with regions we never wanted.
  gr <- gr[as.character(seqnames(gr)) %in% ANALYSIS_CHROMS_UCSC]
  n_in <- length(gr)
  if (!n_in) { cat("  SKIP (no analysis chroms):", r$file_name, "\n"); next }

  lifted_list <- liftOver(gr, chain)
  n_mapped_any <- sum(lengths(lifted_list) > 0)

  # An interval that fragments into several hg19 pieces is NOT cleanly lifted.
  # Keeping fragments would silently inflate region counts and distort widths, so
  # only 1:1 mappings are retained and the split rate is reported.
  n_split <- sum(lengths(lifted_list) > 1)
  keep    <- lengths(lifted_list) == 1
  lifted  <- unlist(lifted_list[keep])

  # Width sanity: a 1:1 lift should not change width much.
  w_in  <- width(gr)[keep]
  w_out <- width(lifted)
  wr    <- w_out / w_in
  n_wdist <- sum(wr < 0.9 | wr > 1.1)
  lifted  <- lifted[wr >= 0.9 & wr <= 1.1]

  n_out <- length(lifted)
  loss  <- 1 - n_out / n_in

  outfile <- file.path(OUT, paste0(r$dataset_id, "__",
                                   sub("\\.gz$", "", r$file_name), ".hg19.rds"))
  saveRDS(lifted, outfile)

  report[[length(report) + 1L]] <- data.frame(
    dataset = r$dataset_id, file = r$file_name,
    n_input = n_in, n_unmapped = n_in - n_mapped_any, n_split = n_split,
    n_width_distorted = n_wdist, n_retained = n_out,
    loss_rate = round(loss, 4),
    out = basename(outfile), stringsAsFactors = FALSE)

  cat(sprintf("  %-46s %7d -> %7d  (loss %5.2f%%, split %d, width-distorted %d)\n",
              substr(r$file_name, 1, 46), n_in, n_out, 100 * loss, n_split, n_wdist))
}

rep <- do.call(rbind, report)
write.csv(rep, "data/metadata/liftover_report.csv", row.names = FALSE)

cat("\n=========== liftOver summary ===========\n")
cat("files lifted   : ", nrow(rep), "\n", sep = "")
cat("intervals in   : ", sum(rep$n_input), "\n", sep = "")
cat("intervals out  : ", sum(rep$n_retained), "\n", sep = "")
cat("overall loss   : ", sprintf("%.2f%%", 100 * (1 - sum(rep$n_retained) / sum(rep$n_input))), "\n", sep = "")
cat("wrote data/metadata/liftover_report.csv\n")

# A lift losing an unreasonable share of intervals means something is wrong with
# the chain, the build assumption, or the input — not something to push through.
MAX_LOSS <- 0.10
bad <- rep[rep$loss_rate > MAX_LOSS, ]
if (nrow(bad)) {
  cat("\nFILES EXCEEDING ", 100 * MAX_LOSS, "% LOSS:\n", sep = "")
  for (i in seq_len(nrow(bad)))
    cat(sprintf("  %s : %.2f%%\n", bad$file[i], 100 * bad$loss_rate[i]))
  cat("\nRESULT: FAIL — investigate before building the region universe.\n")
  quit(status = 1)
}
cat("\nRESULT: PASS — all lifts within the ", 100 * MAX_LOSS, "% loss tolerance.\n", sep = "")
cat("Next: Rscript scripts/02_regulatory/02_consensus_regions.R\n")

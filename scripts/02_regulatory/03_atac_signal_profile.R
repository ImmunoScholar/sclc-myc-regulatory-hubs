# -----------------------------------------------------------------------------
# 03_atac_signal_profile.R — bp-weighted signal histogram for each keystone ATAC track.
#
#   Rscript scripts/02_regulatory/03_atac_signal_profile.R
#
# One streaming pass per file. The histogram is the point: once we know how many
# base pairs sit at each signal level, ANY threshold can be evaluated instantly
# without re-reading 45 million lines. Nine files streamed once here replaces
# nine files streamed once per candidate threshold later.
#
# Memory: nothing is loaded into R. awk accumulates a small associative array and
# emits a few hundred rows per file, which matters on a 10 GB machine (R-13).
#
# Sequence names in these files are ENSEMBL ("20", not "chr20") — R-16. The awk
# filter is written in that convention deliberately; do not "fix" it to chr*.
#
# Output: data/processed/atac_profile/<line>.hist.tsv
#         data/metadata/atac_signal_profile.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml) })
source("R/genome_utils.R")

OUT <- "data/processed/atac_profile"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
man$dest <- file.path("data/raw", man$dataset_id, man$file_name)

# Keystone ATAC = GSE230649 bedGraphs with NO target token in the file name.
ks <- man[man$dataset_id == "GSE230649" & file.exists(man$dest), ]
is_atac <- !grepl("^GSM[0-9]+_(H3K27Ac|MYC|MYCN|MYCL1)_", ks$file_name)
atac <- ks[is_atac, ]
atac$line <- sub("^GSM[0-9]+_([^_]+)_treat_pileup.*$", "\\1", atac$file_name)

cat("keystone ATAC tracks:", nrow(atac), "\n")
cat("lines:", paste(atac$line, collapse = ", "), "\n\n")
stopifnot(nrow(atac) == 9)

# Ensembl-named analysis chromosomes, as a regex for awk.
ENS_CHROMS <- c(as.character(1:22), "X")
awk_chrom_filter <- paste0("($1==\"", paste(ENS_CHROMS, collapse = "\"||$1==\""), "\")")

# Histogram resolution: 0.5 signal units, capped at 500 (values above are lumped
# into the top bin — they are a vanishing fraction of bp and only affect the tail).
AWK <- paste0(
  "awk -F'\\t' '", awk_chrom_filter,
  " { w=$3-$2; v=$4+0; b=int(v*2); if(b>1000) b=1000; h[b]+=w; tot+=w; if(v>0) cov+=w }",
  " END { for (k in h) print k\"\\t\"h[k]; print \"#TOTAL\\t\"tot; print \"#COVERED\\t\"cov }'")

summary_rows <- list()

for (i in seq_len(nrow(atac))) {
  line <- atac$line[i]
  dest <- atac$dest[i]
  outf <- file.path(OUT, paste0(line, ".hist.tsv"))

  if (file.exists(outf)) {
    cat(sprintf("  [cached] %-9s\n", line))
  } else {
    t0 <- Sys.time()
    cmd <- paste("zcat", shQuote(dest), "|", AWK)
    res <- system(cmd, intern = TRUE)
    writeLines(res, outf)
    cat(sprintf("  [done]   %-9s %s  (%.1f min)\n", line, basename(dest),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  raw <- readLines(outf)
  tot <- as.numeric(sub("^#TOTAL\t", "", grep("^#TOTAL", raw, value = TRUE)))
  cov <- as.numeric(sub("^#COVERED\t", "", grep("^#COVERED", raw, value = TRUE)))
  hb  <- raw[!grepl("^#", raw)]
  h   <- utils::read.delim(text = paste(hb, collapse = "\n"), header = FALSE,
                           col.names = c("bin", "bp"))
  h$signal <- h$bin / 2
  h <- h[order(h$signal), ]

  # bp-weighted mean over COVERED bases (signal > 0). Zero-signal bases dominate
  # the genome and would drag any genome-wide mean to near zero.
  covered <- h[h$signal > 0, ]
  mean_cov <- sum(covered$signal * covered$bp) / sum(covered$bp)

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    line = line, total_bp = tot, covered_bp = cov,
    pct_covered = round(100 * cov / tot, 2),
    mean_signal_covered = round(mean_cov, 3),
    max_signal_bin = max(covered$signal),
    stringsAsFactors = FALSE)
}

s <- do.call(rbind, summary_rows)
write.csv(s, "data/metadata/atac_signal_profile.csv", row.names = FALSE)

cat("\n=========== ATAC signal profiles ===========\n")
print(s, row.names = FALSE)

cat("\nwrote data/metadata/atac_signal_profile.csv\n")
cat("histograms in ", OUT, "\n", sep = "")
cat("\nRESULT: profiles built.\n")
cat("Next: Rscript scripts/02_regulatory/04_calibrate_threshold.R\n")

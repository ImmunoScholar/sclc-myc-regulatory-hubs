# -----------------------------------------------------------------------------
# 05_derive_regions.R — the accessible-region universe, from keystone ATAC.
#
#   Rscript scripts/02_regulatory/05_derive_regions.R
#
# Applies the H524-calibrated threshold (04_) to all nine keystone ATAC tracks and
# assembles the universe that every downstream step quantifies over.
#
# PER-LINE SCALING. The threshold is stored as a MULTIPLE of each line's own mean
# signal over covered bases, not as an absolute value. Sequencing depth differs
# several-fold across these nine tracks; an absolute cutoff would let the deeply
# sequenced lines dominate the universe and quietly under-represent the others —
# which would bias a paralog comparison, because paralog groups are confounded
# with which lines they come from.
#
# EXTERNAL PEAKS ARE ANNOTATION, NOT A FILTER (D-023). Each region records how
# many keystone lines support it and whether an independent MACS2 peak set
# corroborates it. Nothing is removed for lacking external support, because the
# external lines are mostly not keystone lines and line-specific regions are the
# ones the hypothesis is about.
#
# Output: data/processed/regions/atac_regions_per_line.rds
#         data/processed/regions/universe.rds          <- the analysis grid
#         data/metadata/universe_report.csv
#         results/tables/m5_universe_final.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(yaml)
})
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
RP  <- CFG$regions
OUT <- "data/processed/regions"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MULT <- as.numeric(readLines("data/metadata/chosen_threshold.txt")[1])
stopifnot(is.finite(MULT), MULT > 0)
cat("calibrated threshold multiple: x", MULT, " of mean covered signal\n\n", sep = "")

ENS_CHROMS <- c(as.character(1:22), "X")

man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
man$dest <- file.path("data/raw", man$dataset_id, man$file_name)
ks <- man[man$dataset_id == "GSE230649" & file.exists(man$dest), ]
atac <- ks[!grepl("^GSM[0-9]+_(H3K27Ac|MYC|MYCN|MYCL1)_", ks$file_name), ]
atac$line <- sub("^GSM[0-9]+_([^_]+)_treat_pileup.*$", "\\1", atac$file_name)
stopifnot(nrow(atac) == 9)

line_threshold <- function(line) {
  hf <- file.path("data/processed/atac_profile", paste0(line, ".hist.tsv"))
  raw <- readLines(hf)
  h <- utils::read.delim(text = paste(raw[!grepl("^#", raw)], collapse = "\n"),
                         header = FALSE, col.names = c("bin", "bp"))
  h$signal <- h$bin / 2
  cov <- h[h$signal > 0, ]
  MULT * sum(cov$signal * cov$bp) / sum(cov$bp)
}

extract_regions <- function(path, thr) {
  filt <- paste0("($1==\"", paste(ENS_CHROMS, collapse = "\"||$1==\""), "\")")
  awk  <- sprintf("awk -F'\\t' '%s && $4+0 >= %f { print $1\"\\t\"$2\"\\t\"$3 }'", filt, thr)
  con  <- pipe(paste("zcat", shQuote(path), "|", awk), "r")
  on.exit(close(con))
  d <- try(utils::read.delim(con, header = FALSE,
                             col.names = c("chrom", "start", "end")), silent = TRUE)
  if (inherits(d, "try-error") || !nrow(d)) return(GRanges())
  gr <- GRanges(to_ucsc_seqnames(d$chrom), IRanges(d$start + 1L, d$end))
  gr <- GenomicRanges::reduce(gr, min.gapwidth = RP$merge_gap)
  gr[width(gr) >= RP$min_width & width(gr) <= RP$max_width]
}

# --- per-line region calling --------------------------------------------------
cache <- file.path(OUT, "atac_regions_per_line.rds")
if (file.exists(cache)) {
  per_line <- readRDS(cache)
  cat("per-line regions loaded from cache\n\n")
} else {
  per_line <- list()
  for (i in seq_len(nrow(atac))) {
    ln  <- atac$line[i]
    thr <- line_threshold(ln)
    t0  <- Sys.time()
    gr  <- extract_regions(atac$dest[i], thr)
    per_line[[ln]] <- gr
    cat(sprintf("  %-9s thr=%8.2f  %7d regions, median %4.0f bp, %6.1f Mb  (%.1f min)\n",
                ln, thr, length(gr), median(width(gr)),
                sum(as.numeric(width(gr))) / 1e6,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  saveRDS(per_line, cache)
}

# --- universe -----------------------------------------------------------------
universe <- GenomicRanges::reduce(do.call(c, unname(per_line)),
                                  min.gapwidth = RP$merge_gap)
universe <- universe[width(universe) >= RP$min_width &
                     width(universe) <= RP$max_width]

sup <- vapply(per_line, function(g) IRanges::overlapsAny(universe, g),
              logical(length(universe)))
mcols(universe)$n_lines_supporting <- rowSums(sup)
for (ln in names(per_line)) mcols(universe)[[paste0("atac_", ln)]] <- sup[, ln]

# --- external corroboration, as annotation ------------------------------------
lift <- list.files("data/processed/liftover", pattern = "\\.hg19\\.rds$", full.names = TRUE)
ext_ds <- unique(sub("__.*$", "", basename(lift)))
for (ds in ext_ds) {
  g <- GenomicRanges::reduce(do.call(c, unname(lapply(
    lift[startsWith(basename(lift), ds)], readRDS))), min.gapwidth = RP$merge_gap)
  mcols(universe)[[paste0("ext_", ds)]] <- IRanges::overlapsAny(universe, g)
}
ext_cols <- grep("^ext_", names(mcols(universe)), value = TRUE)
mcols(universe)$n_external_support <-
  rowSums(as.data.frame(mcols(universe)[, ext_cols, drop = FALSE]))

saveRDS(universe, file.path(OUT, "universe.rds"))

# --- report -------------------------------------------------------------------
genome_bp <- sum(load_hg19_sizes()[ANALYSIS_CHROMS_UCSC])
n <- length(universe)
cat("\n=========== universe ===========\n")
cat("regions        : ", format(n, big.mark = ","), "\n", sep = "")
cat("median width   : ", median(width(universe)), " bp\n", sep = "")
cat("total          : ", round(sum(as.numeric(width(universe))) / 1e6, 1), " Mb\n", sep = "")
cat("genome fraction: ", sprintf("%.2f%%", 100 * sum(as.numeric(width(universe))) / genome_bp), "\n", sep = "")

cat("\nsupport across the 9 keystone ATAC lines:\n")
tb <- table(mcols(universe)$n_lines_supporting)
for (k in names(tb))
  cat(sprintf("  %2s line(s): %7s regions (%5.1f%%)\n", k,
              format(as.integer(tb[[k]]), big.mark = ","), 100 * tb[[k]] / n))

cat("\nexternal corroboration (annotation only, nothing filtered):\n")
for (cc in ext_cols)
  cat(sprintf("  %-18s %7s regions (%5.1f%%)\n", sub("^ext_", "", cc),
              format(sum(mcols(universe)[[cc]]), big.mark = ","),
              100 * mean(mcols(universe)[[cc]])))
cat(sprintf("  %-18s %7s regions (%5.1f%%)\n", "none",
            format(sum(mcols(universe)$n_external_support == 0), big.mark = ","),
            100 * mean(mcols(universe)$n_external_support == 0)))

# Sanity: regions supported by more lines should be more often corroborated.
# If that trend is absent, the calling is picking up noise rather than biology.
cat("\nexternal corroboration by line support (expect increasing):\n")
for (k in sort(unique(mcols(universe)$n_lines_supporting))) {
  idx <- mcols(universe)$n_lines_supporting == k
  cat(sprintf("  %2d line(s): %5.1f%% externally corroborated  (n=%s)\n", k,
              100 * mean(mcols(universe)$n_external_support[idx] > 0),
              format(sum(idx), big.mark = ",")))
}

rep <- data.frame(
  metric = c("regions", "median_width_bp", "total_mb", "pct_genome",
             "pct_multi_line", "pct_externally_corroborated"),
  value = c(n, median(width(universe)),
            round(sum(as.numeric(width(universe))) / 1e6, 1),
            round(100 * sum(as.numeric(width(universe))) / genome_bp, 3),
            round(100 * mean(mcols(universe)$n_lines_supporting >= 2), 1),
            round(100 * mean(mcols(universe)$n_external_support > 0), 1)),
  stringsAsFactors = FALSE)
write.csv(rep, "data/metadata/universe_report.csv", row.names = FALSE)

md <- c("# M5 — accessible-region universe (final)", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        paste0("Derived from the **nine keystone ATAC tracks** (GSE230649, hg19, ",
               "cells matched to the ChIP), thresholded at **x", MULT,
               "** each line's mean signal over covered bases — a multiple ",
               "calibrated against H524's independent MACS2 peak set (D-023)."), "",
        "| metric | value |", "|---|---|")
for (i in seq_len(nrow(rep)))
  md <- c(md, paste0("| ", rep$metric[i], " | ", rep$value[i], " |"))
md <- c(md, "",
        "External peak sets are recorded per region as corroboration and are **not**",
        "used to filter: three of the four external lines are absent from the keystone,",
        "so filtering on them would remove the line-specific regions the paralog",
        "hypothesis concerns.")
writeLines(md, "results/tables/m5_universe_final.md")

cat("\nwrote data/processed/regions/universe.rds\n")
cat("wrote results/tables/m5_universe_final.md\n")
cat("\nRESULT: universe built.\n")
cat("Next: signal quantification of the 19 ChIP tracks over this grid.\n")

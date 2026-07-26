# -----------------------------------------------------------------------------
# 04_qc_report.R — structural and biological QC of the acquired data.
#
#   Rscript scripts/01_data/04_qc_report.R [--deep]
#
# 03_verify.R answers "is this the file we meant to download, intact?".
# This script answers the different and harder question: "is what is inside it
# what we think it is?".
#
# Checks:
#   A. bedGraph / BED / narrowPeak — chromosome naming, interval sanity, and a
#      GENOME BUILD ASSERTION against hg19.chrom.sizes. Any interval extending
#      past its chromosome's length is proof of a build mismatch (risk R-05) and
#      is a hard failure, not a warning.
#   B. GSE60052 expression matrix — dimensions against the stated 79 tumour +
#      7 normal, value range consistent with log2, and the two known parsing
#      traps (leading whitespace in headers, '.normal' suffix encoding class).
#   C. GeoMx count matrices — ROI counts against GEO's stated 175 / 121, and
#      panel size, since ~1,800 genes is the reason the spatial layer is
#      restricted.
#   D. liftOver chain — loadable, and covers the chromosomes we care about.
#   E. Cross-dataset cell-line harmonisation — the same line is written three
#      different ways across these datasets (H1048 / NCIH1048, SHP77 / SHP-77).
#      Unharmonised names would silently break every cross-dataset join.
#
# Default mode samples the head of each large file. --deep streams every
# bedGraph in full for a definitive per-chromosome maximum coordinate; that is
# slower and only needs doing once per download.
#
# Writes: data/metadata/qc_results.csv , results/tables/qc_data_acquisition.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
})

args  <- commandArgs(trailingOnly = TRUE)
DEEP  <- "--deep" %in% args
HEAD_LINES <- 200000L

man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
man$dest <- file.path("data/raw", man$dataset_id, man$file_name)
present <- man[file.exists(man$dest), ]

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

qc <- list()
add_qc <- function(dataset, file, check, value, expected, pass, note = "") {
  qc[[length(qc) + 1L]] <<- data.frame(
    dataset = dataset, file = file, check = check,
    value = as.character(value), expected = as.character(expected),
    pass = pass, note = note, stringsAsFactors = FALSE)
}

cat("=========== QC: data acquisition ===========\n")
cat("mode: ", if (DEEP) "DEEP (full stream)" else "fast (head sample)", "\n", sep = "")
cat("files present: ", nrow(present), " / ", nrow(man), "\n\n", sep = "")

# =============================================================================
# hg19 chromosome sizes — ground truth for build assertions
# =============================================================================
CHROM_FILE <- "data/raw/ucsc_hg19_chrom_sizes/hg19.chrom.sizes"
hg19 <- NULL
if (file.exists(CHROM_FILE)) {
  hg19 <- read.delim(CHROM_FILE, header = FALSE,
                     col.names = c("chrom", "size"), stringsAsFactors = FALSE)
  cat("hg19 reference: ", nrow(hg19), " sequences, chr1 = ",
      format(hg19$size[hg19$chrom == "chr1"], big.mark = ","), " bp\n\n", sep = "")
} else {
  cat("WARNING: hg19.chrom.sizes absent — build assertions cannot run.\n\n")
}

# =============================================================================
# A. interval files
# =============================================================================
read_interval_head <- function(path, n = HEAD_LINES) {
  cmd <- if (grepl("\\.gz$", path)) paste("zcat", shQuote(path)) else paste("cat", shQuote(path))
  con <- pipe(paste(cmd, "| head -n", n), "r")
  on.exit(close(con))
  utils::read.delim(con, header = FALSE, stringsAsFactors = FALSE,
                    comment.char = "#", quote = "")
}

# Definitive per-chromosome maxima without loading the file into memory.
stream_max_by_chrom <- function(path) {
  cmd <- if (grepl("\\.gz$", path)) paste("zcat", shQuote(path)) else paste("cat", shQuote(path))
  awk <- "awk -F'\\t' '$1 ~ /^chr/ { if ($3+0 > m[$1]) m[$1]=$3+0 } END { for (c in m) print c\"\\t\"m[c] }'"
  con <- pipe(paste(cmd, "|", awk), "r")
  on.exit(close(con))
  utils::read.delim(con, header = FALSE, col.names = c("chrom", "max_end"),
                    stringsAsFactors = FALSE)
}

interval_files <- present[grepl("\\.(bedgraph|bed|narrowPeak|broadPeak)(\\.gz)?$",
                                present$file_name, ignore.case = TRUE), ]

cat("--- A. interval files (", nrow(interval_files), ") ---\n", sep = "")

for (i in seq_len(nrow(interval_files))) {
  r <- interval_files[i, ]
  d <- tryCatch(read_interval_head(r$dest), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) {
    add_qc(r$dataset_id, r$file_name, "readable", "no", "yes", FALSE,
           "could not parse as an interval file")
    next
  }

  # chromosome naming convention
  chroms  <- unique(d[[1]])
  ucsc_ok <- all(grepl("^chr", chroms))
  add_qc(r$dataset_id, r$file_name, "chrom_naming",
         if (ucsc_ok) "UCSC (chr*)" else paste(head(chroms, 3), collapse = ","),
         "UCSC (chr*)", ucsc_ok,
         if (ucsc_ok) "" else "Ensembl-style names need renaming before any join")

  # interval sanity: start < end
  if (ncol(d) >= 3 && is.numeric(d[[2]]) && is.numeric(d[[3]])) {
    bad_iv <- sum(d[[3]] <= d[[2]])
    add_qc(r$dataset_id, r$file_name, "start_lt_end", bad_iv, 0, bad_iv == 0)
  }

  # BUILD ASSERTION — only meaningful for files declared hg19
  if (!is.null(hg19) && identical(r$genome_build, "hg19") && ncol(d) >= 3) {
    mx <- if (DEEP) {
      tryCatch(stream_max_by_chrom(r$dest), error = function(e) NULL)
    } else {
      agg <- aggregate(d[[3]], by = list(chrom = d[[1]]), FUN = max)
      names(agg) <- c("chrom", "max_end"); agg
    }
    if (!is.null(mx) && nrow(mx)) {
      mm  <- merge(mx, hg19, by = "chrom")
      ovr <- mm[!is.na(mm$size) & mm$max_end > mm$size, ]
      add_qc(r$dataset_id, r$file_name, "build_assertion_hg19",
             if (nrow(ovr)) paste0(nrow(ovr), " chrom(s) exceed hg19 length") else "within hg19 bounds",
             "within hg19 bounds", nrow(ovr) == 0,
             if (nrow(ovr)) paste0("e.g. ", ovr$chrom[1], " max=", ovr$max_end[1],
                                   " > ", ovr$size[1], " — BUILD MISMATCH") else
               if (DEEP) "definitive (full stream)" else "head sample only")
    }
  }

  # hg38 files: record that they are pending liftOver, do not assert hg19
  if (isTRUE(as.logical(r$liftover_required))) {
    add_qc(r$dataset_id, r$file_name, "liftover_pending", r$genome_build,
           "hg19 after lift", TRUE,
           "intervals lift at M5; signal-only sources have regions called in native build first")
  }
}

# =============================================================================
# B. GSE60052 expression matrix
# =============================================================================
cat("--- B. GSE60052 bulk tumour matrix ---\n")
g60 <- present$dest[present$dataset_id == "GSE60052"]
if (length(g60) && file.exists(g60[1])) {
  con <- gzfile(g60[1], "rt")
  hdr <- readLines(con, n = 1); close(con)
  cols <- strsplit(hdr, "\t")[[1]]
  samples <- cols[-1]

  # Trap 1: leading/trailing whitespace in sample names
  ws <- sum(samples != trimws(samples))
  add_qc("GSE60052", basename(g60[1]), "header_whitespace", ws, "0 (after trimws)",
         TRUE, paste0(ws, " sample name(s) carry whitespace — trimws() is mandatory at load"))

  samples_t <- trimws(samples)
  n_norm <- sum(grepl("\\.normal$", samples_t))
  n_tum  <- length(samples_t) - n_norm
  add_qc("GSE60052", basename(g60[1]), "n_samples", length(samples_t), 86,
         length(samples_t) == 86)
  add_qc("GSE60052", basename(g60[1]), "n_normal", n_norm, 7, n_norm == 7,
         "normals identified by the '.normal' suffix")
  add_qc("GSE60052", basename(g60[1]), "n_tumour", n_tum, 79, n_tum == 79)

  # value range: log2-normalised data should not look like raw counts
  con <- gzfile(g60[1], "rt"); block <- readLines(con, n = 2001); close(con)
  vals <- suppressWarnings(as.numeric(unlist(
    lapply(strsplit(block[-1], "\t"), function(x) x[-1]))))
  vals <- vals[is.finite(vals)]
  mx <- max(vals, na.rm = TRUE)
  looks_log2 <- mx < 100
  add_qc("GSE60052", basename(g60[1]), "value_scale",
         paste0("max=", round(mx, 2)), "log2 scale (max < 100)", looks_log2,
         "confirms raw counts are NOT available -> limma, not DESeq2 (see registry caution)")
  add_qc("GSE60052", basename(g60[1]), "n_genes_sampled", length(block) - 1,
         ">= 2000 readable", (length(block) - 1) >= 2000)
} else {
  cat("  not present yet\n")
}

# =============================================================================
# C. GeoMx count matrices
# =============================================================================
cat("--- C. GeoMx spatial matrices ---\n")
geo_expect <- c(GSE261348 = 175, GSE261345 = 121)
for (ds in names(geo_expect)) {
  f <- present$dest[present$dataset_id == ds & grepl("rawcounts\\.xlsx$", present$file_name)]
  if (!length(f)) { cat("  ", ds, ": raw counts not present yet\n", sep = ""); next }
  ok <- requireNamespace("readxl", quietly = TRUE)
  if (!ok) { cat("  readxl unavailable\n"); next }
  x <- tryCatch(readxl::read_excel(f[1], n_max = 5), error = function(e) NULL)
  if (is.null(x)) {
    add_qc(ds, basename(f[1]), "readable", "no", "yes", FALSE)
    next
  }
  # ROIs are columns (minus the gene-id column)
  n_roi <- ncol(x) - 1
  add_qc(ds, basename(f[1]), "n_roi", n_roi, geo_expect[[ds]],
         n_roi == geo_expect[[ds]],
         "ROI count from GEO series record")
  full <- tryCatch(readxl::read_excel(f[1]), error = function(e) NULL)
  if (!is.null(full)) {
    add_qc(ds, basename(f[1]), "panel_genes", nrow(full), "~1800 (targeted CTA)",
           nrow(full) > 1000 && nrow(full) < 3000,
           "targeted panel — regulon coverage fraction must be reported on every spatial figure")
  }
}

# =============================================================================
# D. liftOver chain
# =============================================================================
cat("--- D. liftOver chain ---\n")
chain <- "data/raw/ucsc_chain_hg38_to_hg19/hg38ToHg19.over.chain.gz"
if (file.exists(chain)) {
  okc <- requireNamespace("rtracklayer", quietly = TRUE)
  if (okc) {
    tmp <- tempfile(fileext = ".chain")
    system2("gunzip", c("-c", shQuote(chain)), stdout = tmp)
    ch <- tryCatch(rtracklayer::import.chain(tmp), error = function(e) NULL)
    if (is.null(ch)) {
      add_qc("ucsc_chain_hg38_to_hg19", basename(chain), "loadable", "no", "yes", FALSE)
    } else {
      nm <- names(ch)
      need <- paste0("chr", c(1:22, "X"))
      miss <- setdiff(need, nm)
      add_qc("ucsc_chain_hg38_to_hg19", basename(chain), "loadable",
             paste0(length(nm), " chains"), "loadable", TRUE)
      add_qc("ucsc_chain_hg38_to_hg19", basename(chain), "covers_analysis_chroms",
             if (length(miss)) paste(miss, collapse = ",") else "all 23",
             "chr1-22,X", length(miss) == 0)
    }
    unlink(tmp)
  }
} else {
  cat("  not present yet\n")
}

# =============================================================================
# E. cell-line name harmonisation across datasets
# =============================================================================
cat("--- E. cell-line harmonisation ---\n")
# Verified spellings encountered in these datasets. Left = canonical.
alias <- list(
  H1048   = c("H1048", "NCIH1048", "NCI-H1048"),
  H211    = c("H211", "NCIH211", "NCI-H211"),
  H526    = c("H526", "NCIH526", "NCI-H526"),
  SHP77   = c("SHP77", "SHP-77"),
  H524    = c("H524", "NCIH524"),
  H69     = c("H69", "NCIH69"),
  H847    = c("H847", "NCIH847"),
  H889    = c("H889", "NCIH889"),
  H196    = c("H196", "NCIH196"),
  COLO668 = c("COLO668", "COLO-668")
)
variants <- sum(vapply(alias, length, integer(1))) - length(alias)
add_qc("cross_dataset", "(all)", "cellline_alias_map",
       paste0(length(alias), " lines, ", variants, " non-canonical spellings"),
       "alias map defined", TRUE,
       "H1048/NCIH1048 and SHP77/SHP-77 both occur; joins must go through this map")

hit <- character(0)
for (canon in names(alias)) {
  for (v in alias[[canon]]) {
    if (any(grepl(v, present$file_name, ignore.case = TRUE))) {
      hit <- c(hit, paste0(v, "->", canon))
    }
  }
}
add_qc("cross_dataset", "(all)", "aliases_observed_in_filenames",
       length(hit), ">0", length(hit) > 0, paste(unique(hit), collapse = "; "))

# =============================================================================
# report
# =============================================================================
res <- do.call(rbind, qc)
write.csv(res, "data/metadata/qc_results.csv", row.names = FALSE)

n_pass <- sum(res$pass); n_fail <- sum(!res$pass)
cat("\n=========== QC summary ===========\n")
cat("checks run   : ", nrow(res), "\n", sep = "")
cat("passed       : ", n_pass, "\n", sep = "")
cat("FAILED       : ", n_fail, "\n", sep = "")

if (n_fail) {
  cat("\n--- failures ---\n")
  f <- res[!res$pass, ]
  for (i in seq_len(nrow(f)))
    cat(sprintf("  [%s] %s :: %s = %s (expected %s) %s\n",
                f$dataset[i], f$file[i], f$check[i], f$value[i],
                f$expected[i], f$note[i]))
}

md <- c("# QC — Data Acquisition (M4)", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste0("Mode: ", if (DEEP) "deep (full file stream)" else "fast (head sample)"),
        paste0("Checks: **", nrow(res), "**, passed **", n_pass,
               "**, failed **", n_fail, "**"), "",
        "Generated by `scripts/01_data/04_qc_report.R`. Regenerate; do not edit.", "",
        "| dataset | file | check | value | expected | pass | note |",
        "|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(res))) {
  md <- c(md, paste0("| ", res$dataset[i], " | `", substr(res$file[i], 1, 44),
                     "` | ", res$check[i], " | ", res$value[i], " | ",
                     res$expected[i], " | ", ifelse(res$pass[i], "yes", "**NO**"),
                     " | ", res$note[i], " |"))
}
writeLines(md, "results/tables/qc_data_acquisition.md")
cat("\nwrote data/metadata/qc_results.csv\n")
cat("wrote results/tables/qc_data_acquisition.md\n")

if (n_fail) {
  cat("\nRESULT: FAIL — resolve the failures above before M5.\n")
  quit(status = 1)
}
cat("\nRESULT: PASS\n")

# -----------------------------------------------------------------------------
# 04_qc_report.R — structural and biological QC of the acquired data.
#
#   Rscript scripts/01_data/04_qc_report.R            # full build assertion
#   Rscript scripts/01_data/04_qc_report.R --fast     # skip the full stream
#
# 03_verify.R answers "is this the file we meant to download, intact?".
# This answers the harder question: "is what is inside it what we think it is?".
#
# TWO DEFECTS IN THE FIRST VERSION OF THIS SCRIPT, both now fixed:
#
#  1. The genome-build assertion SILENTLY PASSED ON NOTHING. It merged bedGraph
#     sequence names against hg19.chrom.sizes, but the bedGraphs use Ensembl
#     names ("20") and the reference uses UCSC names ("chr20"). The join matched
#     zero rows, found zero violations, and recorded a pass. Names are now
#     harmonised first, and assert_hg19_bounds() REFUSES to return a pass when it
#     had fewer than 20 chromosomes to compare. (risk R-16)
#
#  2. Build assertions were computed from a 200,000-line head sample. These
#     bedGraphs are written in BAM-header order (20, 21, 22, 1, 3, 2, ...), so a
#     head sample only ever contains ONE chromosome. The assertion now streams
#     each file in full by default. It is slower and it is the point.
#
#  3. GeoMx counts were read from sheet 1 (SegmentProperties). The count matrix
#     is in the TargetCountMatrix sheet. Read correctly, ROI counts match GEO
#     exactly (175 and 121).
#
# Writes: data/metadata/qc_results.csv , results/tables/qc_data_acquisition.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
})
source("R/genome_utils.R")

args   <- commandArgs(trailingOnly = TRUE)
FAST   <- "--fast" %in% args
RESTREAM <- "--restream" %in% args   # ignore the cache and re-stream every file
HEAD_LINES <- 200000L

CFG <- yaml::read_yaml("config/params.yml")

# --- build-assertion cache ----------------------------------------------------
# Streaming 28 bedGraphs (~1.3 billion lines) takes 30-60 min. QC gets re-run
# whenever any other check changes, so the streamed result is cached against the
# file's size and mtime. A changed file invalidates its own entry automatically;
# --restream forces a full recompute.
CACHE_DIR <- "data/metadata/build_assertion_cache"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

cached_max_by_chrom <- function(path) {
  fi  <- file.info(path)
  sig <- paste0(basename(path), "_", fi$size, "_",
                format(fi$mtime, "%Y%m%d%H%M%S"))
  key <- substr(digest_str(sig), 1, 16)
  cf  <- file.path(CACHE_DIR, paste0(key, ".tsv"))
  if (!RESTREAM && file.exists(cf)) {
    return(utils::read.delim(cf, stringsAsFactors = FALSE))
  }
  out <- stream_max_by_chrom(path, col_end = 3L)
  if (!is.null(out) && nrow(out)) {
    utils::write.table(out, cf, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  out
}

digest_str <- function(s) {
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(s, algo = "md5")
  } else {
    # deterministic fallback, no extra dependency
    paste0(sum(utils::head(utils::strtoi(charToRaw(s), 16L), 64L)), "_", nchar(s))
  }
}

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
cat("build assertion: ", if (FAST) "SKIPPED (--fast)" else "full file stream", "\n", sep = "")
cat("files present  : ", nrow(present), " / ", nrow(man), "\n\n", sep = "")

# =============================================================================
# hg19 ground truth
# =============================================================================
hg19 <- tryCatch(load_hg19_sizes(), error = function(e) NULL)
if (is.null(hg19)) {
  cat("WARNING: hg19.chrom.sizes absent — build assertions cannot run.\n\n")
} else {
  cat("hg19 reference : ", length(hg19), " sequences, chr1 = ",
      format(hg19[["chr1"]], big.mark = ","), " bp\n\n", sep = "")
}

# =============================================================================
# A. interval files
# =============================================================================
read_interval_head <- function(path, n = HEAD_LINES) {
  reader <- if (grepl("\\.gz$", path)) "zcat" else "cat"
  con <- pipe(paste(reader, shQuote(path), "| head -n", n), "r")
  on.exit(close(con))
  utils::read.delim(con, header = FALSE, stringsAsFactors = FALSE,
                    comment.char = "#", quote = "")
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

  chroms <- unique(as.character(d[[1]]))
  style  <- detect_seqname_style(chroms)

  # Naming convention is a PROPERTY to record, not a defect — provided it is one
  # of the two known conventions. "mixed" or "unknown" is a real problem.
  add_qc(r$dataset_id, r$file_name, "seqname_style", style,
         "UCSC or Ensembl", style %in% c("UCSC", "Ensembl"),
         if (identical(style, "Ensembl"))
           "MUST pass through harmonise_seqnames() before any coordinate join"
         else "")

  # interval sanity
  if (ncol(d) >= 3 && is.numeric(d[[2]]) && is.numeric(d[[3]])) {
    bad_iv <- sum(d[[3]] <= d[[2]], na.rm = TRUE)
    add_qc(r$dataset_id, r$file_name, "start_lt_end", bad_iv, 0, bad_iv == 0)
  }

  # --- BUILD ASSERTION, for files declared hg19 ------------------------------
  if (!is.null(hg19) && identical(r$genome_build, "hg19") && ncol(d) >= 3) {
    if (FAST) {
      add_qc(r$dataset_id, r$file_name, "build_assertion_hg19", "not run",
             "within hg19 bounds", FALSE,
             "--fast given; a head sample cannot assert a build on these files")
    } else {
      mx  <- cached_max_by_chrom(r$dest)
      res <- assert_hg19_bounds(mx, hg19)
      add_qc(r$dataset_id, r$file_name, "build_assertion_hg19",
             paste0(res$n_compared, " chroms compared"),
             "within hg19 bounds, >=20 chroms", res$pass, res$detail)

      # Completeness: a truncated file loses whole chromosomes.
      if (!is.null(mx)) {
        seen <- to_ucsc_seqnames(mx$chrom)
        miss <- setdiff(ANALYSIS_CHROMS_UCSC, seen)
        add_qc(r$dataset_id, r$file_name, "analysis_chroms_present",
               paste0(length(intersect(ANALYSIS_CHROMS_UCSC, seen)), "/23"),
               "23 (chr1-22,X)", length(miss) == 0,
               if (length(miss)) paste("missing:", paste(miss, collapse = ",")) else "")
      }
    }
  }

  if (isTRUE(as.logical(r$liftover_required))) {
    add_qc(r$dataset_id, r$file_name, "liftover_pending", r$genome_build,
           "hg19 after lift", TRUE,
           "intervals lift at M5; loss rate reported per dataset (D-014)")
  }
}

# =============================================================================
# B. GSE60052 expression matrix
# =============================================================================
cat("--- B. GSE60052 bulk tumour matrix ---\n")
g60 <- present$dest[present$dataset_id == "GSE60052"]
if (length(g60) && file.exists(g60[1])) {
  con <- gzfile(g60[1], "rt"); hdr <- readLines(con, n = 1); close(con)
  cols    <- strsplit(hdr, "\t")[[1]]
  samples <- cols[-1]

  ws <- sum(samples != trimws(samples))
  add_qc("GSE60052", basename(g60[1]), "header_whitespace", ws,
         "handled by trimws()", TRUE,
         paste0(ws, " sample name(s) carry whitespace — trimws() is mandatory at load"))

  samples_t <- trimws(samples)
  n_norm <- sum(grepl("\\.normal$", samples_t))
  n_tum  <- length(samples_t) - n_norm
  add_qc("GSE60052", basename(g60[1]), "n_samples", length(samples_t), 86,
         length(samples_t) == 86)
  add_qc("GSE60052", basename(g60[1]), "n_normal", n_norm, 7, n_norm == 7,
         "normals identified by the '.normal' suffix")
  add_qc("GSE60052", basename(g60[1]), "n_tumour", n_tum, 79, n_tum == 79)

  con <- gzfile(g60[1], "rt"); block <- readLines(con, n = 2001); close(con)
  vals <- suppressWarnings(as.numeric(unlist(
    lapply(strsplit(block[-1], "\t"), function(x) x[-1]))))
  vals <- vals[is.finite(vals)]
  mx <- max(vals, na.rm = TRUE)
  add_qc("GSE60052", basename(g60[1]), "value_scale",
         paste0("max=", round(mx, 2)), "log2 scale (max < 100)", mx < 100,
         "confirms raw counts are NOT available -> limma, not DESeq2 (R-15)")
} else {
  cat("  not present\n")
}

# =============================================================================
# C. GeoMx count matrices — TargetCountMatrix sheet, not sheet 1
# =============================================================================
cat("--- C. GeoMx spatial matrices ---\n")
geo_expect <- list(GSE261348 = 175L, GSE261345 = 121L)
COUNT_SHEET <- "TargetCountMatrix"
SEG_SHEET   <- "SegmentProperties"

for (ds in names(geo_expect)) {
  f <- present$dest[present$dataset_id == ds & grepl("rawcounts\\.xlsx$", present$file_name)]
  if (!length(f)) { cat("  ", ds, ": raw counts not present\n", sep = ""); next }
  if (!requireNamespace("readxl", quietly = TRUE)) { cat("  readxl unavailable\n"); break }

  sheets <- tryCatch(readxl::excel_sheets(f[1]), error = function(e) character(0))
  add_qc(ds, basename(f[1]), "has_count_sheet",
         if (COUNT_SHEET %in% sheets) COUNT_SHEET else paste(sheets, collapse = ","),
         COUNT_SHEET, COUNT_SHEET %in% sheets,
         "workbook has 5 sheets; sheet 1 is SegmentProperties, not the counts")
  if (!(COUNT_SHEET %in% sheets)) next

  cm <- tryCatch(readxl::read_excel(f[1], sheet = COUNT_SHEET),
                 error = function(e) NULL)
  if (is.null(cm)) {
    add_qc(ds, basename(f[1]), "count_sheet_readable", "no", "yes", FALSE)
    next
  }

  n_roi   <- ncol(cm) - 1L        # first column is TargetName
  n_panel <- nrow(cm)
  exp_roi <- geo_expect[[ds]]

  add_qc(ds, basename(f[1]), "n_roi", n_roi, exp_roi, n_roi == exp_roi,
         "ROI count from the GEO series record")
  add_qc(ds, basename(f[1]), "panel_genes", n_panel, "1000-3000 (targeted CTA)",
         n_panel > 1000 && n_panel < 3000,
         "targeted panel — regulon coverage fraction must be reported on every spatial figure")

  # Annotated segments vs deposited ROIs. Where a gap exists it must match a
  # gap we have already CHARACTERISED in config (spatial.known_segment_gaps) —
  # not merely be tolerated. An undeclared gap, or a declared one whose shape has
  # changed, is a failure. This tests reality against a documented expectation
  # rather than relaxing the check until it goes green.
  if (SEG_SHEET %in% sheets) {
    sp <- tryCatch(readxl::read_excel(f[1], sheet = SEG_SHEET), error = function(e) NULL)
    if (!is.null(sp)) {
      declared <- CFG$spatial$known_segment_gaps[[ds]]
      if (nrow(sp) == n_roi) {
        add_qc(ds, basename(f[1]), "segment_rows_vs_matrix_roi",
               paste0(nrow(sp), " vs ", n_roi), "equal", TRUE)
      } else if (!is.null(declared)) {
        ok <- identical(as.integer(declared$annotated_segments), as.integer(nrow(sp))) &&
              identical(as.integer(declared$deposited_rois),     as.integer(n_roi))
        add_qc(ds, basename(f[1]), "segment_gap_matches_declared",
               paste0(nrow(sp), " vs ", n_roi),
               paste0(declared$annotated_segments, " vs ", declared$deposited_rois),
               ok,
               if (ok) paste0("known deposit gap, characterised (R-17): ",
                              declared$missing_count, " segments from slide ",
                              declared$missing_slide,
                              " — M9 action: exclude that slide")
               else "gap differs from the characterised one — re-characterise before use")
      } else {
        add_qc(ds, basename(f[1]), "segment_rows_vs_matrix_roi",
               paste0(nrow(sp), " vs ", n_roi), "equal, or declared in config",
               FALSE,
               paste0(abs(nrow(sp) - n_roi),
                      " annotated segment(s) absent and NOT declared in ",
                      "config/params.yml spatial.known_segment_gaps — characterise before M9"))
      }
    }
  }
}

# =============================================================================
# D. liftOver chain
# =============================================================================
cat("--- D. liftOver chain ---\n")
chain <- "data/raw/ucsc_chain_hg38_to_hg19/hg38ToHg19.over.chain.gz"
if (file.exists(chain) && requireNamespace("rtracklayer", quietly = TRUE)) {
  tmp <- tempfile(fileext = ".chain")
  system2("gunzip", c("-c", shQuote(chain)), stdout = tmp)
  ch <- tryCatch(rtracklayer::import.chain(tmp), error = function(e) NULL)
  if (is.null(ch)) {
    add_qc("ucsc_chain_hg38_to_hg19", basename(chain), "loadable", "no", "yes", FALSE)
  } else {
    miss <- setdiff(ANALYSIS_CHROMS_UCSC, names(ch))
    add_qc("ucsc_chain_hg38_to_hg19", basename(chain), "loadable",
           paste0(length(names(ch)), " chains"), "loadable", TRUE)
    add_qc("ucsc_chain_hg38_to_hg19", basename(chain), "covers_analysis_chroms",
           if (length(miss)) paste(miss, collapse = ",") else "all 23",
           "chr1-22,X", length(miss) == 0)
  }
  unlink(tmp)
} else {
  cat("  chain absent or rtracklayer unavailable\n")
}

# =============================================================================
# E. cross-dataset naming hazards
# =============================================================================
cat("--- E. cross-dataset harmonisation ---\n")
alias <- list(
  H1048 = c("H1048", "NCIH1048", "NCI-H1048"), H211 = c("H211", "NCIH211"),
  H526  = c("H526", "NCIH526"), SHP77 = c("SHP77", "SHP-77"),
  H524  = c("H524", "NCIH524"), H69 = c("H69"), H847 = c("H847"),
  H889  = c("H889"), H196 = c("H196"), COLO668 = c("COLO668")
)
hit <- character(0)
for (canon in names(alias)) for (v in alias[[canon]])
  if (any(grepl(v, present$file_name, fixed = TRUE))) hit <- c(hit, paste0(v, "->", canon))
add_qc("cross_dataset", "(all)", "cellline_aliases_observed", length(unique(hit)),
       ">0", length(hit) > 0, paste(unique(hit), collapse = "; "))

# Sequence-name conventions actually in play across the project.
styles <- unique(vapply(qc, function(x)
  if (identical(x$check, "seqname_style")) x$value else NA_character_, character(1)))
styles <- styles[!is.na(styles)]
add_qc("cross_dataset", "(all)", "seqname_conventions_in_project",
       paste(sort(styles), collapse = " + "), "documented", TRUE,
       "both conventions occur; harmonise_seqnames() is mandatory before any join (R-16)")

# =============================================================================
# report
# =============================================================================
res <- do.call(rbind, qc)
write.csv(res, "data/metadata/qc_results.csv", row.names = FALSE)

n_pass <- sum(res$pass); n_fail <- sum(!res$pass)
cat("\n=========== QC summary ===========\n")
cat("checks run : ", nrow(res), "\n", sep = "")
cat("passed     : ", n_pass, "\n", sep = "")
cat("FAILED     : ", n_fail, "\n", sep = "")

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
        paste0("Build assertion: ", if (FAST) "SKIPPED (--fast)" else "full file stream"),
        paste0("Checks: **", nrow(res), "**, passed **", n_pass,
               "**, failed **", n_fail, "**"), "",
        "Generated by `scripts/01_data/04_qc_report.R`. Regenerate; do not edit.", "",
        "| dataset | file | check | value | expected | pass | note |",
        "|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(res)))
  md <- c(md, paste0("| ", res$dataset[i], " | `", substr(res$file[i], 1, 44),
                     "` | ", res$check[i], " | ", res$value[i], " | ",
                     res$expected[i], " | ", ifelse(res$pass[i], "yes", "**NO**"),
                     " | ", res$note[i], " |"))
writeLines(md, "results/tables/qc_data_acquisition.md")
cat("\nwrote data/metadata/qc_results.csv\n")
cat("wrote results/tables/qc_data_acquisition.md\n")

if (n_fail) { cat("\nRESULT: FAIL — resolve the failures above before M5.\n"); quit(status = 1) }
cat("\nRESULT: PASS\n")

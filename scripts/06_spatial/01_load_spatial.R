# -----------------------------------------------------------------------------
# 01_load_spatial.R — load the GeoMx cohorts, verify them, and apply the gates.
#
#   Rscript scripts/06_spatial/01_load_spatial.R
#
# The spatial layer is RESTRICTED ORTHOGONAL VALIDATION (datasets.yml D5). It is
# a targeted ~1,800-gene CTA panel, so most regulon members are simply not
# measured, and config/params.yml requires a regulon to clear BOTH an absolute
# floor (>= min_genes_measured) and a coverage fraction before it may be scored.
#
# Three checks run here and each can stop the pipeline:
#   1. deposited ROI counts match the registry (175 / 121)
#   2. the KNOWN deposit gap in GSE261348 still looks exactly as characterised —
#      a DIFFERENT discrepancy on re-download must fail, not be waved through
#   3. per-regulon panel coverage, gated
#
# Output: data/processed/spatial/geomx.rds
#         data/metadata/m9_panel_coverage.csv, m9_segment_qc.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(readxl); library(yaml) })

CFG <- yaml::read_yaml("config/params.yml")
SP  <- CFG$spatial
REG <- readRDS("data/processed/regions/regulons.rds")$regulons
dir.create("data/processed/spatial", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

COHORTS <- list(
  GSE261348 = list(file = "data/raw/GSE261348/GSE261348_IMfirst_DSP_normalizedcounts.xlsx",
                   label = "IMfirst", expect_roi = 175L),
  GSE261345 = list(file = "data/raw/GSE261345/GSE261345_CANTABRICO_DSP_normalizedcounts.xlsx",
                   label = "CANTABRICO", expect_roi = 121L))

read_cohort <- function(cc) {
  # TargetCountMatrix, NOT SegmentProperties. Reading the wrong sheet here yields
  # a table with a plausible shape and entirely wrong contents.
  tcm <- suppressMessages(readxl::read_excel(cc$file, sheet = "TargetCountMatrix"))
  seg <- suppressMessages(readxl::read_excel(cc$file, sheet = "SegmentProperties"))
  tp  <- suppressMessages(readxl::read_excel(cc$file, sheet = "TargetProperties"))
  list(counts = tcm, seg = seg, targets = tp)
}

cat("=========== loading GeoMx cohorts ===========\n")
dat <- list(); cov_rows <- data.frame(); qc_rows <- data.frame()

for (id in names(COHORTS)) {
  cc <- COHORTS[[id]]
  d  <- read_cohort(cc)
  roi_cols <- setdiff(names(d$counts), "TargetName")
  n_roi <- length(roi_cols)

  cat("\n--- ", id, " (", cc$label, ") ---\n", sep = "")
  cat("  panel targets   : ", nrow(d$counts), "\n", sep = "")
  cat("  deposited ROIs  : ", n_roi, "\n", sep = "")
  cat("  annotated segs  : ", nrow(d$seg), "\n", sep = "")

  # --- check 1: ROI count matches the registry -------------------------------
  if (n_roi != cc$expect_roi)
    stop(id, ": expected ", cc$expect_roi, " deposited ROIs, found ", n_roi,
         ". The registry and the deposit disagree; resolve before analysis.")

  # --- check 2: the known deposit gap is UNCHANGED ----------------------------
  # Declared in config precisely so a different discrepancy fails loudly rather
  # than being absorbed as "some segments are missing, as expected".
  gap <- SP$known_segment_gaps[[id]]
  if (!is.null(gap)) {
    if (nrow(d$seg) != gap$annotated_segments || n_roi != gap$deposited_rois)
      stop(id, ": deposit gap differs from the characterised one (config says ",
           gap$annotated_segments, " annotated / ", gap$deposited_rois,
           " deposited; found ", nrow(d$seg), " / ", n_roi, "). Re-characterise it.")
    missing <- setdiff(d$seg$SegmentDisplayName, roi_cols)
    miss_slides <- unique(sub(" \\|.*$", "", missing))
    if (length(missing) != gap$missing_count || !identical(miss_slides, gap$missing_slide))
      stop(id, ": the MISSING segments are not the ones characterised (found ",
           length(missing), " on slide(s) ", paste(miss_slides, collapse = ", "),
           "; config says ", gap$missing_count, " on ", gap$missing_slide, ").")
    cat("  deposit gap     : ", length(missing), " segments on ", gap$missing_slide,
        " — matches characterisation\n", sep = "")
  }

  # --- apply the m9 action: drop the compromised slide ------------------------
  seg <- d$seg[d$seg$SegmentDisplayName %in% roi_cols, ]
  drop_slide <- if (!is.null(gap)) gap$missing_slide else NA_character_
  if (!is.na(drop_slide)) {
    n_before <- nrow(seg)
    seg <- seg[seg$SlideName != drop_slide, ]
    cat("  m9 action       : dropped slide ", drop_slide, " (", n_before - nrow(seg),
        " ROI) — 40x under-sequenced and QC-flagged\n", sep = "")
  }

  mat <- as.matrix(d$counts[, seg$SegmentDisplayName, drop = FALSE])
  rownames(mat) <- d$counts$TargetName
  storage.mode(mat) <- "double"

  # --- check 3: per-regulon panel coverage, gated ----------------------------
  panel_genes <- unique(stats::na.omit(c(d$targets$HUGOSymbol, d$targets$TargetName)))
  cat("  regulon coverage (need >= ", SP$min_genes_measured, " genes AND >= ",
      100 * SP$min_coverage_fraction, "%):\n", sep = "")
  for (p in names(REG)) {
    on <- intersect(REG[[p]], rownames(mat))
    frac <- length(on) / length(REG[[p]])
    pass <- length(on) >= SP$min_genes_measured && frac >= SP$min_coverage_fraction
    cov_rows <- rbind(cov_rows, data.frame(
      cohort = id, paralog = p, regulon_size = length(REG[[p]]),
      measured = length(on), coverage_fraction = round(frac, 4),
      passes_floor = length(on) >= SP$min_genes_measured,
      passes_fraction = frac >= SP$min_coverage_fraction,
      scoreable = pass, stringsAsFactors = FALSE))
    cat(sprintf("    %-6s %3d/%d (%.1f%%) -> %s\n", p, length(on), length(REG[[p]]),
                100 * frac, if (pass) "SCOREABLE" else "NOT SCOREABLE"))
  }

  seg$slide <- seg$SlideName
  # A slide named "IMF-001/002" or "CAN-037 CAN-039" carries material from TWO
  # patients and the deposit gives no per-ROI patient label. Slides whose name
  # holds a single identifier are the only ones where ROIs are unambiguously one
  # tumour. Recorded now so no downstream step can silently treat slide = tumour.
  seg$n_ids <- lengths(regmatches(seg$slide, gregexpr("[0-9]{3}", seg$slide)))
  seg$single_patient <- seg$n_ids == 1L
  cat("  slides          : ", length(unique(seg$slide)), " (",
      length(unique(seg$slide[seg$single_patient])), " unambiguously one patient)\n", sep = "")

  qc_rows <- rbind(qc_rows, data.frame(
    cohort = id, roi = seg$SegmentDisplayName, slide = seg$slide,
    single_patient = seg$single_patient, segment_type = seg$SegmentLabel,
    nuclei = as.numeric(seg$AOINucleiCount), raw_reads = as.numeric(seg$RawReads),
    saturation = as.numeric(seg$SequencingSaturation), stringsAsFactors = FALSE))

  dat[[id]] <- list(mat = mat, seg = seg, panel_genes = panel_genes, label = cc$label)
}

write.csv(cov_rows, "data/metadata/m9_panel_coverage.csv", row.names = FALSE)
write.csv(qc_rows,  "data/metadata/m9_segment_qc.csv", row.names = FALSE)
saveRDS(dat, "data/processed/spatial/geomx.rds")

cat("\n=========== gate summary ===========\n")
sc <- unique(cov_rows$paralog[cov_rows$scoreable])
ns <- setdiff(names(REG), sc)
cat("  scoreable regulons     : ", if (length(sc)) paste(sc, collapse = ", ") else "NONE", "\n", sep = "")
cat("  NOT scoreable          : ", if (length(ns)) paste(ns, collapse = ", ") else "none", "\n", sep = "")
if (!length(sc))
  stop("no regulon clears the coverage gate; spatial scoring is not possible and ",
       "M9 must report that rather than lower the threshold.")

# The margin matters as much as the verdict: a regulon that clears a 10% floor at
# 11.2% is not comfortably measured, and a reader must be told which side of a
# knife edge each verdict sits on.
cat("\n  margins against the ", 100 * SP$min_coverage_fraction, "% fraction floor:\n", sep = "")
for (i in which(cov_rows$cohort == names(COHORTS)[1])) {
  cat(sprintf("    %-6s %.1f%% (%+.1f points)\n", cov_rows$paralog[i],
              100 * cov_rows$coverage_fraction[i],
              100 * (cov_rows$coverage_fraction[i] - SP$min_coverage_fraction)))
}
cat("\nwrote data/processed/spatial/geomx.rds\n")
cat("wrote data/metadata/m9_panel_coverage.csv\n")
cat("wrote data/metadata/m9_segment_qc.csv\n")

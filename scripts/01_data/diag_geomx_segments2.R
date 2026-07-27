# Characterise the 9 GSE261348 segments absent from the count matrix.
# Is it a naming mismatch, a sequencing failure, or an undocumented exclusion?
suppressPackageStartupMessages(library(readxl))

f  <- "data/raw/GSE261348/GSE261348_IMfirst_DSP_rawcounts.xlsx"
sp <- read_excel(f, sheet = "SegmentProperties")
cm <- read_excel(f, sheet = "TargetCountMatrix", n_max = 1)

mat_ids <- names(cm)[-1]
key <- c("SlideName", "ROILabel", "SegmentLabel")
sp_ids <- apply(sp[, key], 1, function(r) paste(trimws(r), collapse = " | "))

cat("=== direction check ===\n")
cat("in SegmentProperties, not in matrix :", length(setdiff(sp_ids, mat_ids)), "\n")
cat("in matrix, not in SegmentProperties :", length(setdiff(mat_ids, sp_ids)), "\n")
cat("  -> if the second is 0, the 9 are genuinely absent, not renamed\n\n")

cat("=== slide IMF-001/002: all its segments ===\n")
sl <- sp$SlideName == "IMF-001/002"
cat("segments on this slide:", sum(sl), "\n")
present_flag <- sp_ids %in% mat_ids
print(data.frame(
  ROI      = sp$ROILabel[sl],
  in_matrix = present_flag[sl],
  RawReads = sp$RawReads[sl],
  Aligned  = sp$AlignedReads[sl],
  Dedup    = sp$DeduplicatedReads[sl],
  Sat      = round(as.numeric(sp$SequencingSaturation[sl]), 1),
  Nuclei   = sp$AOINucleiCount[sl],
  QC       = ifelse(is.na(sp$QCFlags[sl]), "-", substr(sp$QCFlags[sl], 1, 20))
))

cat("\n=== read depth: absent vs present segments (whole cohort) ===\n")
for (v in c("RawReads", "AlignedReads", "DeduplicatedReads", "SequencingSaturation", "AOINucleiCount")) {
  a <- suppressWarnings(as.numeric(sp[[v]][!present_flag]))
  b <- suppressWarnings(as.numeric(sp[[v]][present_flag]))
  cat(sprintf("  %-22s absent: median %-12s | present: median %s\n", v,
              format(median(a, na.rm = TRUE), big.mark = ","),
              format(median(b, na.rm = TRUE), big.mark = ",")))
}

cat("\n=== how many distinct slides / patients? ===\n")
cat("slides in SegmentProperties:", length(unique(sp$SlideName)), "\n")
cat("slides represented in matrix:",
    length(unique(sp$SlideName[present_flag])), "\n")
cat("\nsegments per slide (absent ones marked):\n")
tb <- table(sp$SlideName, ifelse(present_flag, "in_matrix", "ABSENT"))
print(tb[rowSums(tb[, "ABSENT", drop = FALSE]) > 0, , drop = FALSE])

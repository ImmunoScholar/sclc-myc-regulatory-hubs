# Which 9 segments are annotated in SegmentProperties but absent from the
# TargetCountMatrix of GSE261348? Read-only diagnostic.
suppressPackageStartupMessages(library(readxl))

f <- "data/raw/GSE261348/GSE261348_IMfirst_DSP_rawcounts.xlsx"
sp <- read_excel(f, sheet = "SegmentProperties")
cm <- read_excel(f, sheet = "TargetCountMatrix", n_max = 1)

cat("SegmentProperties rows :", nrow(sp), "\n")
cat("TargetCountMatrix ROIs :", ncol(cm) - 1, "\n\n")

cat("SegmentProperties columns:\n")
print(names(sp))

# Matrix column names look like "IMF-030 | 001 | Geometric Segment"
mat_ids <- names(cm)[-1]
cat("\nfirst 3 matrix column names:\n"); print(head(mat_ids, 3))

# Build the same composite key from SegmentProperties.
key_cols <- intersect(c("SlideName", "ROILabel", "SegmentLabel"), names(sp))
cat("\nkey columns used:", paste(key_cols, collapse = " | "), "\n")
sp_ids <- apply(sp[, key_cols], 1, function(r) paste(trimws(r), collapse = " | "))

cat("\nfirst 3 SegmentProperties keys:\n"); print(head(sp_ids, 3))

missing <- setdiff(sp_ids, mat_ids)
cat("\n=== in SegmentProperties but NOT in the count matrix:", length(missing), "===\n")
if (length(missing)) print(missing)

# What distinguishes them? Look for QC flag columns.
idx <- which(sp_ids %in% missing)
qc_cols <- grep("QC|Flag|Status|Pass|Fail|Exclud", names(sp), ignore.case = TRUE, value = TRUE)
cat("\nQC-ish columns present:", paste(qc_cols, collapse = ", "), "\n")
if (length(qc_cols) && length(idx)) {
  cat("\n--- values for the missing segments ---\n")
  print(as.data.frame(sp[idx, c(key_cols, qc_cols)]))
  cat("\n--- same columns, summary across ALL segments ---\n")
  for (cc in qc_cols) { cat("\n", cc, ":\n", sep = ""); print(table(sp[[cc]], useNA = "ifany")) }
}

# Useful counts regardless
for (cc in intersect(c("SegmentLabel", "ScanLabel", "SlideName"), names(sp))) {
  cat("\n", cc, " (missing segments only):\n", sep = "")
  print(table(sp[[cc]][idx], useNA = "ifany"))
}

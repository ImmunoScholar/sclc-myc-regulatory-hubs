# Compact per-dataset view of the manifest. Read-only.
m <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
m$gb <- ifelse(is.na(m$size_bytes_expected), 0, m$size_bytes_expected) / 1024^3

agg <- do.call(rbind, lapply(split(m, m$dataset_id), function(s) data.frame(
  dataset   = s$dataset_id[1],
  build     = paste(unique(s$genome_build), collapse = ","),
  lift      = if (any(s$liftover_required)) "LIFT" else "-",
  files     = nrow(s),
  GB        = round(sum(s$gb), 3),
  status    = paste(unique(s$status), collapse = ","),
  stringsAsFactors = FALSE
)))
agg <- agg[order(-agg$GB), ]
print(agg, row.names = FALSE)

cat("\ntotal files:", nrow(m), "\n")
cat("total GB (automated):", round(sum(m$gb[!m$status %in% c("manual_required","api")]), 2), "\n")
cat("files needing liftOver:", sum(m$liftover_required), "\n")
cat("\nstatus counts:\n"); print(table(m$status))

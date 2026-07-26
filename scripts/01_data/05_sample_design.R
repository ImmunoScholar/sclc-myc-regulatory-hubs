# -----------------------------------------------------------------------------
# 05_sample_design.R — derive the experimental design table from the manifest.
#
#   Rscript scripts/01_data/05_sample_design.R
#
# Parses cell line and assay target out of the deposited FILE NAMES in the
# manifest, so the design table has traceable provenance rather than being a
# hand-typed transcription of what I remember reading.
#
# Three facts cannot be parsed from file names and are declared explicitly with
# their source:
#   1. GSE210113's cell line (file names carry only NTgRNA / ND1-KO). H446 comes
#      from !Sample_source_name_ch1 in the GEO record.
#   2. GSE281523 is primary tissue / PDX, not a cell line at all.
#   3. MYC-paralog amplification status, which is author-declared in Plotnik
#      et al. 2024 and is NOT derivable from any file name.
#
# Writes: data/metadata/sample_design.tsv
# -----------------------------------------------------------------------------

man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)

# --- canonical cell-line names ------------------------------------------------
# The same line is spelled several ways across these deposits. Verified variants:
#   H1048 / NCIH1048 · H211 / NCIH211 · H526 / NCIH526 · SHP77 / SHP-77
ALIAS <- c(
  NCIH1048 = "H1048", `NCI-H1048` = "H1048", H1048 = "H1048",
  NCIH211  = "H211",  `NCI-H211`  = "H211",  H211  = "H211",
  NCIH526  = "H526",  `NCI-H526`  = "H526",  H526  = "H526",
  `SHP-77` = "SHP77", SHP77 = "SHP77",
  NCIH524  = "H524",  H524 = "H524",
  H69 = "H69", H847 = "H847", H889 = "H889", H196 = "H196",
  COLO668 = "COLO668", `COLO-668` = "COLO668",
  H1836 = "H1836", H446 = "H446", H146 = "H146", H82 = "H82",
  Lu139 = "Lu139"
)
canon_line <- function(x) {
  # longest alias first so NCIH1048 wins over H1048
  keys <- names(ALIAS)[order(-nchar(names(ALIAS)))]
  out <- rep(NA_character_, length(x))
  for (k in keys) {
    hit <- is.na(out) & grepl(k, x, fixed = TRUE)
    out[hit] <- ALIAS[[k]]
  }
  out
}

rows <- list()
addr <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

# =============================================================================
# GSE230649 — keystone. Target and line are both in the file name.
#   GSM7230503_MYC_H1048_R1_treat_pileup.bedgraph.gz     -> MYC   / H1048
#   GSM7230512_COLO668_treat_pileup.bedgraph.gz          -> ATAC  / COLO668
# ATAC samples carry no target token; that absence is the marker.
# =============================================================================
ks <- man[man$dataset_id == "GSE230649", ]
TARGETS <- c("H3K27Ac", "MYCL1", "MYCN", "MYC")   # MYCL1/MYCN before MYC
for (i in seq_len(nrow(ks))) {
  fn   <- ks$file_name[i]
  body <- sub("^GSM[0-9]+_", "", fn)
  tgt  <- NA_character_
  for (t in TARGETS) {
    if (startsWith(body, paste0(t, "_"))) { tgt <- t; break }
  }
  assay <- if (is.na(tgt)) "ATAC" else if (tgt == "H3K27Ac") "H3K27ac" else tgt
  addr(dataset = "GSE230649", gsm = ks$gsm[i], assay = assay,
       cell_line = canon_line(body), build = ks$genome_build[i],
       role = "keystone", source = "parsed from file name")
}

# =============================================================================
# Supporting / control datasets
# =============================================================================
for (ds in c("GSE269424", "GSE256345", "GSE281524", "GSE249362")) {
  s <- man[man$dataset_id == ds, ]
  if (!nrow(s)) next
  assay <- switch(ds,
    GSE269424 = "ATAC",      # EGFP control arms only (D-015)
    GSE256345 = "ATAC",
    GSE281524 = "ASCL1",
    GSE249362 = "POU2F3")
  role <- if (ds %in% c("GSE269424", "GSE256345")) "atac_support" else "lineage_tf"
  for (i in seq_len(nrow(s))) {
    addr(dataset = ds, gsm = s$gsm[i], assay = assay,
         cell_line = canon_line(s$file_name[i]), build = s$genome_build[i],
         role = role, source = "parsed from file name")
  }
}

# GSE210113: file names give only NTgRNA / ND1-KO. Line from the GEO record.
s <- man[man$dataset_id == "GSE210113", ]
if (nrow(s)) {
  addr(dataset = "GSE210113", gsm = NA_character_, assay = "NEUROD1",
       cell_line = "H446", build = "hg19", role = "lineage_tf",
       source = "!Sample_source_name_ch1 (not in file name)")
}

# GSE281523: primary tissue and PDX, not a cell line.
s <- man[man$dataset_id == "GSE281523", ]
for (i in seq_len(nrow(s))) {
  addr(dataset = "GSE281523", gsm = s$gsm[i], assay = "ATAC",
       cell_line = "(primary/PDX)", build = s$genome_build[i],
       role = "atac_support", source = "GEO record: primary lung / PDX")
}

design <- do.call(rbind, rows)

# =============================================================================
# MYC-paralog amplification status — author-declared, Plotnik et al. 2024.
# NOT derivable from any file name, and NOT taken from expression: Plotnik's own
# data show MYC-amplified behaves differently from MYC-expressing.
# =============================================================================
AMP <- c(
  H1048 = "MYC-amp", H211 = "MYC-amp", H524 = "MYC-amp",
  H847  = "MYC-amp", SHP77 = "MYC-amp",
  H526  = "MYCN-amp", H69 = "MYCN-amp",
  COLO668 = "MYCL-amp", H889 = "MYCL-amp",
  H196  = "non-amplified"
)
design$paralog_status <- ifelse(design$cell_line %in% names(AMP),
                               AMP[design$cell_line], NA_character_)

write.table(design, "data/metadata/sample_design.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("rows:", nrow(design), "\n\n")
cat("=== keystone GSE230649: assay x cell line ===\n")
print(table(design$assay[design$dataset == "GSE230649"],
            design$cell_line[design$dataset == "GSE230649"]))
cat("\n=== unresolved cell lines (should be none except primary/PDX) ===\n")
print(unique(design$cell_line[is.na(design$cell_line)]))
cat("\n=== lineage-TF coverage ===\n")
lt <- design[design$role == "lineage_tf", ]
print(table(lt$assay, lt$cell_line))
cat("\nwrote data/metadata/sample_design.tsv\n")

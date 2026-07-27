# -----------------------------------------------------------------------------
# 02_resolve_amplification.R — settle MYC-family amplification from real copy number.
#
#   Rscript scripts/04_functional/02_resolve_amplification.R
#
# THE QUESTION THIS ANSWERS, blocked since M5.
# Plotnik profiled "two cell lines harboring the alteration for each amplification
# type", but GSE230649 contains FIVE MYC ChIP samples (H1048, H211, H524, H847,
# SHP77). The paper never names which are amplified, GEO declares no amplification
# field, and our coverage-based proxy FAILED for MYC (D-026): the assays
# contradicted each other (H524 H3K27ac 82x vs ATAC 1.6x) and every MYC-locus
# ratio collapsed over a ±1 Mb window, indicating focal regulatory elements rather
# than copy gain.
#
# DepMap gives direct log2 relative copy number, which settles it.
#
# WHAT IT UNBLOCKS:
#   * M5 gate criterion 3 (distal-fraction contrast, MYC-amplified vs
#     MYC-EXPRESSING) — currently the only criterion marked NOT EVALUABLE
#   * the correctness of the project registry, which asserted all five MYC ChIP
#     lines were amplified — an assertion incompatible with Plotnik's own text
#
# It also tests our coverage proxy against ground truth: MYCN and MYCL1 were
# CONFIRMED by that method, so those calls should be reproduced here. If they are
# not, the proxy was wrong everywhere and D-026's confirmations need revisiting.
#
# Memory: the CN file is 256 MB x 18,620 columns. Only the needed columns are
# extracted with awk before anything enters R (R-13).
#
# Output: data/metadata/depmap_amplification.csv
#         results/tables/m7_amplification.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml); library(data.table) })
CFG <- yaml::read_yaml("config/params.yml")
DM  <- CFG$depmap
CN  <- file.path("data/raw/depmap", DM$files$copy_number_log2)
stopifnot(file.exists(CN))
cat("release: ", DM$release, "   file: ", basename(CN), "\n\n", sep = "")

# keystone line -> DepMap display name. Exact match only: a substring match would
# make "NCIH211" collide with "NCIH2110", a real and different line.
LINE_MAP <- c(COLO668 = "COLO668", H1048 = "NCIH1048", H196 = "NCIH196",
              H211 = "NCIH211", H524 = "NCIH524", H526 = "NCIH526",
              H69 = "NCIH69", H847 = "NCIH847", H889 = "NCIH889", SHP77 = "SHP77")
GENES <- c("MYC", "MYCN", "MYCL", "MYCL1")

# --- extract only the columns we need -----------------------------------------
# data.table::fread(select=), NOT awk -F','.
#
# The first version used awk with a comma field separator. That CANNOT parse this
# file: lineage fields contain QUOTED COMMAS (e.g. "Medulloblastoma | Large Cell
# Lung Carcinoma | ... , Breast Invasive Cancer"), so naive splitting shifted
# columns and mangled rows — 645 of 1,118 lines survived, and five keystone lines
# present in the file were reported as absent. fread is CSV-aware and reads only
# the requested columns.
META <- c("depmap_id","cell_line_display_name","lineage_1","lineage_2","lineage_3")
hdr <- names(data.table::fread(CN, nrows = 0))
gene_cols <- intersect(GENES, hdr)
cat("gene columns found: ", paste(gene_cols, collapse = ", "), "\n", sep = "")
if (!length(gene_cols)) stop("no MYC-family gene columns found in the CN file")

d <- as.data.frame(data.table::fread(
  CN, select = c(intersect(META, hdr), gene_cols), showProgress = FALSE))
cat("extracted: ", nrow(d), " lines x ", ncol(d), " columns\n", sep = "")

n_expected <- length(readLines(CN, n = -1L)) - 1L
if (nrow(d) < n_expected) {
  stop(sprintf("read %d rows but the file has %d — parsing is losing rows",
               nrow(d), n_expected))
}
cat("row count matches the file (", n_expected, ")\n\n", sep = "")

# --- SCLC availability, for later steps ---------------------------------------
lin_col <- intersect(c("lineage_2","lineage_3","lineage_1"), names(d))
sclc <- rep(FALSE, nrow(d))
for (cc in lin_col) sclc <- sclc | grepl("small cell lung", d[[cc]], ignore.case = TRUE)
cat("=========== SCLC lines in DepMap ", DM$release, " ===========\n", sep = "")
cat("lines matching 'small cell lung': ", sum(sclc), "\n", sep = "")
if (sum(sclc)) {
  cat("lineage labels used: ",
      paste(unique(unlist(lapply(lin_col, function(c) unique(d[[c]][sclc])))), collapse = " | "),
      "\n", sep = "")
}

# --- keystone lines ------------------------------------------------------------
cat("\n=========== keystone line lookup ===========\n")
idx <- match(LINE_MAP, d$cell_line_display_name)
found <- data.frame(keystone = names(LINE_MAP), depmap_name = unname(LINE_MAP),
                    row = idx, stringsAsFactors = FALSE)
found$depmap_id <- ifelse(is.na(idx), NA, d$depmap_id[idx])
for (i in seq_len(nrow(found)))
  cat(sprintf("  %-8s -> %-10s %s\n", found$keystone[i], found$depmap_name[i],
              ifelse(is.na(found$row[i]), "NOT IN DEPMAP", found$depmap_id[i])))
n_missing <- sum(is.na(found$row))
if (n_missing) cat("\n", n_missing, " keystone line(s) absent from this release.\n", sep = "")

# --- copy number ---------------------------------------------------------------
gcols <- intersect(GENES, names(d))
myc_l <- if ("MYCL" %in% gcols) "MYCL" else if ("MYCL1" %in% gcols) "MYCL1" else NA

cat("\n=========== log2 relative copy number ===========\n")
cat("DepMap OmicsCNGene is log2(relative CN + 1); 1.0 = neutral diploid.\n")
cat("Amplification is conventionally called above ~1.5 (i.e. >~2x relative CN).\n\n")

ok <- found[!is.na(found$row), ]
res <- data.frame()
for (i in seq_len(nrow(ok))) {
  r <- ok$row[i]
  vals <- vapply(c("MYC","MYCN",myc_l), function(g)
    if (!is.na(g) && g %in% names(d)) as.numeric(d[[g]][r]) else NA_real_, numeric(1))
  names(vals) <- c("MYC","MYCN","MYCL")
  best <- names(vals)[which.max(vals)]
  res <- rbind(res, data.frame(
    line = ok$keystone[i], depmap_id = ok$depmap_id[i],
    MYC = round(vals[["MYC"]], 3), MYCN = round(vals[["MYCN"]], 3),
    MYCL = round(vals[["MYCL"]], 3), highest = best,
    stringsAsFactors = FALSE))
}
res <- res[order(-res$MYC), ]
print(res, row.names = FALSE)

# --- the decisive comparison ---------------------------------------------------
AMP_THRESH <- 1.5
cat("\n=========== amplification calls (log2 > ", AMP_THRESH, ") ===========\n", sep = "")
for (g in c("MYC","MYCN","MYCL")) {
  hit <- res$line[!is.na(res[[g]]) & res[[g]] > AMP_THRESH]
  cat(sprintf("  %-5s amplified: %s\n", g,
              if (length(hit)) paste(hit, collapse = ", ") else "(none)"))
}

cat("\n=========== vs the project registry and D-026 ===========\n")
REG_MYC  <- c("H1048","H211","H524","H847","SHP77")
REG_MYCN <- c("H526","H69"); REG_MYCL <- c("COLO668","H889")
obs_myc  <- res$line[!is.na(res$MYC)  & res$MYC  > AMP_THRESH]
obs_mycn <- res$line[!is.na(res$MYCN) & res$MYCN > AMP_THRESH]
obs_mycl <- res$line[!is.na(res$MYCL) & res$MYCL > AMP_THRESH]

cat("registry claimed MYC-amp : ", paste(REG_MYC, collapse = ", "), "\n", sep = "")
cat("copy number says MYC-amp : ", if (length(obs_myc)) paste(obs_myc, collapse=", ") else "(none)", "\n", sep = "")

# A line absent from this release is UNKNOWN, not "confirmed not amplified".
# Assigning it to the comparator by exclusion would be inference from absence.
measured   <- res$line
unmeasured <- setdiff(REG_MYC, measured)
expressing <- setdiff(intersect(REG_MYC, measured), obs_myc)
cat("  -> MYC-EXPRESSING (measured, not amplified): ",
    if (length(expressing)) paste(expressing, collapse = ", ") else "(none)", "\n", sep = "")
if (length(unmeasured))
  cat("  -> STATUS UNKNOWN (absent from ", DM$release, "): ",
      paste(unmeasured, collapse = ", "),
      "\n     NOT assigned to either group — absence is not evidence.\n", sep = "")
cat("\nD-026 confirmed by coverage proxy — does copy number agree?\n")
cat("  MYCN (H526, H69)      : ", if (setequal(obs_mycn, REG_MYCN)) "AGREES" else
    paste0("DIFFERS -> ", paste(obs_mycn, collapse=", ")), "\n", sep = "")
cat("  MYCL (COLO668, H889)  : ", if (setequal(obs_mycl, REG_MYCL)) "AGREES" else
    paste0("DIFFERS -> ", paste(obs_mycl, collapse=", ")), "\n", sep = "")

n_myc_amp <- length(obs_myc)
cat("\n=========== does this match Plotnik's design? ===========\n")
cat("Plotnik: \"two cell lines harboring the alteration for each amplification type\"\n")
cat("observed MYC-amplified among the 5 MYC-ChIP lines: ", n_myc_amp, "\n", sep = "")
if (n_myc_amp == 2) {
  cat("EXACTLY TWO — consistent with the paper. Gate criterion 3 is now evaluable:\n")
  cat("  MYC-amplified group : ", paste(obs_myc, collapse = ", "), "\n", sep = "")
  cat("  MYC-expressing group: ", paste(setdiff(REG_MYC, obs_myc), collapse = ", "), "\n", sep = "")
} else {
  cat("NOT two. Either the amplification threshold needs revisiting, or Plotnik\n")
  cat("used lines outside this deposit. Report the copy-number values and define\n")
  cat("the groups from them explicitly rather than assuming the paper's split.\n")
}

write.csv(res, "data/metadata/depmap_amplification.csv", row.names = FALSE)
md <- c("# M7 — MYC-family amplification from DepMap copy number", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste0("Release: **", DM$release, "**, log2 relative copy number."), "",
        "Settles the question blocked since M5: our coverage-based proxy could not",
        "resolve MYC amplification (D-026), and the project registry's claim that all",
        "five MYC ChIP lines were amplified is incompatible with Plotnik's own text.", "",
        "| line | DepMap ID | MYC | MYCN | MYCL | highest |", "|---|---|---|---|---|---|")
for (i in seq_len(nrow(res)))
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |", res$line[i], res$depmap_id[i],
                      res$MYC[i], res$MYCN[i], res$MYCL[i], res$highest[i]))
writeLines(md, "results/tables/m7_amplification.md")
cat("\nwrote data/metadata/depmap_amplification.csv\n")
cat("wrote results/tables/m7_amplification.md\n")

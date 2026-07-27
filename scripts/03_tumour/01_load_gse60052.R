# -----------------------------------------------------------------------------
# 01_load_gse60052.R — load and QC the primary tumour cohort.
#
#   Rscript scripts/03_tumour/01_load_gse60052.R
#
# GSE60052: 79 SCLC tumours + 7 normal lung. THREE TRAPS characterised at M4
# (risk R-15) and handled explicitly here rather than discovered again:
#
#   1. ALREADY log2-NORMALISED. No raw counts are deposited, so DESeq2 cannot be
#      used on this cohort. Differential analysis uses limma on the log2 matrix.
#   2. SAMPLE NAMES CARRY LEADING WHITESPACE. Untrimmed, every downstream join
#      fails silently by not matching.
#   3. NORMALS ARE ENCODED IN THE NAME, by a '.normal' suffix. There is no
#      separate class column; miss the suffix and 7 normals are analysed as
#      tumours.
#
# Gene identifiers are SYMBOLS, which must be reconciled against the regulon
# symbols from org.Hs.eg.db. Coverage is reported, because a regulon whose genes
# are largely absent from the cohort cannot be scored in it.
#
# Output: data/processed/tumour/gse60052.rds
#         data/metadata/gse60052_qc.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml) })

CFG <- yaml::read_yaml("config/params.yml")
OUT <- "data/processed/tumour"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

f <- "data/raw/GSE60052/GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv.gz"
stopifnot(file.exists(f))

cat("=========== load ===========\n")
d <- utils::read.delim(gzfile(f), header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE)
cat("raw dimensions: ", nrow(d), " rows x ", ncol(d), " cols\n", sep = "")

# --- trap 2: whitespace in sample names ---------------------------------------
raw_names <- names(d)
names(d)   <- trimws(names(d))
n_ws <- sum(raw_names != names(d))
cat("sample names carrying whitespace: ", n_ws, " (trimmed)\n", sep = "")

gene_col <- names(d)[1]
genes <- trimws(d[[gene_col]])
mat <- as.matrix(d[, -1, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- genes

# --- duplicate symbols --------------------------------------------------------
dup <- sum(duplicated(genes))
if (dup > 0) {
  cat("duplicate gene symbols: ", dup, " — collapsing by max mean expression\n", sep = "")
  o <- order(rowMeans(mat, na.rm = TRUE), decreasing = TRUE)
  mat <- mat[o, , drop = FALSE]
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
}
cat("genes after de-duplication: ", nrow(mat), "\n", sep = "")

# --- trap 3: normals encoded by suffix ----------------------------------------
is_normal <- grepl("\\.normal$", colnames(mat))
n_norm <- sum(is_normal); n_tum <- sum(!is_normal)
cat("\n=========== sample classes ===========\n")
cat("tumour: ", n_tum, "   normal: ", n_norm, "\n", sep = "")
if (n_tum != 79 || n_norm != 7)
  stop(sprintf("expected 79 tumour / 7 normal, got %d / %d — do not proceed", n_tum, n_norm))
cat("matches the deposit's stated 79 + 7\n")

# --- trap 1: confirm log2 scale -----------------------------------------------
mx <- max(mat, na.rm = TRUE); mn <- min(mat, na.rm = TRUE)
cat("\n=========== value scale ===========\n")
cat("range: ", round(mn, 2), " to ", round(mx, 2), "\n", sep = "")
if (mx > 100) stop("values exceed 100 — this does not look log2-transformed")
cat("consistent with log2 (max < 100). DESeq2 is NOT usable; limma only.\n")
pct_zero <- mean(mat == 0, na.rm = TRUE)
cat("zero entries: ", sprintf("%.1f%%", 100 * pct_zero), "\n", sep = "")

# --- expression filter --------------------------------------------------------
minfrac <- CFG$tumour_scoring$min_expressed_fraction
expressed <- rowMeans(mat > 0, na.rm = TRUE) >= minfrac
cat("\ngenes expressed in >=", sprintf("%.0f%%", 100 * minfrac), " of samples: ",
    format(sum(expressed), big.mark = ","), " of ", format(nrow(mat), big.mark = ","),
    "\n", sep = "")
mat_f <- mat[expressed, , drop = FALSE]

# --- regulon coverage ---------------------------------------------------------
cat("\n=========== regulon coverage in this cohort ===========\n")
reg <- readRDS("data/processed/regions/regulons.rds")$regulons
cov <- data.frame()
for (p in names(reg)) {
  g <- reg[[p]]
  n_all <- sum(g %in% rownames(mat))
  n_exp <- sum(g %in% rownames(mat_f))
  cov <- rbind(cov, data.frame(paralog = p, regulon_size = length(g),
                               present = n_all, present_expressed = n_exp,
                               pct_expressed = round(100 * n_exp / length(g), 1)))
  cat(sprintf("  %-6s %3d genes -> %3d present, %3d expressed (%.1f%%)\n",
              p, length(g), n_all, n_exp, 100 * n_exp / length(g)))
}
if (any(cov$present_expressed < 100))
  cat("\nWARNING: a regulon with <100 measurable genes will give noisy scores.\n")

# --- paralog expression, for grouping tumours later --------------------------
cat("\n=========== MYC-family expression in tumours ===========\n")
for (p in c("MYC","MYCN","MYCL1","MYCL")) {
  if (p %in% rownames(mat)) {
    v <- mat[p, !is_normal]
    cat(sprintf("  %-6s median %6.2f  IQR %6.2f-%6.2f  max %6.2f\n",
                p, stats::median(v), stats::quantile(v, .25), stats::quantile(v, .75), max(v)))
  } else cat(sprintf("  %-6s NOT FOUND in the matrix\n", p))
}

cat("\n=========== lineage-TF expression (subtype markers) ===========\n")
for (p in c("ASCL1","NEUROD1","POU2F3","YAP1")) {
  if (p %in% rownames(mat)) {
    v <- mat[p, !is_normal]
    cat(sprintf("  %-8s median %6.2f  range %6.2f-%6.2f\n", p, stats::median(v), min(v), max(v)))
  } else cat(sprintf("  %-8s NOT FOUND\n", p))
}

qc <- data.frame(
  metric = c("genes_raw","genes_dedup","genes_expressed","n_tumour","n_normal",
             "value_min","value_max","pct_zero","whitespace_names_trimmed"),
  value = c(nrow(d), nrow(mat), sum(expressed), n_tum, n_norm,
            round(mn,2), round(mx,2), round(100*pct_zero,1), n_ws),
  stringsAsFactors = FALSE)
write.csv(qc, "data/metadata/gse60052_qc.csv", row.names = FALSE)
write.csv(cov, "data/metadata/gse60052_regulon_coverage.csv", row.names = FALSE)

saveRDS(list(expr = mat, expr_filtered = mat_f, is_normal = is_normal,
             coverage = cov, note = "log2-normalised; limma not DESeq2 (R-15)"),
        file.path(OUT, "gse60052.rds"))

cat("\nwrote data/processed/tumour/gse60052.rds\n")
cat("wrote data/metadata/gse60052_qc.csv\n")
cat("\nRESULT: PASS — cohort loaded, all three traps handled.\n")
cat("Next: Rscript scripts/03_tumour/02_score_regulons.R\n")

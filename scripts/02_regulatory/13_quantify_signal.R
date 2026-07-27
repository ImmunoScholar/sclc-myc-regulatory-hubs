# -----------------------------------------------------------------------------
# 13_quantify_signal.R — signal for all 28 GSE230649 tracks over the universe.
#
#   Rscript scripts/02_regulatory/13_quantify_signal.R
#
# Produces the matrix everything downstream depends on: 102,334 regions x 28
# samples (10 H3K27ac, 5 MYC, 2 MYCN, 2 MYCL1, 9 ATAC).
#
# ALGORITHM. Both inputs are position-sorted within a chromosome and the regions
# are non-overlapping, so this is a linear MERGE-JOIN with a moving pointer per
# chromosome — O(n+m), no interval tree, no genome-sized vector. A bucketed
# lookup was the obvious alternative but would have cost more memory and time for
# no benefit given the data are already sorted.
#
# MEMORY. Nothing genome-scale enters R. awk holds ~102k region records and emits
# one line per region per file. This matters on a 10 GB machine (R-13).
#
# SEQUENCE NAMES. The universe is UCSC-named; the bedGraphs are Ensembl-named
# (R-16). The join file is written in ENSEMBL naming deliberately. Getting this
# wrong yields an empty matrix rather than an error.
#
# RAW VALUES ONLY. Stores sum(signal x overlap_bp) and its per-bp mean, with no
# normalisation. Cross-line normalisation is a separate, explicit step (14_)
# because paralog groups are perfectly confounded with cell lines — MYCN is only
# ever H526+H69, MYCL1 only ever COLO668+H889 — so the normalisation choice can
# manufacture or erase a paralog difference. Baking it into the expensive step
# would hide that choice; keeping it separate makes it auditable and re-runnable.
#
# Output: data/processed/signal/region_signal_raw.rds   (list: sum, mean, meta)
#         data/metadata/signal_quantification_report.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
})
source("R/genome_utils.R")

OUT <- "data/processed/signal"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TMP <- file.path(OUT, "_tmp"); dir.create(TMP, showWarnings = FALSE)

universe <- readRDS("data/processed/regions/universe_final.rds")
n_reg <- length(universe)
cat("universe: ", format(n_reg, big.mark = ","), " regions\n", sep = "")

# --- join file, Ensembl-named, with a stable integer id -----------------------
join_bed <- file.path(TMP, "universe_ens.bed")
if (!file.exists(join_bed)) {
  d <- data.frame(
    chrom = to_ensembl_seqnames(as.character(seqnames(universe))),
    start = start(universe) - 1L,          # back to 0-based half-open
    end   = end(universe),
    id    = seq_len(n_reg))
  d <- d[order(d$chrom, d$start), ]        # merge-join needs sorted input
  utils::write.table(d, join_bed, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
}
cat("join file: ", join_bed, "\n\n", sep = "")

AWK <- '
NR==FNR {
  c=$1; n=++cnt[c]
  rs[c,n]=$2; re[c,n]=$3; rid[c,n]=$4
  next
}
{
  c=$1
  if (!(c in cnt)) next
  v=$4+0
  if (v<=0) next
  s=$2; e=$3
  p = (c in ptr) ? ptr[c] : 1
  while (p <= cnt[c] && re[c,p] <= s) p++
  ptr[c]=p
  q=p
  while (q <= cnt[c] && rs[c,q] < e) {
    os = (s > rs[c,q]) ? s : rs[c,q]
    oe = (e < re[c,q]) ? e : re[c,q]
    if (oe > os) sum[rid[c,q]] += v*(oe-os)
    q++
  }
}
END { for (k in sum) print k"\t"sum[k] }
'

man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
man$dest <- file.path("data/raw", man$dataset_id, man$file_name)
ks <- man[man$dataset_id == "GSE230649" & file.exists(man$dest), ]
stopifnot(nrow(ks) == 28)

parse_sample <- function(fn) {
  body <- sub("^GSM[0-9]+_", "", fn)
  for (t in c("H3K27Ac", "MYCL1", "MYCN", "MYC")) {
    if (startsWith(body, paste0(t, "_"))) {
      ln <- sub(paste0("^", t, "_([^_]+)_.*$"), "\\1", body)
      return(c(assay = if (t == "H3K27Ac") "H3K27ac" else t, line = ln))
    }
  }
  c(assay = "ATAC", line = sub("^([^_]+)_treat_pileup.*$", "\\1", body))
}
meta <- as.data.frame(t(vapply(ks$file_name, parse_sample, character(2))),
                      stringsAsFactors = FALSE)
meta$gsm  <- ks$gsm
meta$file <- ks$file_name
meta$sample <- paste0(meta$assay, "_", meta$line)
rownames(meta) <- NULL
cat("samples by assay:\n"); print(table(meta$assay)); cat("\n")

sum_mat <- matrix(0, nrow = n_reg, ncol = nrow(meta),
                  dimnames = list(NULL, meta$sample))

for (i in seq_len(nrow(meta))) {
  sm  <- meta$sample[i]
  cf  <- file.path(TMP, paste0(sm, ".tsv"))
  if (!file.exists(cf)) {
    t0 <- Sys.time()
    cmd <- sprintf("zcat %s | awk -F'\\t' %s %s - > %s",
                   shQuote(ks$dest[i]), shQuote(AWK), shQuote(join_bed), shQuote(cf))
    st <- system(cmd)
    if (st != 0 || !file.exists(cf)) { cat("  FAILED: ", sm, "\n"); next }
    cat(sprintf("  %-18s %.1f min\n", sm,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  } else {
    cat(sprintf("  %-18s [cached]\n", sm))
  }
  d <- utils::read.delim(cf, header = FALSE, col.names = c("id", "sumsw"))
  sum_mat[d$id, i] <- d$sumsw
}

w <- width(universe)
mean_mat <- sum_mat / w

res <- list(sum = sum_mat, mean = mean_mat, meta = meta,
            regions = universe, width = w,
            note = "RAW signal. No normalisation applied — see 14_.")
saveRDS(res, file.path(OUT, "region_signal_raw.rds"))

# --- report -------------------------------------------------------------------
rep <- data.frame(
  sample = meta$sample, assay = meta$assay, line = meta$line,
  pct_regions_nonzero = round(100 * colMeans(sum_mat > 0), 1),
  median_mean_signal = round(apply(mean_mat, 2, median), 3),
  mean_mean_signal = round(colMeans(mean_mat), 3),
  max_mean_signal = round(apply(mean_mat, 2, max), 1),
  stringsAsFactors = FALSE)
write.csv(rep, "data/metadata/signal_quantification_report.csv", row.names = FALSE)

cat("\n=========== signal per sample ===========\n")
print(rep[order(rep$assay, rep$line), ], row.names = FALSE)

cat("\n=========== the confound to watch ===========\n")
cat("Paralog groups are perfectly confounded with cell lines:\n")
cat("  MYC   : H1048 H211 H524 H847 SHP77\n")
cat("  MYCN  : H526 H69\n")
cat("  MYCL1 : COLO668 H889\n")
cat("A track-quality difference between those pairs is indistinguishable from a\n")
cat("paralog biology difference. Compare median_mean_signal within assay below:\n\n")
for (a in c("MYC", "MYCN", "MYCL1", "H3K27ac", "ATAC")) {
  s <- rep[rep$assay == a, ]
  if (!nrow(s)) next
  cat(sprintf("  %-8s median across lines: %s   (range %.3f - %.3f, %.1fx spread)\n",
              a, paste(sprintf("%.2f", s$median_mean_signal), collapse = " "),
              min(s$median_mean_signal), max(s$median_mean_signal),
              max(s$median_mean_signal) / max(min(s$median_mean_signal), 1e-9)))
}

cat("\nwrote data/processed/signal/region_signal_raw.rds\n")
cat("wrote data/metadata/signal_quantification_report.csv\n")
cat("\nNext: Rscript scripts/02_regulatory/14_normalise_signal.R\n")

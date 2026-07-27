# -----------------------------------------------------------------------------
# 07_local_background.R — re-call ATAC regions using LOCAL background.
#
#   Rscript scripts/02_regulatory/07_local_background.R [--chrom chr1]
#
# WHY THIS EXISTS
#
# 05_ called regions by thresholding raw pileup against a GLOBAL mean. Scoring
# against H524's MACS2 peaks (06_) showed recall pinned at ~0.35 regardless of any
# support rule, with ~60% of our base pairs outside real peaks. Filtering could
# not fix it, which ruled out "we keep too much noise" and pointed at the
# thresholding itself.
#
# The mechanism: these are cancer lines with focal amplifications. An amplicon has
# more DNA, therefore more reads, therefore higher pileup — indistinguishable from
# accessibility under a global threshold. Worse, amplicons drag the genome-wide
# mean upward, raising the threshold everywhere and suppressing genuine peaks in
# normal-copy regions. That produces exactly the observed signature: inflated
# calls in some places, depressed recall everywhere.
#
# WHY IT MATTERS FOR THIS PROJECT SPECIFICALLY. Paralog groups ARE amplification
# groups. If region calling partly tracks copy number, "paralog-specific regions"
# become partly "amplicon-specific regions", and the MYCN-in-MYC overlap at the
# centre of the M5 gate would be measuring aneuploidy. It would look clean.
#
# THE FIX. Fold-enrichment over a LOCAL windowed background, the same principle as
# MACS2's local lambda. Copy number varies at megabase scale, so background inside
# an amplicon is elevated too and largely cancels.
#
#   FE(region) = mean_signal(region) / max(global_mean, local_mean_10kb, local_mean_50kb)
#
# Run on one chromosome first to test the idea cheaply before committing.
#
# Output: data/processed/regions/localbg/<line>_<chrom>.rds
#         data/metadata/local_background_comparison.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(yaml)
})
source("R/genome_utils.R")

args  <- commandArgs(trailingOnly = TRUE)
CHROM <- if (length(args) >= 2 && args[1] == "--chrom") args[2] else "chr1"
CHROM_ENS <- sub("^chr", "", CHROM)

CFG <- yaml::read_yaml("config/params.yml")
RP  <- CFG$regions
OUT <- "data/processed/regions/localbg"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

BIN_SMALL <- 10000L
BIN_LARGE <- 50000L
CAND_MULT <- 2          # permissive candidate cutoff; FE does the real filtering
FE_GRID   <- c(2, 3, 4, 5, 6, 8)

man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
man$dest <- file.path("data/raw", man$dataset_id, man$file_name)
ks   <- man[man$dataset_id == "GSE230649" & file.exists(man$dest), ]
atac <- ks[!grepl("^GSM[0-9]+_(H3K27Ac|MYC|MYCN|MYCL1)_", ks$file_name), ]
atac$line <- sub("^GSM[0-9]+_([^_]+)_treat_pileup.*$", "\\1", atac$file_name)

# Calibrate on H524, the only line with an independent MACS2 reference.
lf  <- list.files("data/processed/liftover", pattern = "H524.*\\.hg19\\.rds$", full.names = TRUE)
ref <- GenomicRanges::reduce(do.call(c, unname(lapply(lf, readRDS))), min.gapwidth = RP$merge_gap)
ref <- ref[as.character(seqnames(ref)) == CHROM]
cat("reference (H524 MACS2, ", CHROM, "): ", length(ref), " peaks, ",
    format(sum(width(ref)), big.mark = ","), " bp\n\n", sep = "")

# --- one pass: binned background AND candidate intervals with mean signal -----
# awk merges contiguous above-cutoff runs itself so R never sees millions of rows.
scan_chrom <- function(path, cutoff) {
  awk <- sprintf(paste0(
    "awk -F'\\t' '$1==\"%s\" {",
    "  w=$3-$2; v=$4+0;",
    "  b1=int($2/%d); s1[b1]+=v*w;",
    "  b2=int($2/%d); s2[b2]+=v*w;",
    "  if (v >= %f) {",
    "    if (open && $2==ce) { ce=$3; sw+=v*w }",
    "    else { if (open) printf \"R\\t%%d\\t%%d\\t%%.4f\\n\", cs, ce, sw;",
    "           cs=$2; ce=$3; sw=v*w; open=1 }",
    "  } else if (open) { printf \"R\\t%%d\\t%%d\\t%%.4f\\n\", cs, ce, sw; open=0 }",
    "}",
    "END { if (open) printf \"R\\t%%d\\t%%d\\t%%.4f\\n\", cs, ce, sw;",
    "      for (k in s1) printf \"S\\t%%d\\t%%.4f\\n\", k, s1[k];",
    "      for (k in s2) printf \"L\\t%%d\\t%%.4f\\n\", k, s2[k] }'"),
    CHROM_ENS, BIN_SMALL, BIN_LARGE, cutoff)
  out <- system(paste("zcat", shQuote(path), "|", awk), intern = TRUE)
  out
}

line_mean <- function(line) {
  raw <- readLines(file.path("data/processed/atac_profile", paste0(line, ".hist.tsv")))
  h <- utils::read.delim(text = paste(raw[!grepl("^#", raw)], collapse = "\n"),
                         header = FALSE, col.names = c("bin", "bp"))
  h$signal <- h$bin / 2
  cov <- h[h$signal > 0, ]
  sum(cov$signal * cov$bp) / sum(cov$bp)
}

res <- list()

for (i in seq_len(nrow(atac))) {
  ln <- atac$line[i]
  gm <- line_mean(ln)
  t0 <- Sys.time()
  out <- scan_chrom(atac$dest[i], CAND_MULT * gm)

  rr <- out[startsWith(out, "R\t")]
  ss <- out[startsWith(out, "S\t")]
  ll <- out[startsWith(out, "L\t")]
  if (!length(rr)) { cat("  ", ln, ": no candidates\n"); next }

  cand <- utils::read.delim(text = paste(rr, collapse = "\n"), header = FALSE,
                            col.names = c("tag","start","end","sumsw"))
  s1 <- utils::read.delim(text = paste(ss, collapse = "\n"), header = FALSE,
                          col.names = c("tag","bin","sumsw"))
  s2 <- utils::read.delim(text = paste(ll, collapse = "\n"), header = FALSE,
                          col.names = c("tag","bin","sumsw"))

  # local background = mean signal per bp within the window
  bg_small <- stats::setNames(s1$sumsw / BIN_SMALL, s1$bin)
  bg_large <- stats::setNames(s2$sumsw / BIN_LARGE, s2$bin)

  cand$width <- cand$end - cand$start
  cand <- cand[cand$width >= RP$min_width, ]
  cand$mean_sig <- cand$sumsw / cand$width
  b1 <- as.character(floor(cand$start / BIN_SMALL))
  b2 <- as.character(floor(cand$start / BIN_LARGE))
  loc <- pmax(gm,
              ifelse(is.na(bg_small[b1]), 0, bg_small[b1]),
              ifelse(is.na(bg_large[b2]), 0, bg_large[b2]), na.rm = TRUE)
  cand$fe <- cand$mean_sig / loc

  for (fe in FE_GRID) {
    sel <- cand[cand$fe >= fe, ]
    if (!nrow(sel)) next
    gr <- GRanges(CHROM, IRanges(sel$start + 1L, sel$end))
    gr <- GenomicRanges::reduce(gr, min.gapwidth = RP$merge_gap)
    gr <- gr[width(gr) >= RP$min_width & width(gr) <= RP$max_width]
    if (ln == "H524" && length(gr)) {
      inter <- sum(width(GenomicRanges::intersect(gr, ref)))
      prec  <- inter / sum(width(gr)); rec <- inter / sum(width(ref))
      f1    <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
      jac   <- inter / sum(width(GenomicRanges::union(gr, ref)))
    } else prec <- rec <- f1 <- jac <- NA_real_
    res[[length(res) + 1L]] <- data.frame(
      line = ln, fe_cutoff = fe, n_regions = length(gr),
      median_width = median(width(gr)),
      mb = round(sum(as.numeric(width(gr))) / 1e6, 2),
      precision = round(prec, 4), recall = round(rec, 4),
      f1 = round(f1, 4), jaccard = round(jac, 4), stringsAsFactors = FALSE)
    if (ln == "H524") saveRDS(gr, file.path(OUT, sprintf("H524_%s_fe%g.rds", CHROM, fe)))
  }
  cat(sprintf("  %-9s candidates=%6d  (%.1f min)\n", ln, nrow(cand),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cmp <- do.call(rbind, res)
write.csv(cmp, "data/metadata/local_background_comparison.csv", row.names = FALSE)

cat("\n=========== H524: local background vs MACS2 reference ===========\n")
h <- cmp[cmp$line == "H524", ]
print(h, row.names = FALSE)
cat("\nglobal-threshold baseline (06_): precision 0.398, recall 0.349, F1 0.372\n")
if (nrow(h)) {
  b <- h[which.max(h$f1), ]
  cat("best local-background F1: ", b$f1, " at FE >= ", b$fe_cutoff,
      "  (precision ", b$precision, ", recall ", b$recall, ")\n", sep = "")
  cat("\nVERDICT: ", if (b$f1 > 0.45) "local background materially improves agreement — adopt."
      else if (b$f1 > 0.40) "modest improvement — worth adopting but not decisive."
      else "no material improvement — the disagreement is NOT copy-number driven; investigate elsewhere.",
      "\n", sep = "")
}

cat("\n=========== region counts by line (FE grid) ===========\n")
print(reshape(cmp[, c("line","fe_cutoff","n_regions")], idvar = "line",
              timevar = "fe_cutoff", direction = "wide"), row.names = FALSE)
cat("\nWatch H526/H69 (MYCN-amplified, highest global-threshold counts) against\n")
cat("H196 (non-amplified). If local background is doing its job, the spread\n")
cat("between them should NARROW relative to the global-threshold counts\n")
cat("(96,200 / 104,017 vs 31,856 genome-wide).\n")

cat("\nwrote data/metadata/local_background_comparison.csv\n")

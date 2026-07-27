# -----------------------------------------------------------------------------
# 18_super_enhancers.R — ROSE-style super-enhancer calling, copy-number aware.
#
#   Rscript scripts/02_regulatory/18_super_enhancers.R
#
# This is genuine white space relative to the prior work: Plotnik performed a 1 kb
# promoter/distal binary split and no super-enhancer analysis at all (D-001).
#
# METHOD (Whyte et al. 2013 / ROSE):
#   1. H3K27ac-positive universe regions become enhancer constituents
#   2. TSS-proximal constituents excluded (config exclude_tss_window), so SEs are
#      not merely large promoters
#   3. constituents stitched within config stitch_distance
#   4. stitched enhancers ranked by total H3K27ac signal
#   5. cut at the inflection point where the ranked-signal curve reaches slope 1
#
# THE COPY-NUMBER PROBLEM, HANDLED UP FRONT RATHER THAN DISCOVERED LATER.
# D-026 measured H3K27ac amplicon coverage at 38-71x for confirmed MYCN and MYCL1
# amplifications. Signal that high will dominate ANY ranking, so amplicon-resident
# stitched enhancers would be called super-enhancers because of DNA copy number,
# not regulatory investment. That is exactly the artefact that broke the
# local-background attempt at the universe stage.
#
# We cannot correct it: D-026 established that our copy-number proxy fails for MYC
# (assays contradict; focal signal collapses over a broad window). Inventing a
# correction we cannot validate would be worse than declaring the limitation. So:
#   * every SE is FLAGGED for overlap with a confirmed amplicon in that same line
#   * amplicon-resident SEs are retained but never reported as regulatory findings
#     without the flag
#   * cross-line reproducibility is reported, since an SE recurring in lines with
#     DIFFERENT amplicons is unlikely to be a copy-number artefact
#
# Output: data/processed/regions/super_enhancers.rds
#         data/metadata/super_enhancer_summary.csv
#         results/tables/m5_super_enhancers.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(yaml)
})
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
SEP <- CFG$super_enhancers
FOLD_K27 <- 2          # H3K27ac-positive, same cutoff as the paralog sets

nrm  <- readRDS("data/processed/signal/region_signal_normalised.rds")
M    <- nrm$mean_fob; meta <- nrm$meta; u <- nrm$regions
w    <- nrm$width

k27_lines <- meta$line[meta$assay == "H3K27ac"]
cat("H3K27ac tracks: ", length(k27_lines), " (", paste(k27_lines, collapse = ", "), ")\n\n", sep = "")

# --- TSS exclusion ------------------------------------------------------------
tx  <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tss <- promoters(tx, upstream = 0, downstream = 1)
tss <- tss[as.character(seqnames(tss)) %in% ANALYSIS_CHROMS_UCSC]
tss_excl <- GenomicRanges::reduce(resize(tss, 2 * SEP$exclude_tss_window, fix = "center"))
is_tss <- IRanges::overlapsAny(u, tss_excl)
cat("constituents excluded as TSS-proximal (+/-", SEP$exclude_tss_window, "): ",
    sum(is_tss), " of ", length(u), "\n\n", sep = "")

# --- confirmed amplicons (D-026) ---------------------------------------------
# MYC deliberately absent: amplification unverified, so no line can be flagged
# for it without asserting something D-026 explicitly could not establish.
AMP <- list(
  list(locus = "MYCN",  gr = GRanges("chr2", IRanges(16080683, 16087129)), lines = c("H526","H69")),
  list(locus = "MYCL1", gr = GRanges("chr1", IRanges(40361098, 40367687)), lines = c("COLO668","H889")))
AMP_WIN <- 1000000
amp_gr <- lapply(AMP, function(a) GenomicRanges::resize(a$gr, width(a$gr) + 2 * AMP_WIN, fix = "center"))
names(amp_gr) <- vapply(AMP, function(a) a$locus, character(1))

# --- ROSE inflection cutoff --------------------------------------------------
# Knee-point (maximum distance from the chord joining the endpoints), NOT the
# first-slope->=1 rule.
#
# The slope rule is numerically fragile: with n stitched regions one step in x is
# 1/n, so any two adjacent low values differing by more than (1/n) of the maximum
# trip it immediately and the cutoff collapses to the minimum. That is exactly
# what happened to H1048 — all 8,238 stitched regions were called
# super-enhancers, median width 986 bp — and it inflated the cross-line union
# from ~5,100 to ~11,000 loci.
#
# The knee point uses the whole curve rather than one local difference, so
# isolated noise near the bottom cannot move it.
rose_cutoff <- function(sig) {
  s <- sort(sig)                      # ascending
  n <- length(s)
  if (n < 20) return(stats::quantile(s, 0.9, names = FALSE))
  rng <- max(s) - min(s)
  if (rng <= 0) return(max(s))
  x <- (seq_len(n) - 1) / (n - 1)
  y <- (s - min(s)) / rng
  # perpendicular distance from the chord (0,0)-(1,1) is proportional to |x - y|;
  # the knee is where the curve sits furthest BELOW the chord
  i <- which.max(x - y)
  s[i]
}

se_all <- list(); summ <- list()

for (ln in k27_lines) {
  j <- which(meta$assay == "H3K27ac" & meta$line == ln)
  if (!length(j)) next

  keep <- (M[, j] >= FOLD_K27) & !is_tss
  if (sum(keep) < 50) { cat("  ", ln, ": too few constituents\n"); next }

  cons <- u[keep]
  mcols(cons)$sig <- M[keep, j] * w[keep]      # signal x bp

  # stitch
  st <- GenomicRanges::reduce(cons, min.gapwidth = SEP$stitch_distance)
  ov <- GenomicRanges::findOverlaps(cons, st)
  agg <- tapply(mcols(cons)$sig[queryHits(ov)], subjectHits(ov), sum)
  mcols(st)$total_signal <- 0
  mcols(st)$total_signal[as.integer(names(agg))] <- as.numeric(agg)
  ncons <- tapply(queryHits(ov), subjectHits(ov), length)
  mcols(st)$n_constituents <- 0L
  mcols(st)$n_constituents[as.integer(names(ncons))] <- as.integer(ncons)

  cut <- rose_cutoff(mcols(st)$total_signal)
  mcols(st)$is_SE <- mcols(st)$total_signal >= cut
  mcols(st)$rank  <- rank(-mcols(st)$total_signal, ties.method = "first")
  mcols(st)$line  <- ln

  # amplicon flag: only for amplicons CONFIRMED in this line
  flag <- rep(NA_character_, length(st))
  for (a in AMP) if (ln %in% a$lines) {
    hit <- IRanges::overlapsAny(st, amp_gr[[a$locus]])
    flag[hit] <- ifelse(is.na(flag[hit]), a$locus, paste(flag[hit], a$locus, sep = ";"))
  }
  mcols(st)$amplicon <- flag

  se <- st[mcols(st)$is_SE]
  se_all[[ln]] <- se

  n_amp <- sum(!is.na(mcols(se)$amplicon))
  top10 <- head(order(-mcols(se)$total_signal), 10)
  n_amp_top10 <- sum(!is.na(mcols(se)$amplicon[top10]))

  summ[[length(summ) + 1L]] <- data.frame(
    line = ln, n_constituents = sum(keep), n_stitched = length(st),
    n_SE = length(se), median_SE_width = median(width(se)),
    max_SE_width = max(width(se)),
    n_SE_in_amplicon = n_amp,
    n_top10_in_amplicon = n_amp_top10,
    stringsAsFactors = FALSE)

  # A cutoff that admits most of the curve has failed, whatever it returned.
  # Published ROSE analyses call roughly 3-15% of stitched enhancers as SEs.
  frac <- length(se) / length(st)
  flagmsg <- if (frac > 0.30) "  <-- CUTOFF FAILED (>30% called)" else ""

  cat(sprintf("  %-9s constituents %6d -> stitched %6d -> SE %5d (%4.1f%%)  (amplicon-resident %d; top10 %d)%s\n",
              ln, sum(keep), length(st), length(se), 100 * frac, n_amp, n_amp_top10, flagmsg))
}

s <- do.call(rbind, summ)
s$pct_of_stitched <- round(100 * s$n_SE / s$n_stitched, 1)
write.csv(s, "data/metadata/super_enhancer_summary.csv", row.names = FALSE)

cat("\n=========== super-enhancers per line ===========\n")
print(s, row.names = FALSE)

bad <- s[s$pct_of_stitched > 30, ]
if (nrow(bad)) {
  cat("\nCUTOFF FAILURE in ", nrow(bad), " line(s): ",
      paste(bad$line, collapse = ", "), "\n", sep = "")
  cat("More than 30% of stitched enhancers called as super-enhancers is not a\n")
  cat("result, it is a broken inflection point. Do not use these SE sets.\n")
  quit(status = 1)
}

# --- cross-line reproducibility ----------------------------------------------
cat("\n=========== cross-line reproducibility ===========\n")
cat("An SE recurring in lines with DIFFERENT amplicons is unlikely to be a\n")
cat("copy-number artefact. This is the main defence available.\n\n")
allse <- GenomicRanges::reduce(do.call(c, unname(se_all)))
supp  <- rowSums(vapply(se_all, function(g) IRanges::overlapsAny(allse, g),
                        logical(length(allse))))
mcols(allse)$n_lines <- supp
tb <- table(supp)
for (k in names(tb))
  cat(sprintf("  in %2s line(s): %6s SE loci (%5.1f%%)\n", k,
              format(as.integer(tb[[k]]), big.mark = ","), 100 * tb[[k]] / length(allse)))
cat(sprintf("\nunion of SE loci: %s; recurrent (>=3 lines): %s\n",
            format(length(allse), big.mark = ","), format(sum(supp >= 3), big.mark = ",")))

# --- annotate recurrent SEs with nearby genes --------------------------------
rec <- allse[supp >= 3]
if (length(rec)) {
  nr <- GenomicRanges::distanceToNearest(rec, tss)
  sym <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = names(tss)[subjectHits(nr)], keytype = "ENTREZID",
    column = "SYMBOL", multiVals = "first"))
  mcols(rec)$nearest_gene <- unname(sym)
  mcols(rec)$distance <- mcols(nr)$distance
  top <- rec[head(order(-mcols(rec)$n_lines, mcols(rec)$distance), 20)]
  cat("\n=========== most recurrent SE loci (top 20) ===========\n")
  print(data.frame(
    locus = paste0(seqnames(top), ":", format(start(top), big.mark = ","), "-",
                   format(end(top), big.mark = ",")),
    width_kb = round(width(top) / 1000, 1),
    n_lines = mcols(top)$n_lines,
    nearest_gene = mcols(top)$nearest_gene,
    dist = mcols(top)$distance, stringsAsFactors = FALSE), row.names = FALSE)
}

saveRDS(list(per_line = se_all, union = allse, summary = s),
        "data/processed/regions/super_enhancers.rds")

md <- c("# M5 — super-enhancers", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "ROSE-style: H3K27ac-positive constituents, TSS-proximal excluded, stitched",
        paste0("within ", SEP$stitch_distance, " bp, ranked by total signal, cut at the"),
        "inflection point. Genuine white space — the prior work performed a 1 kb",
        "promoter/distal split and no SE analysis.", "",
        "**Copy-number caveat.** Confirmed amplicons carry 38-71x H3K27ac coverage",
        "(D-026), so amplicon-resident stitched enhancers rank highly for reasons of",
        "DNA dosage rather than regulatory investment. They are flagged, retained, and",
        "never reported as regulatory findings without the flag. No correction is",
        "applied because D-026 established our copy-number proxy is unreliable, and an",
        "unvalidated correction would be worse than a declared limitation.", "",
        "| line | constituents | stitched | SE | amplicon-resident | of top 10 |",
        "|---|---|---|---|---|---|")
for (i in seq_len(nrow(s)))
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |", s$line[i],
                      format(s$n_constituents[i], big.mark = ","),
                      format(s$n_stitched[i], big.mark = ","),
                      format(s$n_SE[i], big.mark = ","),
                      s$n_SE_in_amplicon[i], s$n_top10_in_amplicon[i]))
writeLines(md, "results/tables/m5_super_enhancers.md")

cat("\nwrote data/processed/regions/super_enhancers.rds\n")
cat("wrote results/tables/m5_super_enhancers.md\n")

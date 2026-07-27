# -----------------------------------------------------------------------------
# 17_infer_amplification.R — determine MYC-family amplification status from the data.
#
#   Rscript scripts/02_regulatory/17_infer_amplification.R
#
# WHY THIS IS NEEDED
#
# Criterion 3 of the M5 gate needs the MYC-amplified vs MYC-expressing split.
# Plotnik states "we profiled two cell lines harboring the alteration for each
# amplification type", but GSE230649 contains FIVE MYC ChIP samples (H1048, H211,
# H524, H847, SHP77). So two of those five are amplified and three are the
# MYC-expressing comparator. The paper never names them — it says only
# "representative SCLC cell lines" — and the GEO records carry no amplification
# field. The project registry lists all five as MYC-amplified, which cannot be
# right.
#
# Rather than guess, measure it. Amplification means extra DNA copies, which
# inflates read pileup across the whole amplicon regardless of assay. So the ratio
#
#     median signal over universe regions within the amplicon window
#     ------------------------------------------------------------
#     median signal over all universe regions, same track
#
# is a copy-number proxy. It needs no new data: the per-region signal matrix
# already covers every track.
#
# RAW signal is used deliberately. The fold-over-background normalisation divides
# by a per-track constant, which is fine here (ratios are unaffected), but using
# raw keeps the quantity interpretable as relative coverage.
#
# CAVEAT, stated up front: pileup depth reflects copy number AND accessibility.
# A locus can be highly accessible without amplification. So this is corroborating
# evidence, strongest when ATAC and H3K27ac agree and when the elevation extends
# across a broad window rather than sitting on a single peak.
#
# Output: data/metadata/inferred_amplification.csv
#         results/tables/m5_amplification_inference.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(GenomicRanges) })

sig  <- readRDS("data/processed/signal/region_signal_raw.rds")
M    <- sig$mean            # raw per-bp mean signal, regions x samples
meta <- sig$meta
u    <- sig$regions

# hg19 coordinates of the MYC-family genes
LOCI <- list(
  MYC   = GRanges("chr8", IRanges(128748315, 128753680)),
  MYCN  = GRanges("chr2", IRanges(16080683,  16087129)),
  MYCL1 = GRanges("chr1", IRanges(40361098,  40367687)))

WINDOWS <- c(tight = 200000, broad = 1000000)

cat("=========== amplicon coverage ratio by line ===========\n")
cat("ratio = median signal in window / median signal genome-wide, same track\n")
cat("a focal amplification should give a ratio well above 1 in BOTH assays\n\n")

rows <- list()
for (wname in names(WINDOWS)) {
  w <- WINDOWS[[wname]]
  idx <- lapply(LOCI, function(g) {
    which(IRanges::overlapsAny(u, GenomicRanges::resize(g, width(g) + 2 * w, fix = "center")))
  })
  for (i in seq_len(ncol(M))) {
    gw <- stats::median(M[, i], na.rm = TRUE)
    if (!is.finite(gw) || gw <= 0) next
    for (locus in names(LOCI)) {
      k <- idx[[locus]]
      if (length(k) < 5) next
      rows[[length(rows) + 1L]] <- data.frame(
        window = wname, locus = locus,
        assay = meta$assay[i], line = meta$line[i],
        n_regions = length(k),
        ratio = round(stats::median(M[k, i], na.rm = TRUE) / gw, 2),
        stringsAsFactors = FALSE)
    }
  }
}
d <- do.call(rbind, rows)
write.csv(d, "data/metadata/inferred_amplification.csv", row.names = FALSE)

# Present ATAC and H3K27ac side by side — the two assays available in every line.
for (wname in names(WINDOWS)) {
  cat("\n---------- window: +/-", format(WINDOWS[[wname]], big.mark = ","), " bp ----------\n", sep = "")
  s <- d[d$window == wname & d$assay %in% c("ATAC", "H3K27ac"), ]
  if (!nrow(s)) next
  tab <- reshape(s[, c("line", "assay", "locus", "ratio")],
                 idvar = c("line", "locus"), timevar = "assay", direction = "wide")
  names(tab) <- sub("^ratio\\.", "", names(tab))
  tab <- tab[order(tab$locus, -rowMeans(tab[, intersect(c("ATAC","H3K27ac"), names(tab)), drop = FALSE],
                                        na.rm = TRUE)), ]
  for (locus in names(LOCI)) {
    ss <- tab[tab$locus == locus, ]
    if (!nrow(ss)) next
    cat("\n", locus, " locus:\n", sep = "")
    print(ss[, c("line", intersect(c("ATAC","H3K27ac"), names(ss)))], row.names = FALSE)
  }
}

# --- verdict per line ---------------------------------------------------------
cat("\n\n=========== inferred status ===========\n")
broad <- d[d$window == "broad" & d$assay %in% c("ATAC", "H3K27ac"), ]
agg <- aggregate(ratio ~ line + locus, data = broad, FUN = mean)
lines_all <- sort(unique(agg$line))
out <- data.frame()
for (l in lines_all) {
  s <- agg[agg$line == l, ]
  best <- s[which.max(s$ratio), ]
  out <- rbind(out, data.frame(
    line = l,
    MYC   = round(s$ratio[s$locus == "MYC"],   2),
    MYCN  = round(s$ratio[s$locus == "MYCN"],  2),
    MYCL1 = round(s$ratio[s$locus == "MYCL1"], 2),
    strongest = best$locus, max_ratio = round(best$ratio, 2),
    stringsAsFactors = FALSE))
}
print(out, row.names = FALSE)

cat("\nRegistry currently claims:\n")
cat("  MYC-amp   H1048 H211 H524 H847 SHP77   <- all five, which conflicts with\n")
cat("                                            Plotnik's 'two cell lines per type'\n")
cat("  MYCN-amp  H526 H69\n")
cat("  MYCL-amp  COLO668 H889\n")
cat("  non-amp   H196\n")
cat("\nCompare against the measured ratios above. A line whose MYC ratio is near 1\n")
cat("is NOT MYC-amplified and belongs in the MYC-expressing comparator group.\n")
cat("\nCAVEAT: pileup reflects copy number AND accessibility. Treat a call as solid\n")
cat("only where ATAC and H3K27ac agree and the elevation holds in the broad window.\n")

md <- c("# M5 — MYC-family amplification inferred from coverage", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "Plotnik profiled \"two cell lines harboring the alteration for each",
        "amplification type\" but GSE230649 has five MYC ChIP samples, and neither",
        "the paper nor the GEO records name which are amplified. Measured here",
        "instead, as amplicon coverage relative to each track's genome-wide median.", "",
        "| line | MYC | MYCN | MYCL1 | strongest | ratio |", "|---|---|---|---|---|---|")
for (i in seq_len(nrow(out)))
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |", out$line[i], out$MYC[i],
                      out$MYCN[i], out$MYCL1[i], out$strongest[i], out$max_ratio[i]))
md <- c(md, "",
        "Coverage reflects copy number and accessibility together, so a call is",
        "treated as solid only where ATAC and H3K27ac agree and the elevation holds",
        "across the broad (+/-1 Mb) window.")
writeLines(md, "results/tables/m5_amplification_inference.md")
cat("\nwrote data/metadata/inferred_amplification.csv\n")
cat("wrote results/tables/m5_amplification_inference.md\n")

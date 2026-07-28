# The previous comparison used identical() on data frames, which is sensitive to
# ROW ORDER and to Hits metadata. Counts were equal, so test SET equality of the
# actual (region, tss) pairs — that is what the pipeline consumes.
suppressPackageStartupMessages({
  library(GenomicRanges); library(TxDb.Hsapiens.UCSC.hg19.knownGene)
})
nrm <- readRDS("data/processed/signal/region_signal_normalised.rds")
u   <- nrm$regions
tx   <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tssg <- promoters(tx, upstream = 0, downstream = 1)
tssg <- tssg[as.character(seqnames(tssg)) %in% paste0("chr", c(1:22, "X"))]

maxd <- 200000
w_raw <- GenomicRanges::resize(tssg, 2 * maxd, fix = "center")
w_trm <- suppressWarnings(GenomicRanges::trim(w_raw))

cat("windows trim() actually changed:",
    sum(start(w_raw) != start(w_trm) | end(w_raw) != end(w_trm)), "\n\n")

o1 <- suppressWarnings(GenomicRanges::findOverlaps(u, w_raw))
o2 <- suppressWarnings(GenomicRanges::findOverlaps(u, w_trm))

k1 <- paste(queryHits(o1), subjectHits(o1))
k2 <- paste(queryHits(o2), subjectHits(o2))

cat("pairs untrimmed :", length(k1), "\n")
cat("pairs trimmed   :", length(k2), "\n")
cat("only untrimmed  :", length(setdiff(k1, k2)), "\n")
cat("only trimmed    :", length(setdiff(k2, k1)), "\n")
cat("set-equal       :", setequal(k1, k2), "\n")
cat("order-identical :", identical(k1, k2), "\n")

cat("\nVERDICT: ",
    if (setequal(k1, k2))
      "pair SETS are identical — the missing trim() does not change any linkage.\n         The 54 warnings are hygiene only; no downstream rerun is required."
    else
      "pair sets DIFFER — downstream results are affected and must be rerun.",
    "\n", sep = "")

# -----------------------------------------------------------------------------
# test_resize_trim.R — is the missing trim() after resize() cosmetic or a bug?
#
#   Rscript scripts/00_setup/test_resize_trim.R
#
# GenomicRanges warns that widening TSS windows produces out-of-bound ranges. The
# claim this project makes is that the warning is HYGIENE: the out-of-bound
# territory contains no universe regions, so intersecting with it adds nothing and
# no linkage changes.
#
# That claim decides whether a pipeline rerun is required, so it is TESTED here
# rather than argued. The test compares overlap pair SETS — not data frames, whose
# comparison is sensitive to row order and to Hits metadata and which reported a
# spurious difference on the first attempt.
#
# Exits non-zero if the sets ever diverge, which would mean downstream results are
# affected and must be regenerated.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(yaml)
})

nrm <- readRDS("data/processed/signal/region_signal_normalised.rds")
u   <- nrm$regions
P2G <- yaml::read_yaml("config/params.yml")$peak_to_gene
maxd <- if (!is.null(P2G$max_distance)) P2G$max_distance else 200000

tx   <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tssg <- promoters(tx, upstream = 0, downstream = 1)
tssg <- tssg[as.character(seqnames(tssg)) %in% paste0("chr", c(1:22, "X"))]

check <- function(label, w) {
  raw <- suppressWarnings(GenomicRanges::resize(tssg, w, fix = "center"))
  trm <- suppressWarnings(GenomicRanges::trim(raw))
  n_changed <- sum(start(raw) != start(trm) | end(raw) != end(trm))
  o1 <- suppressWarnings(GenomicRanges::findOverlaps(u, raw))
  o2 <- suppressWarnings(GenomicRanges::findOverlaps(u, trm))
  k1 <- paste(queryHits(o1), subjectHits(o1))
  k2 <- paste(queryHits(o2), subjectHits(o2))
  ok <- setequal(k1, k2)
  cat(sprintf("  %-28s width %8s | out-of-bound %4d | pairs %s vs %s | %s\n",
              label, format(w, big.mark = ","), n_changed,
              format(length(k1), big.mark = ","), format(length(k2), big.mark = ","),
              if (ok) "SET-IDENTICAL" else "DIFFERENT"))
  ok
}

cat("=========== resize() with and without trim() ===========\n")
results <- c(
  check("promoter window", 2000L),
  check("distal prefilter window", as.integer(2 * maxd))
)

cat("\n")
if (all(results)) {
  cat("PASS — trimming changes no overlap pair. The out-of-bound warnings are\n")
  cat("hygiene, and no downstream regeneration is required on their account.\n")
} else {
  cat("FAIL — trimming CHANGES overlap pairs. Downstream results depend on this\n")
  cat("and must be regenerated before anything else is trusted.\n")
  quit(status = 1)
}

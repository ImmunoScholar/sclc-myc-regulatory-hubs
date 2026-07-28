# Does the missing trim() after resize() change any result, or is it cosmetic?
#
# Reasoning to be TESTED, not trusted: resize() can push a window past a
# chromosome end or below coordinate 1. Overlap is interval intersection, and the
# universe contains no regions out there, so the out-of-bound territory should be
# empty and the overlap unchanged. If that reasoning is right the warnings are
# hygiene; if it is wrong, every downstream result is affected and must be rerun.
suppressPackageStartupMessages({
  library(GenomicRanges); library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db)
})

nrm <- readRDS("data/processed/signal/region_signal_normalised.rds")
u   <- nrm$regions

tx   <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tssg <- promoters(tx, upstream = 0, downstream = 1)
tssg <- tssg[as.character(seqnames(tssg)) %in% paste0("chr", c(1:22, "X"))]

raw     <- GenomicRanges::resize(tssg, 2000, fix = "center")
trimmed <- GenomicRanges::trim(raw)

oob <- sum(start(raw) < 1 | end(raw) > seqlengths(raw)[as.character(seqnames(raw))], na.rm = TRUE)
cat("out-of-bound windows after resize():", oob, "of", length(raw), "\n")
cat("windows whose coordinates trim() changed:", sum(start(raw) != start(trimmed) | end(raw) != end(trimmed)), "\n\n")

link_for <- function(gr) {
  ov <- GenomicRanges::findOverlaps(u, gr)
  sym <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = names(tssg)[subjectHits(ov)], keytype = "ENTREZID",
    column = "SYMBOL", multiVals = "first"))
  d <- data.frame(region = queryHits(ov), gene = unname(sym), stringsAsFactors = FALSE)
  d[!is.na(d$gene), ]
}

a <- link_for(raw); b <- link_for(trimmed)
cat("promoter links, untrimmed:", nrow(a), "\n")
cat("promoter links, trimmed  :", nrow(b), "\n")

ka <- paste(a$region, a$gene); kb <- paste(b$region, b$gene)
cat("links only in untrimmed  :", length(setdiff(ka, kb)), "\n")
cat("links only in trimmed    :", length(setdiff(kb, ka)), "\n")

identical_links <- setequal(ka, kb) && nrow(a) == nrow(b)
cat("\nVERDICT: promoter link sets are",
    if (identical_links) "IDENTICAL — the missing trim() is cosmetic"
    else "DIFFERENT — downstream results are affected and must be rerun", "\n")

# Second check: does it alter the distal window used for peak-to-gene linking?
P2G <- yaml::read_yaml("config/params.yml")$peak_to_gene
maxd <- if (!is.null(P2G$max_distance)) P2G$max_distance else 250000
wide_raw <- GenomicRanges::resize(tssg, 2 * maxd, fix = "center")
wide_trm <- GenomicRanges::trim(wide_raw)
ov1 <- GenomicRanges::findOverlaps(u, wide_raw)
ov2 <- GenomicRanges::findOverlaps(u, wide_trm)
cat("\ndistal window (", 2 * maxd, " bp) overlaps: untrimmed ", length(ov1),
    " vs trimmed ", length(ov2), " -> ",
    if (length(ov1) == length(ov2) && identical(as.data.frame(ov1), as.data.frame(ov2)))
      "IDENTICAL" else "DIFFERENT", "\n", sep = "")

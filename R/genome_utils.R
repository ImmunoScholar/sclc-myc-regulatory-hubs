# -----------------------------------------------------------------------------
# genome_utils.R — sequence-name harmonisation and genome-build assertions.
#
#   source("R/genome_utils.R")
#
# WHY THIS EXISTS
#
# The keystone deposit (GSE230649) uses ENSEMBL sequence names — bare "1", "20",
# "X" — while the project's annotation resources (TxDb.Hsapiens.UCSC.hg19.knownGene,
# hg19.chrom.sizes, UCSC liftOver chains) all use UCSC names, "chr1", "chr20",
# "chrX". Both describe identical hg19 coordinates.
#
# This mismatch does not raise an error. GenomicRanges operations between objects
# with disjoint seqlevels return ZERO overlaps, which looks like a biological
# result: "MYC binds none of these regions." Every coordinate join in this project
# must therefore pass through harmonise_seqnames() first.
#
# It also caused a real defect at M4: the build assertion in 04_qc_report.R
# merged bedGraph seqnames against hg19.chrom.sizes, matched nothing, found no
# violations, and reported a PASS. See risk R-16. assert_hg19_bounds() below
# refuses to return a pass when it had nothing to compare.
# -----------------------------------------------------------------------------

# Chromosomes the analysis uses. chrY excluded: these lines are of mixed sex and
# Y is unreliable in cancer lines. chrM excluded: no regulatory interpretation.
ANALYSIS_CHROMS_UCSC <- paste0("chr", c(1:22, "X"))

#' Normalise any sequence name to UCSC style.
#' Idempotent: "chr1" -> "chr1", "1" -> "chr1", "MT" -> "chrM".
to_ucsc_seqnames <- function(x) {
  x <- as.character(x)
  bare <- sub("^chr", "", x, ignore.case = TRUE)
  bare[bare %in% c("MT", "M", "mt")] <- "M"
  paste0("chr", bare)
}

#' Normalise any sequence name to Ensembl style. "chr1" -> "1", "chrM" -> "MT".
to_ensembl_seqnames <- function(x) {
  x <- as.character(x)
  bare <- sub("^chr", "", x, ignore.case = TRUE)
  bare[bare %in% c("M", "MT", "mt")] <- "MT"
  bare
}

#' Identify the naming convention of a set of sequence names.
#' Returns "UCSC", "Ensembl", "mixed", or "unknown".
detect_seqname_style <- function(chroms) {
  chroms <- unique(as.character(chroms))
  chroms <- chroms[!is.na(chroms) & nzchar(chroms)]
  if (!length(chroms)) return("unknown")

  ens_core  <- c(as.character(1:22), "X", "Y", "MT")
  n_ucsc    <- sum(grepl("^chr", chroms, ignore.case = TRUE))
  n_ens     <- sum(chroms %in% ens_core)

  if (n_ucsc == length(chroms))                return("UCSC")
  if (n_ens  == length(chroms))                return("Ensembl")
  if (n_ucsc > 0 && n_ens > 0)                 return("mixed")
  if (n_ucsc > 0 || n_ens > 0)                 return("mixed")
  "unknown"
}

#' Harmonise a GRanges (or seqname vector) to UCSC style.
#' Use before ANY cross-object coordinate operation.
harmonise_seqnames <- function(x, style = c("UCSC", "Ensembl")) {
  style <- match.arg(style)
  conv  <- if (style == "UCSC") to_ucsc_seqnames else to_ensembl_seqnames
  if (inherits(x, "GRanges")) {
    GenomeInfoDb::seqlevels(x) <- conv(GenomeInfoDb::seqlevels(x))
    return(x)
  }
  conv(x)
}

#' Load hg19 chromosome sizes as a named integer vector, UCSC names.
load_hg19_sizes <- function(
    path = "data/raw/ucsc_hg19_chrom_sizes/hg19.chrom.sizes") {
  if (!file.exists(path)) stop("hg19.chrom.sizes not found at ", path)
  d <- utils::read.delim(path, header = FALSE,
                         col.names = c("chrom", "size"),
                         stringsAsFactors = FALSE)
  stats::setNames(as.numeric(d$size), to_ucsc_seqnames(d$chrom))
}

#' Assert that per-chromosome maximum coordinates fall within hg19 bounds.
#'
#' @param max_by_chrom data.frame(chrom, max_end) — names in EITHER convention.
#' @param sizes named vector from load_hg19_sizes().
#' @param min_chroms refuse to pass unless at least this many chromosomes were
#'   actually compared. This is the guard against the vacuous pass (R-16).
#' @return list(pass, n_compared, n_violations, detail)
assert_hg19_bounds <- function(max_by_chrom, sizes, min_chroms = 20L) {
  if (is.null(max_by_chrom) || !nrow(max_by_chrom)) {
    return(list(pass = FALSE, n_compared = 0L, n_violations = NA_integer_,
                detail = "no coordinates supplied — nothing to assert"))
  }

  key <- to_ucsc_seqnames(max_by_chrom$chrom)
  # Compare only against real hg19 sequences; scaffolds are out of scope.
  keep <- key %in% names(sizes)
  n_cmp <- sum(keep)

  if (n_cmp == 0L) {
    return(list(pass = FALSE, n_compared = 0L, n_violations = NA_integer_,
                detail = paste0(
                  "NO sequence name matched hg19 after harmonisation (saw: ",
                  paste(utils::head(unique(max_by_chrom$chrom), 4), collapse = ","),
                  "). This is the vacuous-pass failure mode, not a pass.")))
  }

  lim  <- sizes[key[keep]]
  obs  <- as.numeric(max_by_chrom$max_end[keep])
  over <- which(obs > lim)

  if (n_cmp < min_chroms) {
    return(list(pass = FALSE, n_compared = n_cmp, n_violations = length(over),
                detail = paste0("only ", n_cmp, " chromosome(s) compared; ",
                                "expected at least ", min_chroms,
                                " — sample is too partial to assert a build")))
  }

  if (length(over)) {
    i <- over[1]
    return(list(pass = FALSE, n_compared = n_cmp, n_violations = length(over),
                detail = paste0(length(over), " chromosome(s) exceed hg19 length, e.g. ",
                                key[keep][i], " max=", format(obs[i], scientific = FALSE),
                                " > ", format(lim[i], scientific = FALSE),
                                " — BUILD MISMATCH")))
  }

  list(pass = TRUE, n_compared = n_cmp, n_violations = 0L,
       detail = paste0("all ", n_cmp, " chromosomes within hg19 bounds"))
}

#' Stream an interval file and return the maximum end coordinate per chromosome.
#'
#' Streams the WHOLE file with awk rather than sampling. Necessary, not merely
#' thorough: these bedGraphs are written in BAM-header order
#' (20, 21, 22, 1, 3, 2, ...), so a head sample sees only one chromosome and any
#' build assertion built on it is meaningless.
#'
#' @param col_end which column holds the interval end (3 for bedGraph/BED).
stream_max_by_chrom <- function(path, col_end = 3L) {
  reader <- if (grepl("\\.gz$", path)) "zcat" else "cat"
  awk_prog <- sprintf(
    "awk -F'\\t' '$1 !~ /^(track|#|browser)/ && NF >= %d { if ($%d+0 > m[$1]) m[$1]=$%d+0 } END { for (c in m) print c\"\\t\"m[c] }'",
    col_end, col_end, col_end)
  con <- pipe(paste(reader, shQuote(path), "|", awk_prog), "r")
  on.exit(close(con))
  out <- try(utils::read.delim(con, header = FALSE,
                               col.names = c("chrom", "max_end"),
                               stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

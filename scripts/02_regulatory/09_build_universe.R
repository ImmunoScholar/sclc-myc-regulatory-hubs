# -----------------------------------------------------------------------------
# 09_build_universe.R — the final accessible-region universe, genome-wide.
#
#   Rscript scripts/02_regulatory/09_build_universe.R
#
# Supersedes 05_. Two changes, both from D-024:
#   * threshold x2.5 (was x3, calibrated against the wrong target)
#   * universes built at x2, x2.5, x3 and x4 so every gate metric can be reported
#     across thresholds. The primary is x2.5; the others exist so the conclusion
#     can be shown independent of the choice rather than resting on it.
#
# Each line is thresholded at a multiple of ITS OWN mean signal over covered bases,
# so differing depth does not let deep tracks dominate the universe — which would
# matter here because paralog groups are confounded with cell lines.
#
# Output: data/processed/regions/universe_x<mult>.rds   (4 files)
#         data/processed/regions/universe.rds           (primary, x2.5)
#         data/metadata/universe_multi_threshold.csv
#         results/tables/m5_universe_final.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(yaml)
})
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
RP  <- CFG$regions
OUT <- "data/processed/regions"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MULTIPLES <- c(2, 2.5, 3, 4)
PRIMARY   <- 2.5
ENS_CHROMS <- c(as.character(1:22), "X")
TSS_WIN <- 1000L

man <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)
man$dest <- file.path("data/raw", man$dataset_id, man$file_name)
ks   <- man[man$dataset_id == "GSE230649" & file.exists(man$dest), ]
atac <- ks[!grepl("^GSM[0-9]+_(H3K27Ac|MYC|MYCN|MYCL1)_", ks$file_name), ]
atac$line <- sub("^GSM[0-9]+_([^_]+)_treat_pileup.*$", "\\1", atac$file_name)
stopifnot(nrow(atac) == 9)

line_mean <- function(line) {
  raw <- readLines(file.path("data/processed/atac_profile", paste0(line, ".hist.tsv")))
  h <- utils::read.delim(text = paste(raw[!grepl("^#", raw)], collapse = "\n"),
                         header = FALSE, col.names = c("bin","bp"))
  h$signal <- h$bin / 2; cov <- h[h$signal > 0, ]
  sum(cov$signal * cov$bp) / sum(cov$bp)
}

# One pass per file emits merged runs above the LOWEST multiple; stricter sets are
# derived from the same scan by re-thresholding on the run's mean signal. Nine
# passes instead of thirty-six.
scan_line <- function(path, cutoff) {
  filt <- paste0("($1==\"", paste(ENS_CHROMS, collapse = "\"||$1==\""), "\")")
  awk <- sprintf(paste0(
    "awk -F'\\t' '%s {",
    " w=$3-$2; v=$4+0;",
    " if (v >= %f) {",
    "   if (open && $1==cc && $2==ce) { ce=$3; sw+=v*w }",
    "   else { if (open) printf \"%%s\\t%%d\\t%%d\\t%%.4f\\n\", cc, cs, ce, sw;",
    "          cc=$1; cs=$2; ce=$3; sw=v*w; open=1 }",
    " } else if (open) { printf \"%%s\\t%%d\\t%%d\\t%%.4f\\n\", cc, cs, ce, sw; open=0 }",
    "}",
    "END { if (open) printf \"%%s\\t%%d\\t%%d\\t%%.4f\\n\", cc, cs, ce, sw }'"),
    filt, cutoff)
  con <- pipe(paste("zcat", shQuote(path), "|", awk), "r"); on.exit(close(con))
  d <- try(utils::read.delim(con, header = FALSE,
                             col.names = c("chrom","start","end","sumsw")), silent = TRUE)
  if (inherits(d, "try-error") || !nrow(d)) return(NULL)
  d
}

cache <- file.path(OUT, "per_line_runs.rds")
if (file.exists(cache)) {
  runs <- readRDS(cache); cat("per-line scans loaded from cache\n\n")
} else {
  runs <- list()
  for (i in seq_len(nrow(atac))) {
    ln <- atac$line[i]; gm <- line_mean(ln)
    t0 <- Sys.time()
    d  <- scan_line(atac$dest[i], min(MULTIPLES) * gm)
    d$mean_sig <- d$sumsw / (d$end - d$start)
    d$gm <- gm
    runs[[ln]] <- d
    cat(sprintf("  %-9s mean=%.3f  runs=%8d  (%.1f min)\n", ln, gm, nrow(d),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  saveRDS(runs, cache)
}

to_gr <- function(d, mult) {
  sel <- d[d$mean_sig >= mult * d$gm, ]
  if (!nrow(sel)) return(GRanges())
  gr <- GenomicRanges::reduce(
    GRanges(to_ucsc_seqnames(sel$chrom), IRanges(sel$start + 1L, sel$end)),
    min.gapwidth = RP$merge_gap)
  gr[width(gr) >= RP$min_width & width(gr) <= RP$max_width]
}

# TSS reference for distal fraction
tx  <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tss <- promoters(tx, upstream = 0, downstream = 1)
tss <- tss[as.character(seqnames(tss)) %in% ANALYSIS_CHROMS_UCSC]
tss_win <- GenomicRanges::reduce(resize(tss, width = 2 * TSS_WIN, fix = "center"))

lift <- list.files("data/processed/liftover", pattern = "\\.hg19\\.rds$", full.names = TRUE)
ext_ds <- unique(sub("__.*$", "", basename(lift)))
ext_gr <- lapply(ext_ds, function(ds)
  GenomicRanges::reduce(do.call(c, unname(lapply(lift[startsWith(basename(lift), ds)], readRDS))),
                        min.gapwidth = RP$merge_gap))
names(ext_gr) <- ext_ds

genome_bp <- sum(load_hg19_sizes()[ANALYSIS_CHROMS_UCSC])
summ <- list()

for (m in MULTIPLES) {
  per_line <- lapply(runs, to_gr, mult = m)
  u <- GenomicRanges::reduce(do.call(c, unname(per_line)), min.gapwidth = RP$merge_gap)
  u <- u[width(u) >= RP$min_width & width(u) <= RP$max_width]

  sup <- vapply(per_line, function(g) IRanges::overlapsAny(u, g), logical(length(u)))
  mcols(u)$n_lines_supporting <- rowSums(sup)
  for (ln in names(per_line)) mcols(u)[[paste0("atac_", ln)]] <- sup[, ln]
  for (ds in ext_ds) mcols(u)[[paste0("ext_", ds)]] <- IRanges::overlapsAny(u, ext_gr[[ds]])
  ecols <- grep("^ext_", names(mcols(u)), value = TRUE)
  mcols(u)$n_external_support <- rowSums(as.data.frame(mcols(u)[, ecols, drop = FALSE]))
  mcols(u)$is_promoter <- IRanges::overlapsAny(u, tss_win)

  saveRDS(u, file.path(OUT, sprintf("universe_x%g.rds", m)))
  if (m == PRIMARY) saveRDS(u, file.path(OUT, "universe.rds"))

  summ[[length(summ) + 1L]] <- data.frame(
    multiple = m, primary = (m == PRIMARY),
    n_regions = length(u), median_width = median(width(u)),
    total_mb = round(sum(as.numeric(width(u))) / 1e6, 1),
    pct_genome = round(100 * sum(as.numeric(width(u))) / genome_bp, 3),
    pct_promoter = round(100 * mean(mcols(u)$is_promoter), 1),
    pct_distal = round(100 * (1 - mean(mcols(u)$is_promoter)), 1),
    pct_multiline = round(100 * mean(mcols(u)$n_lines_supporting >= 2), 1),
    pct_ext_corrob = round(100 * mean(mcols(u)$n_external_support > 0), 1),
    stringsAsFactors = FALSE)

  cat(sprintf("\nx%-4g  %8s regions  %5.2f%% genome  distal %4.1f%%  multi-line %4.1f%%  ext %4.1f%%\n",
              m, format(length(u), big.mark = ","),
              100 * sum(as.numeric(width(u))) / genome_bp,
              100 * (1 - mean(mcols(u)$is_promoter)),
              100 * mean(mcols(u)$n_lines_supporting >= 2),
              100 * mean(mcols(u)$n_external_support > 0)))
}

s <- do.call(rbind, summ)
write.csv(s, "data/metadata/universe_multi_threshold.csv", row.names = FALSE)

cat("\n=========== universes across thresholds ===========\n")
print(s, row.names = FALSE)

md <- c("# M5 — accessible-region universe", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        paste0("Derived from the nine keystone ATAC tracks (GSE230649, hg19, cells ",
               "matched to the ChIP). Primary threshold **x", PRIMARY,
               "** each line's mean signal over covered bases (D-024)."), "",
        "Four universes are retained so every M5 gate metric can be reported across",
        "thresholds. The conclusion should not depend on the choice; if it does,",
        "that fragility is the finding.", "",
        "| x mean | regions | % genome | distal % | multi-line % | ext. corrob % |",
        "|---|---|---|---|---|---|")
for (i in seq_len(nrow(s)))
  md <- c(md, sprintf("| %s%s | %s | %s | %s | %s | %s |",
                      s$multiple[i], ifelse(s$primary[i], " **(primary)**", ""),
                      format(s$n_regions[i], big.mark = ","), s$pct_genome[i],
                      s$pct_distal[i], s$pct_multiline[i], s$pct_ext_corrob[i]))
md <- c(md, "",
        "Validated intrinsically: H524/chr1 TSS enrichment 8.9x and matched-cell",
        "H3K27ac fold 2.83x at the primary threshold (D-024). Agreement with an",
        "external laboratory's ATAC peak calls was abandoned as a target — it",
        "measured inter-laboratory reproducibility, not fitness for purpose.")
writeLines(md, "results/tables/m5_universe_final.md")

cat("\nwrote universe.rds (primary, x", PRIMARY, ")\n", sep = "")
cat("wrote results/tables/m5_universe_final.md\n")
cat("\nNext: quantify the 19 ChIP tracks over this grid.\n")

# -----------------------------------------------------------------------------
# 12_finalise_universe.R — select the threshold and write the final universe.
#
#   Rscript scripts/02_regulatory/12_finalise_universe.R
#
# Uses the EXACT per-line scans from 11_ (not the approximations from 09_).
#
# Support rule is fixed at >=2 keystone ATAC lines and is NOT a free choice:
# MYCN-amplified and MYCL-amplified lines number two each, so a region accessible
# only in H526+H69 — exactly a MYCN-specific enhancer — has support level 2.
# Requiring >=3 would structurally delete every MYCN- and MYCL1-specific region and
# with them the MYCN-in-MYC gate criterion. Step 10's metric-optimal suggestion
# (x2.5 / >=3 lines) would have done precisely that.
#
# Threshold is chosen against published targets — TSS enrichment >=8x, promoter
# fraction 15-30%, genome fraction 1-3% — maximising regions subject to all three.
#
# ALSO ADDED HERE, and absent from every earlier version: per-PARALOG-GROUP support
# annotation. The analysis is about paralog specificity, so each region records how
# many MYC-, MYCN- and MYCL-amplified lines support it, and is classified
# accordingly. Without this, "paralog-specific" would have to be inferred later
# from ChIP signal alone, with no accessibility evidence behind it.
#
# Output: data/processed/regions/universe_final.rds
#         data/metadata/universe_threshold_selection.csv
#         results/tables/m5_universe_selected.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(yaml)
})
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
RP  <- CFG$regions
EX  <- "data/processed/regions/exact"
OUT <- "data/processed/regions"
MIN_LINES <- 2L
MULTS <- c(2.5, 3, 4, 5)

# Paralog groups among the nine ATAC lines (author-declared amplification status).
GROUPS <- list(
  MYC   = c("H1048", "H524", "H847", "SHP77"),
  MYCN  = c("H526", "H69"),
  MYCL1 = c("COLO668", "H889"),
  none  = c("H196")
)
ALL_LINES <- unlist(GROUPS, use.names = FALSE)

sizes     <- load_hg19_sizes()
genome_bp <- sum(sizes[ANALYSIS_CHROMS_UCSC])

tx  <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tss <- promoters(tx, upstream = 0, downstream = 1)
tss <- tss[as.character(seqnames(tss)) %in% ANALYSIS_CHROMS_UCSC]
tss_win  <- GenomicRanges::reduce(resize(tss, width = 2000L, fix = "center"))
tss_frac <- sum(as.numeric(width(tss_win))) / genome_bp

# --- refuse to run against scans that are still being written -----------------
# Reading a BED while awk is mid-write yields truncated lines, which surface as
# NA coordinates deep inside GRanges construction rather than as a clear error.
#
# NOT done by process matching: `pgrep -f 11_rescan_exact.sh` also matches the
# calling shell's own command line if that command mentions the script, which
# produced a false positive on a run that had actually completed.
#
# Instead: require every expected output to exist and to have been untouched for
# a quiet period. That cannot self-match and it tests the actual precondition —
# are the files finished — rather than a proxy for it.
QUIET_SECS <- 25
expected <- as.vector(outer(ALL_LINES, MULTS,
                            function(l, m) sprintf("%s_x%g.bed", l, m)))
paths <- file.path(EX, expected)
missing <- expected[!file.exists(paths)]
if (length(missing)) {
  stop(sprintf("%d expected scan file(s) missing, e.g. %s\n  Run: bash scripts/02_regulatory/11_rescan_exact.sh",
               length(missing), missing[1]), call. = FALSE)
}
age <- as.numeric(difftime(Sys.time(), file.info(paths)$mtime, units = "secs"))
if (any(age < QUIET_SECS, na.rm = TRUE)) {
  i <- which.min(age)
  stop(sprintf("scan files are still being written (%s modified %.0fs ago).\n  Wait for 11_rescan_exact.sh to finish, then re-run.",
               expected[i], age[i]), call. = FALSE)
}
cat("scan files: ", length(paths), " present, all quiet for >", QUIET_SECS, "s\n\n", sep = "")

read_bed <- function(path) {
  if (!file.exists(path) || file.size(path) == 0) return(GRanges())
  d <- utils::read.delim(path, header = FALSE, col.names = c("chrom","start","end"),
                         stringsAsFactors = FALSE)
  # A truncated final line, or any malformed row, must fail loudly and name the
  # file — not propagate as an opaque NA error from IRanges.
  bad <- is.na(d$start) | is.na(d$end) | is.na(d$chrom) | !nzchar(d$chrom)
  if (any(bad)) {
    stop(sprintf("%s: %d malformed row(s) (first at line %d). The scan is incomplete or corrupt.\n  Delete it and re-run 11_rescan_exact.sh.",
                 basename(path), sum(bad), which(bad)[1]), call. = FALSE)
  }
  if (!nrow(d)) return(GRanges())
  gr <- GenomicRanges::reduce(
    GRanges(to_ucsc_seqnames(d$chrom), IRanges(d$start + 1L, d$end)),
    min.gapwidth = RP$merge_gap)
  gr[width(gr) >= RP$min_width & width(gr) <= RP$max_width]
}

lift <- list.files("data/processed/liftover", pattern = "\\.hg19\\.rds$", full.names = TRUE)
ext_ds <- unique(sub("__.*$", "", basename(lift)))
ext_gr <- lapply(ext_ds, function(ds)
  GenomicRanges::reduce(do.call(c, unname(lapply(
    lift[startsWith(basename(lift), ds)], readRDS))), min.gapwidth = RP$merge_gap))
names(ext_gr) <- ext_ds

annotate <- function(u, per_line) {
  sup <- vapply(per_line, function(g) IRanges::overlapsAny(u, g), logical(length(u)))
  if (is.null(dim(sup))) sup <- matrix(sup, nrow = length(u),
                                       dimnames = list(NULL, names(per_line)))
  mcols(u)$n_lines_supporting <- rowSums(sup)
  for (ln in colnames(sup)) mcols(u)[[paste0("atac_", ln)]] <- sup[, ln]

  # per-paralog-group support
  for (g in names(GROUPS)) {
    cols <- intersect(GROUPS[[g]], colnames(sup))
    mcols(u)[[paste0("n_", g, "_lines")]] <-
      if (length(cols)) rowSums(sup[, cols, drop = FALSE]) else 0L
  }
  for (ds in ext_ds) mcols(u)[[paste0("ext_", ds)]] <- IRanges::overlapsAny(u, ext_gr[[ds]])
  ecols <- grep("^ext_", names(mcols(u)), value = TRUE)
  mcols(u)$n_external_support <- rowSums(as.data.frame(mcols(u)[, ecols, drop = FALSE]))
  mcols(u)$is_promoter <- IRanges::overlapsAny(u, tss_win)

  # Accessibility class. Descriptive only — it says where a region is OPEN, not
  # where a paralog BINDS. Binding evidence comes from the ChIP quantification.
  nm  <- mcols(u)$n_MYC_lines;   nn <- mcols(u)$n_MYCN_lines
  nl  <- mcols(u)$n_MYCL1_lines
  grp <- (nm > 0) + (nn > 0) + (nl > 0)
  cls <- rep("mixed", length(u))
  cls[grp >= 3]            <- "shared_all_paralogs"
  cls[grp == 1 & nm > 0]   <- "MYC_lines_only"
  cls[grp == 1 & nn > 0]   <- "MYCN_lines_only"
  cls[grp == 1 & nl > 0]   <- "MYCL1_lines_only"
  cls[grp == 0]            <- "non_amplified_only"
  mcols(u)$accessibility_class <- cls
  u
}

sel <- list()
universes <- list()

for (m in MULTS) {
  per_line <- lapply(ALL_LINES, function(ln) read_bed(file.path(EX, sprintf("%s_x%g.bed", ln, m))))
  names(per_line) <- ALL_LINES
  if (all(lengths(per_line) == 0)) { cat("x", m, ": no scans found — run 11_ first\n", sep = ""); next }

  u <- GenomicRanges::reduce(do.call(c, unname(per_line)), min.gapwidth = RP$merge_gap)
  u <- u[width(u) >= RP$min_width & width(u) <= RP$max_width]
  u <- annotate(u, per_line)
  u <- u[mcols(u)$n_lines_supporting >= MIN_LINES]

  prom <- mean(mcols(u)$is_promoter)
  bp   <- sum(as.numeric(width(u)))
  universes[[as.character(m)]] <- u

  sel[[length(sel) + 1L]] <- data.frame(
    multiple = m, n_regions = length(u), median_width = median(width(u)),
    pct_genome = round(100 * bp / genome_bp, 3),
    pct_promoter = round(100 * prom, 1),
    tss_enrichment = round(prom / tss_frac, 2),
    pct_ext_corrob = round(100 * mean(mcols(u)$n_external_support > 0), 1),
    n_MYCN_only = sum(mcols(u)$accessibility_class == "MYCN_lines_only"),
    n_MYCL1_only = sum(mcols(u)$accessibility_class == "MYCL1_lines_only"),
    n_MYC_only = sum(mcols(u)$accessibility_class == "MYC_lines_only"),
    n_shared = sum(mcols(u)$accessibility_class == "shared_all_paralogs"),
    stringsAsFactors = FALSE)

  cat(sprintf("x%-4g %8s regions  %5.2f%% genome  prom %4.1f%%  TSS %5.2fx  ext %4.1f%%\n",
              m, format(length(u), big.mark = ","), 100 * bp / genome_bp,
              100 * prom, prom / tss_frac, 100 * mean(mcols(u)$n_external_support > 0)))
}

s <- do.call(rbind, sel)
write.csv(s, "data/metadata/universe_threshold_selection.csv", row.names = FALSE)

cat("\n=========== exact scans, >=", MIN_LINES, " lines ===========\n", sep = "")
print(s, row.names = FALSE)

cat("\n=========== paralog-specific accessibility ===========\n")
print(s[, c("multiple","n_MYC_only","n_MYCN_only","n_MYCL1_only","n_shared")], row.names = FALSE)
cat("\nThese counts matter: if MYCN_only or MYCL1_only regions are scarce, the\n")
cat("paralog-specific arms have little accessibility evidence to build on.\n")

ok <- s[s$tss_enrichment >= 8 & s$pct_promoter >= 15 & s$pct_promoter <= 30 &
        s$pct_genome >= 1 & s$pct_genome <= 3, ]
cat("\n=========== selection ===========\n")
if (nrow(ok)) {
  pick <- ok[which.max(ok$n_regions), ]
  cat("meets all targets (TSS>=8x, promoter 15-30%, genome 1-3%):\n")
  print(ok[, c("multiple","n_regions","pct_genome","pct_promoter","tss_enrichment")],
        row.names = FALSE)
  cat("\nSELECTED: x", pick$multiple, " -> ", format(pick$n_regions, big.mark = ","),
      " regions, ", pick$pct_genome, "% genome, ", pick$pct_promoter,
      "% promoter, TSS ", pick$tss_enrichment, "x\n", sep = "")
  final <- universes[[as.character(pick$multiple)]]
  mcols(final)$threshold_multiple <- pick$multiple
  saveRDS(final, file.path(OUT, "universe_final.rds"))
  saveRDS(universes, file.path(OUT, "universes_all_thresholds.rds"))
  cat("\nwrote data/processed/regions/universe_final.rds\n")
  cat("wrote data/processed/regions/universes_all_thresholds.rds (sensitivity analysis)\n")

  md <- c("# M5 — selected accessible-region universe", "",
          paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
          paste0("**Threshold x", pick$multiple, ", support >= ", MIN_LINES,
                 " keystone ATAC lines.**"), "",
          "Support level is fixed by the experimental design, not optimised: MYCN- and",
          "MYCL-amplified lines number two each, so requiring >=3 lines would delete",
          "every paralog-specific region for two of three paralogs.", "",
          "| x mean | regions | % genome | promoter % | TSS enrich | ext corrob % |",
          "|---|---|---|---|---|---|")
  for (i in seq_len(nrow(s)))
    md <- c(md, sprintf("| %s%s | %s | %s | %s | %s | %s |",
                        s$multiple[i], ifelse(s$multiple[i] == pick$multiple, " **(selected)**", ""),
                        format(s$n_regions[i], big.mark = ","), s$pct_genome[i],
                        s$pct_promoter[i], s$tss_enrichment[i], s$pct_ext_corrob[i]))
  md <- c(md, "", "## Accessibility class (descriptive — openness, not binding)", "",
          "| x mean | MYC lines only | MYCN only | MYCL1 only | shared |",
          "|---|---|---|---|---|")
  for (i in seq_len(nrow(s)))
    md <- c(md, sprintf("| %s | %s | %s | %s | %s |", s$multiple[i],
                        format(s$n_MYC_only[i], big.mark = ","),
                        format(s$n_MYCN_only[i], big.mark = ","),
                        format(s$n_MYCL1_only[i], big.mark = ","),
                        format(s$n_shared[i], big.mark = ",")))
  writeLines(md, "results/tables/m5_universe_selected.md")
  cat("wrote results/tables/m5_universe_selected.md\n")
} else {
  cat("NO threshold meets all targets at >=", MIN_LINES, " lines.\n", sep = "")
  cat("Do not proceed. Options: accept the external MACS2 peaks with the\n")
  cat("cell-line mismatch declared, or narrow Aim 1's scope. Both go in the log.\n")
  quit(status = 1)
}

# -----------------------------------------------------------------------------
# 02_consensus_regions.R — candidate accessible-region universes, with numbers.
#
#   Rscript scripts/02_regulatory/02_consensus_regions.R
#
# This script DELIBERATELY DOES NOT COMMIT to a universe. It builds the candidate
# options, measures each, and writes a report, because the choice of support rule
# is a scientific decision that should be made against real numbers rather than
# assumed from the frozen spec.
#
# The problem it exposes:
#
#   The frozen spec says a region needs support from >=2 independent ATAC datasets
#   (config: regions.min_datasets_supporting). After the M4 findings there are only
#   TWO usable external ATAC peak sources, and between them they cover four cell
#   lines of which only ONE (H524) is in the keystone MYC-family ChIP set:
#
#     GSE269424  H524, Lu139   (EGFP control arms only — D-015)
#     GSE256345  H146, H82
#
#   The keystone's own ATAC covers 9 lines and is the right cellular match, but it
#   ships as bedGraph SIGNAL with no peak calls, so it cannot contribute intervals
#   without a signal-thresholding step (handled in 03_).
#
#   Requiring >=2 datasets therefore means requiring a region to be accessible in
#   two different, largely non-keystone cell lines — which selects for SHARED
#   regulatory elements and against line-specific ones. For a project about
#   paralog-SPECIFIC regulation that filter runs against the hypothesis.
#
# Outputs: data/processed/regions/candidate_universes.rds
#          data/metadata/consensus_region_report.csv
#          results/tables/m5_universe_options.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(yaml)
})
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
RP  <- CFG$regions
OUT <- "data/processed/regions"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

LIFT <- "data/processed/liftover"
files <- list.files(LIFT, pattern = "\\.hg19\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)

# --- provenance: which dataset and cell line does each file represent? --------
meta <- data.frame(path = files, stringsAsFactors = FALSE)
meta$dataset <- sub("__.*$", "", basename(meta$path))
meta$line <- NA_character_
for (ln in c("H524", "Lu139", "H146", "H82")) {
  meta$line[grepl(ln, basename(meta$path), fixed = TRUE)] <- ln
}
KEYSTONE_LINES <- c("H1048","H211","H524","H847","SHP77","COLO668","H889","H526","H69","H196")
meta$in_keystone <- meta$line %in% KEYSTONE_LINES

cat("=== inputs ===\n")
for (i in seq_len(nrow(meta)))
  cat(sprintf("  %-12s %-7s keystone:%-5s %s\n", meta$dataset[i], meta$line[i],
              meta$in_keystone[i], basename(meta$path[i])))

grl <- lapply(meta$path, readRDS)
names(grl) <- basename(meta$path)
cat("\nintervals per file:", paste(lengths(grl), collapse = ", "), "\n")

# --- normalise each interval set ----------------------------------------------
# Merge within a file first: overlapping/adjacent calls in one sample are one
# element, and counting them separately would inflate apparent support.
norm_one <- function(gr) {
  gr <- gr[as.character(seqnames(gr)) %in% ANALYSIS_CHROMS_UCSC]
  gr <- GenomicRanges::reduce(gr, min.gapwidth = RP$merge_gap)
  gr[width(gr) >= RP$min_width]
}
grl <- lapply(grl, norm_one)
cat("after per-file merge:", paste(lengths(grl), collapse = ", "), "\n\n")

# --- per-DATASET consensus ----------------------------------------------------
# Within a dataset, take the union across its replicates/lines. Datasets are the
# unit of independence in the frozen spec.
by_ds <- split(seq_len(nrow(meta)), meta$dataset)
ds_gr <- lapply(by_ds, function(idx) {
  GenomicRanges::reduce(do.call(c, unname(grl[idx])), min.gapwidth = RP$merge_gap)
})
cat("=== per-dataset consensus ===\n")
for (d in names(ds_gr))
  cat(sprintf("  %-12s %8d regions, median width %5.0f bp\n",
              d, length(ds_gr[[d]]), median(width(ds_gr[[d]]))))

# --- per-CELL-LINE consensus (the more biologically meaningful unit) ----------
by_line <- split(seq_len(nrow(meta)), meta$line)
line_gr <- lapply(by_line, function(idx) {
  GenomicRanges::reduce(do.call(c, unname(grl[idx])), min.gapwidth = RP$merge_gap)
})
cat("\n=== per-cell-line consensus ===\n")
for (l in names(line_gr))
  cat(sprintf("  %-7s %8d regions  (keystone line: %s)\n",
              l, length(line_gr[[l]]), l %in% KEYSTONE_LINES))

# --- candidate universes ------------------------------------------------------
all_regions <- GenomicRanges::reduce(do.call(c, unname(ds_gr)),
                                     min.gapwidth = RP$merge_gap)
all_regions <- all_regions[width(all_regions) >= RP$min_width &
                           width(all_regions) <= RP$max_width]

support_ds   <- rowSums(vapply(ds_gr,   function(g) IRanges::overlapsAny(all_regions, g), logical(length(all_regions))))
support_line <- rowSums(vapply(line_gr, function(g) IRanges::overlapsAny(all_regions, g), logical(length(all_regions))))
mcols(all_regions)$support_datasets   <- support_ds
mcols(all_regions)$support_celllines  <- support_line

universes <- list(
  union_all        = all_regions,
  ds_support_ge2   = all_regions[support_ds   >= 2],
  line_support_ge2 = all_regions[support_line >= 2],
  line_support_ge3 = all_regions[support_line >= 3]
)

genome_bp <- sum(load_hg19_sizes()[ANALYSIS_CHROMS_UCSC])

rep <- do.call(rbind, lapply(names(universes), function(nm) {
  g <- universes[[nm]]
  data.frame(
    universe     = nm,
    n_regions    = length(g),
    median_width = if (length(g)) median(width(g)) else NA,
    total_mb     = round(sum(as.numeric(width(g))) / 1e6, 1),
    pct_genome   = round(100 * sum(as.numeric(width(g))) / genome_bp, 3),
    stringsAsFactors = FALSE)
}))

cat("\n=========== candidate universes ===========\n")
print(rep, row.names = FALSE)

# --- the decisive number: keystone-line representation ------------------------
# H524 is the only external ATAC line that is also a keystone MYC-family ChIP
# line. How much of each universe is supported by it?
h524 <- line_gr[["H524"]]
cat("\n=========== keystone-line coverage ===========\n")
cat("H524 is the ONLY external ATAC line present in the keystone ChIP set.\n\n")
for (nm in names(universes)) {
  g <- universes[[nm]]
  if (!length(g)) next
  f <- mean(IRanges::overlapsAny(g, h524))
  cat(sprintf("  %-17s %7d regions, %5.1f%% supported by H524\n",
              nm, length(g), 100 * f))
}

# --- how much does a >=2 rule cost? -------------------------------------------
cat("\n=========== cost of the >=2 support rule ===========\n")
u  <- universes$union_all
d2 <- universes$ds_support_ge2
cat(sprintf("  union                    : %7d regions\n", length(u)))
cat(sprintf("  >=2 datasets             : %7d regions  (%.1f%% of union retained)\n",
            length(d2), 100 * length(d2) / length(u)))
cat(sprintf("  regions unique to ONE dataset: %7d  (%.1f%% discarded)\n",
            length(u) - length(d2), 100 * (1 - length(d2) / length(u))))
cat(sprintf("  of those discarded, %.1f%% are H524-supported\n",
            100 * mean(IRanges::overlapsAny(u[support_ds < 2], h524))))

saveRDS(universes, file.path(OUT, "candidate_universes.rds"))
write.csv(rep, "data/metadata/consensus_region_report.csv", row.names = FALSE)

md <- c("# M5 — candidate accessible-region universes", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "**No universe is committed to by this script.** These are options measured",
        "against real data so the support rule can be chosen deliberately.", "",
        "| universe | regions | median width | total Mb | % genome |",
        "|---|---|---|---|---|")
for (i in seq_len(nrow(rep)))
  md <- c(md, sprintf("| `%s` | %s | %s | %s | %s |", rep$universe[i],
                      format(rep$n_regions[i], big.mark = ","), rep$median_width[i],
                      rep$total_mb[i], rep$pct_genome[i]))
md <- c(md, "",
        "## Cell-line composition problem",
        "",
        "External ATAC peak sources cover H524, Lu139, H146, H82.",
        "**Only H524 is also a keystone MYC-family ChIP line.**",
        "The keystone's own ATAC covers 9 lines but ships as bedGraph signal with",
        "no peak calls, so it contributes no intervals until 03_ derives regions",
        "from signal.",
        "",
        "A `>=2 independent datasets` rule therefore requires accessibility in two",
        "largely non-keystone cell lines, selecting for shared elements and against",
        "line-specific ones — the opposite of what a paralog-specificity question needs.")
writeLines(md, "results/tables/m5_universe_options.md")

cat("\nwrote data/processed/regions/candidate_universes.rds\n")
cat("wrote results/tables/m5_universe_options.md\n")
cat("\nRESULT: candidates built. The support rule is NOT yet decided.\n")

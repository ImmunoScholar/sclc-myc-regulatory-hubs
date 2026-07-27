# -----------------------------------------------------------------------------
# 15_paralog_regions.R — paralog-occupied and ACTIVE regions, under every
# defensible replicate convention.
#
#   Rscript scripts/02_regulatory/15_paralog_regions.R
#
# WHY THIS IS NOT A SINGLE CALCULATION
#
# The paralogs have unequal replicate counts: MYC 5 lines, MYCN 2, MYCL1 2. The
# rule for combining lines therefore changes the comparison in a way that has
# nothing to do with biology:
#
#   union (>=1 line)  — set size GROWS with replicate count, favouring MYC
#   >=2 lines         — intermediate
#   all lines         — set size SHRINKS with replicate count, penalising MYC
#
# A first attempt used "all lines" and produced MYC as the SMALLEST paralog set
# (4,977 vs MYCL1 7,956), inverting the published ordering (18,823 / 4,017 /
# 5,688). That was the convention, not the data.
#
# Plotnik's 18,823 MYC regions almost certainly come from pooling across their
# five MYC lines, i.e. a union — but their exact convention is not stated, and
# the choice largely determines whether the published MYCN-in-MYC value of 0.84
# is reproduced. So every convention is computed and the RANGE reported. Picking
# the convention that best matches the published number would be fitting the
# method to the expected answer.
#
# ALSO: a replicate-matched comparison. MYC is subsampled to 2 lines over all
# choose(5,2)=10 pairs, giving a like-for-like contrast against MYCN and MYCL1
# that is free of the replicate-count asymmetry entirely.
#
# ACTIVE = paralog-occupied AND H3K27ac-positive in the same line(s), mirroring
# Plotnik's binding-intersect-H3K27ac definition.
#
# Output: data/processed/regions/paralog_regions.rds
#         data/metadata/paralog_region_conventions.csv
#         results/tables/m5_paralog_regions.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(GenomicRanges); library(yaml) })

CFG <- yaml::read_yaml("config/params.yml")
# Read from config, not hardcoded. These were literals here until the D-028
# parameter sweep, which contradicted this project's own rule that a threshold
# hard-coded in a script is a bug.
FOLD_GRID    <- unlist(CFG$active_regions$fold_grid)
PRIMARY_FOLD <- CFG$active_regions$primary_fold
stopifnot(length(FOLD_GRID) > 0, is.finite(PRIMARY_FOLD))

nrm  <- readRDS("data/processed/signal/region_signal_normalised.rds")
M    <- nrm$mean_fob                  # fold over own background
meta <- nrm$meta
u    <- nrm$regions
n    <- nrow(M)

PARALOG_LINES <- list(
  MYC   = c("H1048","H211","H524","H847","SHP77"),
  MYCN  = c("H526","H69"),
  MYCL1 = c("COLO668","H889"))

col_for <- function(assay, line) {
  j <- which(meta$assay == assay & meta$line == line)
  if (length(j) == 1) j else NA_integer_
}

# H3K27ac-positive per line, at the same fold cutoff
k27_pos <- function(line, fold) {
  j <- col_for("H3K27ac", line)
  if (is.na(j)) return(rep(NA, n))
  M[, j] >= fold
}

# paralog-bound per line
bound <- function(paralog, line, fold) {
  j <- col_for(paralog, line)
  if (is.na(j)) return(rep(NA, n))
  M[, j] >= fold
}

# ACTIVE in a line = bound by the paralog AND H3K27ac-positive, same line
active_in_line <- function(paralog, line, fold) {
  b <- bound(paralog, line, fold); k <- k27_pos(line, fold)
  if (all(is.na(b)) || all(is.na(k))) return(NULL)
  b & k
}

combine <- function(mats, rule) {
  if (!length(mats)) return(rep(FALSE, n))
  m <- do.call(cbind, mats)
  switch(rule,
         union   = rowSums(m) >= 1,
         ge2     = rowSums(m) >= 2,
         all     = rowSums(m) == ncol(m))
}

paralog_active <- function(paralog, fold, rule, lines = NULL) {
  lns <- lines %||% PARALOG_LINES[[paralog]]
  mats <- Filter(Negate(is.null), lapply(lns, function(l) active_in_line(paralog, l, fold)))
  combine(mats, rule)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

rows <- list()
for (fold in FOLD_GRID) {
  for (rule in c("union", "ge2", "all")) {
    sets <- lapply(names(PARALOG_LINES), paralog_active, fold = fold, rule = rule)
    names(sets) <- names(PARALOG_LINES)
    nn <- vapply(sets, sum, integer(1))
    mycn_in_myc  <- if (nn[["MYCN"]]  > 0) sum(sets$MYCN  & sets$MYC) / nn[["MYCN"]]  else NA
    mycl1_in_myc <- if (nn[["MYCL1"]] > 0) sum(sets$MYCL1 & sets$MYC) / nn[["MYCL1"]] else NA
    rows[[length(rows) + 1L]] <- data.frame(
      fold = fold, rule = rule,
      n_MYC = nn[["MYC"]], n_MYCN = nn[["MYCN"]], n_MYCL1 = nn[["MYCL1"]],
      mycn_in_myc = round(mycn_in_myc, 3), mycl1_in_myc = round(mycl1_in_myc, 3),
      stringsAsFactors = FALSE)
  }
}
conv <- do.call(rbind, rows)

cat("=========== active-region counts and overlaps by convention ===========\n")
print(conv, row.names = FALSE)
cat("\nPublished (Plotnik 2024): MYC 18,823  MYCN 4,017  MYCL1 5,688;\n")
cat("MYCN-in-MYC ~0.84; MYCL1 'largely non-overlapping'.\n")

# --- replicate-matched: MYC subsampled to 2 lines -----------------------------
cat("\n=========== replicate-matched (MYC subsampled to 2 lines) ===========\n")
cat("Removes the replicate-count asymmetry entirely: every paralog gets 2 lines.\n\n")
pairs <- utils::combn(PARALOG_LINES$MYC, 2, simplify = FALSE)
mm <- list()
for (fold in FOLD_GRID) {
  for (rule in c("union", "all")) {
    vals <- vapply(pairs, function(p) {
      myc <- paralog_active("MYC", fold, rule, lines = p)
      mycn <- paralog_active("MYCN", fold, rule)
      if (sum(mycn) == 0) return(NA_real_)
      sum(mycn & myc) / sum(mycn)
    }, numeric(1))
    nvals <- vapply(pairs, function(p) sum(paralog_active("MYC", fold, rule, lines = p)), integer(1))
    mm[[length(mm) + 1L]] <- data.frame(
      fold = fold, rule = rule,
      median_n_MYC_2line = stats::median(nvals),
      mycn_in_myc_median = round(stats::median(vals, na.rm = TRUE), 3),
      mycn_in_myc_min = round(min(vals, na.rm = TRUE), 3),
      mycn_in_myc_max = round(max(vals, na.rm = TRUE), 3),
      stringsAsFactors = FALSE)
  }
}
matched <- do.call(rbind, mm)
print(matched, row.names = FALSE)

write.csv(conv, "data/metadata/paralog_region_conventions.csv", row.names = FALSE)

prim <- conv[conv$fold == PRIMARY_FOLD, ]
cat("\n=========== at the primary fold cutoff (", PRIMARY_FOLD, "x) ===========\n", sep = "")
print(prim, row.names = FALSE)
cat("\nThe spread across conventions is the point. If MYCN-in-MYC varies widely,\n")
cat("the published 0.84 cannot be reproduced or refuted without knowing which\n")
cat("convention generated it, and the gate must say so rather than pick one.\n")

sets_primary <- lapply(names(PARALOG_LINES), paralog_active,
                       fold = PRIMARY_FOLD, rule = "union")
names(sets_primary) <- names(PARALOG_LINES)
saveRDS(list(conventions = conv, matched = matched,
             sets_union_primary = sets_primary,
             regions = u, fold = PRIMARY_FOLD),
        "data/processed/regions/paralog_regions.rds")

md <- c("# M5 — paralog active regions by replicate convention", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "Paralogs have unequal replicate counts (MYC 5 lines, MYCN 2, MYCL1 2), so",
        "the line-combining rule changes set sizes for reasons unrelated to biology.",
        "Every convention is reported; none is selected for agreeing with the",
        "published value.", "",
        "| fold | rule | MYC | MYCN | MYCL1 | MYCN-in-MYC | MYCL1-in-MYC |",
        "|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(conv)))
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | **%s** | %s |",
                      conv$fold[i], conv$rule[i],
                      format(conv$n_MYC[i], big.mark = ","),
                      format(conv$n_MYCN[i], big.mark = ","),
                      format(conv$n_MYCL1[i], big.mark = ","),
                      conv$mycn_in_myc[i], conv$mycl1_in_myc[i]))
md <- c(md, "", "## Replicate-matched (MYC subsampled to 2 of 5 lines, all 10 pairs)", "",
        "| fold | rule | median MYC n | MYCN-in-MYC median | min | max |",
        "|---|---|---|---|---|---|")
for (i in seq_len(nrow(matched)))
  md <- c(md, sprintf("| %s | %s | %s | **%s** | %s | %s |",
                      matched$fold[i], matched$rule[i],
                      format(matched$median_n_MYC_2line[i], big.mark = ","),
                      matched$mycn_in_myc_median[i],
                      matched$mycn_in_myc_min[i], matched$mycn_in_myc_max[i]))
writeLines(md, "results/tables/m5_paralog_regions.md")

cat("\nwrote data/processed/regions/paralog_regions.rds\n")
cat("wrote results/tables/m5_paralog_regions.md\n")

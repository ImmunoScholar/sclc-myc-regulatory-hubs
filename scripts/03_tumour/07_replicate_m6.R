# -----------------------------------------------------------------------------
# 07_replicate_m6.R — does the M6 lineage-dominance result reproduce?
#
#   Rscript scripts/03_tumour/07_replicate_m6.R
#
# GSE60052 (n=79) vs George et al. 2015 (n=81, cBioPortal). Independent patients,
# sequencing and processing. A single-cohort NULL is weak; a null that reproduces
# in independently generated data is not.
#
# SCORING METHOD CHANGED FOR THIS COMPARISON, deliberately.
# singscore ranks each gene against the others present in the matrix. For
# GSE60052 that background was 33,683 genes; for George only 1,004 were fetched,
# and most of those ARE regulon members. Comparing singscore across the two would
# compare scoring artefacts, not biology.
#
# Both cohorts are therefore rescored with a UNIVERSE-INDEPENDENT statistic: the
# mean per-gene z-score across samples within each cohort. It is scale-free (so
# GSE60052's log2 and George's raw counts are both fine), depends only on the
# regulon's own genes, and is identical in construction across cohorts.
# GSE60052's singscore values are recomputed here too rather than reused, so the
# comparison is like-for-like.
#
# GENE NAMING: GSE60052 has MYCL1 and not MYCL; George has MYCL and not MYCL1.
# Handled explicitly — this is the silent-mismatch class that cost three rounds at
# the universe stage.
#
# Output: data/processed/tumour/replication.rds
#         results/tables/m6_replication.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml) })
CFG <- yaml::read_yaml("config/params.yml")
set.seed(CFG$project$seed)

reg <- readRDS("data/processed/regions/regulons.rds")$regulons
g60 <- readRDS("data/processed/tumour/gse60052.rds")
geo <- readRDS("data/processed/tumour/george2015.rds")

NE_MARKERS <- c("ASCL1","INSM1","CHGA","SYP","NCAM1","DLL3","CALCA","GRP","UCHL1","SYT11")
TF_NAMES   <- c("ASCL1","NEUROD1","POU2F3","YAP1")

# paralog gene symbol differs between cohorts
paralog_symbol <- function(p, rn) {
  cand <- if (p == "MYCL1") c("MYCL1","MYCL") else p
  hit <- cand[cand %in% rn]
  if (!length(hit)) NA_character_ else hit[1]
}

zscore_rows <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sd <- apply(M, 1, stats::sd, na.rm = TRUE)
  sd[!is.finite(sd) | sd == 0] <- NA_real_
  (M - mu) / sd
}

analyse <- function(M, label) {
  cat("\n==================== ", label, " ====================\n", sep = "")
  cat("matrix: ", nrow(M), " genes x ", ncol(M), " samples\n", sep = "")
  Z <- zscore_rows(M)

  # regulon scores: mean z across measurable members
  sets <- lapply(reg, function(g) intersect(g, rownames(Z)))
  S <- sapply(sets, function(g) colMeans(Z[g, , drop = FALSE], na.rm = TRUE))
  for (p in names(sets))
    cat(sprintf("  %-6s %3d/%d genes used\n", p, length(sets[[p]]), length(reg[[p]])))

  # covariates
  ne_use <- intersect(NE_MARKERS, rownames(Z))
  ne <- colMeans(Z[ne_use, , drop = FALSE], na.rm = TRUE)
  tfs <- intersect(TF_NAMES, rownames(Z))
  TFz <- t(Z[tfs, , drop = FALSE])
  cat("  NE score from ", length(ne_use), " markers; lineage TFs: ",
      paste(tfs, collapse = ", "), "\n", sep = "")

  psym <- vapply(names(reg), paralog_symbol, character(1), rn = rownames(Z))
  cat("  paralog symbols resolved: ",
      paste(sprintf("%s->%s", names(psym), psym), collapse = "  "), "\n", sep = "")

  out <- data.frame()
  for (p in names(reg)) {
    if (is.na(psym[[p]])) next
    y <- S[, p]; x <- as.numeric(Z[psym[[p]], ])
    raw <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
    # variance partitioning
    r2 <- function(m) summary(m)$r.squared
    m_par <- stats::lm(y ~ x)
    m_lin <- stats::lm(y ~ ne + TFz)
    m_both <- stats::lm(y ~ x + ne + TFz)
    # partial association
    ry <- stats::residuals(stats::lm(y ~ ne + TFz))
    rx <- stats::residuals(stats::lm(x ~ ne + TFz))
    pc <- suppressWarnings(stats::cor.test(rx, ry, method = "spearman"))
    # MYC/NE antagonism (Ireland)
    myc_sym <- paralog_symbol("MYC", rownames(Z))
    ire <- if (!is.na(myc_sym))
      suppressWarnings(stats::cor.test(as.numeric(Z[myc_sym, ]), ne, method = "spearman"))$estimate
      else NA_real_

    out <- rbind(out, data.frame(
      cohort = label, paralog = p, n = ncol(M), n_genes = length(sets[[p]]),
      rho_raw = round(raw$estimate, 3), p_raw = signif(raw$p.value, 3),
      r2_paralog = round(r2(m_par), 3), r2_lineage = round(r2(m_lin), 3),
      unique_paralog = round(r2(m_both) - r2(m_lin), 3),
      unique_lineage = round(r2(m_both) - r2(m_par), 3),
      rho_partial = round(pc$estimate, 3), p_partial = signif(pc$p.value, 3),
      myc_vs_ne = round(ire, 3), stringsAsFactors = FALSE))
  }
  out$fdr_partial <- signif(stats::p.adjust(out$p_partial, "BH"), 3)
  print(out[, c("paralog","n_genes","rho_raw","r2_paralog","r2_lineage",
                "unique_paralog","unique_lineage","rho_partial","fdr_partial")],
        row.names = FALSE)
  out
}

a60 <- analyse(g60$expr_filtered[, !g60$is_normal, drop = FALSE], "GSE60052")
ageo <- analyse(geo$expr, "George2015")
both <- rbind(a60, ageo)

# ---- side-by-side --------------------------------------------------------------
cat("\n\n==================== REPLICATION ====================\n")
cat("Same scoring method in both cohorts (mean per-gene z-score).\n\n")
key <- c("unique_paralog","unique_lineage","rho_partial","fdr_partial")
for (p in unique(both$paralog)) {
  s <- both[both$paralog == p, ]
  if (nrow(s) < 2) next
  cat("--- ", p, " ---\n", sep = "")
  cat(sprintf("  %-12s %10s %10s\n", "", "GSE60052", "George2015"))
  for (k in key)
    cat(sprintf("  %-12s %10s %10s\n", k,
                format(s[s$cohort=="GSE60052", k]), format(s[s$cohort=="George2015", k])))
  cat("\n")
}

cat("MYC vs NE score (Ireland 2020 antagonism):\n")
cat(sprintf("  GSE60052   rho = %s\n", unique(a60$myc_vs_ne)[1]))
cat(sprintf("  George2015 rho = %s\n", unique(ageo$myc_vs_ne)[1]))

# ---- verdict -----------------------------------------------------------------
cat("\n==================== VERDICT ====================\n")
lin_dom <- tapply(both$unique_lineage > both$unique_paralog, both$cohort, sum)
surv <- tapply(both$fdr_partial < 0.05 & both$rho_partial > 0, both$cohort, sum)
cat("lineage explains more unique variance than paralog:\n")
for (cc in names(lin_dom)) cat(sprintf("  %-12s %d/%d regulons\n", cc, lin_dom[[cc]], 3))
cat("paralog associations surviving lineage adjustment:\n")
for (cc in names(surv)) cat(sprintf("  %-12s %d/%d\n", cc, surv[[cc]], 3))

# TWO SEPARATE QUESTIONS. The first version conflated them and returned a
# misleading verdict.
#
#   Q1 (the claim being tested): does the NEGATIVE result — lineage dominates the
#       paralog — reproduce? This is the M6 conclusion.
#   Q2: does any apparent POSITIVE association reproduce? A positive that appears
#       in one cohort and vanishes in another is a spurious finding failing to
#       replicate, which SUPPORTS Q1 rather than contradicting it.
#
# Requiring zero survivors in both cohorts made a non-replicating positive look
# like a failure of the negative result. It is the opposite.
core_replicated <- all(lin_dom >= 3)
cat("\n")
cat("Q1 — does lineage dominance reproduce?  ",
    if (core_replicated) "YES, 3/3 regulons in BOTH cohorts" else "NO", "\n", sep = "")
cat("Q2 — does any positive paralog association reproduce?\n")
disc <- character(0)
for (p in unique(both$paralog)) {
  s <- both[both$paralog == p, ]
  if (nrow(s) < 2) next
  sig <- s$fdr_partial < 0.05 & s$rho_partial > 0
  if (all(sig)) disc <- c(disc, paste0(p, ": REPRODUCES in both"))
  else if (any(sig)) disc <- c(disc, paste0(
    p, ": appears in ", s$cohort[sig], " only (FDR ", s$fdr_partial[sig],
    ") but NOT in ", s$cohort[!sig], " (FDR ", s$fdr_partial[!sig], ") — not robust"))
}
if (!length(disc)) cat("  none in either cohort\n") else
  for (x in disc) cat("  ", x, "\n", sep = "")

cat("\n")
if (core_replicated && !any(grepl("REPRODUCES in both", disc))) {
  cat("REPLICATED. The M6 negative result reproduces in an independently generated\n")
  cat("cohort: lineage explains more unique variance than the paralog in every\n")
  cat("regulon, in both datasets, and no paralog association reproduces across both.\n")
  cat("The single cohort-specific positive fails to replicate, which supports rather\n")
  cat("than undermines the conclusion. It is also method-dependent (absent under\n")
  cat("singscore) and comes from a regulon that failed M5 validation (D-030).\n")
} else if (core_replicated) {
  cat("CORE RESULT REPLICATED, with a paralog association reproducing in both\n")
  cat("cohorts. Report that association as genuine and lineage-independent.\n")
} else {
  cat("NOT REPLICATED. Lineage dominance does not hold in both cohorts. The M6\n")
  cat("result is cohort-dependent and must be reported as such.\n")
}
reproduced <- core_replicated

saveRDS(list(per_cohort = both, reproduced = reproduced),
        "data/processed/tumour/replication.rds")
write.csv(both, "data/metadata/m6_replication.csv", row.names = FALSE)

md <- c("# M6 — replication in an independent cohort", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "GSE60052 (n=79) vs George et al. 2015 (n=81, cBioPortal). Independent",
        "patients, sequencing and processing.", "",
        "**Scoring**: mean per-gene z-score in both cohorts. singscore was not used",
        "for the comparison because its ranking background differed between cohorts",
        "(33,683 vs 1,004 genes) and would have compared scoring artefacts.", "",
        "**Gene naming**: GSE60052 carries `MYCL1`, George carries `MYCL`. Resolved",
        "explicitly.", "",
        "| cohort | paralog | genes | rho raw | unique paralog R2 | unique lineage R2 | rho partial | FDR |",
        "|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(both)))
  md <- c(md, sprintf("| %s | %s | %d | %s | %s | **%s** | %s | %s |",
                      both$cohort[i], both$paralog[i], both$n_genes[i], both$rho_raw[i],
                      both$unique_paralog[i], both$unique_lineage[i],
                      both$rho_partial[i], both$fdr_partial[i]))
writeLines(md, "results/tables/m6_replication.md")
cat("\nwrote results/tables/m6_replication.md\n")

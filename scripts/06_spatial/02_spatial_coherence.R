# -----------------------------------------------------------------------------
# 02_spatial_coherence.R — restricted orthogonal validation in spatial data.
#
#   Rscript scripts/06_spatial/02_spatial_coherence.R
#
# Two questions, both narrow on purpose (datasets.yml D5: restricted orthogonal
# validation only):
#
#   1. COHERENCE. Is the MYC regulon score a property of the tumour, or does it
#      vary from region to region within one? Variance partitioned between and
#      within slides.
#   2. REPLICATION. Does the M6/M7 result — regulon scores are explained by
#      neuroendocrine lineage state, not by the paralog's own expression — hold in
#      a third modality with independent tissue, chemistry and normalisation?
#
# THREE constraints bound every number below and are attached to every output:
#   * only the MYC regulon clears the coverage gate, on 56 of 500 members. This is
#     a 56-gene proxy, not the regulon. MYCN and MYCL1 are NOT scoreable.
#   * "slide" is not "tumour". Most slides carry two patients with no per-ROI
#     patient label, so within-slide variance is an UPPER bound on within-tumour
#     variance. Only 3 slides (IMfirst) are unambiguously one patient.
#   * the NE score rests on 3 markers of 10, because the rest are off-panel.
#
# Output: data/metadata/m9_coherence.csv, m9_variance.csv, m9_associations.csv
#         results/tables/m9_spatial.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml) })

CFG <- yaml::read_yaml("config/params.yml")
SP  <- CFG$spatial
set.seed(CFG$project$seed)

dat <- readRDS("data/processed/spatial/geomx.rds")
REG <- readRDS("data/processed/regions/regulons.rds")$regulons
cov <- read.csv("data/metadata/m9_panel_coverage.csv", stringsAsFactors = FALSE)

SCOREABLE <- unique(cov$paralog[cov$scoreable])
if (!length(SCOREABLE)) stop("nothing scoreable; 01_load_spatial.R should have stopped first")
cat("scoreable regulons: ", paste(SCOREABLE, collapse = ", "), "\n", sep = "")
cat("NOT scoreable     : ", paste(setdiff(names(REG), SCOREABLE), collapse = ", "),
    " — excluded from every result below\n\n", sep = "")

LIN <- c("ASCL1", "NEUROD1", "POU2F3", "YAP1")
NE_MARKERS <- c("CHGA","SYP","INSM1","NCAM1","ASCL1","DLL3","UCHL1","CALCA","GRP","CEACAM5")

# Mean per-gene z-score: universe-independent, and identical to the scoring used
# in M6 so the spatial and bulk results are comparable rather than merely similar.
zscore_mean <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (!length(g)) return(rep(NA_real_, ncol(mat)))
  m <- log2(mat[g, , drop = FALSE] + 1)
  z <- t(scale(t(m)))
  z <- z[stats::complete.cases(z), , drop = FALSE]
  colMeans(z)
}

# Unique R^2: the variance a term explains that no other term can. Same
# decomposition as the M6 lineage-confounding analysis.
unique_r2 <- function(y, full_terms, target, df) {
  ok <- stats::complete.cases(df[, full_terms, drop = FALSE]) & !is.na(y)
  if (sum(ok) < 10) return(NA_real_)
  f_full <- stats::lm(y[ok] ~ ., data = df[ok, full_terms, drop = FALSE])
  rest <- setdiff(full_terms, target)
  r2_full <- summary(f_full)$r.squared
  r2_rest <- if (!length(rest)) 0 else
    summary(stats::lm(y[ok] ~ ., data = df[ok, rest, drop = FALSE]))$r.squared
  max(0, r2_full - r2_rest)
}

coh <- data.frame(); vpart <- data.frame(); assoc <- data.frame()

for (id in names(dat)) {
  D <- dat[[id]]; mat <- D$mat; seg <- D$seg
  cat("=========== ", id, " (", D$label, ") ===========\n", sep = "")
  cat("  ROIs: ", ncol(mat), " across ", length(unique(seg$slide)), " slides\n", sep = "")

  score <- zscore_mean(mat, REG[["MYC"]])
  n_used <- length(intersect(REG[["MYC"]], rownames(mat)))
  ne_on <- intersect(NE_MARKERS, rownames(mat))
  ne <- zscore_mean(mat, ne_on)
  cat("  MYC regulon scored on ", n_used, "/500 members | NE score on ",
      length(ne_on), "/", length(NE_MARKERS), " markers (", paste(ne_on, collapse = ", "),
      ")\n", sep = "")

  expr <- function(g) if (g %in% rownames(mat)) as.numeric(log2(mat[g, ] + 1)) else rep(NA_real_, ncol(mat))
  df <- data.frame(MYC = expr("MYC"), NE = ne,
                   ASCL1 = expr("ASCL1"), NEUROD1 = expr("NEUROD1"),
                   POU2F3 = expr("POU2F3"), YAP1 = expr("YAP1"))

  # --- 1. coherence: between-slide vs within-slide --------------------------
  fit <- stats::lm(score ~ seg$slide)
  r2_between <- summary(fit)$r.squared
  cat("\n  COHERENCE\n")
  cat(sprintf("    between-slide R2 : %.3f   (within-slide %.3f)\n",
              r2_between, 1 - r2_between))
  coh <- rbind(coh, data.frame(
    cohort = id, subset = "all slides", n_roi = ncol(mat),
    n_slides = length(unique(seg$slide)),
    between_slide_r2 = round(r2_between, 4),
    within_slide_r2 = round(1 - r2_between, 4),
    slide_equals_tumour = FALSE, stringsAsFactors = FALSE))

  # The interpretable version: slides carrying a single patient identifier, where
  # within-slide really is within-tumour. Reported separately rather than pooled,
  # because pooling would let two-patient slides inflate "within-tumour" variance.
  if (any(seg$single_patient)) {
    k <- seg$single_patient
    fit1 <- stats::lm(score[k] ~ seg$slide[k])
    r2b1 <- summary(fit1)$r.squared
    cat(sprintf("    single-patient slides only (%d ROI, %d slides): between %.3f, within %.3f\n",
                sum(k), length(unique(seg$slide[k])), r2b1, 1 - r2b1))
    coh <- rbind(coh, data.frame(
      cohort = id, subset = "single-patient slides", n_roi = sum(k),
      n_slides = length(unique(seg$slide[k])),
      between_slide_r2 = round(r2b1, 4), within_slide_r2 = round(1 - r2b1, 4),
      slide_equals_tumour = TRUE, stringsAsFactors = FALSE))
  } else {
    cat("    single-patient slides: NONE — no unambiguous within-tumour estimate\n")
    coh <- rbind(coh, data.frame(
      cohort = id, subset = "single-patient slides", n_roi = 0L, n_slides = 0L,
      between_slide_r2 = NA_real_, within_slide_r2 = NA_real_,
      slide_equals_tumour = TRUE, stringsAsFactors = FALSE))
  }

  # --- 2. replication: lineage vs the paralog's own expression ---------------
  terms <- c("MYC", "NE", "ASCL1", "NEUROD1", "POU2F3", "YAP1")
  terms <- terms[vapply(df[terms], function(v) sum(!is.na(v)) > 10, logical(1))]
  u_myc <- unique_r2(score, terms, "MYC", df)
  lin_terms <- setdiff(terms, "MYC")
  ok <- stats::complete.cases(df[, terms, drop = FALSE])
  r2_full <- summary(stats::lm(score[ok] ~ ., data = df[ok, terms, drop = FALSE]))$r.squared
  r2_nolin <- summary(stats::lm(score[ok] ~ ., data = df[ok, "MYC", drop = FALSE]))$r.squared
  u_lin <- max(0, r2_full - r2_nolin)

  cat("\n  REPLICATION (does lineage or MYC explain the score?)\n")
  cat(sprintf("    unique R2, MYC expression : %.4f\n", u_myc))
  cat(sprintf("    unique R2, lineage block  : %.4f\n", u_lin))
  cat(sprintf("    full model R2             : %.4f\n", r2_full))
  vpart <- rbind(vpart, data.frame(
    cohort = id, n_roi = sum(ok), unique_myc = round(u_myc, 4),
    unique_lineage = round(u_lin, 4), full_r2 = round(r2_full, 4),
    lineage_dominates = u_lin > u_myc, stringsAsFactors = FALSE))

  # --- pairwise associations, raw and lineage-adjusted -----------------------
  for (v in terms) {
    rho <- suppressWarnings(stats::cor(score, df[[v]], method = "spearman",
                                       use = "complete.obs"))
    p <- suppressWarnings(stats::cor.test(score, df[[v]], method = "spearman",
                                          exact = FALSE))$p.value
    assoc <- rbind(assoc, data.frame(cohort = id, variable = v,
                                     rho = round(rho, 4), p = signif(p, 3),
                                     stringsAsFactors = FALSE))
  }
  # MYC association after removing lineage — the direct spatial analogue of the
  # partial correlation that collapsed in both bulk cohorts.
  if (all(c("MYC", "NE") %in% terms)) {
    rs <- stats::residuals(stats::lm(score[ok] ~ ., data = df[ok, lin_terms, drop = FALSE]))
    rm_ <- stats::residuals(stats::lm(df$MYC[ok] ~ ., data = df[ok, lin_terms, drop = FALSE]))
    pr <- suppressWarnings(stats::cor(rs, rm_, method = "spearman"))
    pp <- suppressWarnings(stats::cor.test(rs, rm_, method = "spearman", exact = FALSE))$p.value
    cat(sprintf("    MYC rho raw %+.3f -> lineage-adjusted %+.3f (p %.3g)\n",
                assoc$rho[assoc$cohort == id & assoc$variable == "MYC"], pr, pp))
    assoc <- rbind(assoc, data.frame(cohort = id, variable = "MYC (lineage-adjusted)",
                                     rho = round(pr, 4), p = signif(pp, 3),
                                     stringsAsFactors = FALSE))
  }
  cat("\n")
}

write.csv(coh,   "data/metadata/m9_coherence.csv", row.names = FALSE)
write.csv(vpart, "data/metadata/m9_variance.csv", row.names = FALSE)
write.csv(assoc, "data/metadata/m9_associations.csv", row.names = FALSE)

cat("=========== verdict ===========\n")
cat("  lineage dominates in ", sum(vpart$lineage_dominates), "/", nrow(vpart),
    " cohorts\n", sep = "")

mycov <- cov[cov$paralog == "MYC", ]
md <- c("# M9 — spatial coherence and orthogonal validation", "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "## Scope and what limits it", "",
  paste0("Restricted orthogonal validation only (datasets.yml D5). Two GeoMx DSP ",
         "cohorts on a targeted CTA panel: GSE261348 (IMfirst) and GSE261345 ",
         "(CANTABRICO)."), "",
  "**Only the MYC regulon can be scored.** Panel coverage against the gate",
  paste0("(>= ", SP$min_genes_measured, " genes AND >= ",
         100 * SP$min_coverage_fraction, "% of members):"), "",
  "| paralog | measured / size | coverage | scoreable |", "|---|---|---|---|",
  vapply(seq_len(nrow(mycov)), function(i) "", character(1))[0],
  {
    u <- cov[!duplicated(cov$paralog), ]
    sprintf("| %s | %d / %d | %.1f%% | %s |", u$paralog, u$measured, u$regulon_size,
            100 * u$coverage_fraction, ifelse(u$scoreable, "**yes**", "no"))
  }, "",
  paste0("MYCL1 fails by 0.4 percentage points and MYC clears by 1.2. These are ",
         "knife-edge verdicts against a threshold that was itself raised during ",
         "the project, and should be read as such rather than as clean ",
         "separations. The MYC result below is a **56-gene proxy**, not the ",
         "500-gene regulon."), "",
  "Two further limits apply to every number here:", "",
  paste0("1. **Slide is not tumour.** Most slides carry two patient identifiers ",
         "with no per-ROI patient label, so within-slide variance is an upper ",
         "bound on within-tumour variance. Only 3 IMfirst slides are ",
         "unambiguously one patient; CANTABRICO has none."),
  paste0("2. **The NE score rests on 3 of 10 markers** (NCAM1, ASCL1, DLL3); the ",
         "rest are off-panel."), "",
  "## Coherence — is the score a property of the tumour or the region?", "",
  "| cohort | subset | ROIs | slides | between-slide R² | within-slide R² |",
  "|---|---|---|---|---|---|",
  sprintf("| %s | %s | %d | %d | %s | %s |", coh$cohort, coh$subset, coh$n_roi,
          coh$n_slides,
          ifelse(is.na(coh$between_slide_r2), "—", sprintf("%.3f", coh$between_slide_r2)),
          ifelse(is.na(coh$within_slide_r2), "—", sprintf("%.3f", coh$within_slide_r2))), "",
  "## Replication — lineage versus the paralog's own expression", "",
  "| cohort | ROIs | unique R² MYC | unique R² lineage | full R² | lineage dominates |",
  "|---|---|---|---|---|---|",
  sprintf("| %s | %d | %.4f | %.4f | %.4f | %s |", vpart$cohort, vpart$n_roi,
          vpart$unique_myc, vpart$unique_lineage, vpart$full_r2,
          ifelse(vpart$lineage_dominates, "**yes**", "no")), "",
  "## Associations", "",
  "| cohort | variable | Spearman rho | p |", "|---|---|---|---|",
  sprintf("| %s | %s | %+.3f | %.3g |", assoc$cohort, assoc$variable, assoc$rho, assoc$p))
writeLines(md, "results/tables/m9_spatial.md")

cat("\nwrote data/metadata/m9_coherence.csv\n")
cat("wrote data/metadata/m9_variance.csv\n")
cat("wrote data/metadata/m9_associations.csv\n")
cat("wrote results/tables/m9_spatial.md\n")

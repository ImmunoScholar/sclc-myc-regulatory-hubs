# -----------------------------------------------------------------------------
# 02_score_regulons.R — score cell-line-derived regulons in patient tumours.
#
#   Rscript scripts/03_tumour/02_score_regulons.R
#
# THIS IS THE CENTRAL TEST OF AIM 2. Everything in M5 established that
# paralog-resolved regulatory programmes can be built from cell-line chromatin.
# The question the project exists to answer is whether they OPERATE IN PATIENTS.
#
# The decisive analysis is the SPECIFICITY MATRIX: correlate each regulon's score
# against each paralog's own expression across 79 tumours. If the MYC regulon
# tracks MYC expression more closely than MYCN or MYCL1 expression — i.e. the
# diagonal dominates — the programmes transfer AND retain paralog identity. If all
# regulons correlate with everything equally, they carry a generic
# proliferation/activity signal and paralog specificity does not survive the move
# from lines to tumours.
#
# A NULL RESULT HERE IS A REAL FINDING, not a failure (R-01, gap statement §5).
# It is reported with a power analysis rather than reframed.
#
# METHODS. singscore is primary: rank-based, single-sample, independent of cohort
# composition — which matters because a second cohort of different design follows.
# GSVA is the secondary sensitivity check (D-012 rationale: it is kernel-based and
# so methodologically distinct, not a restatement).
#
# Output: data/processed/tumour/regulon_scores.rds
#         data/metadata/m6_specificity_matrix.csv
#         results/tables/m6_regulon_scores.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(singscore); library(GSVA); library(yaml)
})

CFG <- yaml::read_yaml("config/params.yml")
set.seed(CFG$project$seed)

tum <- readRDS("data/processed/tumour/gse60052.rds")
reg <- readRDS("data/processed/regions/regulons.rds")
regulons <- reg$regulons
prog <- utils::read.csv("data/metadata/regulon_programme.csv", stringsAsFactors = FALSE)

E   <- tum$expr_filtered
isN <- tum$is_normal
Et  <- E[, !isN, drop = FALSE]          # 79 tumours
cat("tumour matrix: ", nrow(Et), " genes x ", ncol(Et), " samples\n\n", sep = "")

# restrict each regulon to measurable genes, and carry the coverage forward
sets <- lapply(regulons, function(g) intersect(g, rownames(Et)))
for (p in names(sets))
  cat(sprintf("  %-6s %3d of %3d genes measurable (%.1f%%)  programme: %s\n",
              p, length(sets[[p]]), length(regulons[[p]]),
              100 * length(sets[[p]]) / length(regulons[[p]]),
              ifelse(prog$pass[prog$paralog == p], prog$programme[prog$paralog == p],
                     "NONE (unvalidated, D-030)")))

# ---- scoring -----------------------------------------------------------------
cat("\n=========== scoring ===========\n")
rk <- singscore::rankGenes(Et)
ss <- sapply(sets, function(g)
  singscore::simpleScore(rk, upSet = g)$TotalScore)
rownames(ss) <- colnames(Et)
cat("singscore done\n")

gp <- GSVA::gsvaParam(exprData = as.matrix(Et), geneSets = sets)
gv <- t(GSVA::gsva(gp, verbose = FALSE))
cat("GSVA done\n")

agree <- diag(stats::cor(ss, gv[, colnames(ss)], method = "spearman"))
cat("\nsingscore vs GSVA agreement (Spearman, per regulon):\n")
for (p in names(agree)) cat(sprintf("  %-6s %.3f\n", p, agree[[p]]))
if (any(agree < 0.5))
  cat("  WARNING: a regulon where the two methods disagree is not robustly scored.\n")

# ---- the specificity matrix --------------------------------------------------
PARA_GENE <- c(MYC = "MYC", MYCN = "MYCN", MYCL1 = "MYCL1")
avail <- PARA_GENE[PARA_GENE %in% rownames(Et)]
cat("\n=========== specificity matrix (singscore) ===========\n")
cat("rows = regulon, cols = paralog expression. Diagonal should dominate.\n\n")

mk_matrix <- function(S) {
  M <- matrix(NA_real_, nrow = ncol(S), ncol = length(avail),
              dimnames = list(colnames(S), names(avail)))
  P <- M
  for (r in colnames(S)) for (cc in names(avail)) {
    ct <- suppressWarnings(stats::cor.test(S[, r], Et[avail[[cc]], ],
                                          method = "spearman"))
    M[r, cc] <- ct$estimate; P[r, cc] <- ct$p.value
  }
  list(rho = M, p = P)
}
sm <- mk_matrix(ss)
print(round(sm$rho, 3))
cat("\np-values:\n"); print(signif(sm$p, 3))

diag_wins <- vapply(rownames(sm$rho), function(r) {
  if (!(r %in% colnames(sm$rho))) return(NA)
  which.max(sm$rho[r, ]) == which(colnames(sm$rho) == r)
}, logical(1))
cat("\nregulon most strongly correlated with its OWN paralog:\n")
for (r in names(diag_wins))
  cat(sprintf("  %-6s %s\n", r, if (isTRUE(diag_wins[[r]])) "yes" else "NO"))
n_diag <- sum(diag_wins, na.rm = TRUE)
p_diag <- stats::pbinom(n_diag - 1, length(diag_wins), 1/length(avail), lower.tail = FALSE)
cat(sprintf("\n%d/%d diagonal; P(>=%d by chance) = %.3f\n",
            n_diag, length(diag_wins), n_diag, p_diag))

# ---- own-paralog association, with effect size and FDR ----------------------
cat("\n=========== own-paralog association ===========\n")
res <- data.frame()
for (p in intersect(colnames(ss), names(avail))) {
  x <- ss[, p]; y <- Et[avail[[p]], ]
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  # Fisher z CI for Spearman
  n <- length(x); z <- atanh(ct$estimate); se <- 1 / sqrt(n - 3)
  ci <- tanh(c(z - 1.96 * se, z + 1.96 * se))
  res <- rbind(res, data.frame(
    paralog = p, n = n, rho = round(ct$estimate, 3),
    ci_low = round(ci[1], 3), ci_high = round(ci[2], 3),
    p = ct$p.value, n_genes_used = length(sets[[p]]),
    programme_validated = prog$pass[prog$paralog == p], stringsAsFactors = FALSE))
}
res$fdr <- stats::p.adjust(res$p, method = "BH")
res$p <- signif(res$p, 3); res$fdr <- signif(res$fdr, 3)
print(res, row.names = FALSE)

# power: rho detectable at n=79, alpha 0.05, 80% power
cat(sprintf("\npower note: at n=%d, |rho| >= %.2f is detectable at alpha 0.05 with 80%% power.\n",
            ncol(Et), tanh(2.80 / sqrt(ncol(Et) - 3))))

# ---- tumour vs normal --------------------------------------------------------
cat("\n=========== tumour vs normal ===========\n")
rk_all <- singscore::rankGenes(E)
ss_all <- sapply(sets, function(g) singscore::simpleScore(rk_all, upSet = g)$TotalScore)
tn <- data.frame()
for (p in colnames(ss_all)) {
  a <- ss_all[!isN, p]; b <- ss_all[isN, p]
  wt <- suppressWarnings(stats::wilcox.test(a, b))
  tn <- rbind(tn, data.frame(paralog = p, tumour_median = round(median(a), 4),
                             normal_median = round(median(b), 4),
                             p = signif(wt$p.value, 3), stringsAsFactors = FALSE))
}
print(tn, row.names = FALSE)
cat("(n=7 normals — underpowered; reported for orientation, not inference)\n")

saveRDS(list(singscore = ss, gsva = gv, sets = sets, specificity = sm,
             own_association = res, tumour_vs_normal = tn, agreement = agree),
        "data/processed/tumour/regulon_scores.rds")
write.csv(as.data.frame(sm$rho), "data/metadata/m6_specificity_matrix.csv")
write.csv(res, "data/metadata/m6_own_association.csv", row.names = FALSE)

md <- c("# M6 — regulon scores in patient tumours (GSE60052, n=79)", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "Primary scorer singscore (rank-based, cohort-independent); GSVA as the",
        "methodologically distinct sensitivity check.", "",
        "| paralog | genes used | rho with own expression | 95% CI | FDR | programme |",
        "|---|---|---|---|---|---|")
for (i in seq_len(nrow(res)))
  md <- c(md, sprintf("| %s | %d/500 | %s | %s to %s | %s | %s |",
                      res$paralog[i], res$n_genes_used[i], res$rho[i],
                      res$ci_low[i], res$ci_high[i], res$fdr[i],
                      ifelse(res$programme_validated[i], "validated", "**unvalidated**")))
md <- c(md, "", "## Specificity matrix (Spearman rho)", "",
        paste0("| regulon | ", paste(colnames(sm$rho), collapse = " | "), " |"),
        paste0("|", paste(rep("---", ncol(sm$rho) + 1), collapse = "|"), "|"))
for (r in rownames(sm$rho))
  md <- c(md, paste0("| ", r, " | ", paste(round(sm$rho[r, ], 3), collapse = " | "), " |"))
md <- c(md, "",
        paste0("Diagonal dominance: ", n_diag, "/", length(diag_wins),
               " regulons correlate most strongly with their own paralog ",
               "(P = ", signif(p_diag, 3), " by chance)."), "",
        "A null result here is a reportable finding, not a failure: it would mean",
        "paralog-resolved programmes do not retain paralog identity in patient",
        "tumours (gap statement section 5).")
writeLines(md, "results/tables/m6_regulon_scores.md")

cat("\nwrote data/processed/tumour/regulon_scores.rds\n")
cat("wrote results/tables/m6_regulon_scores.md\n")

# -----------------------------------------------------------------------------
# 03_lineage_confounding.R — is the regulon signal paralog-specific or lineage?
#
#   Rscript scripts/03_tumour/03_lineage_confounding.R
#
# A PRIMARY ANALYSIS, not a robustness check (risk R-01, logged at M1 as the
# project's largest scientific threat). It became the decisive analysis when
# 02_score_regulons.R returned a null: no regulon tracked its own paralog's
# expression (MYC rho -0.202 with the WRONG SIGN, MYCN 0.135 ns, MYCL1 0.268
# borderline), diagonal dominance P = 0.259 — while every regulon was NEGATIVELY
# correlated with MYC expression, strongly for two (-0.416, -0.468).
#
# THE HYPOTHESIS UNDER TEST. Ireland et al. 2020 showed MYC drives SCLC AWAY from
# the neuroendocrine state (ASCL1 -> NEUROD1 -> YAP1). Our MYC regulon is
# neurogenesis-enriched (p = 7e-6, D-029). If MYC-high tumours are the least
# neuroendocrine, a neurogenesis-weighted regulon should score LOW in them — which
# is what we observe. Under that reading the regulons index NE/subtype identity
# rather than paralog activity, and R-01 has materialised.
#
# THE DECISIVE TEST is variance partitioning. If lineage/NE state explains far
# more of the regulon score than the paralog's own expression does, the regulons
# are subtype readouts and paralog specificity does not survive translation. That
# is a NEGATIVE RESULT the gap statement committed in advance to reporting (§5),
# with the framework and audit trail as the contribution.
#
# RESOLUTION IS REPORTED PER TF, never uniformly (R-14, config forbids a blanket
# claim). Tumour-level subtyping uses all four TFs because it is expression-based
# and well powered at n=79. Occupancy-level confounding is a separate analysis and
# is limited to POU2F3 (3 keystone lines), marginally ASCL1 (1), and is impossible
# for NEUROD1 (0).
#
# Output: data/processed/tumour/lineage_confounding.rds
#         data/metadata/m6_variance_partition.csv
#         results/tables/m6_lineage_confounding.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml) })

CFG <- yaml::read_yaml("config/params.yml")
LC  <- CFG$lineage_confounding
set.seed(CFG$project$seed)

tum <- readRDS("data/processed/tumour/gse60052.rds")
sc  <- readRDS("data/processed/tumour/regulon_scores.rds")
E   <- tum$expr_filtered[, !tum$is_normal, drop = FALSE]
S   <- sc$singscore
stopifnot(identical(rownames(S), colnames(E)))
n <- ncol(E)
cat("tumours: ", n, "\n\n", sep = "")

TFS <- c(unlist(LC$lineage_tfs_with_chip), unlist(LC$lineage_tfs_expression_only))
TFS <- TFS[TFS %in% rownames(E)]
cat("lineage TFs measurable: ", paste(TFS, collapse = ", "), "\n", sep = "")
cat("  with ChIP data: ", paste(intersect(unlist(LC$lineage_tfs_with_chip), TFS), collapse = ", "),
    "\n  expression only: ", paste(intersect(unlist(LC$lineage_tfs_expression_only), TFS), collapse = ", "),
    "\n\n", sep = "")

# ---- subtype assignment ------------------------------------------------------
# Standard SCLC convention: the highest-expressing lineage TF defines the subtype.
tfm <- E[TFS, , drop = FALSE]
tfz <- t(scale(t(tfm)))                       # z within TF, so scales are comparable
subtype <- rownames(tfz)[apply(tfz, 2, which.max)]
subtype <- factor(paste0("SCLC-", substr(subtype, 1, 1)))
cat("=========== subtype distribution ===========\n")
print(table(subtype))

# ---- neuroendocrine score ----------------------------------------------------
NE_MARKERS <- c("ASCL1","INSM1","CHGA","SYP","NCAM1","DLL3","CALCA","GRP","UCHL1","SYT11")
ne_use <- intersect(NE_MARKERS, rownames(E))
cat("\nNE score from ", length(ne_use), " canonical markers: ",
    paste(ne_use, collapse = ", "), "\n", sep = "")
ne <- colMeans(t(scale(t(E[ne_use, , drop = FALSE]))))
cat("NE score range: ", round(min(ne), 2), " to ", round(max(ne), 2), "\n", sep = "")

# ---- test the Ireland prediction --------------------------------------------
cat("\n=========== is MYC anti-correlated with the NE programme? ===========\n")
cat("Ireland 2020: MYC drives SCLC away from the NE state. If so, MYC expression\n")
cat("should be NEGATIVELY correlated with NE score and with ASCL1.\n\n")
PARA <- intersect(c("MYC","MYCN","MYCL1"), rownames(E))
ire <- data.frame()
for (p in PARA) {
  for (target in c("NE_score", TFS)) {
    y <- if (target == "NE_score") ne else E[target, ]
    ct <- suppressWarnings(stats::cor.test(E[p, ], y, method = "spearman"))
    ire <- rbind(ire, data.frame(paralog = p, target = target,
                                 rho = round(ct$estimate, 3),
                                 p = signif(ct$p.value, 3), stringsAsFactors = FALSE))
  }
}
ire$fdr <- signif(stats::p.adjust(ire$p, "BH"), 3)
print(ire, row.names = FALSE)
myc_ne <- ire$rho[ire$paralog == "MYC" & ire$target == "NE_score"]
cat(sprintf("\nMYC vs NE score: rho = %.3f -> %s\n", myc_ne,
            if (myc_ne < -0.2) "SUPPORTS the Ireland reading"
            else if (myc_ne > 0.2) "CONTRADICTS it"
            else "inconclusive"))

# ---- regulon scores by subtype ----------------------------------------------
cat("\n=========== regulon score by subtype ===========\n")
sub_tab <- data.frame()
for (r in colnames(S)) {
  kw <- stats::kruskal.test(S[, r] ~ subtype)
  meds <- tapply(S[, r], subtype, stats::median)
  sub_tab <- rbind(sub_tab, data.frame(
    regulon = r, kruskal_p = signif(kw$p.value, 3),
    t(round(meds, 4)), stringsAsFactors = FALSE, check.names = FALSE))
}
print(sub_tab, row.names = FALSE)

# ---- THE DECISIVE TEST: variance partitioning -------------------------------
cat("\n=========== variance partitioning (the decisive test) ===========\n")
cat("R2 of regulon score explained by: own paralog expression alone, lineage/NE\n")
cat("alone, and both. If lineage dominates, the regulons index subtype not paralog.\n\n")
vp <- data.frame()
for (r in intersect(colnames(S), PARA)) {
  y  <- S[, r]
  m_para <- stats::lm(y ~ E[r, ])
  m_lin  <- stats::lm(y ~ ne + t(tfz))
  m_both <- stats::lm(y ~ E[r, ] + ne + t(tfz))
  r2 <- function(m) summary(m)$r.squared
  # unique contribution of paralog = R2(both) - R2(lineage only)
  vp <- rbind(vp, data.frame(
    regulon = r,
    r2_paralog_only = round(r2(m_para), 3),
    r2_lineage_only = round(r2(m_lin), 3),
    r2_both = round(r2(m_both), 3),
    unique_paralog = round(r2(m_both) - r2(m_lin), 3),
    unique_lineage = round(r2(m_both) - r2(m_para), 3),
    stringsAsFactors = FALSE))
}
print(vp, row.names = FALSE)

dom <- vp$unique_lineage > vp$unique_paralog
cat("\nlineage explains more unique variance than the paralog in ",
    sum(dom), "/", nrow(vp), " regulons\n", sep = "")

# ---- partial correlation ----------------------------------------------------
cat("\n=========== paralog association, adjusted for lineage ===========\n")
cat("Does any paralog-regulon association survive removing lineage/NE state?\n\n")
pc <- data.frame()
for (r in intersect(colnames(S), PARA)) {
  ry <- stats::residuals(stats::lm(S[, r] ~ ne + t(tfz)))
  rx <- stats::residuals(stats::lm(E[r, ] ~ ne + t(tfz)))
  ct <- suppressWarnings(stats::cor.test(rx, ry, method = "spearman"))
  raw <- suppressWarnings(stats::cor.test(E[r, ], S[, r], method = "spearman"))
  pc <- rbind(pc, data.frame(
    regulon = r, rho_raw = round(raw$estimate, 3),
    rho_partial = round(ct$estimate, 3), p_partial = signif(ct$p.value, 3),
    stringsAsFactors = FALSE))
}
pc$fdr_partial <- signif(stats::p.adjust(pc$p_partial, "BH"), 3)
print(pc, row.names = FALSE)
survives <- sum(pc$fdr_partial < 0.05 & pc$rho_partial > 0)
cat("\nassociations surviving adjustment (positive, FDR<0.05): ", survives, "/", nrow(pc), "\n", sep = "")

# ---- verdict ----------------------------------------------------------------
cat("\n=========== VERDICT ===========\n")
if (survives == 0 && sum(dom) >= 2) {
  cat("R-01 HAS MATERIALISED. Regulon scores are explained substantially better by\n")
  cat("lineage/NE state than by the paralog's own expression, and no paralog\n")
  cat("association survives adjustment. Paralog-resolved programmes do NOT retain\n")
  cat("paralog identity in patient tumours independently of lineage state.\n")
  cat("This is the negative result the gap statement committed to reporting (S5).\n")
} else if (survives > 0) {
  cat("PARTIAL SURVIVAL: ", survives, " paralog association(s) persist after\n", sep = "")
  cat("adjusting for lineage. Report those as lineage-independent and the rest as\n")
  cat("confounded. Per-TF resolution caveats still apply (R-14).\n")
} else {
  cat("MIXED / INCONCLUSIVE — report descriptively; do not force a conclusion.\n")
}
cat("\nPER-TF RESOLUTION (R-14) — occupancy-level confounding is NOT uniform:\n")
cat("  POU2F3  within-line in 3 keystone lines (hg19)          strongest arm\n")
cat("  ASCL1   within-line in SHP-77 only, n=1 (hg38)          descriptive only\n")
cat("  NEUROD1 H446 only, ZERO keystone overlap                subtype-level only\n")
cat("  YAP1    NO ChIP data anywhere                           expression only\n")

write.csv(vp, "data/metadata/m6_variance_partition.csv", row.names = FALSE)
write.csv(pc, "data/metadata/m6_partial_correlation.csv", row.names = FALSE)
write.csv(ire, "data/metadata/m6_paralog_vs_lineage.csv", row.names = FALSE)
saveRDS(list(subtype = subtype, ne_score = ne, ireland = ire,
             by_subtype = sub_tab, variance = vp, partial = pc),
        "data/processed/tumour/lineage_confounding.rds")

md <- c("# M6 — lineage confounding (PRIMARY analysis, risk R-01)", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste0("GSE60052, n = ", n, " tumours."), "",
        "## Subtype distribution", "",
        paste0("`", paste(names(table(subtype)), table(subtype), sep = "=", collapse = "  "), "`"), "",
        "## Variance partitioning", "",
        "| regulon | R2 paralog only | R2 lineage only | R2 both | unique paralog | unique lineage |",
        "|---|---|---|---|---|---|")
for (i in seq_len(nrow(vp)))
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | **%s** |", vp$regulon[i],
                      vp$r2_paralog_only[i], vp$r2_lineage_only[i], vp$r2_both[i],
                      vp$unique_paralog[i], vp$unique_lineage[i]))
md <- c(md, "", "## Paralog association before and after adjusting for lineage", "",
        "| regulon | rho raw | rho partial | FDR |", "|---|---|---|---|")
for (i in seq_len(nrow(pc)))
  md <- c(md, sprintf("| %s | %s | %s | %s |", pc$regulon[i], pc$rho_raw[i],
                      pc$rho_partial[i], pc$fdr_partial[i]))
md <- c(md, "", "## Per-TF resolution (R-14)", "",
        "Occupancy-level confounding cannot be tested uniformly and no blanket",
        "\"we controlled for lineage TFs\" claim is permitted:", "",
        "- **POU2F3** — within-line in 3 keystone lines, hg19. Strongest arm.",
        "- **ASCL1** — within-line in SHP-77 only (n=1), hg38. Descriptive.",
        "- **NEUROD1** — H446 only, zero keystone overlap. Subtype-level only.",
        "- **YAP1** — no ChIP data in any acquired dataset. Expression only.")
writeLines(md, "results/tables/m6_lineage_confounding.md")

cat("\nwrote results/tables/m6_lineage_confounding.md\n")

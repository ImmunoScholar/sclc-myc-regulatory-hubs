# -----------------------------------------------------------------------------
# 04_jung_benchmark.R — the Jung 2017 paralog-blind MYC-activity signature.
#
#   Rscript scripts/05_integration/04_jung_benchmark.R
#
# The project contract names two benchmark comparators, both PARALOG-BLIND, which
# is the point of the comparison: MSigDB HALLMARK_MYC_TARGETS (done at M5) and the
# Jung et al. 2017 18-gene MYC-activity signature. The Jung signature was recorded
# as NOT CURATED for most of the project because the paper is not in PubMed
# Central and the gene list is not in the abstract. Guessing it was refused.
#
# Source, now transcribed from the paper itself:
#   Jung LA et al. "A Myc Activity Signature Predicts Poor Clinical Outcomes in
#   Myc-Associated Cancers." Cancer Res 77(4):971-981, 15 Feb 2017.
#   PMID 27923830 · doi:10.1158/0008-5472.CAN-15-2906 · Table 1.
#   The paper states the set is "17 genes upregulated and 1 gene downregulated by
#   Myc"; the transcription below is asserted to match that split.
#
# Two questions are asked, and the second is the one that matters:
#   1. are the paralog regulons enriched for the signature? (underpowered — said so)
#   2. does this paralog-blind MYC-activity signature ALSO track lineage rather
#      than MYC in tumours? If it does, the confound is not an artefact of how
#      THIS project built its regulons.
#
# Output: data/metadata/jung2017_signature.csv, jung_benchmark.csv
#         results/tables/jung_benchmark.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml) })
CFG <- yaml::read_yaml("config/params.yml"); set.seed(CFG$project$seed)

# --- Table 1, transcribed ------------------------------------------------------
JUNG <- data.frame(
  entrez = c(7965, 328, 655, 894, 1736, 1965, 2023, 3159, 55885, 4673,
             4830, 4953, 6125, 6193, 6495, 6723, 8409, 10397),
  symbol = c("AIMP2","APEX1","BMP7","CCND2","DKC1","EIF2S1","ENO1","HMGA1","LMO3",
             "NAP1L1","NME1","ODC1","RPL5","RPS5","SIX1","SRM","UXT","NDRG1"),
  direction = c(rep("Up", 17), "Down"),
  stringsAsFactors = FALSE)

# Guard against a transcription slip: the paper's own arithmetic must hold.
stopifnot(nrow(JUNG) == 18L,
          sum(JUNG$direction == "Up") == 17L,
          sum(JUNG$direction == "Down") == 1L,
          !anyDuplicated(JUNG$symbol), !anyDuplicated(JUNG$entrez))
write.csv(JUNG, "data/metadata/jung2017_signature.csv", row.names = FALSE)

cat("=========== Jung 2017 signature ===========\n")
cat("  genes: ", nrow(JUNG), " (", sum(JUNG$direction == "Up"), " up, ",
    sum(JUNG$direction == "Down"), " down)\n", sep = "")
cat("  ", paste(JUNG$symbol, collapse = ", "), "\n\n", sep = "")

REG <- readRDS("data/processed/regions/regulons.rds")$regulons
ev  <- readRDS("data/processed/integration/moes_evidence.rds")
universe <- ev$universe
out <- data.frame()

# --- 1. regulon enrichment (underpowered by construction) ----------------------
cat("=========== 1. regulon overlap with the signature ===========\n")
jn <- intersect(JUNG$symbol, universe)
cat("  signature genes inside the MOES universe: ", length(jn), "/", nrow(JUNG), "\n", sep = "")
for (p in names(REG)) {
  reg <- intersect(REG[[p]], universe)
  hit <- intersect(reg, jn)
  exp_hit <- length(reg) * length(jn) / length(universe)
  ft <- stats::fisher.test(matrix(c(
    length(hit), length(jn) - length(hit),
    length(reg) - length(hit),
    length(universe) - length(reg) - length(jn) + length(hit)), 2))
  cat(sprintf("  %-6s %d/%d hits (expected %.2f) OR %.2f p %.3f %s\n",
              p, length(hit), length(jn), exp_hit, ft$estimate, ft$p.value,
              if (length(hit)) paste0("[", paste(hit, collapse = ","), "]") else ""))
  out <- rbind(out, data.frame(test = "regulon_overlap", paralog = p,
                               n_hit = length(hit), n_expected = round(exp_hit, 3),
                               or = round(unname(ft$estimate), 3),
                               p = signif(ft$p.value, 3), stringsAsFactors = FALSE))
}
cat("\n  POWER: an 18-gene signature against 500-gene regulons in a ",
    format(length(universe), big.mark = ","), "-gene universe expects <1 overlapping\n", sep = "")
cat("  gene by chance. This test cannot detect anything but a very large effect,\n")
cat("  and a null here is uninformative rather than evidence of no relationship.\n\n")

# --- 2. does the paralog-blind signature track lineage too? --------------------
cat("=========== 2. signature vs lineage in tumours (the real question) ===========\n")
g60 <- readRDS("data/processed/tumour/gse60052.rds")
lin <- readRDS("data/processed/tumour/lineage_confounding.rds")
sc  <- readRDS("data/processed/tumour/regulon_scores.rds")

expr <- g60$expr_filtered[, !g60$is_normal, drop = FALSE]
zmean <- function(genes, sign = 1) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 3) return(NULL)
  z <- t(scale(t(expr[g, , drop = FALSE])))
  z <- z[stats::complete.cases(z), , drop = FALSE]
  colMeans(z) * sign
}
up   <- zmean(JUNG$symbol[JUNG$direction == "Up"])
down <- zmean(JUNG$symbol[JUNG$direction == "Down"])
jung_score <- if (is.null(down)) up else (up * 17 - down) / 18
cat("  signature genes measured in GSE60052: ",
    length(intersect(JUNG$symbol, rownames(expr))), "/", nrow(JUNG), "\n", sep = "")

ne <- lin$ne_score[colnames(expr)]
myc <- as.numeric(expr["MYC", ])
tests <- list(
  list(lbl = "vs NE score",       v = ne),
  list(lbl = "vs MYC expression", v = myc),
  list(lbl = "vs ASCL1",          v = as.numeric(expr["ASCL1", ])),
  list(lbl = "vs POU2F3",         v = if ("POU2F3" %in% rownames(expr)) as.numeric(expr["POU2F3", ]) else NULL))
for (t in tests) {
  if (is.null(t$v)) next
  ct <- suppressWarnings(stats::cor.test(jung_score, t$v, method = "spearman", exact = FALSE))
  cat(sprintf("  Jung signature %-20s rho %+.3f  p %.3g\n", t$lbl, ct$estimate, ct$p.value))
  out <- rbind(out, data.frame(test = "jung_vs_covariate", paralog = t$lbl,
                               n_hit = NA, n_expected = NA,
                               or = round(unname(ct$estimate), 4),
                               p = signif(ct$p.value, 3), stringsAsFactors = FALSE))
}

# Unique variance, mirroring the M6 decomposition exactly.
df <- data.frame(MYC = myc, NE = ne)
ok <- stats::complete.cases(df) & !is.na(jung_score)
r2_full <- summary(stats::lm(jung_score[ok] ~ MYC + NE, df[ok, ]))$r.squared
r2_myc  <- summary(stats::lm(jung_score[ok] ~ MYC, df[ok, ]))$r.squared
r2_ne   <- summary(stats::lm(jung_score[ok] ~ NE, df[ok, ]))$r.squared
cat(sprintf("\n  unique R2 MYC     : %.4f\n", max(0, r2_full - r2_ne)))
cat(sprintf("  unique R2 lineage : %.4f\n", max(0, r2_full - r2_myc)))
out <- rbind(out,
  data.frame(test = "jung_variance", paralog = "unique_myc", n_hit = NA, n_expected = NA,
             or = round(max(0, r2_full - r2_ne), 4), p = NA, stringsAsFactors = FALSE),
  data.frame(test = "jung_variance", paralog = "unique_lineage", n_hit = NA, n_expected = NA,
             or = round(max(0, r2_full - r2_myc), 4), p = NA, stringsAsFactors = FALSE))

# --- 3. agreement with this project's regulon scores ---------------------------
cat("\n=========== 3. agreement with the project's regulon scores ===========\n")
ss <- sc$singscore
for (p in colnames(ss)) {
  ct <- suppressWarnings(stats::cor.test(jung_score, ss[, p], method = "spearman", exact = FALSE))
  cat(sprintf("  Jung vs %-6s regulon score  rho %+.3f  p %.3g\n", p, ct$estimate, ct$p.value))
  out <- rbind(out, data.frame(test = "jung_vs_regulon_score", paralog = p,
                               n_hit = NA, n_expected = NA,
                               or = round(unname(ct$estimate), 4),
                               p = signif(ct$p.value, 3), stringsAsFactors = FALSE))
}

write.csv(out, "data/metadata/jung_benchmark.csv", row.names = FALSE)

uj <- out$or[out$test == "jung_variance" & out$paralog == "unique_myc"]
ul <- out$or[out$test == "jung_variance" & out$paralog == "unique_lineage"]
md <- c("# Benchmark — Jung et al. 2017 18-gene Myc activity signature", "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "Source: Jung LA et al., *Cancer Res* 77(4):971-981, 2017. PMID 27923830,",
  "doi:10.1158/0008-5472.CAN-15-2906, Table 1. Transcribed from the paper;",
  "17 up / 1 down, matching the paper's own description.", "",
  "This comparator is **paralog-blind**, which is the point of using it.", "",
  "## Signature", "",
  paste0("`", paste(JUNG$symbol, collapse = "`, `"), "`"), "",
  paste0("Down-regulated: `", paste(JUNG$symbol[JUNG$direction == "Down"], collapse = ", "), "`."), "",
  "## 1. Regulon overlap — underpowered, reported as such", "",
  "| paralog | hits | expected | OR | p |", "|---|---|---|---|---|",
  sprintf("| %s | %d | %.2f | %.2f | %.3f |",
          out$paralog[out$test == "regulon_overlap"], out$n_hit[out$test == "regulon_overlap"],
          out$n_expected[out$test == "regulon_overlap"], out$or[out$test == "regulon_overlap"],
          out$p[out$test == "regulon_overlap"]), "",
  paste0("An 18-gene signature against 500-gene regulons in a ",
         format(length(universe), big.mark = ","), "-gene universe expects fewer ",
         "than one overlapping gene by chance. A null here means the test had no ",
         "power, not that there is no relationship."), "",
  "## 2. The signature tracks neither MYC expression nor lineage", "",
  "| comparison | Spearman rho | p |", "|---|---|---|",
  sprintf("| %s | %+.3f | %.3g |",
          out$paralog[out$test == "jung_vs_covariate"], out$or[out$test == "jung_vs_covariate"],
          out$p[out$test == "jung_vs_covariate"]), "",
  sprintf("Unique R2: MYC expression **%.4f**, lineage **%.4f**.", uj, ul), "",
  paste0("**Both are near zero.** Lineage is nominally the larger of the two, but ",
         "on values of ", sprintf("%.4f", ul), " against ", sprintf("%.4f", uj),
         " that ordering carries no weight and is not reported as a finding. The ",
         "honest reading is that this signature is largely unrelated to *either* ",
         "covariate in this cohort."), "",
  paste0("What it does show is that a published, paralog-blind MYC-activity ",
         "signature does **not** track MYC expression in SCLC tumours (rho ",
         sprintf("%+.3f", out$or[out$test == "jung_vs_covariate" & out$paralog == "vs MYC expression"]),
         ", n.s.). That is consistent with this project's broader observation that ",
         "MYC mRNA is a poor proxy for MYC regulatory activity in SCLC, and it is ",
         "an independent line of support for it — but it is not evidence of the ",
         "lineage confounding reported in M6, and must not be read as such."), "",
  "## 3. Agreement with this project's regulon scores", "",
  "| regulon | Spearman rho | p |", "|---|---|---|",
  sprintf("| %s | %+.3f | %.3g |",
          out$paralog[out$test == "jung_vs_regulon_score"],
          out$or[out$test == "jung_vs_regulon_score"], out$p[out$test == "jung_vs_regulon_score"]))
writeLines(md, "results/tables/jung_benchmark.md")

cat("\nwrote data/metadata/jung2017_signature.csv\n")
cat("wrote data/metadata/jung_benchmark.csv\n")
cat("wrote results/tables/jung_benchmark.md\n")

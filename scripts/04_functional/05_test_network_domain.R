# -----------------------------------------------------------------------------
# 05_test_network_domain.R — is CCLE SCLC expression usable for the network domain?
#
#   Rscript scripts/04_functional/05_test_network_domain.R
#
# TEST BEFORE BUILD. The tumour-based network domain was barred (D-033) because
# GENIE3 recovers regulators from expression, and expression in SCLC tumours is
# dominated by neuroendocrine lineage state. CCLE SCLC cell lines were proposed as
# the replacement — but SCLC lines also span the SCLC-A/N/P/Y subtypes, so the same
# confound may apply. Building GENIE3 first and checking afterwards would risk
# admitting a confounded domain that LOOKS independent.
#
# TWO PRE-CONDITIONS, both of which can fail:
#   1. POWER. GENIE3 conventionally wants hundreds of samples. How many SCLC lines
#      actually have expression?
#   2. CONFOUNDING. Does paralog expression across SCLC lines track NE/lineage
#      state, as it did across tumours (MYC vs NE rho -0.59 and -0.47)?
#
# If paralog expression is lineage-driven within cell lines too, a GENIE3 network
# built here inherits the confound and the domain must be excluded — on evidence,
# the same standard applied to the transcriptional domain.
#
# Output: data/metadata/m7_network_feasibility.csv
#         results/tables/m7_network_feasibility.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(yaml) })
CFG <- yaml::read_yaml("config/params.yml")
DM  <- CFG$depmap
f <- file.path("data/raw/depmap", DM$files$expression)
stopifnot(file.exists(f))
cat("release: ", DM$release, "\nfile   : ", basename(f), "\n\n", sep = "")

cat("reading expression matrix...\n")
d <- data.table::fread(f, showProgress = FALSE)
cat("dimensions: ", nrow(d), " lines x ", ncol(d), " columns\n\n", sep = "")

META <- c("depmap_id","cell_line_display_name","lineage_1","lineage_2",
          "lineage_3","lineage_4","lineage_6")
meta_present <- intersect(META, names(d))

# exact subtype match — no fallback (the fallback swept in the NSCLC panel, D-036)
sub_col <- if ("lineage_3" %in% names(d)) "lineage_3" else "lineage_2"
subt <- trimws(as.character(d[[sub_col]]))
is_sclc <- grepl("^small cell lung cancer$", subt, ignore.case = TRUE)

cat("=========== pre-condition 1: POWER ===========\n")
cat("SCLC lines with expression: ", sum(is_sclc), "\n", sep = "")
cat("GENIE3 convention: hundreds of samples for stable importance estimates.\n")
power_ok <- sum(is_sclc) >= 100
cat("power adequate (>=100): ", if (power_ok) "YES" else "NO", "\n", sep = "")
if (!power_ok)
  cat("  At n=", sum(is_sclc), ", tree-ensemble importances are unstable and a\n",
      "  network built here would carry uncertainty no bootstrap can repair.\n", sep = "")

E <- d[is_sclc, ]
lines_sclc <- E$cell_line_display_name
gene_cols <- setdiff(names(E), meta_present)

get_gene <- function(g) if (g %in% gene_cols) as.numeric(E[[g]]) else NULL

cat("\n=========== pre-condition 2: CONFOUNDING ===========\n")
NE_MARKERS <- c("ASCL1","INSM1","CHGA","SYP","NCAM1","DLL3","CALCA","GRP","UCHL1","SYT11")
ne_use <- intersect(NE_MARKERS, gene_cols)
cat("NE markers available: ", length(ne_use), " (", paste(ne_use, collapse=", "), ")\n", sep = "")
NEm <- as.matrix(E[, ..ne_use]); storage.mode(NEm) <- "double"
ne <- rowMeans(scale(NEm), na.rm = TRUE)

TFS <- intersect(c("ASCL1","NEUROD1","POU2F3","YAP1"), gene_cols)
PARA <- intersect(c("MYC","MYCN","MYCL","MYCL1"), gene_cols)
cat("lineage TFs: ", paste(TFS, collapse=", "), "\n", sep = "")
cat("paralogs   : ", paste(PARA, collapse=", "), "\n\n", sep = "")

cat("Does paralog expression track NE/lineage state ACROSS SCLC CELL LINES?\n")
cat("(In tumours it did: MYC vs NE rho -0.590 and -0.466 — the reason the\n")
cat(" tumour-based network domain was barred.)\n\n")
res <- data.frame()
for (p in PARA) {
  v <- get_gene(p); if (is.null(v)) next
  for (tgt in c("NE_score", TFS)) {
    y <- if (tgt == "NE_score") ne else get_gene(tgt)
    ct <- suppressWarnings(stats::cor.test(v, y, method = "spearman"))
    res <- rbind(res, data.frame(paralog = p, target = tgt,
                                 rho = round(ct$estimate, 3),
                                 p = signif(ct$p.value, 3), stringsAsFactors = FALSE))
  }
}
res$fdr <- signif(stats::p.adjust(res$p, "BH"), 3)
print(res, row.names = FALSE)

myc_ne <- res$rho[res$paralog == "MYC" & res$target == "NE_score"]
strong <- res[abs(res$rho) > 0.4 & res$fdr < 0.05, ]
cat("\nMYC vs NE score in SCLC cell lines: rho = ", myc_ne, "\n", sep = "")
cat("paralog-lineage associations with |rho|>0.4 and FDR<0.05: ", nrow(strong), "\n", sep = "")
if (nrow(strong)) {
  for (i in seq_len(nrow(strong)))
    cat(sprintf("  %s vs %s: rho %.3f (FDR %.3g)\n", strong$paralog[i],
                strong$target[i], strong$rho[i], strong$fdr[i]))
}
confounded <- abs(myc_ne) > 0.4 || nrow(strong) >= 2

cat("\n=========== VERDICT ===========\n")
if (!power_ok && confounded) {
  cat("EXCLUDE the network domain. It fails BOTH pre-conditions: n=", sum(is_sclc),
      " is far\n", sep = "")
  cat("below what tree-ensemble inference needs, AND paralog expression tracks\n")
  cat("lineage state in cell lines just as it did in tumours. A GENIE3 network\n")
  cat("built here would be underpowered and would carry the same confound that\n")
  cat("removed the transcriptional domain, while appearing to be independent\n")
  cat("evidence. MOES runs on TWO domains: cis-regulatory and functional.\n")
} else if (!power_ok) {
  cat("EXCLUDE on power. n=", sum(is_sclc), " is far below the hundreds of samples\n", sep = "")
  cat("GENIE3 needs. Not confounded, but not estimable either.\n")
} else if (confounded) {
  cat("EXCLUDE on confounding. Adequately powered, but paralog expression tracks\n")
  cat("lineage state, so the network would inherit the D-033 confound.\n")
} else {
  cat("ADMIT the network domain. Adequately powered and paralog expression does\n")
  cat("NOT track lineage state in cell lines — unlike in tumours. Build GENIE3.\n")
}

out <- data.frame(n_sclc_lines = sum(is_sclc), power_ok = power_ok,
                  myc_vs_ne_rho = myc_ne, n_strong_lineage_assoc = nrow(strong),
                  confounded = confounded,
                  admit_network_domain = power_ok && !confounded,
                  stringsAsFactors = FALSE)
write.csv(out, "data/metadata/m7_network_feasibility.csv", row.names = FALSE)
write.csv(res, "data/metadata/m7_paralog_lineage_celllines.csv", row.names = FALSE)

md <- c("# M7 — network domain feasibility", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste0("DepMap ", DM$release, " expression, ", sum(is_sclc), " SCLC cell lines."), "",
        "Tested BEFORE building, so a confounded domain cannot be admitted while",
        "appearing independent.", "",
        paste0("- **Power**: n = ", sum(is_sclc), " SCLC lines. GENIE3 conventionally needs hundreds. ",
               if (power_ok) "Adequate." else "**Inadequate.**"),
        paste0("- **Confounding**: MYC vs NE score rho = ", myc_ne,
               " (tumours gave -0.590 and -0.466). ",
               if (confounded) "**Confounded.**" else "Not confounded."), "",
        paste0("**Verdict: ", if (power_ok && !confounded) "ADMIT" else "EXCLUDE",
               "** the network domain."), "",
        "| paralog | target | rho | FDR |", "|---|---|---|---|")
for (i in seq_len(nrow(res)))
  md <- c(md, sprintf("| %s | %s | %s | %s |", res$paralog[i], res$target[i],
                      res$rho[i], res$fdr[i]))
writeLines(md, "results/tables/m7_network_feasibility.md")
cat("\nwrote results/tables/m7_network_feasibility.md\n")

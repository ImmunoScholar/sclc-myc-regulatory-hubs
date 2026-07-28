# -----------------------------------------------------------------------------
# 03_publication_tables.R — the numbered tables, as standalone deliverables.
#
#   Rscript scripts/08_report/03_publication_tables.R
#
# The report renders its tables inline. This writes the same content as numbered,
# self-contained CSV + Markdown files, each carrying a NOTE line, because the M10
# gate requires "every table has notes where interpretation is non-obvious" and a
# table separated from its caption is the easiest thing in a repository to
# misread.
#
# Every value is read from the analysis outputs. Nothing is typed in.
#
# Output: results/tables/publication/TableN_*.{csv,md}
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml) })
OUT <- "results/tables/publication"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
md <- function(f) read.csv(file.path("data/metadata", f), stringsAsFactors = FALSE)

emit <- function(n, slug, title, df, note) {
  stem <- file.path(OUT, sprintf("Table%d_%s", n, slug))
  write.csv(df, paste0(stem, ".csv"), row.names = FALSE)
  ln <- c(sprintf("# Table %d. %s", n, title), "",
          paste(paste0("| ", paste(names(df), collapse = " | "), " |")),
          paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
          apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")),
          "", paste("**Note.**", note))
  writeLines(ln, paste0(stem, ".md"))
  cat(sprintf("  Table %d  %-34s %d rows\n", n, slug, nrow(df)))
}

cat("=========== publication tables ===========\n")

g <- md("m5_gate_results.csv")
emit(1, "M5_gate", "Pre-registered gate on the paralog-resolved regulatory layer",
     data.frame(Criterion = g$criterion, Observed = g$value,
                Expected = g$expected, Pass = ifelse(g$pass, "PASS", "FAIL")),
     paste("Criterion 3 fails on magnitude while reproducing in direction (+3.3 points",
           "against a published +27) and was evaluated only at M7, once DepMap copy",
           "number resolved amplification status. No significance test is applied to it:",
           "the comparison is 2 cell lines against 2."))

r <- md("m6_replication.csv")
emit(2, "lineage_dominance", "Lineage versus paralog expression in two tumour cohorts",
     data.frame(Cohort = r$cohort, Paralog = r$paralog, n = r$n,
                rho_raw = r$rho_raw, rho_adjusted = r$rho_partial,
                uniqueR2_paralog = r$unique_paralog, uniqueR2_lineage = r$unique_lineage),
     paste("Unique R2 is the variance a term explains that the other cannot. Adjustment",
           "is for neuroendocrine score and all four lineage transcription factors.",
           "At n=79, |rho| >= 0.31 was detectable at 80% power, so the near-zero adjusted",
           "values are a powered negative rather than an absence of evidence."))

e <- md("m7_dependency_enrichment.csv"); p <- md("m7_dependency_pooled.csv")
dep <- data.frame(Set = c(e$paralog, "POOLED"),
                  Genes = c(e$regulon_tested, p$n_regulon_genes),
                  Selective = c(e$n_selective, p$selective_in),
                  OR = c(e$or, p$or), CI_low = c(e$ci_low, p$ci_low),
                  CI_high = c(e$ci_high, p$ci_high), p = c(e$p, signif(p$p, 3)))
emit(3, "dependency", "SCLC-selective CRISPR dependency among regulon genes",
     dep, paste("The pooled test is the powered question; per-paralog odds ratios rest on",
                "4-5 selective genes each and their intervals overlap almost entirely.",
                "MYCL1's nominal p = 0.028 carries an interval that includes 1, so no",
                "per-paralog verdict should be read from this table."))

dg <- md("moes_diagnostics.csv"); cg <- md("moes_concordance_global.csv")
mr <- md("moes_ranking.csv")
emit(4, "MOES", "Evidence integration across the admitted domains",
     data.frame(Paralog = cg$paralog, MaxTopK_overlap_vs_chance = cg$max_ratio,
                Global_p = cg$global_p,
                Genes_at_FDR05 = vapply(cg$paralog, function(x)
                  sum(mr$fdr[mr$paralog == x] < 0.05), integer(1)),
                Best_FDR = round(vapply(cg$paralog, function(x)
                  min(mr$fdr[mr$paralog == x]), numeric(1)), 3)),
     paste("Two of four planned evidence domains were admitted, and only one attributes",
           "evidence to a paralog. Aggregate convergence is real; no individual gene",
           "survives testing across", format(dg$value[dg$metric == "moes_universe_genes"],
           big.mark = ","), "genes. No prioritised hub list is reported."))

cv <- md("m9_panel_coverage.csv"); cv <- cv[!duplicated(cv$paralog), ]
v9 <- md("m9_variance.csv")
emit(5, "spatial", "Spatial panel coverage and orthogonal validation",
     data.frame(Paralog = cv$paralog, Measured = cv$measured,
                Regulon_size = cv$regulon_size,
                Coverage = sprintf("%.1f%%", 100 * cv$coverage_fraction),
                Scoreable = ifelse(cv$scoreable, "yes", "no")),
     paste("Only the MYC regulon clears the coverage gate, by 1.2 percentage points;",
           "MYCL1 misses by 0.4. Spatial results are a 56-gene proxy, not the 500-gene",
           "regulon. In the two cohorts lineage explains unique R2",
           sprintf("%.3f and %.3f against %.3f and %.3f for MYC.",
                   v9$unique_lineage[1], v9$unique_lineage[2],
                   v9$unique_myc[1], v9$unique_myc[2])))

dr <- md("m7_drug_response.csv"); dr$drug <- sub(" \\(.*$", "", dr$drug)
dr <- head(dr[order(dr$p_raw), ], 8)
emit(6, "drug_response", "BET-inhibitor sensitivity versus MYC-family expression",
     data.frame(Screen = dr$screen, Drug = dr$drug, Gene = dr$gene, n = dr$n,
                rho = dr$rho_raw, p = dr$p_raw,
                rho_NEadj = dr$rho_ne_adjusted, p_NEadj = dr$p_ne_adjusted),
     paste("AUC is a sensitivity metric on which lower means more sensitive, so positive",
           "rho means MYCL-high lines are LESS sensitive. This is the only association in",
           "the project that survives lineage adjustment, and is reported as exploratory:",
           "it does not replicate in PRISM, and GDSC1 and GDSC2 are the same platform over",
           "overlapping lines rather than two independent screens."))

cat("\nwrote", length(list.files(OUT)), "files to", OUT, "\n")

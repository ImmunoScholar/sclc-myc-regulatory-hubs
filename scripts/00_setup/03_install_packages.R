# -----------------------------------------------------------------------------
# 03_install_packages.R — install the frozen dependency set into the renv library.
#
# Run from the project root:
#   Rscript scripts/00_setup/03_install_packages.R
#
# Design notes:
#  * Installs group by group and records every failure rather than aborting on
#    the first one, so a single broken package does not hide the others.
#  * Does NOT snapshot. Snapshot is a separate, deliberate step (04_snapshot.R)
#    run only after this reports a clean result.
#  * The single-cell stack (Seurat etc.) is deliberately absent — layer dropped,
#    risk R-03. Do not reinstate without reopening that decision.
# -----------------------------------------------------------------------------

options(warn = 1)

groups <- list(
  infrastructure = c("here", "sessioninfo", "renv"),

  # Genomic ranges & signal quantification over bedGraph.
  # BSgenome.Hsapiens.UCSC.hg19 is DEFERRED: ~700 MB and only needed if a step
  # requires genomic sequence. No planned step does. See decision log D-008.
  genomics = c("GenomicRanges", "rtracklayer", "GenomeInfoDb", "IRanges",
               "S4Vectors", "SummarizedExperiment",
               "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db", "annotatr"),

  expression = c("DESeq2", "limma", "edgeR", "sva"),

  scoring = c("GSVA", "singscore", "AUCell", "msigdbr"),

  network = c("GENIE3", "igraph"),

  integration = c("RobustRankAggreg"),

  data_access = c("GEOquery", "readxl", "data.table", "httr2", "jsonlite"),

  viz = c("ggplot2", "ComplexHeatmap", "circlize", "patchwork", "ggrepel",
          "scales", "ragg", "viridisLite", "colorspace"),

  reporting = c("rmarkdown", "knitr", "testthat")
)

cat("=========================================================\n")
cat("R           :", as.character(getRversion()), "\n")
cat("started     :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("library     :", renv::paths$library(), "\n")
cat("CRAN repo   :", getOption("repos")[["CRAN"]], "\n")
cat("BioC repo   :", getOption("repos")[["BioCsoft"]], "\n")
cat("=========================================================\n\n")

failures <- character(0)

for (g in names(groups)) {
  pkgs <- groups[[g]]
  cat("\n#### GROUP:", g, "(", length(pkgs), "packages )", "####\n")
  for (p in pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      cat(sprintf("  [have]    %-38s %s\n", p, as.character(packageVersion(p))))
      next
    }
    cat(sprintf("  [install] %s ...\n", p))
    ok <- tryCatch({
      renv::install(p, prompt = FALSE)
      TRUE
    }, error = function(e) {
      cat(sprintf("  [FAIL]    %s : %s\n", p, conditionMessage(e)))
      FALSE
    })
    if (!ok || !requireNamespace(p, quietly = TRUE)) {
      failures <<- c(failures, p)
    } else {
      cat(sprintf("  [ok]      %-38s %s\n", p, as.character(packageVersion(p))))
    }
  }
}

cat("\n=========================================================\n")
cat("finished    :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
if (length(failures)) {
  cat("FAILED (", length(failures), "):", paste(failures, collapse = ", "), "\n")
  cat("RESULT: FAIL — do not snapshot. Investigate the failures above.\n")
  cat("Typical cause on Ubuntu: a missing system .so. The error text names it;\n")
  cat("find the providing apt package and install, then re-run this script.\n")
} else {
  cat("RESULT: PASS — all packages present. Next: Rscript scripts/00_setup/04_snapshot.R\n")
}
cat("=========================================================\n")

# -----------------------------------------------------------------------------
# test_environment.R — assertions that must hold before any analysis runs.
#
#   Rscript tests/test_environment.R
#
# Fails loudly rather than letting a broken environment produce plausible-looking
# but wrong results. Run this after renv::restore() on any new machine.
# -----------------------------------------------------------------------------

library(testthat)

cat("Environment checks\n")
cat("R:", as.character(getRversion()), "|", R.version$platform, "\n\n")

test_that("R is new enough", {
  expect_gte(getRversion(), package_version("4.4.0"))
})

test_that("project structure is intact", {
  for (d in c("config", "data/raw", "data/processed", "data/metadata", "docs",
              "figures", "logs", "results", "scripts", "tests")) {
    expect_true(dir.exists(d), info = paste("missing directory:", d))
  }
  for (f in c(".gitignore", "README.md", "LICENSE", "renv.lock",
              "config/params.yml", "docs/project_contract.md",
              "docs/decision_log.md", "docs/research_journal.md")) {
    expect_true(file.exists(f), info = paste("missing file:", f))
  }
})

test_that("core packages load", {
  core <- c("here", "sessioninfo", "GenomicRanges", "rtracklayer", "GenomeInfoDb",
            "IRanges", "S4Vectors", "SummarizedExperiment", "org.Hs.eg.db",
            "TxDb.Hsapiens.UCSC.hg19.knownGene", "annotatr", "DESeq2", "limma",
            "edgeR", "sva", "singscore", "GSVA", "AUCell", "msigdbr", "GENIE3", "igraph",
            "RobustRankAggreg", "GEOquery", "readxl", "data.table", "ggplot2",
            "ComplexHeatmap", "circlize", "patchwork", "ggrepel", "scales",
            "ragg", "viridisLite", "colorspace")
  for (p in core) {
    expect_true(requireNamespace(p, quietly = TRUE), info = paste("cannot load:", p))
  }
})

test_that("config parses and its critical invariants hold", {
  skip_if_not_installed("yaml")
  cfg <- yaml::read_yaml("config/params.yml")

  # Genome build must be hg19 — mixing builds is the project's most likely
  # silent failure mode (risk R-05).
  expect_identical(cfg$project$genome_build, "hg19")

  # Memory discipline is not optional on a 10 GB machine (risk R-13).
  expect_true(isTRUE(cfg$compute$chromosome_wise))

  # MOES must remain weight-free (decision D-003).
  expect_null(cfg$moes$weights)
  expect_identical(cfg$moes$method, "robust_rank_aggregation")
  expect_identical(cfg$moes$stages, 2L)

  # Pre-registered controls must be present before any result exists.
  expect_true("RPS26" %in% cfg$controls$negative_control_genes)

  # The secondary scoring method must be methodologically distinct from the
  # primary. singscore and AUCell are both rank-based, so pairing them would
  # make the sensitivity analysis a restatement rather than a cross-check.
  expect_identical(cfg$tumour_scoring$method, "singscore")
  expect_identical(cfg$tumour_scoring$secondary_method, "gsva")

  # DepMap release must be pinned before use — null is a deliberate tripwire
  # (risk R-08). This test documents the tripwire; it does not enforce a value.
  if (is.null(cfg$depmap$release)) {
    message("NOTE: depmap.release is unpinned. Pin it before running Aim 3.")
  }
})

test_that("genome annotation is the expected build", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  expect_identical(unname(GenomeInfoDb::genome(txdb)[1]), "hg19")
})

test_that("no data files are tracked by git", {
  skip_if(!dir.exists(".git"), "not a git repository")
  tracked <- system2("git", c("ls-files"), stdout = TRUE)
  bad <- grep("^data/(raw|processed)/(?!.*\\.gitkeep$)", tracked, perl = TRUE, value = TRUE)
  expect_length(bad, 0)

  big_ext <- "\\.(bedGraph|bw|bam|fastq\\.gz|rds|RData|h5|tar|zip)$"
  expect_length(grep(big_ext, tracked, ignore.case = TRUE, value = TRUE), 0)
})

cat("\nAll environment checks passed.\n")

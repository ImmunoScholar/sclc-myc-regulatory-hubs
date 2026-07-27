# -----------------------------------------------------------------------------
# verify_config.R — assert config/params.yml parses and its invariants hold.
#
#   Rscript scripts/00_setup/verify_config.R
#
# Run after any config edit and at every milestone boundary, alongside
# audit_decisions.sh. The D-028 sweep found 11 stale or dead parameters, three of
# which were code/config divergences; this catches that class automatically.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(yaml))
cfg <- yaml::read_yaml("config/params.yml")
cat("config parses OK\n\n")

ok <- TRUE
chk <- function(label, cond, detail = "") {
  cat(sprintf("  %-46s %s %s\n", label, if (isTRUE(cond)) "PASS" else "FAIL", detail))
  if (!isTRUE(cond)) ok <<- FALSE
}

cat("=== live parameters ===\n")
cat(sprintf("  %-46s %s\n", "fold_grid",
            paste(unlist(cfg$active_regions$fold_grid), collapse = ", ")))
cat(sprintf("  %-46s %s\n", "primary_fold", cfg$active_regions$primary_fold))
cat(sprintf("  %-46s %s\n", "regions.min_lines_supporting", cfg$regions$min_lines_supporting))
cat(sprintf("  %-46s %s\n", "regions.atac_threshold_multiple", cfg$regions$atac_threshold_multiple))
cat(sprintf("  %-46s %s\n", "regulon_validity.method", cfg$regulon_validity$method))
cat(sprintf("  %-46s %s / %s\n", "spatial min_genes / coverage",
            cfg$spatial$min_genes_measured, cfg$spatial$min_coverage_fraction))

cat("\n=== dead parameters removed (D-028) ===\n")
chk("active_regions.signal_quantile removed",      is.null(cfg$active_regions$signal_quantile))
chk("active_regions.sensitivity_quantiles removed", is.null(cfg$active_regions$sensitivity_quantiles))
chk("motif_validation.n_shuffles removed",         is.null(cfg$m5_gate$motif_validation$n_shuffles))
chk("peak_to_gene.min_abs_correlation removed",    is.null(cfg$peak_to_gene$min_abs_correlation))

cat("\n=== frozen invariants ===\n")
chk("genome build is hg19",            identical(cfg$project$genome_build, "hg19"))
chk("chromosome-wise processing on",   isTRUE(cfg$compute$chromosome_wise))
chk("MOES has NO weights",             is.null(cfg$moes$weights))
chk("MOES is rank aggregation",        identical(cfg$moes$method, "robust_rank_aggregation"))
chk("MOES bootstrap enabled",          isTRUE(cfg$moes$bootstrap$enabled))
chk("negative control RPS26 present",  "RPS26" %in% unlist(cfg$controls$negative_control_genes))
chk("DepMap release still unpinned",   is.null(cfg$depmap$release),
    "(tripwire — pin at M7)")
chk("secondary scorer is GSVA",        identical(cfg$tumour_scoring$secondary_method, "gsva"))

cat("\n=== data-availability honesty ===\n")
chk("YAP1 not in ChIP-based lineage list",
    !("YAP1" %in% unlist(cfg$lineage_confounding$lineage_tfs_with_chip)),
    "(no YAP1 ChIP exists)")
chk("Jung 2017 flagged NOT_YET_CURATED",
    identical(cfg$controls$benchmark_signatures_manual$JUNG2017_MYC_ACTIVITY$status,
              "NOT_YET_CURATED"))
chk("DepMap requires CCLE expression",
    any(grepl("Expression", unlist(cfg$depmap$required_files))))
chk("peak_to_gene labelled as activity proxy",
    identical(cfg$peak_to_gene$activity_proxy, "promoter_h3k27ac"))
chk("regulon validity is leave-one-line-out",
    identical(cfg$regulon_validity$method, "leave_one_line_out"))

cat("\n")
if (!ok) { cat("RESULT: FAIL — config invariants violated.\n"); quit(status = 1) }
cat("RESULT: PASS — config is internally consistent.\n")

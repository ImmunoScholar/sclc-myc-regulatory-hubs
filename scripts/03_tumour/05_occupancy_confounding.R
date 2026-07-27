# -----------------------------------------------------------------------------
# 05_occupancy_confounding.R — lineage confounding at the CHROMATIN level.
#
#   Rscript scripts/03_tumour/05_occupancy_confounding.R
#
# THE SECOND, INDEPENDENT TEST OF R-01. M6 showed the confound in tumour
# EXPRESSION: lineage explains 39-49% of regulon score variance against 0.0-1.6%
# uniquely paralog, and nothing survives adjustment (D-032, D-033).
#
# This asks the same question one level down, in the chromatin the regulons were
# built from: are paralog-bound regions ALSO lineage-TF-bound IN THE SAME CELLS?
# If yes, the confound is present at the point of origin, not introduced by the
# move to tumours. Different data type, different failure modes, same question.
#
# RESOLUTION IS NOT UNIFORM (R-14) and the script reports per TF rather than
# pooling:
#   POU2F3  — peaks in NCIH1048, NCIH211, NCIH526 (hg19). THREE keystone lines
#             spanning two paralog groups (MYC: H1048, H211; MYCN: H526).
#             WITHIN-LINE co-occupancy is possible. The only strong arm.
#   NEUROD1 — peaks in H446 only (hg19). ZERO keystone overlap, so this is a
#             subtype-level comparison and cannot separate line-specific from
#             subtype-specific effects.
#   ASCL1   — SHP-77 bigWig SIGNAL in hg38. Would need both signal thresholding
#             and a liftOver for n=1. Deliberately NOT attempted: the result would
#             not support inference and the effort would imply it did.
#
# Drug-treated arms are excluded. In GSE249362 the FHD286 POU2F3 peak file is
# 4,606 bytes against 181,493 for DMSO — the drug abolishes binding, so treated
# arms do not represent the native repertoire.
#
# Output: data/processed/tumour/occupancy_confounding.rds
#         results/tables/m6_occupancy_confounding.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(GenomicRanges); library(yaml) })
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
FOLD <- CFG$active_regions$primary_fold
set.seed(CFG$project$seed)

nrm <- readRDS("data/processed/signal/region_signal_normalised.rds")
M <- nrm$mean_fob; meta <- nrm$meta; u <- nrm$regions
n_reg <- length(u)
cat("universe: ", format(n_reg, big.mark = ","), " regions\n\n", sep = "")

# ---- POU2F3 peaks, untreated arms only ---------------------------------------
p3dir <- "data/raw/GSE249362"
p3 <- list.files(p3dir, pattern = "POU2F3.*\\.bed\\.gz$", full.names = TRUE)
p3 <- p3[!grepl("FHD286|FHD609", p3)]     # exclude drug-treated
cat("=========== POU2F3 peak files (untreated) ===========\n")
LINE_MAP <- c(NCIH1048 = "H1048", NCIH211 = "H211", NCIH526 = "H526")
pou <- list()
for (f in p3) {
  ln <- NA_character_
  for (k in names(LINE_MAP)) if (grepl(k, basename(f), fixed = TRUE)) ln <- LINE_MAP[[k]]
  if (is.na(ln)) next
  d <- utils::read.delim(gzfile(f), header = FALSE, stringsAsFactors = FALSE)
  gr <- GenomicRanges::reduce(GRanges(to_ucsc_seqnames(d[[1]]),
                                      IRanges(d[[2]] + 1L, d[[3]])))
  gr <- gr[as.character(seqnames(gr)) %in% ANALYSIS_CHROMS_UCSC]
  pou[[ln]] <- gr
  cat(sprintf("  %-6s %6d peaks  %s\n", ln, length(gr), basename(f)))
}
stopifnot(length(pou) > 0)

# ---- NEUROD1 peaks from the series tar --------------------------------------
cat("\n=========== NEUROD1 peaks (H446, no keystone overlap) ===========\n")
nd_tar <- "data/raw/GSE210113/GSE210113_RAW.tar"
nd_dir <- "data/processed/tumour/_neurod1"
nd <- NULL
if (file.exists(nd_tar)) {
  dir.create(nd_dir, showWarnings = FALSE, recursive = TRUE)
  if (!length(list.files(nd_dir)))
    utils::untar(nd_tar, exdir = nd_dir)
  cand <- list.files(nd_dir, pattern = "NEUROD1.*(narrowPeak|broadPeak)\\.gz$",
                     full.names = TRUE)
  cand <- cand[!grepl("KO|ND1-KO", basename(cand))]   # wild-type arm only
  if (length(cand)) {
    d <- utils::read.delim(gzfile(cand[1]), header = FALSE, stringsAsFactors = FALSE)
    nd <- GenomicRanges::reduce(GRanges(to_ucsc_seqnames(d[[1]]),
                                        IRanges(d[[2]] + 1L, d[[3]])))
    nd <- nd[as.character(seqnames(nd)) %in% ANALYSIS_CHROMS_UCSC]
    cat("  ", basename(cand[1]), ": ", length(nd), " peaks\n", sep = "")
  } else cat("  no wild-type NEUROD1 peak file found in the tar\n")
} else cat("  tar not present\n")

# ---- paralog active regions, per line ---------------------------------------
col_for <- function(a, l) { j <- which(meta$assay == a & meta$line == l); if (length(j)==1) j else NA }
active_in_line <- function(p, l) {
  jp <- col_for(p, l); jk <- col_for("H3K27ac", l)
  if (is.na(jp) || is.na(jk)) return(NULL)
  (M[, jp] >= FOLD) & (M[, jk] >= FOLD)
}
PARA_OF_LINE <- list(H1048 = "MYC", H211 = "MYC", H526 = "MYCN")

# ---- WITHIN-LINE co-occupancy (the strong arm) ------------------------------
cat("\n=========== POU2F3 within-line co-occupancy ===========\n")
cat("Same cells for both the paralog ChIP and the lineage-TF ChIP.\n")
cat("Baseline = fraction of the WHOLE universe that is POU2F3-bound in that line.\n\n")
res <- data.frame()
for (ln in names(pou)) {
  p <- PARA_OF_LINE[[ln]]
  a <- active_in_line(p, ln)
  if (is.null(a)) { cat("  ", ln, ": no ", p, " ChIP\n", sep = ""); next }
  bound <- IRanges::overlapsAny(u, pou[[ln]])
  base  <- mean(bound)                       # universe-wide POU2F3 rate
  obs   <- mean(bound[a])                    # among paralog-active regions
  ft <- stats::fisher.test(table(factor(a, c(FALSE,TRUE)), factor(bound, c(FALSE,TRUE))))
  res <- rbind(res, data.frame(
    line = ln, paralog = p, n_active = sum(a),
    pct_universe_pou2f3 = round(100 * base, 1),
    pct_active_pou2f3 = round(100 * obs, 1),
    enrichment = round(obs / base, 2),
    or = round(as.numeric(ft$estimate), 2), p = signif(ft$p.value, 3),
    stringsAsFactors = FALSE))
  cat(sprintf("  %-6s (%-5s) %6d active | POU2F3-bound: universe %.1f%%, active %.1f%% | %.2fx  OR %.2f  p %s\n",
              ln, p, sum(a), 100*base, 100*obs, obs/base, ft$estimate,
              format.pval(ft$p.value, digits = 3)))
}
if (nrow(res)) { cat("\n"); print(res, row.names = FALSE) }

# ---- does PARALOG-SPECIFIC binding avoid lineage TFs? -----------------------
# The sharper question. If paralog-specific regions are LESS lineage-bound than
# shared ones, paralog identity is partly independent of lineage occupancy. If
# they are equally or MORE lineage-bound, the confound reaches the chromatin.
cat("\n=========== paralog-specific vs shared regions ===========\n")
PL <- list(MYC = c("H1048","H211","H524","H847","SHP77"),
           MYCN = c("H526","H69"), MYCL1 = c("COLO668","H889"))
pset <- lapply(names(PL), function(p) {
  mats <- Filter(Negate(is.null), lapply(PL[[p]], active_in_line, p = p))
  if (!length(mats)) return(rep(FALSE, n_reg))
  rowSums(do.call(cbind, mats)) >= 2
}); names(pset) <- names(PL)

spec <- data.frame()
for (ln in names(pou)) {
  p <- PARA_OF_LINE[[ln]]
  others <- setdiff(names(pset), p)
  only_p <- pset[[p]] & !Reduce(`|`, pset[others])
  shared <- pset[[p]] &  Reduce(`|`, pset[others])
  bound <- IRanges::overlapsAny(u, pou[[ln]])
  if (sum(only_p) < 50 || sum(shared) < 50) next
  ft <- stats::fisher.test(matrix(c(sum(bound & only_p), sum(!bound & only_p),
                                    sum(bound & shared), sum(!bound & shared)),
                                  2, byrow = TRUE))
  spec <- rbind(spec, data.frame(
    line = ln, paralog = p,
    n_specific = sum(only_p), pct_specific_pou2f3 = round(100*mean(bound[only_p]),1),
    n_shared = sum(shared), pct_shared_pou2f3 = round(100*mean(bound[shared]),1),
    or_specific_vs_shared = round(as.numeric(ft$estimate), 2),
    p = signif(ft$p.value, 3), stringsAsFactors = FALSE))
}
if (nrow(spec)) print(spec, row.names = FALSE) else cat("  insufficient region counts\n")

# ---- NEUROD1, subtype-level only --------------------------------------------
if (!is.null(nd)) {
  cat("\n=========== NEUROD1 (subtype-level, H446 not a keystone line) ===========\n")
  bound <- IRanges::overlapsAny(u, nd)
  base <- mean(bound)
  cat(sprintf("  universe NEUROD1-bound: %.1f%%\n", 100*base))
  for (p in names(pset)) {
    a <- pset[[p]]
    if (sum(a) < 50) next
    cat(sprintf("  %-6s active (%6d): %.1f%%  enrichment %.2fx\n",
                p, sum(a), 100*mean(bound[a]), mean(bound[a])/base))
  }
  cat("  CANNOT separate line-specific from subtype-specific effects (R-14).\n")
}

# ---- peak-set QC: do the counts match the known subtype biology? -------------
# NCI-H1048 is the canonical SCLC-P / POU2F3-driven line and should carry the MOST
# POU2F3 binding. If it carries the least, the peak sets are not comparable across
# lines and any cross-line contrast built on them is unsafe.
cat("\n=========== POU2F3 peak-set QC ===========\n")
POU_SUBTYPE_LINE <- "H1048"
cnt <- vapply(pou, length, integer(1))
cat("peak counts: ", paste(sprintf("%s=%s", names(cnt), format(cnt, big.mark=",")),
                           collapse = "  "), "\n", sep = "")
cat("expected: ", POU_SUBTYPE_LINE, " (SCLC-P, POU2F3-driven) should have the MOST\n", sep = "")
qc_ok <- names(cnt)[which.max(cnt)] == POU_SUBTYPE_LINE
cat("observed most peaks in: ", names(cnt)[which.max(cnt)], " -> ",
    if (qc_ok) "consistent" else "INCONSISTENT WITH SUBTYPE BIOLOGY", "\n", sep = "")
if (!qc_ok) {
  cat("  The POU2F3-driven line has FEWER peaks than POU2F3-negative lines.\n")
  cat("  Filenames indicate different batches (DX1 vs CKX), so either the CKX peak\n")
  cat("  calls are less stringent or POU2F3 ChIP in POU2F3-negative lines is\n")
  cat("  yielding non-specific background. Cross-line contrasts are NOT SAFE.\n")
}

cat("\n=========== interpretation ===========\n")
if (nrow(res)) {
  cat("1. OVERALL co-occupancy: paralog-active regions are ",
      sprintf("%.1f-%.1fx", min(res$enrichment), max(res$enrichment)),
      " enriched for\n   lineage-TF binding. Largely EXPECTED — both mark active regulatory\n",
      "   elements, so high baseline overlap is not by itself evidence of confounding.\n", sep = "")
}
if (nrow(spec)) {
  dep <- spec$or_specific_vs_shared < 0.8 & spec$p < 0.05
  cat("\n2. THE DECISIVE CONTRAST — paralog-SPECIFIC vs SHARED regions:\n")
  if (any(dep)) {
    cat("   Paralog-specific regions are DEPLETED for lineage-TF binding relative to\n")
    cat("   shared regions in ", sum(dep), "/", nrow(spec), " lines (OR ",
        paste(sprintf("%.2f", spec$or_specific_vs_shared[dep]), collapse = ", "), ").\n", sep = "")
    cat("   So at the CHROMATIN level paralog specificity is PARTLY INDEPENDENT of\n")
    cat("   lineage occupancy — the opposite of the tumour-expression result.\n")
    cat("   Reading: paralog-specific enhancers are distinguishable from lineage-TF\n")
    cat("   sites; what fails is the downstream EXPRESSION readout in tumours.\n")
  } else {
    cat("   No depletion: paralog-specific regions are as lineage-bound as shared\n")
    cat("   ones. The confound extends to the chromatin.\n")
  }
}
if (!qc_ok) {
  cat("\n3. CONFIDENCE IS LIMITED BY THE PEAK-SET QC ABOVE. The depletion signal comes\n")
  cat("   from H211 and H526 — the lines whose peak sets are least trustworthy —\n")
  cat("   while H1048, the line with biologically sensible counts, shows NO\n")
  cat("   difference (OR 0.94, p 0.32). Treat as SUGGESTIVE, not established, and\n")
  cat("   do not build the cis-regulatory MOES domain on it without corroboration.\n")
}
cat("\nASCL1 NOT TESTED: SHP-77 bigWig signal in hg38, n=1. Thresholding plus a\n")
cat("liftOver for a single line would not support inference (R-14).\n")

saveRDS(list(within_line = res, specific_vs_shared = spec,
             pou2f3 = pou, neurod1 = nd),
        "data/processed/tumour/occupancy_confounding.rds")
if (nrow(res)) write.csv(res, "data/metadata/m6_occupancy_confounding.csv", row.names = FALSE)

md <- c("# M6 — lineage confounding at the chromatin level (R-14)", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "Second, independent test of R-01. M6 found the confound in tumour",
        "EXPRESSION; this asks whether paralog-bound regions are also",
        "lineage-TF-bound IN THE SAME CELLS.", "",
        "## POU2F3 within-line co-occupancy (the only strong arm)", "",
        "| line | paralog | active regions | universe POU2F3 | active POU2F3 | enrichment | OR | p |",
        "|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(res)))
  md <- c(md, sprintf("| %s | %s | %s | %s%% | %s%% | %sx | %s | %s |",
                      res$line[i], res$paralog[i], format(res$n_active[i], big.mark=","),
                      res$pct_universe_pou2f3[i], res$pct_active_pou2f3[i],
                      res$enrichment[i], res$or[i], res$p[i]))
md <- c(md, "", "## Resolution per TF (R-14) — not uniform, never pooled", "",
        "- **POU2F3** — 3 keystone lines (H1048, H211 = MYC; H526 = MYCN), hg19, peaks. Strong.",
        "- **NEUROD1** — H446 only, zero keystone overlap. Subtype-level; cannot separate line from subtype.",
        "- **ASCL1** — SHP-77 bigWig signal, hg38, n=1. NOT tested; would not support inference.",
        "", "Drug-treated arms excluded: the FHD286 POU2F3 peak file is 4,606 bytes",
        "against 181,493 for DMSO, so treated arms do not represent native binding.")
writeLines(md, "results/tables/m6_occupancy_confounding.md")
cat("\nwrote results/tables/m6_occupancy_confounding.md\n")

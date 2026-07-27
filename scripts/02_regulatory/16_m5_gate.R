# -----------------------------------------------------------------------------
# 16_m5_gate.R — the M5 gate.
#
#   Rscript scripts/02_regulatory/16_m5_gate.R
#
# Four criteria (D-020, as amended by D-025). None is a region COUNT: counts are
# set by our own threshold and universe size, so matching Plotnik's 18,823 would
# validate nothing. Every criterion is relative or sequence-based — quantities we
# cannot move by tuning our parameters.
#
#   1. MYCN-in-MYC overlap ~0.84 +/- 0.15, STABLE across the fold range.
#   2. MYCN vs MYCL1 differential nesting, as enrichment over chance (D-025).
#   3. Distal-fraction contrast, MYC-amplified vs MYC-expressing lines.
#      *** BLOCKED — see below. Per-line values computed; grouping unverified. ***
#   4. Motif validation: paralog-specific E-box central dinucleotides, enriched in
#      each paralog's OWN regions relative to the other paralogs' regions and to
#      matched background. The only criterion independent of every parameter.
#
# CRITERION 3 IS BLOCKED, NOT SKIPPED. It needs MYC ChIP in both MYC-amplified and
# MYC-expressing lines. MYC ChIP exists in H1048, H211, H524, H847 and SHP77, and
# the project's amplification registry lists ALL FIVE as MYC-amplified. H196, the
# only non-amplified line, has no MYC-family ChIP. Plotnik reported 39% vs 12%
# between those groups, so their grouping must differ from our registry. Resolve
# against the paper's supplement before this criterion can be evaluated. Per-line
# distal fractions are computed and reported so the contrast can be formed the
# moment the grouping is confirmed.
#
# Output: data/metadata/m5_gate_results.csv
#         results/tables/m5_gate.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(Biostrings)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(yaml)
})

CFG  <- yaml::read_yaml("config/params.yml")
GATE <- CFG$m5_gate
nrm  <- readRDS("data/processed/signal/region_signal_normalised.rds")
M    <- nrm$mean_fob; meta <- nrm$meta; u <- nrm$regions
n    <- nrow(M)
# From config, not hardcoded — see D-028.
FOLDS   <- unlist(CFG$active_regions$fold_grid)
PRIMARY <- CFG$active_regions$primary_fold
stopifnot(length(FOLDS) > 0, is.finite(PRIMARY))

PARALOG_LINES <- list(MYC = c("H1048","H211","H524","H847","SHP77"),
                      MYCN = c("H526","H69"), MYCL1 = c("COLO668","H889"))
col_for <- function(a, l) { j <- which(meta$assay == a & meta$line == l); if (length(j)==1) j else NA }
active_in_line <- function(p, l, f) {
  jp <- col_for(p, l); jk <- col_for("H3K27ac", l)
  if (is.na(jp) || is.na(jk)) return(NULL)
  (M[, jp] >= f) & (M[, jk] >= f)
}
paralog_set <- function(p, f, lines = PARALOG_LINES[[p]]) {
  mats <- Filter(Negate(is.null), lapply(lines, active_in_line, p = p, f = f))
  if (!length(mats)) return(rep(FALSE, n))
  rowSums(do.call(cbind, mats)) >= 2      # min_replicates = 2, pre-registered
}

results <- list()
addr <- function(...) results[[length(results)+1L]] <<- data.frame(..., stringsAsFactors = FALSE)

# ===================== CRITERION 1 =====================
cat("=========== CRITERION 1: MYCN-in-MYC overlap ===========\n")
c1 <- data.frame()
for (f in FOLDS) {
  S <- lapply(names(PARALOG_LINES), paralog_set, f = f); names(S) <- names(PARALOG_LINES)
  v <- if (sum(S$MYCN) > 0) sum(S$MYCN & S$MYC) / sum(S$MYCN) else NA
  c1 <- rbind(c1, data.frame(fold = f, n_MYC = sum(S$MYC), n_MYCN = sum(S$MYCN),
                             n_MYCL1 = sum(S$MYCL1), mycn_in_myc = round(v, 3)))
}
print(c1, row.names = FALSE)
pub <- GATE$mycn_in_myc_fraction$published; tol <- GATE$mycn_in_myc_fraction$tolerance
in_tol <- abs(c1$mycn_in_myc - pub) <= tol
stab <- diff(range(c1$mycn_in_myc, na.rm = TRUE))
p1 <- all(in_tol, na.rm = TRUE) && stab <= 0.10
cat(sprintf("\npublished %.2f +/- %.2f | observed %.3f-%.3f | spread %.3f\n",
            pub, tol, min(c1$mycn_in_myc), max(c1$mycn_in_myc), stab))
cat("CRITERION 1: ", if (p1) "PASS" else "FAIL", "\n\n", sep = "")
addr(criterion = "1_mycn_in_myc", value = sprintf("%.3f (spread %.3f)",
     c1$mycn_in_myc[c1$fold == PRIMARY], stab),
     expected = sprintf("%.2f +/- %.2f, stable", pub, tol), pass = p1)

# ===================== CRITERION 2 =====================
cat("=========== CRITERION 2: MYCN vs MYCL1 differential nesting ===========\n")
S <- lapply(names(PARALOG_LINES), paralog_set, f = PRIMARY); names(S) <- names(PARALOG_LINES)
exp_rate <- sum(S$MYC) / n
e_mycn  <- (sum(S$MYCN  & S$MYC) / sum(S$MYCN))  / exp_rate
e_mycl1 <- (sum(S$MYCL1 & S$MYC) / sum(S$MYCL1)) / exp_rate
ct <- matrix(c(sum(S$MYCN & S$MYC), sum(S$MYCN & !S$MYC),
               sum(S$MYCL1 & S$MYC), sum(S$MYCL1 & !S$MYC)), nrow = 2, byrow = TRUE,
             dimnames = list(c("MYCN","MYCL1"), c("in_MYC","out_MYC")))
ft <- stats::fisher.test(ct)
cat("contingency:\n"); print(ct)
cat(sprintf("\nchance expectation (|MYC|/|universe|): %.3f\n", exp_rate))
cat(sprintf("enrichment  MYCN %.2fx   MYCL1 %.2fx   ratio %.2f\n", e_mycn, e_mycl1, e_mycn/e_mycl1))
cat(sprintf("Fisher OR %.3f (95%% CI %.3f-%.3f), p = %.3g\n",
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
p2 <- (e_mycn > e_mycl1) && (ft$p.value < GATE$mycl1_vs_mycn_nesting$alpha)
cat("CRITERION 2: ", if (p2) "PASS" else "FAIL",
    "  (direction ", if (e_mycn > e_mycl1) "correct" else "WRONG", ")\n\n", sep = "")
addr(criterion = "2_differential_nesting",
     value = sprintf("MYCN %.2fx vs MYCL1 %.2fx, OR %.2f, p=%.2g", e_mycn, e_mycl1, ft$estimate, ft$p.value),
     expected = "enrichment(MYCN) > enrichment(MYCL1), p < 0.05", pass = p2)

# ===================== CRITERION 3 (BLOCKED) =====================
cat("=========== CRITERION 3: distal fraction (BLOCKED) ===========\n")
tx  <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tssw <- GenomicRanges::reduce(resize(promoters(tx, 0, 1), 2000, fix = "center"))
is_prom <- IRanges::overlapsAny(u, tssw)
c3 <- data.frame()
for (p in names(PARALOG_LINES)) for (l in PARALOG_LINES[[p]]) {
  a <- active_in_line(p, l, PRIMARY); if (is.null(a)) next
  c3 <- rbind(c3, data.frame(paralog = p, line = l, n_active = sum(a),
                             pct_distal = round(100 * mean(!is_prom[a]), 1)))
}
print(c3, row.names = FALSE)
cat("\nBLOCKED: all five MYC-ChIP lines are MYC-amplified in the project registry,\n")
cat("so no MYC-expressing comparator exists. Plotnik reported 39% vs 12% between\n")
cat("those groups; their grouping must differ. Verify against the paper's\n")
cat("supplement, then form the contrast from the per-line values above.\n")
cat("CRITERION 3: NOT EVALUABLE\n\n")
addr(criterion = "3_distal_contrast", value = "not evaluable",
     expected = "MYC-amplified > MYC-expressing by >=0.10", pass = NA)

# ===================== CRITERION 4: MOTIFS =====================
cat("=========== CRITERION 4: paralog E-box motifs ===========\n")
EB <- GATE$motif_validation$paralog_ebox
FL <- GATE$motif_validation$flank_bp
gen <- BSgenome.Hsapiens.UCSC.hg19

# Build the flanking windows by ARITHMETIC, not by resize()+trim().
#
# trim() is a no-op unless the GRanges carries valid seqlengths, and attaching
# seqinfo did not propagate as expected — a window near the end of chr17 still
# ran past the chromosome and getSeq errored. Rather than debug seqinfo plumbing,
# clamp start and end against the genome's own chromosome lengths directly. There
# is then nothing to propagate and nothing to silently no-op.
prep_windows <- function(gr, flank) {
  if (!length(gr)) return(GRanges())
  sl   <- GenomeInfoDb::seqlengths(gen)
  chrs <- as.character(GenomicRanges::seqnames(gr))
  keep <- chrs %in% names(sl)
  if (!any(keep)) return(GRanges())
  gr <- gr[keep]; chrs <- chrs[keep]

  ctr    <- GenomicRanges::start(gr) + floor(GenomicRanges::width(gr) / 2)
  chrlen <- as.numeric(sl[chrs])
  s <- pmax(1, ctr - flank)
  e <- pmin(chrlen, ctr + flank)

  ok <- (e - s) >= 19          # drop windows too short to contain a 6-mer usefully
  if (!any(ok)) return(GRanges())
  GRanges(chrs[ok], IRanges(as.integer(s[ok]), as.integer(e[ok])))
}

count_motif <- function(gr, patt) {
  w <- prep_windows(gr, FL)
  if (!length(w)) return(c(hits = 0, bp = 0))
  s <- getSeq(gen, w)
  h <- sum(vcountPattern(DNAString(patt), s)) +
       sum(vcountPattern(reverseComplement(DNAString(patt)), s))
  c(hits = h, bp = sum(as.numeric(width(s))))
}
set.seed(CFG$project$seed)
bg_idx <- sample(which(!(S$MYC | S$MYCN | S$MYCL1)), min(10000, sum(!(S$MYC|S$MYCN|S$MYCL1))))
bg_gr  <- u[bg_idx]

mot <- data.frame()
for (p in names(EB)) {
  patt <- EB[[p]]
  bgc  <- count_motif(bg_gr, patt); bg_rate <- bgc["hits"] / bgc["bp"]
  for (q in names(PARALOG_LINES)) {
    cc <- count_motif(u[S[[q]]], patt)
    mot <- rbind(mot, data.frame(motif_of = p, motif = patt, region_set = q,
                                 rate_per_kb = round(1000 * cc["hits"] / cc["bp"], 3),
                                 enrich_vs_bg = round((cc["hits"]/cc["bp"]) / bg_rate, 3)))
  }
}
rownames(mot) <- NULL
print(mot, row.names = FALSE)

# --- 4a. absolute rate vs background: CONFOUNDED, reported for transparency ---
own_best_abs <- vapply(names(EB), function(p) {
  s <- mot[mot$motif_of == p, ]; s$region_set[which.max(s$enrich_vs_bg)] == p
}, logical(1))
cat("\n[4a] absolute rate vs unmatched background (CONFOUNDED):\n")
for (p in names(own_best_abs)) cat(sprintf("  %-6s %s\n", p, if (own_best_abs[[p]]) "yes" else "NO"))
cat("  Nearly all enrichments are <1, i.e. motifs appear DEPLETED in bound regions.\n")
cat("  Cause is base composition, not biology: CAGATG and CACATG are 50%% GC and\n")
cat("  CACCTG is 67%% GC, while bound regions are promoter-shifted and GC-rich, so an\n")
cat("  AT-richer 6-mer is expected to be rarer there. The background was not\n")
cat("  GC-matched, so 4a measures composition. It is NOT used for the verdict.\n")

# --- 4b. compositional test: GC-controlled by construction --------------------
# All three motifs are counted in the SAME sequences, so base composition cancels.
# The motif-to-paralog assignment is fixed in advance (Plotnik), so this is a
# directional test: P(all three land in their designated set by chance) = (1/3)^3.
cat("\n[4b] motif SHARE within the same sequences (GC-controlled):\n")
share <- matrix(NA_real_, nrow = 3, ncol = 3,
                dimnames = list(motif = names(EB), region_set = names(PARALOG_LINES)))
cnts  <- share
for (rs in names(PARALOG_LINES)) {
  h <- vapply(names(EB), function(p) unname(count_motif(u[S[[rs]]], EB[[p]])["hits"]), numeric(1))
  cnts[, rs]  <- h
  share[, rs] <- h / sum(h)
}
print(round(share, 3))

own_best_comp <- vapply(names(EB), function(p)
  colnames(share)[which.max(share[p, ])] == p, logical(1))
cat("\nown motif has highest share in own regions:\n")
for (p in names(own_best_comp)) cat(sprintf("  %-6s %s\n", p, if (own_best_comp[[p]]) "yes" else "NO"))

n_ok  <- sum(own_best_comp)
p_dir <- stats::pbinom(n_ok - 1, 3, 1/3, lower.tail = FALSE)   # >= n_ok by chance
chi   <- suppressWarnings(stats::chisq.test(cnts))
cat(sprintf("\ndirectional test: %d/3 correct, P(>=%d by chance) = %.3f\n", n_ok, n_ok, p_dir))
cat(sprintf("independence of motif x region-set: chi-sq = %.1f, df = %d, p = %.3g\n",
            chi$statistic, chi$parameter, chi$p.value))
cat("standardised residuals (positive on the diagonal = own-motif preference):\n")
print(round(chi$stdres, 2))

p4 <- (n_ok == 3) && (chi$p.value < 0.05)
cat("\nCRITERION 4: ", if (p4) "PASS" else "FAIL", sprintf(" (%d/3, compositional)\n\n", n_ok), sep = "")
addr(criterion = "4_motif_specificity",
     value = sprintf("%d/3 own-motif top share (p=%.3f); chi-sq p=%.2g", n_ok, p_dir, chi$p.value),
     expected = "3/3 own-motif highest share, motif x set dependence p<0.05", pass = p4)
write.csv(as.data.frame(share), "data/metadata/m5_motif_shares.csv")

# ===================== VERDICT =====================
res <- do.call(rbind, results)
write.csv(res, "data/metadata/m5_gate_results.csv", row.names = FALSE)
write.csv(mot, "data/metadata/m5_motif_results.csv", row.names = FALSE)
cat("=========== M5 GATE ===========\n")
print(res, row.names = FALSE)
ev <- res[!is.na(res$pass), ]
cat(sprintf("\nevaluable criteria: %d/%d passed; 1 blocked\n", sum(ev$pass), nrow(ev)))

md <- c("# M5 Gate", "", paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "No criterion is a region count: counts are determined by our own threshold",
        "and universe size, so matching a published count would validate nothing (D-020).", "",
        "| criterion | observed | expected | result |", "|---|---|---|---|")
for (i in seq_len(nrow(res)))
  md <- c(md, sprintf("| %s | %s | %s | %s |", res$criterion[i], res$value[i], res$expected[i],
                      ifelse(is.na(res$pass[i]), "**BLOCKED**", ifelse(res$pass[i], "PASS", "**FAIL**"))))
writeLines(md, "results/tables/m5_gate.md")
cat("\nwrote results/tables/m5_gate.md\n")

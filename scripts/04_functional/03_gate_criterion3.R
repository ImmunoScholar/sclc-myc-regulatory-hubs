# -----------------------------------------------------------------------------
# 03_gate_criterion3.R — M5 gate criterion 3, finally evaluable.
#
#   Rscript scripts/04_functional/03_gate_criterion3.R
#
# Blocked since M5 because it needs a MYC-amplified vs MYC-EXPRESSING split that
# neither the paper, GEO, nor our coverage proxy could supply (D-026). DepMap
# copy number settled it:
#
#   MYC-AMPLIFIED  : H524 (log2 6.73), H211 (2.35)
#   MYC-EXPRESSING : H1048 (0.98), SHP77 (1.30)
#   UNKNOWN        : H847 — absent from DepMap 26Q1, excluded from BOTH groups
#                    because absence is not evidence
#
# THE TEST (config m5_gate$distal_fraction_contrast). Plotnik reported 39% of
# MYC-bound active regions distal in MYC-amplified lines against 12% in
# MYC-expressing lines — a 27-point contrast. Required here: correct DIRECTION
# (amplified > expressing) and a difference of at least 0.10.
#
# WHAT A FAILURE WOULD MEAN. Direction without magnitude would say the effect
# exists but is far weaker in our hands, which is a real discordance worth
# reporting — the gate was rebuilt (D-020) precisely so criteria could fail on
# evidence rather than on our own free parameters.
#
# POWER IS THE OBVIOUS LIMIT: 2 lines vs 2 lines. Per-line values are reported
# individually and no significance test is applied to n=2 vs n=2, because none
# would be honest.
#
# Output: data/metadata/m5_criterion3.csv
#         results/tables/m5_criterion3.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(yaml)
})
source("R/genome_utils.R")

CFG  <- yaml::read_yaml("config/params.yml")
C3   <- CFG$m5_gate$distal_fraction_contrast
FOLD <- CFG$active_regions$primary_fold

AMP  <- c("H524", "H211")
EXPR <- c("H1048", "SHP77")
UNKNOWN <- c("H847")
cat("MYC-amplified  : ", paste(AMP, collapse = ", "), "\n", sep = "")
cat("MYC-expressing : ", paste(EXPR, collapse = ", "), "\n", sep = "")
cat("excluded (status unknown): ", paste(UNKNOWN, collapse = ", "), "\n\n", sep = "")

nrm <- readRDS("data/processed/signal/region_signal_normalised.rds")
M <- nrm$mean_fob; meta <- nrm$meta; u <- nrm$regions

tx  <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tss <- promoters(tx, upstream = 0, downstream = 1)
tss <- tss[as.character(seqnames(tss)) %in% ANALYSIS_CHROMS_UCSC]
# Plotnik used a 1 kb cutoff; matched here rather than reusing our 2.5 kb SE window.
tss_win <- GenomicRanges::reduce(resize(tss, 2000, fix = "center"))
is_prom <- IRanges::overlapsAny(u, tss_win)
cat("promoter definition: TSS +/-1 kb (matching Plotnik)\n")
cat("universe promoter-proximal: ", sprintf("%.1f%%", 100 * mean(is_prom)), "\n\n", sep = "")

col_for <- function(a, l) { j <- which(meta$assay == a & meta$line == l); if (length(j)==1) j else NA }
active_in_line <- function(l) {
  jp <- col_for("MYC", l); jk <- col_for("H3K27ac", l)
  if (is.na(jp) || is.na(jk)) return(NULL)
  (M[, jp] >= FOLD) & (M[, jk] >= FOLD)
}

cat("=========== per-line distal fraction of MYC-active regions ===========\n")
per <- data.frame()
for (grp in list(list(g = "amplified", ls = AMP), list(g = "expressing", ls = EXPR))) {
  for (l in grp$ls) {
    a <- active_in_line(l)
    if (is.null(a)) { cat("  ", l, ": no MYC ChIP\n", sep = ""); next }
    per <- rbind(per, data.frame(group = grp$g, line = l, n_active = sum(a),
                                 pct_distal = round(100 * mean(!is_prom[a]), 1),
                                 stringsAsFactors = FALSE))
  }
}
print(per, row.names = FALSE)

amp_mean  <- mean(per$pct_distal[per$group == "amplified"]) / 100
expr_mean <- mean(per$pct_distal[per$group == "expressing"]) / 100
diff <- amp_mean - expr_mean

cat("\n=========== contrast ===========\n")
cat(sprintf("  amplified  mean distal: %.1f%%  (%s)\n", 100 * amp_mean,
            paste(per$pct_distal[per$group=="amplified"], collapse = ", ")))
cat(sprintf("  expressing mean distal: %.1f%%  (%s)\n", 100 * expr_mean,
            paste(per$pct_distal[per$group=="expressing"], collapse = ", ")))
cat(sprintf("  difference            : %+.1f points\n", 100 * diff))
cat(sprintf("\n  published (Plotnik)   : %.0f%% vs %.0f%%  = %+.0f points\n",
            100 * C3$myc_amplified_published, 100 * C3$myc_expressing_published,
            100 * (C3$myc_amplified_published - C3$myc_expressing_published)))

dir_ok <- diff > 0
mag_ok <- diff >= C3$min_difference
pass <- dir_ok && mag_ok

cat("\n=========== CRITERION 3 ===========\n")
cat("  direction (amplified > expressing): ", if (dir_ok) "PASS" else "FAIL", "\n", sep = "")
cat("  magnitude (>= ", C3$min_difference, "): ", if (mag_ok) "PASS" else "FAIL",
    sprintf("  (observed %.3f)", diff), "\n", sep = "")
cat("  CRITERION 3: ", if (pass) "PASS" else "FAIL", "\n", sep = "")

if (dir_ok && !mag_ok) {
  cat("\nDirection reproduces, magnitude does not. The effect is present but far\n")
  cat("weaker than published. Plausible contributors, none of which can be\n")
  cat("separated at n=2 vs n=2:\n")
  cat("  * our regions come from a shared ATAC-defined universe that is ~84%\n")
  cat("    distal by construction, compressing the achievable range\n")
  cat("  * Plotnik called peaks de novo per line with no shared grid\n")
  cat("  * H847 is excluded (status unknown), so the expressing group is n=2\n")
  cat("Report as a partial reproduction: direction yes, magnitude no.\n")
}
cat("\nPOWER: 2 lines vs 2 lines. No significance test is applied because none\n")
cat("would be honest at this n. Per-line values are reported above.\n")

out <- data.frame(criterion = "3_distal_contrast",
                  amplified_pct = round(100*amp_mean,1),
                  expressing_pct = round(100*expr_mean,1),
                  difference = round(diff,3),
                  published_difference = C3$myc_amplified_published - C3$myc_expressing_published,
                  direction_pass = dir_ok, magnitude_pass = mag_ok, pass = pass,
                  stringsAsFactors = FALSE)
write.csv(rbind(per[,c("group","line","n_active","pct_distal")],
                data.frame(group="SUMMARY", line="", n_active=NA,
                           pct_distal=round(100*diff,1))),
          "data/metadata/m5_criterion3.csv", row.names = FALSE)
write.csv(out, "data/metadata/m5_criterion3_summary.csv", row.names = FALSE)

# --- close the loop on the canonical gate table --------------------------------
# This script resolves a criterion that 16_m5_gate.R could only record as "not
# evaluable". Without writing the verdict back, m5_gate_results.csv — the table
# the report and the audit both read — keeps saying the criterion was never
# evaluated, while the decision log says it FAILED. Two sources of truth, and the
# stale one is the canonical one.
gate_path <- "data/metadata/m5_gate_results.csv"
gate <- read.csv(gate_path, stringsAsFactors = FALSE)
i <- which(gate$criterion == "3_distal_contrast")
if (length(i) != 1L)
  stop("expected exactly one 3_distal_contrast row in ", gate_path, "; found ", length(i))

gate$value[i] <- sprintf(
  "amplified %.1f%% vs expressing %.1f%%, %+.1f pts (published %+.0f); direction %s, magnitude %s",
  100*amp_mean, 100*expr_mean, 100*diff,
  100*(C3$myc_amplified_published - C3$myc_expressing_published),
  if (dir_ok) "reproduces" else "fails", if (mag_ok) "reproduces" else "fails")
gate$pass[i] <- pass
write.csv(gate, gate_path, row.names = FALSE)

cat("\nupdated ", gate_path, " — criterion 3 now recorded as ",
    if (pass) "PASS" else "FAIL", " (was \"not evaluable\")\n", sep = "")
n_pass <- sum(gate$pass %in% TRUE); n_fail <- sum(gate$pass %in% FALSE)
n_open <- sum(is.na(gate$pass))
cat("M5 gate now stands at ", n_pass, " PASS / ", n_fail, " FAIL",
    if (n_open) paste0(" / ", n_open, " open") else "", "\n", sep = "")

md <- c("# M5 gate criterion 3 — distal-fraction contrast", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "Blocked since M5; unblocked by DepMap Public 26Q1 copy number.", "",
        "| group | line | MYC-active regions | % distal |", "|---|---|---|---|")
for (i in seq_len(nrow(per)))
  md <- c(md, sprintf("| %s | %s | %s | %s |", per$group[i], per$line[i],
                      format(per$n_active[i], big.mark=","), per$pct_distal[i]))
md <- c(md, "",
        sprintf("**Observed contrast: %+.1f points** (%.1f%% vs %.1f%%).",
                100*diff, 100*amp_mean, 100*expr_mean),
        sprintf("**Published (Plotnik): %+.0f points** (39%% vs 12%%).",
                100*(C3$myc_amplified_published - C3$myc_expressing_published)), "",
        paste0("**Result: ", if (pass) "PASS" else "FAIL",
               "** — direction ", if (dir_ok) "reproduces" else "does not reproduce",
               ", magnitude ", if (mag_ok) "reproduces" else "does not reproduce", "."), "",
        "H847 excluded: absent from DepMap 26Q1, so its amplification status is",
        "unknown and it belongs to neither group. Power is 2 lines vs 2 lines and",
        "no significance test is applied.")
writeLines(md, "results/tables/m5_criterion3.md")
cat("\nwrote results/tables/m5_criterion3.md\n")

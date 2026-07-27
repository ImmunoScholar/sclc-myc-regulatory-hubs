# -----------------------------------------------------------------------------
# 14_normalise_signal.R — make tracks comparable without inventing the answer.
#
#   Rscript scripts/02_regulatory/14_normalise_signal.R
#
# THE PROBLEM THIS HAS TO SOLVE WITHOUT CHEATING
#
# Paralog groups are perfectly confounded with cell lines:
#     MYC    H1048 H211 H524 H847 SHP77      MYCN   H526 H69
#     MYCL1  COLO668 H889
# So a difference in ChIP efficiency or depth between the MYCN pair and the MYCL1
# pair is indistinguishable from a difference in paralog biology. The M5 gate
# turns on exactly such comparisons (MYCN-in-MYC ~0.84; MYCL1 NOT a MYC subset),
# which makes the normalisation choice load-bearing for the headline result.
#
# WHY THE OBVIOUS METHODS ARE WRONG HERE
#
#  * Quantile normalisation forces identical signal distributions across samples.
#    If MYCN genuinely binds fewer regions than MYC, quantile normalisation
#    ERASES that — it would manufacture agreement with the published 0.84 whether
#    or not the data support it. Unusable for the primary comparison.
#  * Median-ratio (DESeq-style) assumes most features are non-differential. The
#    paralogs are expected to differ across a large share of regions, so the
#    assumption fails by construction.
#  * No normalisation lets sequencing depth masquerade as binding breadth.
#
# WHAT IS USED INSTEAD: fold-over-own-background. Each track is divided by its own
# background level, estimated as a low quantile of its signal across the universe
# (a given paralog binds a minority of accessible regions, so the bulk of the
# distribution IS background). This corrects technical depth while leaving the
# SHAPE of each distribution free — a track with genuinely more strong binding
# keeps more high-fold regions. It cannot force two paralogs to agree.
#
# All methods are computed and compared so the choice is auditable, but only
# fold-over-background is used downstream.
#
# Output: data/processed/signal/region_signal_normalised.rds
#         data/metadata/normalisation_comparison.csv
#         results/tables/m5_normalisation.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(GenomicRanges) })

sig <- readRDS("data/processed/signal/region_signal_raw.rds")
M   <- sig$mean                      # regions x samples, per-bp mean signal
meta <- sig$meta
n_reg <- nrow(M)
cat("matrix: ", format(n_reg, big.mark = ","), " regions x ", ncol(M), " samples\n\n", sep = "")

BG_Q <- 0.50    # background estimate: median across the universe

# --- method 1: raw ------------------------------------------------------------
raw <- M

# --- method 2: fold over own background (ADOPTED) -----------------------------
bg <- apply(M, 2, function(x) {
  b <- stats::quantile(x, BG_Q, na.rm = TRUE)
  if (b <= 0) b <- mean(x[x > 0], na.rm = TRUE) / 4   # fallback for sparse tracks
  b
})
fob <- sweep(M, 2, bg, "/")

# --- method 3: quantile, WITHIN assay (computed for comparison only) ----------
qnorm_within <- M
for (a in unique(meta$assay)) {
  j <- which(meta$assay == a)
  if (length(j) < 2) next
  sub <- M[, j, drop = FALSE]
  ranks <- apply(sub, 2, rank, ties.method = "average")
  target <- rowMeans(apply(sub, 2, sort))
  qnorm_within[, j] <- apply(ranks, 2, function(r) target[round(r)])
}

# --- method 4: scale to common median -----------------------------------------
med <- apply(M, 2, stats::median, na.rm = TRUE)
med[med <= 0] <- 1
medscale <- sweep(M, 2, med / mean(med), "/")

methods <- list(raw = raw, fold_over_background = fob,
                quantile_within_assay = qnorm_within, median_scaled = medscale)

# --- diagnostic: what does each method do to the paralog comparison? ----------
# Proxy for the gate: at a fixed fold cutoff, how many regions does each paralog
# "occupy", and what is the MYCN-in-MYC overlap? If a method drives these to
# agreement regardless of the data, it is deciding the result rather than
# enabling the measurement.
FOLD <- 2
paralogs <- c("MYC", "MYCN", "MYCL1")

occupied <- function(mat, assay) {
  j <- which(meta$assay == assay)
  if (!length(j)) return(rep(FALSE, n_reg))
  # a region counts as occupied if it exceeds the cutoff in EVERY line assayed
  # for that paralog (min_replicates = 2 in config)
  rowSums(mat[, j, drop = FALSE] >= FOLD) == length(j)
}

rows <- list()
for (nm in names(methods)) {
  mat <- methods[[nm]]
  # put every method on a comparable footing: express as fold over own median
  mm <- sweep(mat, 2, apply(mat, 2, function(x) max(stats::median(x), 1e-9)), "/")
  occ <- lapply(paralogs, function(p) occupied(mm, p))
  names(occ) <- paralogs
  n <- vapply(occ, sum, integer(1))
  mycn_in_myc  <- if (n[["MYCN"]]  > 0) sum(occ$MYCN  & occ$MYC) / n[["MYCN"]]  else NA
  mycl1_in_myc <- if (n[["MYCL1"]] > 0) sum(occ$MYCL1 & occ$MYC) / n[["MYCL1"]] else NA
  rows[[length(rows) + 1L]] <- data.frame(
    method = nm, n_MYC = n[["MYC"]], n_MYCN = n[["MYCN"]], n_MYCL1 = n[["MYCL1"]],
    mycn_in_myc = round(mycn_in_myc, 3), mycl1_in_myc = round(mycl1_in_myc, 3),
    stringsAsFactors = FALSE)
}
cmp <- do.call(rbind, rows)

cat("=========== per-sample background estimates ===========\n")
bgt <- data.frame(sample = meta$sample, assay = meta$assay, line = meta$line,
                  background = round(bg, 4), stringsAsFactors = FALSE)
for (a in c("MYC","MYCN","MYCL1","H3K27ac","ATAC")) {
  s <- bgt[bgt$assay == a, ]
  if (!nrow(s)) next
  cat(sprintf("  %-8s %s\n", a,
              paste(sprintf("%s=%.3f", s$line, s$background), collapse = "  ")))
}

cat("\n=========== effect of each method on the paralog comparison ===========\n")
print(cmp, row.names = FALSE)
cat("\nPublished reference (Plotnik 2024): MYCN-in-MYC ~= 0.84;\n")
cat("MYCL1 'largely non-overlapping' with MYC.\n")
cat("\nREAD THIS CAREFULLY: if quantile_within_assay lands close to 0.84 while the\n")
cat("others do not, that is NOT corroboration. Quantile normalisation forces the\n")
cat("distributions to match, so agreement is an artefact of the method. Only\n")
cat("fold_over_background leaves the paralogs free to disagree.\n")

write.csv(cmp, "data/metadata/normalisation_comparison.csv", row.names = FALSE)

out <- list(mean_raw = raw, mean_fob = fob, background = bg,
            meta = meta, regions = sig$regions, width = sig$width,
            adopted = "fold_over_background",
            rationale = paste(
              "Corrects technical depth without constraining distribution shape.",
              "Quantile normalisation was rejected: it forces identical",
              "distributions and would manufacture agreement with the published",
              "MYCN-in-MYC value regardless of the data. Median-ratio was rejected:",
              "it assumes most features are non-differential, which fails when the",
              "paralogs are expected to differ across many regions."))
saveRDS(out, "data/processed/signal/region_signal_normalised.rds")

md <- c("# M5 — signal normalisation", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "**Adopted: fold-over-own-background.**", "",
        "Paralog groups are perfectly confounded with cell lines (MYCN is only",
        "H526+H69; MYCL1 only COLO668+H889), so the normalisation choice is",
        "load-bearing for the headline comparison and is recorded explicitly.", "",
        "| method | MYC regions | MYCN | MYCL1 | MYCN-in-MYC | MYCL1-in-MYC |",
        "|---|---|---|---|---|---|")
for (i in seq_len(nrow(cmp)))
  md <- c(md, sprintf("| %s%s | %s | %s | %s | **%s** | %s |",
                      cmp$method[i],
                      ifelse(cmp$method[i] == "fold_over_background", " **(adopted)**", ""),
                      format(cmp$n_MYC[i], big.mark = ","),
                      format(cmp$n_MYCN[i], big.mark = ","),
                      format(cmp$n_MYCL1[i], big.mark = ","),
                      cmp$mycn_in_myc[i], cmp$mycl1_in_myc[i]))
md <- c(md, "",
        "Quantile normalisation is shown for comparison and is **not** used.",
        "It forces identical signal distributions, so if MYCN genuinely binds",
        "fewer regions than MYC it erases that difference — it would reproduce the",
        "published 0.84 whether or not the data support it. Agreement produced that",
        "way is an artefact of the method, not evidence.")
writeLines(md, "results/tables/m5_normalisation.md")

cat("\nwrote data/processed/signal/region_signal_normalised.rds\n")
cat("wrote results/tables/m5_normalisation.md\n")

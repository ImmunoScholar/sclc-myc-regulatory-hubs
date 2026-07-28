# -----------------------------------------------------------------------------
# fig04_moes_null.R — the evidence-integration framework returns nothing, and why.
#
#   Rscript scripts/07_figures/fig04_moes_null.R
#
# MOES was the methodological centrepiece specified at M1: four domains, two
# stages, paralog-specific prioritisation. It runs on two domains and produces no
# gene at FDR < 0.05. This figure exists so that null is legible as a RESULT with
# a mechanism and a power statement, rather than as a missing table.
#
#   A  the two admitted domains are independent — the mechanism
#   B  no FDR curve approaches the threshold
#   C  the method detects convergence when it exists — the power statement
#   D  between-paralog agreement is inflated by the shared functional layer
#
# Panel C is the panel that makes A and B mean something. Without it this figure
# would show an absence and claim it was a finding.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg); library(colorspace)
})
source("R/theme_project.R")

rank_tab <- read.csv("data/metadata/moes_ranking.csv", stringsAsFactors = FALSE)
pc       <- read.csv("data/metadata/moes_positive_control.csv", stringsAsFactors = FALSE)
FDR_T    <- yaml::read_yaml("config/params.yml")$moes$fdr_threshold

# All inputs are committed CSVs. This figure previously read
# data/processed/integration/moes_results.rds, which is correctly gitignored, so
# it could not be rebuilt from a fresh clone — the M11 clean-clone check caught
# it. MR is reconstructed from the diagnostics and attribution tables the
# analysis now writes.
diagn    <- read.csv("data/metadata/moes_diagnostics.csv", stringsAsFactors = FALSE)
dgv      <- function(k) diagn$value[diagn$metric == k]
MR <- list(
  layer_rho   = setNames(c(dgv("layer_rho_MYC"), dgv("layer_rho_MYCN"), dgv("layer_rho_MYCL1")),
                         c("MYC", "MYCN", "MYCL1")),
  attribution = read.csv("data/metadata/moes_attribution.csv", stringsAsFactors = FALSE),
  n_perm      = dgv("n_permutations"),
  universe    = seq_len(dgv("moes_universe_genes")))

cat("accessibility checks\n")
check_palette(PAL_PARALOG, label = "paralog palette")
check_text_size()
cat("\n")

rank_tab$paralog <- factor(rank_tab$paralog, levels = c("MYC", "MYCN", "MYCL1"))
N <- length(MR$universe)

# Guard: this figure's title asserts a null. If a rerun ever produces a
# discovery, the figure must not keep claiming otherwise.
n_sig <- sum(rank_tab$fdr < FDR_T)
if (n_sig > 0)
  stop("MOES now yields ", n_sig, " genes at FDR < ", FDR_T,
       "; this figure asserts a null and must be rewritten before it is drawn.")
cat("confirmed: 0 genes at FDR < ", FDR_T, " across ", format(N, big.mark = ","),
    " genes x 3 paralogs\n\n", sep = "")

# ---- A. the two domains are independent --------------------------------------
# Restricted to genes that actually CARRY cis evidence for that paralog. Genes
# with no peak-to-gene link share one tied rank, and plotting them produces a
# dense vertical band that reads as structure while being nothing but the tie
# block. The tie mass is large and differs sharply by paralog — MYCN has cis
# evidence for only ~3.4k of the 10.4k universe genes — so it is reported in the
# subtitle rather than drawn as if it were data.
# Per-paralog cis coverage from the committed summary rather than from
# moes_evidence.rds. n_in_universe is defined as the count of genes with cis
# evidence inside the MOES universe, which is exactly what this figure reports.
# The .rds read was the second one in this script and was missed on the first
# pass — the clean-clone gate caught it.
evs <- read.csv("data/metadata/moes_evidence_summary.csv", stringsAsFactors = FALSE)
n_scored <- setNames(evs$n_in_universe, evs$paralog)

# Concordance-at-K with its permutation envelope, both computed in 02_moes.R.
# This replaced a rank-rank scatter, which is the wrong instrument: the cis scores
# are heavily TIED (a gene with a single promoter link scores exactly 1), so the
# scatter bands into vertical stripes that show the tie structure rather than the
# association. Only the K range where >= 5 genes are expected to overlap by chance
# is drawn — below that the ratio is integer noise.
conc <- read.csv("data/metadata/moes_concordance.csv", stringsAsFactors = FALSE)
cglob <- read.csv("data/metadata/moes_concordance_global.csv", stringsAsFactors = FALSE)
conc <- conc[conc$testable, ]
conc$paralog <- factor(conc$paralog, levels = levels(rank_tab$paralog))
cglob$paralog <- factor(cglob$paralog, levels = levels(rank_tab$paralog))

pA <- ggplot(conc, aes(K, ratio, colour = paralog)) +
  geom_ribbon(aes(ymin = null_lo, ymax = null_hi), fill = "#DDDDDD",
              colour = NA, alpha = 0.55, show.legend = FALSE) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.4, colour = "#666666") +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = PAL_PARALOG, name = NULL) +
  scale_x_log10(breaks = c(300, 1000, 3000, 10000),
                labels = c("300", "1,000", "3,000", "10,000")) +
  # Limits must contain the WHOLE null envelope (its lower bound reaches 0.31 at
  # small K). Clipping a null band makes it read narrower than it is, and an
  # observed curve then looks like it clears a threshold that was cropped.
  scale_y_continuous(limits = c(0.25, 2.4)) +
  labs(title = "A · Weak convergence, above the permutation null",
       subtitle = sprintf("top-K overlap vs chance; grey = 95%% null envelope\nglobal p: MYC %.3f · MYCN %.3f · MYCL1 %.3f",
                          cglob$global_p[cglob$paralog == "MYC"],
                          cglob$global_p[cglob$paralog == "MYCN"],
                          cglob$global_p[cglob$paralog == "MYCL1"]),
       x = "K (top-ranked genes from each domain)", y = "observed / expected overlap") +
  theme_project() +
  theme(legend.position = "bottom", plot.margin = margin(4, 14, 4, 4))

# ---- B. no FDR curve reaches the threshold -----------------------------------
fc <- do.call(rbind, lapply(levels(rank_tab$paralog), function(p) {
  f <- sort(rank_tab$fdr[rank_tab$paralog == p])
  data.frame(paralog = p, idx = seq_along(f), fdr = f)
}))
fc$paralog <- factor(fc$paralog, levels = levels(rank_tab$paralog))
fc <- fc[fc$idx <= 2000, ]

pB <- ggplot(fc, aes(idx, fdr, colour = paralog)) +
  geom_hline(yintercept = FDR_T, linetype = "22", linewidth = 0.45, colour = "#A63603") +
  geom_line(linewidth = 0.7) +
  annotate("text", x = 2000, y = FDR_T, label = sprintf("FDR = %.2f", FDR_T),
           hjust = 1, vjust = -0.6, size = 2.1, colour = "#A63603") +
  scale_colour_manual(values = PAL_PARALOG, name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "B · No gene comes close to significance",
       subtitle = sprintf("empirical FDR by MOES rank; best achieved is %.3f",
                          min(rank_tab$fdr)),
       x = "MOES rank (best 2,000 shown)", y = "empirical FDR") +
  theme_project() +
  theme(legend.position = "bottom", plot.margin = margin(4, 14, 4, 4))

# ---- C. the power statement --------------------------------------------------
obs_rho <- mean(MR$layer_rho)
pc$disc <- pc$n_fdr05 + 1     # +1 so zeros are visible on a log axis

pC <- ggplot(pc, aes(layer_rho, disc)) +
  geom_vline(xintercept = obs_rho, linewidth = 0.5, colour = "#A63603") +
  geom_line(linewidth = 0.6, colour = "#2C5985") +
  geom_point(size = 1.9, colour = "#2C5985") +
  annotate("text", x = obs_rho, y = 220, label = "observed", hjust = -0.12,
           size = 2.1, colour = "#A63603", fontface = "bold") +
  scale_y_log10(breaks = c(1, 2, 11, 101, 289),
                labels = c("0", "1", "10", "100", "288")) +
  scale_x_continuous(limits = c(-0.12, 1.03)) +
  labs(title = "C · The method finds signal when it is there",
       subtitle = "synthetic signal blended into the functional layer,\nsame FDR machinery; discoveries appear from rho ≈ 0.2",
       x = "correlation between domains", y = "genes at FDR < 0.05") +
  theme_project() +
  theme(plot.margin = margin(4, 14, 4, 4))

# ---- D. the shared layer inflates agreement ----------------------------------
at <- MR$attribution
ad <- rbind(
  data.frame(pair = at$pair, kind = "MOES (shared functional layer)", rho = at$moes_rho),
  data.frame(pair = at$pair, kind = "cis-regulatory only", rho = at$cis_rho))
ad$kind <- factor(ad$kind, levels = c("cis-regulatory only", "MOES (shared functional layer)"))
ad$pair <- factor(ad$pair, levels = rev(at$pair))

pD <- ggplot(ad, aes(rho, pair, fill = kind)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", rho)),
            position = position_dodge(width = 0.72), hjust = -0.18, size = 2.0,
            colour = "#1A1A1A") +
  scale_fill_manual(values = c("cis-regulatory only" = "#8C8C8C",
                               "MOES (shared functional layer)" = "#2C5985"),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 0.78), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "D · Paralogs look alike only via the shared layer",
       subtitle = "agreement between paralog rankings; the shared functional\nlayer pulls all three toward the same genes",
       x = expression(Spearman~rho), y = NULL) +
  theme_project() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Evidence convergence is real but too diffuse to attribute to any gene",
    subtitle = paste0(
      "MOES was specified as four domains and two stages; three were excluded on their own evidence. The two that remain agree weakly about which\n",
      "genes matter most (~2x top-K overlap, global p 0.008–0.065) but no single gene survives testing across 10,387, so no hub list is reported."),
    caption = "Generated by scripts/07_figures/fig04_moes_null.R · 10,387 genes · 10,000 permutations · RRA verified against RobustRankAggreg (D-041)",
    theme = theme_project() +
      theme(plot.title = element_text(size = FIG_BASE + 2, face = "bold"),
            plot.subtitle = element_text(size = FIG_BASE - 0.5, lineheight = 1.25)))

save_fig(fig, "fig04_moes_null", width = 200, height = 170,
         script = "scripts/07_figures/fig04_moes_null.R",
         caption = "The reduced two-domain MOES returns no prioritised hubs; layer independence, FDR curves, power, and shared-layer inflation")

writeLines(c(
  "# Figure 4 — Convergence is real, weak, and not localisable to genes", "",
  "**A.** Overlap between the top K genes of each domain, divided by the K²/N",
  "overlap expected by chance, with the 95% permutation envelope in grey. This is",
  "Robust Rank Aggregation's own question: do genes rank highly in both lists more",
  "often than chance? They do, by up to ~2x, and the excursion clears the envelope",
  sprintf("for two of three paralogs — global p %.3f (MYC), %.3f (MYCN), %.3f (MYCL1).",
          cglob$global_p[cglob$paralog == "MYC"],
          cglob$global_p[cglob$paralog == "MYCN"],
          cglob$global_p[cglob$paralog == "MYCL1"]), "",
  "The global p compares the maximum ratio over K against the permutation",
  "distribution of that maximum, controlling family-wise error across the curve.",
  "Counting pointwise exceedances would badly overstate this: the K values are",
  "nested, so one excursion produces a run of them.", "",
  "This coexists with near-zero Spearman correlation between the domains (-0.023,",
  "+0.012, -0.008) without contradiction. Spearman measures monotone association",
  "across all 10,387 genes; concordance-at-K measures agreement at the TOP of both",
  "lists. The domains agree weakly about which genes matter most and not at all",
  "about the ordering of the rest.", "",
  "Plotted rather than a rank-rank scatter because the cis scores are heavily",
  "tied — a gene with a single promoter link scores exactly 1 — so a scatter bands",
  "into vertical stripes showing tie structure rather than association. Only the K",
  "range where at least 5 genes are expected to overlap by chance is drawn; below",
  "that the ratio is integer noise.", "",
  "Coverage of the cis layer is also uneven. It is",
  sprintf("uneven — %s genes for MYC but %s for MYCN and %s for MYCL1, of %s in the",
          format(n_scored[["MYC"]], big.mark = ","),
          format(n_scored[["MYCN"]], big.mark = ","),
          format(n_scored[["MYCL1"]], big.mark = ","), format(N, big.mark = ",")),
  "MOES universe. Genes with no peak-to-gene link share a single tied rank, and",
  "including them would draw a dense band that looks like structure but is only",
  "the tie block. This unevenness is a real limit on what MOES could have found",
  "for the two smaller paralogs, and it follows directly from MYCN and MYCL1",
  "having two contributing ChIP lines each against MYC's five.", "",
  sprintf("**B.** Empirical FDR by MOES rank, from %s permutations. No curve approaches",
          format(MR$n_perm, big.mark = ",")),
  sprintf("the %.2f threshold; the best FDR achieved across all three paralogs is %.3f.",
          FDR_T, min(rank_tab$fdr)),
  "This is not a marginal miss that a softer threshold would rescue.", "",
  "Panels A and B answer different questions and the pair is the finding: the",
  "aggregate overlap in A is detectable, and the per-gene evidence in B is not.",
  "A signal of ~2x spread across the top of both rankings leaves no individual",
  "gene strong enough to survive multiple testing over 10,387. Aggregate",
  "detectable, per-gene unattributable — which is why no hub list is reported",
  "even though the domains are not independent.", "",
  "**C.** The power statement, without which panel B would be an absence",
  "presented as a finding. A synthetic functional layer was blended with the real",
  "cis ranking at increasing strength and passed through the identical FDR",
  "machinery. The response is monotone: nothing on pure noise, discoveries",
  "appearing once the domains correlate at rho ≈ 0.2, and 288 genes when the two",
  "layers are identical. The observed correlation sits at the far left of this",
  "axis, well below where the method demonstrably detects convergence.", "",
  "**D.** Why any apparent paralog-specificity in a MOES ranking would be an",
  "artefact. The functional layer is identical for all three paralogs, so it pulls",
  "every ranking toward the same genes: paralogs agree at rho +0.50 to +0.59 in",
  "MOES but only +0.17 to +0.33 in the chromatin that actually distinguishes them.",
  "A top-N table drawn from this would present chromatin evidence under a",
  "multi-omics label.", "",
  "No prioritised hub list is reported. Doing so would be wrong twice over: it",
  "would present genes whose individual FDR is ~0.36 as prioritised hits, and it",
  "would imply a paralog attribution that panel D shows comes from chromatin",
  "alone. The full ranking with its FDR column is in",
  "`data/metadata/moes_ranking.csv` so this is inspectable rather than hidden.", "",
  "MOES is a heuristic prioritisation, not a predictive or clinically validated",
  "model (project contract, section 7)."
), "figures/fig04_moes_null_caption.md")
cat("wrote figures/fig04_moes_null_caption.md\n")

# -----------------------------------------------------------------------------
# fig02_m5_gate.R — is the cell-line regulatory layer sound enough to build on?
#
#   Rscript scripts/07_figures/fig02_m5_gate.R
#
# This figure exists to make Figure 1 interpretable. Figure 1 reports a NULL in
# tumours: paralog-resolved programmes lose paralog identity to lineage state. A
# null is only worth reporting if the thing that failed to transfer was real in
# the first place. This figure is that evidence — and it is deliberately shown
# including the criterion that FAILED.
#
#   A  the four-criterion M5 gate, verdicts as they stand (3 PASS / 1 FAIL)
#   B  criterion 4 — each motif is most enriched in its OWN paralog's regions
#   C  criterion 3 — the one that failed: direction yes, magnitude no
#   D  internal validity — leave-one-line-out AUC, 8/9 above threshold
#
# Panel C is not buried. A gate reported only where it passed is not a gate.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg); library(colorspace)
})
source("R/theme_project.R")

gate  <- read.csv("data/metadata/m5_gate_results.csv", stringsAsFactors = FALSE)
shares<- as.matrix(read.csv("data/metadata/m5_motif_shares.csv", row.names = 1))
c3    <- read.csv("data/metadata/m5_criterion3.csv", stringsAsFactors = FALSE)
val   <- read.csv("data/metadata/regulon_validity.csv", stringsAsFactors = FALSE)

cat("accessibility checks\n")
check_palette(PAL_PARALOG, label = "paralog palette")
PAL_VERDICT <- c("PASS" = "#2C5985", "FAIL" = "#A63603")
check_palette(PAL_VERDICT, label = "verdict palette")
check_text_size()
cat("\n")

# Guard: this figure asserts a specific tally in its own subtitle. If the gate
# table moves, the subtitle must not silently keep claiming the old score —
# exactly the staleness D-039 found in the gate table itself.
n_pass <- sum(gate$pass %in% TRUE); n_fail <- sum(gate$pass %in% FALSE)
if (any(is.na(gate$pass)))
  stop("gate table still has an unevaluated criterion; refusing to draw a verdict figure")
cat("gate tally read from disk: ", n_pass, " PASS / ", n_fail, " FAIL\n\n", sep = "")

# ---- A. the gate scorecard ---------------------------------------------------
# Short labels written here rather than taken from the CSV: the stored `value`
# strings are full provenance sentences, correct for a table and unreadable in a
# panel.
sc <- data.frame(
  criterion = c("1 · MYCN sites nest inside MYC",
                "2 · nesting is differential",
                "3 · distal contrast by amplification",
                "4 · motif specificity"),
  observed  = c("0.912 (spread 0.037 across thresholds)",
                "MYCN 5.50x vs MYCL1 4.63x  ·  OR 3.14",
                "+3.3 points  (published +27)",
                "3/3 own-motif top share  ·  p = 1.6e-101"),
  verdict   = ifelse(gate$pass, "PASS", "FAIL"),
  stringsAsFactors = FALSE)
sc$criterion <- factor(sc$criterion, levels = rev(sc$criterion))

pA <- ggplot(sc, aes(y = criterion)) +
  geom_tile(aes(x = 1, fill = verdict), width = 0.9, height = 0.82) +
  geom_text(aes(x = 1, label = verdict), colour = "white",
            size = 2.5, fontface = "bold") +
  geom_text(aes(x = 1.62, label = observed), hjust = 0, size = 2.3,
            colour = "#1A1A1A") +
  scale_fill_manual(values = PAL_VERDICT, guide = "none") +
  scale_x_continuous(limits = c(0.5, 4.7)) +
  labs(title = "A · The M5 gate, including where it failed",
       subtitle = sprintf("cell-line chromatin, %d of 4 criteria pass", n_pass)) +
  theme_project() +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank(), plot.margin = margin(4, 14, 4, 4))

# ---- B. criterion 4, compositional motif share -------------------------------
# Columns of `shares` sum to 1: each column is one paralog's region set, split
# into the three E-box variants. The claim is NOT that a paralog's own motif is
# the commonest in its regions — CACCTG dominates everywhere, which is why the
# earlier enrichment formulation failed on GC content. The claim is that each
# motif reaches its HIGHEST share in its own paralog's regions, i.e. the row-wise
# maximum sits on the diagonal, 3/3.
sm <- as.data.frame(as.table(shares))
names(sm) <- c("motif", "region_set", "share")
ebox <- c(MYC = "CAGATG", MYCN = "CACATG", MYCL1 = "CACCTG")
sm$motif_lab <- factor(paste0(sm$motif, "\n", ebox[as.character(sm$motif)]),
                       levels = paste0(names(ebox), "\n", ebox))
sm$region_set <- factor(sm$region_set, levels = c("MYC", "MYCN", "MYCL1"))
sm$own <- as.character(sm$motif) == as.character(sm$region_set)

# The own-set bar is marked by OUTLINING it in the bar layer itself, not by a
# separate marker layer. A marker layer drawn from subset(sm, own) has one row
# per x group, leaving position_dodge nothing to dodge against, so every marker
# lands on the group centre and points at the wrong bar — the same defect as
# fig01 panel A. An aesthetic carried inside geom_col cannot decouple from it.
pB <- ggplot(sm, aes(motif_lab, share, fill = region_set)) +
  # group = region_set is REQUIRED, not decorative. position_dodge orders bars by
  # the interaction of every discrete aesthetic, so adding `colour` silently
  # re-sorted each group by own-ness: purple sat 3rd in the MYC group and 1st in
  # the others. Pinning the group keeps the order MYC, MYCN, MYCL1 everywhere,
  # which is the only way the legend means the same thing in all three groups.
  geom_col(aes(colour = own, group = region_set),
           position = position_dodge(width = 0.78), width = 0.72, linewidth = 0.5) +
  scale_colour_manual(values = c("TRUE" = "#1A1A1A", "FALSE" = NA), guide = "none") +
  scale_fill_manual(values = PAL_PARALOG, name = "regions from") +
  scale_y_continuous(limits = c(0, 0.72), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "B · Each motif peaks in its own regions",
       subtitle = "compositional share within each region set;\noutlined bar is the own-set bar — highest in 3/3",
       y = "share of E-box occurrences", x = NULL) +
  theme_project() +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

# ---- C. criterion 3, the failure ---------------------------------------------
per <- c3[c3$group != "SUMMARY", ]
per$group <- factor(per$group, levels = c("amplified", "expressing"),
                    labels = c("MYC-amplified", "MYC-expressing"))
obs <- data.frame(
  source = factor(c("this analysis", "Plotnik et al. 2024"),
                  levels = c("this analysis", "Plotnik et al. 2024")),
  contrast = c(c3$pct_distal[c3$group == "SUMMARY"], 27))

pC <- ggplot(obs, aes(source, contrast, fill = source)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%+.1f", contrast)), vjust = -0.45, size = 2.5,
            colour = "#1A1A1A") +
  scale_fill_manual(values = c("this analysis" = PAL_VERDICT[["FAIL"]],
                               "Plotnik et al. 2024" = "#8C8C8C"), guide = "none") +
  scale_y_continuous(limits = c(0, 32), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "C · Criterion 3 fails on magnitude, not direction",
       # Per-line values are grouped, not listed flat: without the group prefix the
       # reader cannot tell which two lines form the amplified arm, which is the
       # whole comparison.
       subtitle = paste0("distal-fraction contrast, MYC-amplified minus MYC-expressing\n",
                         paste(vapply(split(per, per$group), function(d)
                           sprintf("%s: %s", tolower(sub("MYC-", "", d$group[1])),
                                   paste(sprintf("%s %.1f%%", d$line, d$pct_distal),
                                         collapse = ", ")),
                           character(1)), collapse = "   ·   ")),
       y = "contrast (percentage points)", x = NULL) +
  theme_project() +
  theme(panel.grid.major.x = element_blank(), plot.margin = margin(4, 14, 4, 4))

# ---- D. internal validity ----------------------------------------------------
loo <- val[val$type == "loo", ]
loo$paralog <- factor(loo$paralog, levels = c("MYC", "MYCN", "MYCL1"))
loo <- loo[order(loo$paralog, loo$auc), ]
loo$held_out <- factor(loo$held_out, levels = loo$held_out)
LOO_MIN <- 0.70

pD <- ggplot(loo, aes(auc, held_out, colour = paralog)) +
  geom_vline(xintercept = LOO_MIN, linetype = "22", linewidth = 0.4,
             colour = "#666666") +
  geom_segment(aes(x = 0.5, xend = auc, yend = held_out), linewidth = 0.4,
               alpha = 0.5) +
  geom_point(size = 2) +
  scale_colour_manual(values = PAL_PARALOG, name = NULL) +
  scale_x_continuous(limits = c(0.5, 0.88)) +
  labs(title = "D · Regulons survive holding out a line",
       # The threshold is named here rather than annotated inside the panel: an
       # in-panel label anchored to the first category is clipped by the panel edge.
       subtitle = sprintf("leave-one-line-out AUC, %d of %d at or above the\ndashed 0.70 threshold",
                          sum(loo$auc >= LOO_MIN), nrow(loo)),
       x = "AUC (held-out line)", y = NULL) +
  theme_project() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(0.8, 1)) +
  plot_annotation(
    title = "The paralog-resolved regulatory layer is internally sound, with one criterion failed",
    subtitle = paste0(
      "Cell-line chromatin, 10 SCLC lines. This is the evidence that the programmes tested in Figure 1 were real before they were\n",
      "tested in tumours; the criterion that failed is shown rather than omitted."),
    caption = "Generated by scripts/07_figures/fig02_m5_gate.R · criterion 3 evaluated at M7 from DepMap Public 26Q1 copy number (D-035, D-039)",
    theme = theme_project() +
      theme(plot.title = element_text(size = FIG_BASE + 2, face = "bold"),
            plot.subtitle = element_text(size = FIG_BASE - 0.5, lineheight = 1.25)))

save_fig(fig, "fig02_m5_gate", width = 200, height = 165,
         script = "scripts/07_figures/fig02_m5_gate.R",
         caption = "M5 gate verdicts, motif specificity, the failed distal-contrast criterion, and leave-one-line-out validity")

writeLines(c(
  "# Figure 2 — The M5 gate on the paralog-resolved regulatory layer", "",
  "**A.** The four pre-registered gate criteria and their verdicts. Three pass,",
  "one fails. Criterion 3 was blocked at M5 because amplification status could not",
  "be verified, and was evaluated at M7 once DepMap Public 26Q1 copy number",
  "resolved it (D-035); the verdict is recorded in the gate table itself (D-039).", "",
  "**B.** Criterion 4, compositional motif share. Each column of the underlying",
  "matrix is one paralog's region set split across the three E-box variants, so",
  "shares sum to 1 within a set. CACCTG is the commonest variant in every set,",
  "which is why an enrichment-over-background formulation failed on GC content.",
  "The test asked instead whether each motif reaches its highest share in its own",
  "paralog's regions. It does, in 3/3 (chi-square p = 1.6e-101).", "",
  "**C.** Criterion 3, the failure. The direction reproduces — MYC-amplified lines",
  "carry a higher distal fraction than MYC-expressing lines — but the magnitude is",
  "+3.3 points against a published +27. Contributing factors, none separable at",
  "2 lines vs 2: the shared ATAC-defined universe is ~84% distal by construction,",
  "which compresses the achievable range; Plotnik called peaks de novo per line",
  "with no shared grid; and H847 is absent from DepMap so the expressing group is",
  "n=2. No significance test is applied, because none would be honest at this n.", "",
  "**D.** Leave-one-line-out AUC. Each regulon is rebuilt with one contributing",
  "line withheld and scored on that line. 8 of 9 reach the 0.70 threshold; the",
  "exception is MYC/H1048 at 0.69. This is internal validity only — it shows the",
  "programmes are not driven by a single line, not that they generalise to tumours.",
  "Figure 1 shows that they do not."
), "figures/fig02_m5_gate_caption.md")
cat("wrote figures/fig02_m5_gate_caption.md\n")

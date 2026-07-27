# -----------------------------------------------------------------------------
# fig01_lineage_dominance.R — the project's headline result.
#
#   Rscript scripts/07_figures/fig01_lineage_dominance.R
#
# Paralog-resolved regulatory programmes do NOT retain paralog identity in patient
# tumours independently of neuroendocrine lineage state — replicated across two
# independent cohorts (160 patients), with the mechanism reproduced in three
# independent datasets.
#
# Four panels, each carrying one claim:
#   A  variance partitioning — lineage vs paralog, both cohorts. THE result.
#   B  paralog association collapses to zero after lineage adjustment.
#   C  MYC/NE antagonism reproduces in 3 datasets (the mechanism).
#   D  specificity matrix — the diagonal does not dominate.
#
# Panels A and B carry the negative finding; C shows it has a cause rather than
# being an unexplained null.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg); library(colorspace)
})
source("R/theme_project.R")

rep  <- readRDS("data/processed/tumour/replication.rds")$per_cohort
spec <- as.matrix(utils::read.csv("data/metadata/m6_specificity_matrix.csv", row.names = 1))
cell <- utils::read.csv("data/metadata/m7_paralog_lineage_celllines.csv", stringsAsFactors = FALSE)

# Panel A encodes VARIANCE SOURCE, not gene, so it must not borrow a paralog
# colour: doing so implies "this bar is MYC" in a panel where paralog is the
# x-axis. Blue carries lineage (the term that explains the variance); neutral
# dark grey carries paralog (the term that does not), which also keeps the
# reader from reading a paralog identity into it.
PAL_SRC <- c("paralog" = "#4D4D4D", "lineage" = "#2C5985")
cat("accessibility checks\n")
check_palette(PAL_SRC, label = "source palette")
check_palette(PAL_PARALOG, label = "paralog palette")
check_text_size()
cat("\n")

rep$cohort <- factor(rep$cohort, levels = c("GSE60052", "George2015"),
                     labels = c("GSE60052 (n=79)", "George 2015 (n=81)"))
rep$paralog <- factor(rep$paralog, levels = c("MYC", "MYCN", "MYCL1"))

# ---- A. variance partitioning ------------------------------------------------
vp <- rbind(
  data.frame(cohort = rep$cohort, paralog = rep$paralog,
             source = "paralog", r2 = rep$unique_paralog),
  data.frame(cohort = rep$cohort, paralog = rep$paralog,
             source = "lineage", r2 = rep$unique_lineage))
vp$source <- factor(vp$source, levels = c("lineage", "paralog"))

# The label layer keeps BOTH fill levels and blanks the lineage ones. Subsetting
# to source == "paralog" leaves the layer with a single fill level, so
# position_dodge has nothing to dodge against and drops every label onto the
# group centre instead of over its own bar.
pA <- ggplot(vp, aes(paralog, r2, fill = source)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  # Rotated: the values being labelled are near zero, so a horizontal label is
  # several times wider than its own bar and collides with the neighbouring bar
  # or the panel edge. Vertical text has a narrow footprint and clears both.
  geom_text(aes(label = ifelse(source == "paralog", sprintf("%.3f", r2), "")),
            position = position_dodge(width = 0.72), angle = 90,
            hjust = -0.15, vjust = 0.5, size = 2.0, colour = PAL_SRC[["paralog"]]) +
  facet_wrap(~cohort) +
  scale_fill_manual(values = PAL_SRC, name = NULL,
                    labels = c("lineage / NE state", "paralog expression")) +
  scale_y_continuous(limits = c(0, 0.55), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "A · Lineage explains the scores, not the paralog",
       subtitle = "unique variance explained, after removing the other term",
       y = expression(unique~R^2), x = NULL) +
  theme_project() +
  # Extra right margin: titles are plot-positioned, so the left column's title
  # would otherwise butt straight into panel B's.
  theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
        plot.margin = margin(4, 14, 4, 4))

# ---- B. association before vs after adjustment -------------------------------
pc <- rbind(
  data.frame(cohort = rep$cohort, paralog = rep$paralog, stage = "raw", rho = rep$rho_raw),
  data.frame(cohort = rep$cohort, paralog = rep$paralog, stage = "lineage-adjusted",
             rho = rep$rho_partial))
pc$stage <- factor(pc$stage, levels = c("raw", "lineage-adjusted"))

pB <- ggplot(pc, aes(stage, rho, group = interaction(cohort, paralog), colour = paralog)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#999999") +
  geom_line(linewidth = 0.5, alpha = 0.85) +
  geom_point(size = 1.7) +
  facet_wrap(~cohort) +
  scale_colour_manual(values = PAL_PARALOG, name = NULL) +
  labs(title = "B · The association collapses after adjustment",
       subtitle = "Spearman rho, regulon score vs its own paralog's expression",
       y = expression(rho), x = NULL) +
  theme_project() +
  theme(legend.position = "bottom")

# ---- C. the mechanism, three datasets ----------------------------------------
myc_ne <- data.frame(
  dataset = c("GSE60052\ntumours", "George 2015\ntumours", "CCLE\ncell lines"),
  n = c(79, 81, 59),
  rho = c(unique(rep$myc_vs_ne[rep$cohort == "GSE60052 (n=79)"])[1],
          unique(rep$myc_vs_ne[rep$cohort == "George 2015 (n=81)"])[1],
          cell$rho[cell$paralog == "MYC" & cell$target == "NE_score"][1]))
myc_ne$dataset <- factor(myc_ne$dataset, levels = myc_ne$dataset)

pC <- ggplot(myc_ne, aes(dataset, rho)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#999999") +
  geom_col(width = 0.55, fill = PAL_PARALOG[["MYC"]]) +
  # Bars run DOWNWARD from zero, so vjust must be negative to seat the label
  # inside the bar. A positive vjust puts white text below the bar end, on white
  # background, where it is invisible rather than merely misplaced.
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -0.8, size = 2.3, colour = "white") +
  scale_y_continuous(limits = c(-0.7, 0.05)) +
  labs(title = "C · The mechanism, reproduced three times",
       subtitle = "MYC expression vs neuroendocrine score (Ireland et al. 2020)",
       y = expression(Spearman~rho), x = NULL) +
  theme_project() +
  theme(panel.grid.major.x = element_blank(), plot.margin = margin(4, 14, 4, 4))

# ---- D. specificity matrix ---------------------------------------------------
sm <- as.data.frame(as.table(spec))
names(sm) <- c("regulon", "paralog_expression", "rho")
sm$regulon <- factor(sm$regulon, levels = rev(c("MYC","MYCN","MYCL1")))
sm$paralog_expression <- factor(sm$paralog_expression, levels = c("MYC","MYCN","MYCL1"))
sm$diag <- as.character(sm$regulon) == as.character(sm$paralog_expression)

pD <- ggplot(sm, aes(paralog_expression, regulon, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_tile(data = subset(sm, diag), colour = "#1A1A1A", linewidth = 0.7, fill = NA) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.3,
            colour = ifelse(abs(sm$rho) > 0.3, "white", "#1A1A1A")) +
  scale_fill_diverging(limits = c(-0.5, 0.5), name = expression(rho)) +
  labs(title = "D · The diagonal does not dominate",
       subtitle = "regulon score vs each paralog's expression (GSE60052);\nboxed cells are the expected self-associations") +
  theme_matrix() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 0.95)) +
  plot_annotation(
    # Kept short enough to render inside 200 mm at this weight — the longer
    # phrasing clipped mid-word at the right edge.
    title = "Regulatory programmes lose paralog identity to lineage state in SCLC tumours",
    subtitle = paste0(
      "Replicated across two independent patient cohorts (160 tumours). Regulons were built from paralog-resolved cell-line chromatin and are\n",
      "internally valid there; in tumours their scores are explained by neuroendocrine lineage state, not by the paralog's own expression."),
    caption = "Generated by scripts/07_figures/fig01_lineage_dominance.R · regulon scores are mean per-gene z-scores, identical in both cohorts",
    theme = theme_project() +
      theme(plot.title = element_text(size = FIG_BASE + 2.5, face = "bold"),
            plot.subtitle = element_text(size = FIG_BASE - 0.5, lineheight = 1.25)))

save_fig(fig, "fig01_lineage_dominance", width = 200, height = 175,
         script = "scripts/07_figures/fig01_lineage_dominance.R",
         caption = "Lineage dominance over paralog identity, replicated in two cohorts")

writeLines(c(
  "# Figure 1 — Lineage dominance over paralog identity", "",
  "**A.** Unique variance in regulon score attributable to the paralog's own expression",
  "versus lineage/NE state, after removing the other term. Lineage explains 0.24-0.43;",
  "the paralog explains 0.000-0.068. Holds in 3/3 regulons in both cohorts.", "",
  "**B.** Spearman correlation between regulon score and its own paralog's expression,",
  "before and after regressing out NE score and all four lineage TFs. Every association",
  "collapses toward zero. 0/3 survive in George 2015; the single GSE60052 survivor",
  "(MYCL1) is method-dependent and comes from the regulon that failed M5 validation.", "",
  "**C.** MYC expression versus neuroendocrine score in three independent datasets,",
  "reproducing Ireland et al. 2020. This is why the effect runs in the direction it does:",
  "MYC-high SCLC is the least neuroendocrine, and the MYC regulon is neurogenesis-enriched",
  "(p = 7e-6), so MYC-high tumours score LOW on it.", "",
  "**D.** Regulon score against each paralog's expression. If programmes retained paralog",
  "identity the boxed diagonal would dominate. It does not (2/3 nominal, P = 0.259).", "",
  "A null result here was pre-committed in the gap statement (section 5) as a reportable",
  "finding rather than a failure, and is reported with power analysis: at n=79,",
  "|rho| >= 0.31 was detectable at 80% power."
), "figures/fig01_lineage_dominance_caption.md")
cat("wrote figures/fig01_lineage_dominance_caption.md\n")

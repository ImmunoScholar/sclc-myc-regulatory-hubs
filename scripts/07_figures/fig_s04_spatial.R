# -----------------------------------------------------------------------------
# fig_s04_spatial.R — the same question, a third modality.
#
#   Rscript scripts/07_figures/fig_s04_spatial.R
#
# Panels C and D deliberately mirror Figure 1A and 1B. If lineage dominance were
# an artefact of bulk RNA-seq scoring, it should not survive a change of tissue,
# chemistry, normalisation and gene panel. It does.
#
#   A  the coverage gate — only MYC is scoreable, and only just
#   B  coherence: about a third of the score varies WITHIN one tumour
#   C  unique variance, lineage vs MYC (mirrors Fig 1A)
#   D  the association collapses under adjustment (mirrors Fig 1B)
#
# Panel A leads because it bounds everything after it: this is a 56-of-500-gene
# proxy, and a figure that showed the result without the limitation would be
# claiming more than the panel can measure.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg); library(colorspace); library(yaml)
})
source("R/theme_project.R")

cov   <- read.csv("data/metadata/m9_panel_coverage.csv", stringsAsFactors = FALSE)
coh   <- read.csv("data/metadata/m9_coherence.csv", stringsAsFactors = FALSE)
vp    <- read.csv("data/metadata/m9_variance.csv", stringsAsFactors = FALSE)
assoc <- read.csv("data/metadata/m9_associations.csv", stringsAsFactors = FALSE)
SP    <- yaml::read_yaml("config/params.yml")$spatial

cat("accessibility checks\n")
check_palette(PAL_PARALOG, label = "paralog palette")
PAL_SRC <- c("paralog" = "#4D4D4D", "lineage" = "#2C5985")
check_palette(PAL_SRC, label = "source palette")
check_text_size()
cat("\n")

LAB <- c(GSE261348 = "IMfirst\n(174 ROI, 18 slides)",
         GSE261345 = "CANTABRICO\n(121 ROI, 14 slides)")
relab <- function(x) factor(LAB[x], levels = unname(LAB))

# Guard: the figure states that exactly one regulon is scoreable.
u <- cov[!duplicated(cov$paralog), ]
if (sum(u$scoreable) != 1L)
  stop("panel coverage now makes ", sum(u$scoreable),
       " regulons scoreable; this figure asserts exactly one and must be rewritten.")

# ---- A. the coverage gate ----------------------------------------------------
u$paralog <- factor(u$paralog, levels = c("MYC", "MYCN", "MYCL1"))
FLOOR <- SP$min_coverage_fraction

pA <- ggplot(u, aes(paralog, 100 * coverage_fraction, fill = paralog)) +
  geom_hline(yintercept = 100 * FLOOR, linetype = "22", linewidth = 0.5,
             colour = "#A63603") +
  geom_col(aes(colour = scoreable), width = 0.62, linewidth = 0.5) +
  # Labels sit INSIDE the bars. Above them, MYCL1's label (bar at 9.6%) lands
  # exactly on the 10% floor line and its annotation.
  geom_text(aes(label = sprintf("%d/%d", measured, regulon_size)),
            vjust = 1.6, size = 2.1, colour = "white") +
  annotate("text", x = 3.4, y = 100 * FLOOR, label = sprintf("%.0f%% floor", 100 * FLOOR),
           hjust = 1, vjust = -0.6, size = 2.0, colour = "#A63603") +
  scale_fill_manual(values = PAL_PARALOG, guide = "none") +
  scale_colour_manual(values = c("TRUE" = "#1A1A1A", "FALSE" = NA), guide = "none") +
  scale_y_continuous(limits = c(0, 14), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A · Only the MYC regulon is measurable",
       subtitle = "fraction of regulon members on the CTA panel;\noutlined bar passes the gate",
       y = "% of regulon on panel", x = NULL) +
  theme_project() +
  theme(panel.grid.major.x = element_blank(), plot.margin = margin(4, 14, 4, 4))

# ---- B. coherence ------------------------------------------------------------
cb <- coh[!is.na(coh$between_slide_r2), ]
cb$lab <- ifelse(cb$subset == "single-patient slides",
                 sprintf("%s\nsingle-patient (%d ROI)", sub("\n.*", "", LAB[cb$cohort]), cb$n_roi),
                 sprintf("%s\nall slides (%d ROI)", sub("\n.*", "", LAB[cb$cohort]), cb$n_roi))
cb$lab <- factor(cb$lab, levels = rev(cb$lab))
cbl <- rbind(
  data.frame(lab = cb$lab, part = "between slides", r2 = cb$between_slide_r2),
  data.frame(lab = cb$lab, part = "within slides", r2 = cb$within_slide_r2))
cbl$part <- factor(cbl$part, levels = c("within slides", "between slides"))

pB <- ggplot(cbl, aes(r2, lab, fill = part)) +
  geom_col(width = 0.6) +
  geom_text(data = cb, aes(x = between_slide_r2 / 2, y = lab,
                           label = sprintf("%.2f", between_slide_r2)),
            inherit.aes = FALSE, size = 2.0, colour = "white") +
  scale_fill_manual(values = c("between slides" = "#2C5985",
                               "within slides" = "#C6C6C6"), name = NULL) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  labs(title = "B · A third of the score varies within one tumour",
       subtitle = "variance in MYC regulon score, partitioned by slide",
       x = expression(R^2), y = NULL) +
  theme_project() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

# ---- C. unique variance (mirrors Fig 1A) -------------------------------------
vl <- rbind(
  data.frame(cohort = vp$cohort, source = "paralog", r2 = vp$unique_myc),
  data.frame(cohort = vp$cohort, source = "lineage", r2 = vp$unique_lineage))
vl$cohort <- relab(vl$cohort)
vl$source <- factor(vl$source, levels = c("lineage", "paralog"))

pC <- ggplot(vl, aes(cohort, r2, fill = source)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = ifelse(source == "paralog", sprintf("%.3f", r2), "")),
            position = position_dodge(width = 0.72), angle = 90,
            hjust = -0.15, vjust = 0.5, size = 2.0, colour = PAL_SRC[["paralog"]]) +
  scale_fill_manual(values = PAL_SRC, name = NULL,
                    labels = c("lineage / NE state", "MYC expression")) +
  scale_y_continuous(limits = c(0, 0.62), expand = expansion(mult = c(0, 0.03))) +
  labs(title = "C · Lineage explains the score here too",
       subtitle = "unique variance explained, after removing the other term",
       y = expression(unique~R^2), x = NULL) +
  theme_project() +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
        plot.margin = margin(4, 14, 4, 4))

# ---- D. adjustment collapse (mirrors Fig 1B) ---------------------------------
ar <- assoc[assoc$variable %in% c("MYC", "MYC (lineage-adjusted)"), ]
ar$stage <- factor(ifelse(ar$variable == "MYC", "raw", "lineage-adjusted"),
                   levels = c("raw", "lineage-adjusted"))
ar$cohort_l <- relab(ar$cohort)

pD <- ggplot(ar, aes(stage, rho, group = cohort_l, colour = cohort_l)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#999999") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  geom_text(data = subset(ar, stage == "lineage-adjusted"),
            aes(label = sprintf("%+.3f", rho)), hjust = -0.25, size = 2.0,
            show.legend = FALSE) +
  scale_colour_manual(values = c("#2C5985", "#A63603"), name = NULL) +
  scale_y_continuous(limits = c(-0.28, 0.28)) +
  labs(title = "D · The association does not survive adjustment",
       subtitle = "Spearman rho, MYC regulon score vs MYC expression",
       y = expression(rho), x = NULL) +
  theme_project() +
  theme(legend.position = "bottom")

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Lineage dominance reproduces in spatial data, on the one regulon the panel can measure",
    subtitle = paste0(
      "GeoMx DSP, two cohorts, 295 ROIs. A targeted ~1,800-gene panel measures only 56 of 500 MYC regulon members and fewer for MYCN and\n",
      "MYCL1, which are not scoreable at all. Within that limit, the M6/M7 result holds in tissue with different chemistry and normalisation."),
    caption = "Generated by scripts/07_figures/fig_s04_spatial.R · GSE261348 + GSE261345 · slide IMF-001/002 excluded (R-17) · NE score from 3 of 10 markers",
    theme = theme_project() +
      theme(plot.title = element_text(size = FIG_BASE + 2, face = "bold"),
            plot.subtitle = element_text(size = FIG_BASE - 0.5, lineheight = 1.25)))

save_fig(fig, "fig_s04_spatial", width = 200, height = 175,
         script = "scripts/07_figures/fig_s04_spatial.R",
         caption = "Spatial coverage gate, within-tumour coherence, and orthogonal replication of lineage dominance")

writeLines(c(
  "# Figure S4 — Spatial coherence and orthogonal validation", "",
  "**A.** The coverage gate. The GeoMx CTA panel measures ~1,800 targets, so most",
  "regulon members are simply absent. Against a floor of 20 measured genes AND 10%",
  "of members, only the MYC regulon qualifies, at 56 of 500 (11.2%). MYCN reaches",
  "8.4% and MYCL1 9.6%, and neither is scored anywhere in this figure. These are",
  "narrow verdicts — MYC clears by 1.2 percentage points and MYCL1 misses by 0.4 —",
  "against a threshold raised mid-project, and are reported as narrow rather than",
  "as clean separations. Everything below is a 56-gene proxy, not the regulon.", "",
  "**B.** Variance in the MYC regulon score partitioned by slide. Roughly half sits",
  "between slides and half within them. Restricted to the three IMfirst slides that",
  "carry a single patient identifier — the only ones where within-slide is",
  "unambiguously within-tumour — 0.70 is between and 0.30 within. So about a third",
  "of the score is regional variation inside one tumour, not a tumour-level",
  "property. CANTABRICO contributes no such estimate: every one of its slides",
  "carries two patient identifiers with no per-ROI label.", "",
  "**C.** The Figure 1A test in a third modality. Unique variance in the regulon",
  "score attributable to MYC's own expression versus lineage state. Lineage",
  "explains 0.474 and 0.508; MYC explains 0.019 and 0.010. The direction and the",
  "magnitude of the gap both match the bulk cohorts, in tissue with different",
  "chemistry, normalisation and gene panel.", "",
  "**D.** The Figure 1B test. The raw association between regulon score and MYC",
  "expression is weak to begin with (+0.11, +0.14) and does not survive adjustment",
  "for lineage, falling to -0.02 (p = 0.77) and -0.17 (p = 0.06). Note the sign",
  "change in CANTABRICO: after adjustment the residual association runs negative,",
  "consistent with the MYC/NE antagonism reported in Figure 1C rather than with any",
  "residual paralog-specific signal.", "",
  "This is the strongest available evidence that the M6/M7 null is a property of",
  "the biology rather than of bulk deconvolution or of one scoring choice. It is",
  "also bounded: it speaks only for MYC, only through 56 genes, and only where the",
  "panel looks."
), "figures/fig_s04_spatial_caption.md")
cat("wrote figures/fig_s04_spatial_caption.md\n")

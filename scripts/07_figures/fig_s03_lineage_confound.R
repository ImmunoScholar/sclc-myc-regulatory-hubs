# -----------------------------------------------------------------------------
# fig_s03_lineage_confound.R — the confound, at every level it was tested.
#
#   Rscript scripts/07_figures/fig_s03_lineage_confound.R
#
# Figure 1 reports that lineage explains the regulon scores and the paralog does
# not. The obvious objection is that this is an artefact of tumour-level scoring.
# It is not: the same confound is present in the chromatin the regulons were
# built from, in the cell lines, and gene-by-gene genome-wide.
#
#   A  chromatin — MYC/MYCN-active regions are POU2F3-bound far above background
#   B  cell lines — MYC expression tracks NE state and ASCL1 (R-01, D-037)
#   C  genes — almost nothing survives lineage adjustment, and what does is not
#      enriched in the regulons
#
# This is the figure that says the negative result is a property of the biology
# as measured here, not of one scoring choice.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg); library(colorspace)
})
source("R/theme_project.R")

occ  <- read.csv("data/metadata/m6_occupancy_confounding.csv", stringsAsFactors = FALSE)
cell <- read.csv("data/metadata/m7_paralog_lineage_celllines.csv", stringsAsFactors = FALSE)
gl   <- read.csv("data/metadata/m6_gene_level_summary.csv", stringsAsFactors = FALSE)

cat("accessibility checks\n")
check_palette(PAL_PARALOG, label = "paralog palette")
PAL_SET <- c("universe" = "#BBBBBB", "paralog-active" = "#2C5985")
check_palette(PAL_SET, label = "region-set palette")
check_text_size()
cat("\n")

# ---- A. chromatin-level co-occupancy -----------------------------------------
occ$lab <- sprintf("%s\n(%s)", occ$line, occ$paralog)
oc <- rbind(
  data.frame(lab = occ$lab, set = "universe",       pct = occ$pct_universe_pou2f3),
  data.frame(lab = occ$lab, set = "paralog-active", pct = occ$pct_active_pou2f3))
oc$set <- factor(oc$set, levels = c("universe", "paralog-active"))

pA <- ggplot(oc, aes(lab, pct, fill = set)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(data = occ, aes(x = lab, y = pct_active_pou2f3 + 4,
                            label = sprintf("%.1fx", enrichment)),
            inherit.aes = FALSE, size = 2.1, colour = "#1A1A1A") +
  scale_fill_manual(values = PAL_SET, name = NULL) +
  scale_y_continuous(limits = c(0, 85), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A · The confound is already in the chromatin",
       subtitle = "POU2F3-bound fraction: paralog-active regions vs the universe\nlabels give the enrichment; all three p < 1e-15",
       y = "% of regions POU2F3-bound", x = NULL) +
  theme_project() +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank(),
        plot.margin = margin(4, 14, 4, 4))

# ---- B. cell-line level ------------------------------------------------------
cell$paralog <- factor(c(MYC = "MYC", MYCN = "MYCN", MYCL = "MYCL1")[cell$paralog],
                       levels = c("MYC", "MYCN", "MYCL1"))
cell$target <- factor(cell$target,
                      levels = rev(c("NE_score", "ASCL1", "NEUROD1", "POU2F3", "YAP1")),
                      labels = rev(c("NE score", "ASCL1", "NEUROD1", "POU2F3", "YAP1")))
cell$sig <- ifelse(cell$fdr < 0.05, "*", "")

pB <- ggplot(cell, aes(paralog, target, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f%s", rho, sig)), size = 2.2,
            colour = ifelse(abs(cell$rho) > 0.28, "white", "#1A1A1A")) +
  scale_fill_diverging(limits = c(-0.5, 0.5), name = expression(rho)) +
  labs(title = "B · MYC tracks lineage in cell lines too",
       subtitle = "Spearman rho, CCLE SCLC lines (n=59); * is FDR < 0.05") +
  theme_matrix() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# ---- C. gene level, genome-wide ----------------------------------------------
gl$paralog <- factor(gl$paralog, levels = c("MYC", "MYCN", "MYCL1"))
gl$lab <- ifelse(gl$n_fdr05 == 0, "0 of 33,682",
                 sprintf("%s of %s", format(gl$n_fdr05, big.mark = ","),
                         format(gl$n_tested, big.mark = ",")))
gl$enr_lab <- ifelse(is.na(gl$enrich_or), "no survivors to test",
                     sprintf("OR %.2f, p = %.2f", gl$enrich_or, gl$enrich_p))

pC <- ggplot(gl, aes(paralog, n_fdr05, fill = paralog)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = lab), vjust = -0.55, size = 2.1, colour = "#1A1A1A") +
  geom_text(aes(y = -14, label = enr_lab), size = 1.9, colour = "#666666") +
  scale_fill_manual(values = PAL_PARALOG, guide = "none") +
  scale_y_continuous(limits = c(-26, 250), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "C · Survivors are not enriched in the regulons",
       subtitle = "genes still associated with the paralog after lineage adjustment\ngrey text: enrichment of survivors inside vs outside the regulon",
       y = "genes surviving (FDR < 0.05)", x = NULL) +
  theme_project() +
  theme(panel.grid.major.x = element_blank())

fig <- pA / (pB | pC) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = "The lineage confound is present at chromatin, cell-line and gene level — not only in tumour scores",
    subtitle = paste0(
      "The negative result in Figure 1 does not depend on how regulon scores were computed. Neuroendocrine lineage state is entangled with\n",
      "MYC-family regulatory occupancy in the cell lines the regulons were built from, and survives every level at which it was tested."),
    caption = "Generated by scripts/07_figures/fig_s03_lineage_confound.R · POU2F3 ChIP from GSE210113 · CCLE expression · GSE60052 (D-033, R-01, R-14)",
    theme = theme_project() +
      theme(plot.title = element_text(size = FIG_BASE + 2, face = "bold"),
            plot.subtitle = element_text(size = FIG_BASE - 0.5, lineheight = 1.25)))

save_fig(fig, "fig_s03_lineage_confound", width = 200, height = 175,
         script = "scripts/07_figures/fig_s03_lineage_confound.R",
         caption = "The lineage confound at chromatin, cell-line and gene level")

writeLines(c(
  "# Figure S3 — The lineage confound at every level tested", "",
  "**A.** Chromatin-level co-occupancy. In the three keystone lines with POU2F3",
  "ChIP, regions active for MYC or MYCN are POU2F3-bound at 2.5–3.1x the rate of",
  "the ATAC universe they were drawn from (all p < 1e-15). The entanglement",
  "between MYC-family occupancy and lineage-factor occupancy is therefore present",
  "in the very data the regulons were built from, before any tumour is scored.",
  "Resolution is not uniform across lineage TFs and was never pooled: POU2F3 has",
  "three keystone lines, ASCL1 one, NEUROD1 none (risk R-14).", "",
  "**B.** Cell-line level. Across 59 CCLE SCLC lines, MYC expression correlates",
  "negatively with neuroendocrine score (rho -0.441, FDR 0.008) and with ASCL1",
  "(rho -0.409, FDR 0.010). MYCN and MYCL show no comparable structure. This is",
  "the independent reproduction of Ireland et al. 2020 in a third dataset, and it",
  "is also one of the two grounds on which the network domain was excluded: a",
  "cell-line network would inherit exactly this confound (D-037).", "",
  "**C.** Gene level, genome-wide. Of 33,682 genes tested per paralog, the number",
  "still associated with the paralog after adjusting for NE score and all four",
  "lineage TFs is 206 for MYC, 62 for MYCL1, and zero for MYCN. Crucially, the",
  "survivors are not concentrated in the regulons: enrichment inside versus",
  "outside is OR 1.21 (p = 0.57) for MYC and OR 1.57 (p = 0.48) for MYCL1 —",
  "background rate. Survivors exist, but they are not the regulon genes, so they",
  "do not rescue the paralog-specific claim.", "",
  "Taken together these three levels rule out the most obvious objection to",
  "Figure 1 — that the null is an artefact of tumour-level regulon scoring."
), "figures/fig_s03_lineage_confound_caption.md")
cat("wrote figures/fig_s03_lineage_confound_caption.md\n")

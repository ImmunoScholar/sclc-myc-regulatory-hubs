# -----------------------------------------------------------------------------
# fig_s02_regulatory_pipeline.R — how signal became regulons, and what guards it.
#
#   Rscript scripts/07_figures/fig_s02_regulatory_pipeline.R
#
# Supplementary methods. Every panel shows a step where a defensible-looking
# choice could have produced a wrong answer, together with the check that caught
# it. Several of these guards exist because the first attempt failed:
#
#   A  threshold calibration — regions are real signal, and the cutoff is not
#      load-bearing (D-024)
#   B  criterion 1 is stable within its convention, and the convention matters
#   C  super-enhancer calling, with the guard added after the first rule
#      returned every stitched enhancer as a super-enhancer
#   D  peak-to-gene linking, which is not nearest-gene assignment
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg); library(colorspace); library(yaml)
})
source("R/theme_project.R")

iv   <- read.csv("data/metadata/intrinsic_validation.csv", stringsAsFactors = FALSE)
conv <- read.csv("data/metadata/paralog_region_conventions.csv", stringsAsFactors = FALSE)
se   <- read.csv("data/metadata/super_enhancer_summary.csv", stringsAsFactors = FALSE)
p2g  <- read.csv("data/metadata/peak_to_gene_summary.csv", stringsAsFactors = FALSE)
bl   <- read.csv("data/metadata/blacklist_filter_report.csv", stringsAsFactors = FALSE)
lo   <- read.csv("data/metadata/liftover_report.csv", stringsAsFactors = FALSE)
cfg  <- yaml::read_yaml("config/params.yml")

CHOSEN_FOLD <- cfg$regions$atac_threshold_multiple
PRIMARY     <- cfg$active_regions$primary_fold
SE_GUARD    <- 30   # hard guard, % of stitched enhancers (see panel C)

cat("accessibility checks\n")
PAL_RULE <- c("ge2" = "#2C5985", "union" = "#E08214", "all" = "#8C8C8C")
check_palette(PAL_RULE, label = "convention palette")
check_text_size()
cat("\n")

getval <- function(m) p2g$value[p2g$metric == m]

# ---- A. threshold calibration ------------------------------------------------
# Two independent marks of real regulatory signal — TSS enrichment and H3K27ac
# fold — both rise monotonically with the ATAC threshold. That is the evidence
# the regions are signal rather than background, and it holds across the whole
# range, which is why the chosen cutoff is not load-bearing (D-024).
ivl <- rbind(
  data.frame(multiple = iv$multiple, metric = "TSS enrichment", value = iv$tss_enrich),
  data.frame(multiple = iv$multiple, metric = "H3K27ac fold",   value = iv$k27_fold))

pA <- ggplot(ivl, aes(multiple, value, colour = metric)) +
  geom_vline(xintercept = CHOSEN_FOLD, linetype = "22", linewidth = 0.4,
             colour = "#666666") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = c("TSS enrichment" = "#2C5985",
                                 "H3K27ac fold" = "#A63603"), name = NULL) +
  scale_x_continuous(breaks = iv$multiple) +
  labs(title = "A · Regions carry real regulatory signal",
       subtitle = sprintf("both marks rise monotonically; dashed line is the chosen %sx cutoff",
                          CHOSEN_FOLD),
       x = "ATAC threshold (multiple of genome median)", y = "fold enrichment") +
  theme_project() +
  theme(legend.position = "bottom", plot.margin = margin(4, 14, 4, 4))

# ---- B. nesting stability vs convention --------------------------------------
# Criterion 1 asks what fraction of MYCN-bound regions sit inside MYC-bound
# regions. That fraction depends on how a "paralog region" is defined, and the
# figure says so: the pre-specified >=2-replicate rule is stable across a 2.7-fold
# change in threshold, but the union and all-replicate conventions sit elsewhere.
# Reporting only the chosen convention would hide that the choice does work.
conv$rule <- factor(conv$rule, levels = c("ge2", "union", "all"))
spread_ge2 <- diff(range(conv$mycn_in_myc[conv$rule == "ge2"]))

pB <- ggplot(conv, aes(fold, mycn_in_myc, colour = rule, group = rule)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.6) +
  geom_point(data = subset(conv, rule == "ge2" & fold == PRIMARY),
             size = 3.4, shape = 21, fill = NA, colour = "#1A1A1A", stroke = 0.6) +
  scale_colour_manual(values = PAL_RULE, name = "replicate rule",
                      labels = c("ge2" = "≥2 replicates (used)", "union" = "union",
                                 "all" = "all replicates")) +
  scale_y_continuous(limits = c(0.3, 1)) +
  labs(title = "B · Criterion 1 is stable within its convention",
       subtitle = sprintf("MYCN-in-MYC nesting; circled point is the reported value\n≥2-replicate rule varies by only %.3f across thresholds",
                          spread_ge2),
       x = "fold threshold", y = "fraction of MYCN regions inside MYC") +
  theme_project() +
  theme(legend.position = "bottom")

# ---- C. super-enhancer calling ----------------------------------------------
# The first stitching rule (first slope >= 1) returned 8,238 of 8,238 stitched
# enhancers as super-enhancers for H1048 — every one. Replaced by a knee-point
# cutoff with a hard guard at 30% of stitched enhancers. Every line now lands
# far below the guard, which is the check that the collapse has not recurred.
se <- se[order(se$pct_of_stitched), ]
se$line <- factor(se$line, levels = se$line)

pC <- ggplot(se, aes(pct_of_stitched, line)) +
  geom_segment(aes(x = 0, xend = pct_of_stitched, yend = line),
               linewidth = 0.4, colour = "#BBBBBB") +
  geom_point(size = 2, colour = "#2C5985") +
  geom_text(aes(label = format(n_SE, big.mark = ",")), hjust = -0.45, size = 1.9,
            colour = "#666666") +
  scale_x_continuous(limits = c(0, SE_GUARD), breaks = seq(0, 30, 10)) +
  labs(title = "C · Super-enhancer calls stay far below the guard",
       subtitle = sprintf("%% of stitched enhancers called SE; hard guard is %d%%\nlabels give the SE count per line",
                          SE_GUARD),
       x = "% of stitched enhancers called SE", y = NULL) +
  theme_project() +
  theme(panel.grid.major.y = element_blank())

# ---- D. peak-to-gene linking -------------------------------------------------
funnel <- data.frame(
  stage = c("candidate pairs", "retained links", "unique regions", "unique genes"),
  n = c(getval("candidate_pairs"), getval("retained_links"),
        getval("unique_regions"), getval("unique_genes")))
funnel$stage <- factor(funnel$stage, levels = rev(funnel$stage))

pD <- ggplot(funnel, aes(n, stage)) +
  geom_col(fill = "#2C5985", width = 0.62) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.12, size = 2.1,
            colour = "#1A1A1A") +
  scale_x_continuous(limits = c(0, max(funnel$n) * 1.30),
                     breaks = seq(0, 120000, 40000),
                     labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "D · Links are activity-correlated, not nearest-gene",
       subtitle = sprintf("only %.1f%% of links agree with the nearest gene;\n%.1fx enrichment over a distance-matched null",
                          getval("pct_agree_nearest_gene"), getval("null_enrichment")),
       x = "count", y = NULL) +
  theme_project() +
  theme(panel.grid.major.y = element_blank())

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Constructing the paralog-resolved regulatory layer, and the checks on each step",
    subtitle = sprintf(paste0(
      "Supplementary methods. hg38 deposits were lifted to hg19 with %.2f–%.2f%% loss across %d files; the ENCODE v2 blacklist removed\n",
      "%s of %s universe regions (%.1f%%). Panels A and C show guards added after an initial rule failed."),
      100*min(lo$loss_rate), 100*max(lo$loss_rate), nrow(lo),
      format(bl$value[bl$metric == "universe_removed"], big.mark = ","),
      format(bl$value[bl$metric == "universe_before"], big.mark = ","),
      bl$value[bl$metric == "pct_removed"]),
    caption = "Generated by scripts/07_figures/fig_s02_regulatory_pipeline.R · hg19 throughout · thresholds from config/params.yml",
    theme = theme_project() +
      theme(plot.title = element_text(size = FIG_BASE + 2, face = "bold"),
            plot.subtitle = element_text(size = FIG_BASE - 0.5, lineheight = 1.25)))

save_fig(fig, "fig_s02_regulatory_pipeline", width = 200, height = 175,
         script = "scripts/07_figures/fig_s02_regulatory_pipeline.R",
         caption = "Region construction: threshold calibration, convention stability, super-enhancer guard, peak-to-gene linking")

writeLines(c(
  "# Figure S2 — Constructing the regulatory layer", "",
  "**A.** Intrinsic validation of the ATAC threshold. TSS enrichment and H3K27ac",
  "fold both rise monotonically with the threshold, which is the evidence that the",
  "retained regions are regulatory signal rather than background. Because the",
  "relationship holds across the whole range rather than at one point, the chosen",
  sprintf("%sx cutoff is not load-bearing (D-024): conclusions do not depend on it.", CHOSEN_FOLD), "",
  "**B.** The MYCN-in-MYC nesting fraction behind gate criterion 1, across four",
  "thresholds and three replicate conventions. Within the pre-specified",
  sprintf("≥2-replicate rule the value varies by only %.3f across a 2.7-fold change in", spread_ge2),
  "threshold. The union and all-replicate conventions sit at materially different",
  "values, so the convention is a real choice and is shown rather than hidden.",
  "The circled point is the value reported in the gate.", "",
  "**C.** Super-enhancer calling. The first stitching rule — first slope ≥ 1 —",
  "returned 8,238 of 8,238 stitched enhancers as super-enhancers for H1048, i.e.",
  "all of them, which is not a super-enhancer call at all. It was replaced by a",
  sprintf("knee-point cutoff with a hard guard at %d%% of stitched enhancers. Every line", SE_GUARD),
  sprintf("now lands between %.1f%% and %.1f%%, far below the guard.",
          min(se$pct_of_stitched), max(se$pct_of_stitched)), "",
  "**D.** Peak-to-gene linking. Distance-weighted, activity-correlated linking",
  sprintf("retains %s of %s candidate pairs. Only %.1f%% of retained links agree with the",
          format(getval("retained_links"), big.mark = ","),
          format(getval("candidate_pairs"), big.mark = ","),
          getval("pct_agree_nearest_gene")),
  "nearest gene, so this is not nearest-gene assignment under another name, and the",
  sprintf("link set is enriched %.1fx over a distance-matched null.", getval("null_enrichment")),
  "An earlier version of this step assigned distal regions only, which left",
  "promoter-bound genes out of every regulon and understated the enrichment",
  "(D-029)."
), "figures/fig_s02_regulatory_pipeline_caption.md")
cat("wrote figures/fig_s02_regulatory_pipeline_caption.md\n")

# -----------------------------------------------------------------------------
# fig03_functional_domain.R — what functional evidence survived, and what it can
# and cannot attribute.
#
#   Rscript scripts/07_figures/fig03_functional_domain.R
#
# The project set out to do multi-layer paralog-specific prioritisation. Three of
# the four planned evidence domains did not survive contact with the data. This
# figure reports that outcome directly, because a prioritisation framework
# presented without its own attrition is a claim about evidence that the evidence
# does not support.
#
#   A  dependency enrichment — real when pooled, indistinguishable per paralog
#   B  none of MYC/MYCN/MYCL is itself an SCLC-selective dependency
#   C  the regulon programmes, including where they disagree with the source paper
#   D  MOES domain accounting: four planned, two admitted, one paralog-resolved
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ragg); library(colorspace); library(yaml)
})
source("R/theme_project.R")

enr    <- read.csv("data/metadata/m7_dependency_enrichment.csv", stringsAsFactors = FALSE)
pooled <- read.csv("data/metadata/m7_dependency_pooled.csv", stringsAsFactors = FALSE)
own    <- read.csv("data/metadata/m7_paralog_own_dependency.csv", stringsAsFactors = FALSE)
onto   <- read.csv("data/metadata/regulon_ontology.csv", stringsAsFactors = FALSE)
prog   <- read.csv("data/metadata/regulon_programme.csv", stringsAsFactors = FALSE)
moes   <- yaml::read_yaml("config/params.yml")$moes
DEP_THRESH <- yaml::read_yaml("config/params.yml")$depmap$selective_dependency_threshold

cat("accessibility checks\n")
check_palette(PAL_PARALOG, label = "paralog palette")
PAL_STATUS <- c("admitted" = "#2C5985", "excluded" = "#A63603")
check_palette(PAL_STATUS, label = "domain status palette")
check_text_size()
cat("\n")

# The figure's headline is the domain count. Read it from config rather than
# writing "two" into a title that config could later contradict.
n_admit <- length(moes$domains)
n_drop  <- length(moes$domains_dropped)
n_resolved <- sum(unlist(moes$domain_attribution) == "paralog_resolved")
cat("MOES domains: ", n_admit, " admitted, ", n_drop, " dropped, ",
    n_resolved, " paralog-resolved\n\n", sep = "")

# ---- A. dependency enrichment, per paralog vs pooled -------------------------
fa <- rbind(
  data.frame(label = enr$paralog, or = enr$or, lo = enr$ci_low, hi = enr$ci_high,
             n = enr$n_selective, kind = "per paralog", stringsAsFactors = FALSE),
  data.frame(label = "POOLED", or = pooled$or, lo = pooled$ci_low, hi = pooled$ci_high,
             n = pooled$selective_in, kind = "pooled", stringsAsFactors = FALSE))
fa$label <- factor(fa$label, levels = rev(c("MYC", "MYCN", "MYCL1", "POOLED")))
# as.character() is REQUIRED. Indexing a named vector with a factor uses the
# integer level codes, not the names: with levels reversed for the plot, MYC
# became index 4 (out of range -> NA -> the point was silently dropped) and
# MYCN/MYCL1 became indices 3/2, i.e. each drew the other's colour. Wrong
# colours, no error, and only a dropped-row warning to hint at it.
fa$col <- ifelse(fa$kind == "pooled", "#1A1A1A",
                 PAL_PARALOG[as.character(fa$label)])
stopifnot(!any(is.na(fa$col)))

pA <- ggplot(fa, aes(or, label)) +
  geom_vline(xintercept = 1, linewidth = 0.35, colour = "#999999") +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = I(col)), orientation = "y",
                width = 0.18, linewidth = 0.5) +
  geom_point(aes(colour = I(col), shape = kind), size = 2.2, fill = "white") +
  geom_text(aes(label = sprintf("%d gene%s", n, ifelse(n == 1, "", "s"))),
            y = as.numeric(fa$label) + 0.30, size = 1.9, colour = "#666666") +
  scale_shape_manual(values = c("per paralog" = 16, "pooled" = 23), guide = "none") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10)) +
  labs(title = "A · Real pooled, not separable per paralog",
       subtitle = paste0("SCLC-selective dependency among regulon genes\n",
                         "odds ratio, 95% CI, log scale"),
       x = "odds ratio", y = NULL) +
  theme_project() +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(4, 14, 4, 4))

# ---- B. the paralogs' own dependency -----------------------------------------
# DepMap names the third paralog MYCL; the project uses MYCL1 throughout.
own$paralog <- c(MYC = "MYC", MYCN = "MYCN", MYCL = "MYCL1")[own$gene]
own$paralog <- factor(own$paralog, levels = rev(c("MYC", "MYCN", "MYCL1")))
db <- rbind(
  data.frame(paralog = own$paralog, grp = "SCLC lines", eff = own$mean_sclc),
  data.frame(paralog = own$paralog, grp = "other lineages", eff = own$mean_other))
db$grp <- factor(db$grp, levels = c("other lineages", "SCLC lines"))

pB <- ggplot(db, aes(eff, paralog)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#999999") +
  geom_vline(xintercept = DEP_THRESH, linetype = "22", linewidth = 0.4, colour = "#666666") +
  geom_line(aes(group = paralog), linewidth = 0.5, colour = "#BBBBBB") +
  geom_point(aes(fill = grp), shape = 21, size = 2.4, colour = "#1A1A1A", stroke = 0.3) +
  # Anchored to the panel's left edge, not to the points. Positioning relative to
  # the leftmost point puts the label on top of the points whenever both sit near
  # zero, which is exactly the case for MYCN and MYCL1.
  geom_text(data = own, aes(x = -2.6, y = paralog, label = sprintf("%+.2f", delta)),
            hjust = 0, size = 2.1, colour = "#1A1A1A", inherit.aes = FALSE) +
  scale_fill_manual(values = c("other lineages" = "#CCCCCC", "SCLC lines" = "#2C5985"),
                    name = NULL) +
  scale_x_continuous(limits = c(-2.65, 0.35)) +
  labs(title = "B · No paralog is selectively essential",
       subtitle = paste0("mean CRISPR gene effect; dashed line = ", DEP_THRESH,
                         "\nlabels give SCLC minus other lineages"),
       x = "gene effect (Chronos)", y = NULL) +
  theme_project() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

# ---- C. regulon programmes ---------------------------------------------------
ob <- onto[onto$regulon == "both", ]
pc <- rbind(
  data.frame(paralog = ob$paralog, set = "Hallmark", or = ob$hallmark_or),
  data.frame(paralog = ob$paralog, set = "neurogenesis", or = ob$neuro_or))
pc$paralog <- factor(pc$paralog, levels = c("MYC", "MYCN", "MYCL1"))
# Mark the programme each regulon was ASSIGNED, from the gate that assigned it.
assigned <- setNames(ifelse(prog$pass, prog$programme, NA), prog$paralog)
pc$is_assigned <- mapply(function(p, s)
  !is.na(assigned[p]) && grepl(substr(assigned[p], 1, 4), s, ignore.case = TRUE),
  as.character(pc$paralog), as.character(pc$set))

pC <- ggplot(pc, aes(paralog, or, fill = paralog)) +
  geom_hline(yintercept = 1, linewidth = 0.35, colour = "#999999") +
  geom_col(aes(colour = is_assigned, group = paralog), width = 0.7, linewidth = 0.5) +
  facet_wrap(~set) +
  scale_colour_manual(values = c("TRUE" = "#1A1A1A", "FALSE" = NA), guide = "none") +
  scale_fill_manual(values = PAL_PARALOG, guide = "none") +
  scale_y_continuous(limits = c(0, 3.3), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "C · MYC is neurogenesis-weighted; MYCN is not",
       subtitle = "regulon-gene enrichment; outlined bar = assigned programme",
       y = "odds ratio", x = NULL) +
  theme_project() +
  theme(panel.grid.major.x = element_blank(), plot.margin = margin(4, 14, 4, 4))

# ---- D. domain accounting ----------------------------------------------------
dom <- data.frame(
  domain = c("cis-regulatory", "functional", "transcriptional", "network"),
  status = c("admitted", "admitted", "excluded", "excluded"),
  attribution = c("paralog-resolved", "gene-level only", "—", "—"),
  why = c("chromatin, per paralog",
          "dependency, ORs indistinguishable",
          "lineage-confounded (D-033)",
          "underpowered + confounded (D-037)"),
  stringsAsFactors = FALSE)
dom$domain <- factor(dom$domain, levels = rev(dom$domain))

# Column x-positions are spaced against the widest string each column carries.
# The first attempt sized the status tile to 0.85 data units while "ADMITTED" at
# 2 pt needs more than that, so the word overflowed its own tile and the next
# column started underneath it.
pD <- ggplot(dom, aes(y = domain)) +
  geom_tile(aes(x = 1.2, fill = status), width = 1.95, height = 0.78) +
  geom_text(aes(x = 1.2, label = toupper(status)), colour = "white", size = 1.75,
            fontface = "bold") +
  geom_text(aes(x = 2.45, label = attribution), hjust = 0, size = 2.0,
            fontface = ifelse(dom$attribution == "paralog-resolved", "bold", "plain"),
            colour = "#1A1A1A") +
  geom_text(aes(x = 5.9, label = why), hjust = 0, size = 1.9, colour = "#666666") +
  scale_fill_manual(values = PAL_STATUS, guide = "none") +
  scale_x_continuous(limits = c(0.15, 12.1)) +
  labs(title = sprintf("D · %d of 4 domains admitted, %d paralog-resolved",
                       n_admit, n_resolved),
       subtitle = "paralog attribution rests on chromatin alone") +
  theme_project() +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank())

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "Functional evidence supports the regulon genes collectively, not any one paralog",
    subtitle = paste0(
      "Dependency is not expression, so it escapes the lineage confound that removed the transcriptional domain — but on 4–5 selective genes\n",
      "per regulon it cannot separate the paralogs. Of four planned evidence domains, two are admitted and one attributes evidence to a paralog."),
    caption = "Generated by scripts/07_figures/fig03_functional_domain.R · DepMap Public 26Q1, 25 SCLC lines vs 1,183 other lineages (D-036, D-037)",
    theme = theme_project() +
      theme(plot.title = element_text(size = FIG_BASE + 2, face = "bold"),
            plot.subtitle = element_text(size = FIG_BASE - 0.5, lineheight = 1.25)))

save_fig(fig, "fig03_functional_domain", width = 200, height = 170,
         script = "scripts/07_figures/fig03_functional_domain.R",
         caption = "Dependency enrichment, paralog non-essentiality, regulon programmes, and MOES domain accounting")

writeLines(c(
  "# Figure 3 — Functional evidence and domain accounting", "",
  "**A.** SCLC-selective CRISPR dependency among paralog regulon genes, against",
  "the rest of the tested genome. Pooled across all three regulons the enrichment",
  sprintf("is real: OR %.2f (95%% CI %.2f–%.2f), p = %.4f, on %d selective genes in %s tested.",
          pooled$or, pooled$ci_low, pooled$ci_high, pooled$p,
          pooled$selective_in, format(pooled$n_regulon_genes, big.mark = ",")),
  "Per paralog it is not separable: ORs 2.45 / 2.21 / 3.09 on 4, 4 and 5 genes,",
  "with confidence intervals that overlap almost completely. MYCL1's nominal",
  "p = 0.028 comes with a 95% CI of 0.97–7.59 that includes 1 — Fisher's exact",
  "p and its conditional interval disagree at this count, which is itself a reason",
  "not to read a per-paralog result here. The pooled test is the powered question.", "",
  "**B.** The paralogs' own dependency. None is SCLC-selective. MYC is strongly",
  sprintf("essential in both groups but LESS so in SCLC (%+.3f, FDR %.2f); MYCN and MYCL",
          own$delta[own$gene == "MYC"], own$fdr[own$gene == "MYC"]),
  "sit near zero in every lineage. Nothing in this analysis supports targeting the",
  "paralogs themselves, and the figure is drawn so that cannot be misread.", "",
  "**C.** Regulon programme enrichment. The MYC regulon is neurogenesis-weighted",
  "(OR 1.85, p = 7e-6), independently reproducing the enhancer-to-neurogenesis link",
  "reported by Plotnik et al. The MYCN regulon is not (OR 1.09, p = 0.33); it is",
  "Hallmark/housekeeping-weighted instead (OR 2.99, p = 3.6e-4). This is a",
  "discordance with the source paper, which grouped MYCN with MYC. The MYCL1",
  "regulon reaches neither programme and was flagged rather than confirmed (D-030).", "",
  "**D.** MOES domain accounting. Four evidence domains were specified at M1; two",
  "are admitted. Transcriptional was dropped as lineage-confounded at both the",
  "aggregate and per-gene level (D-033). Network was excluded on both of its",
  "pre-conditions — 59 SCLC lines is below what tree-ensemble inference needs, and",
  "MYC expression tracks NE state across those lines (rho -0.441), the same confound",
  "that barred a tumour-based network (D-037). Of the two admitted domains only",
  "cis-regulatory attributes evidence to a paralog, so the multi-layer",
  "paralog-specific prioritisation specified at M1 is not achievable on this",
  "evidence, and is reported as such rather than presented at reduced strength."
), "figures/fig03_functional_domain_caption.md")
cat("wrote figures/fig03_functional_domain_caption.md\n")

# -----------------------------------------------------------------------------
# fig_s01_data_landscape.R — the data landscape, from verified metadata only.
#
#   Rscript scripts/07_figures/fig_s01_data_landscape.R
#
# THIS FIGURE CONTAINS NO RESULTS. Every value is metadata verified against GEO
# records at M4: which assay exists in which cell line, which genome build each
# deposit uses, and where the lineage-TF controls do and do not overlap the
# keystone. It is an inventory, and its caption says so.
#
# It exists because three findings from M4 are much easier to see than to read:
#   A. the keystone design, including its two gaps
#   B. the hg19/hg38 split that makes liftOver load-bearing (R-05)
#   C. the lineage-TF overlap shortfall that limits the R-01 analysis (R-14)
#
# Inputs : data/metadata/sample_design.tsv , data/metadata/dataset_manifest.csv
# Output : figures/fig_s01_data_landscape.png
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
})
source("R/theme_project.R")

# This figure was written at M4, before check_palette() existed, and so shipped
# without ever running the gate — which is how the failing M1 paralog palette
# survived into a committed figure (D-038). PAL_AMP is checked here specifically
# because it is NOT the paralog trio: it adds a fourth colour (grey) that the
# trio's own score says nothing about.
cat("accessibility checks\n")
check_palette(PAL_AMP, label = "amplification palette")
check_palette(PAL_PRESENT, label = "presence palette")
check_text_size()
cat("\n")

design <- read.delim("data/metadata/sample_design.tsv", stringsAsFactors = FALSE)
man    <- read.csv("data/metadata/dataset_manifest.csv", stringsAsFactors = FALSE)

# Keystone lines ordered by paralog group, so the block structure is visible.
AMP_ORDER <- c("MYC-amp", "MYCN-amp", "MYCL-amp", "non-amplified")
ks <- design[design$dataset == "GSE230649", ]
line_amp <- unique(ks[, c("cell_line", "paralog_status")])
line_amp$paralog_status <- factor(line_amp$paralog_status, levels = AMP_ORDER)
line_amp <- line_amp[order(line_amp$paralog_status, line_amp$cell_line), ]
LINE_LEVELS <- line_amp$cell_line

# =============================================================================
# Panel A — keystone assay x cell line
# =============================================================================
ASSAY_LEVELS <- c("MYC", "MYCN", "MYCL1", "H3K27ac", "ATAC")
gridA <- expand.grid(cell_line = LINE_LEVELS, assay = ASSAY_LEVELS,
                     stringsAsFactors = FALSE)
have  <- unique(ks[, c("cell_line", "assay")])
gridA$state <- ifelse(paste(gridA$cell_line, gridA$assay) %in%
                        paste(have$cell_line, have$assay), "present", "absent")
gridA$assay     <- factor(gridA$assay, levels = rev(ASSAY_LEVELS))
gridA$cell_line <- factor(gridA$cell_line, levels = LINE_LEVELS)

# Amplification status strip beneath the matrix.
stripA <- line_amp
stripA$cell_line <- factor(stripA$cell_line, levels = LINE_LEVELS)

# Cell-line labels sit on the amplification strip at the BOTTOM of the panel.
# Putting them on top collided with the neighbouring panel's title.
pA_mat <- ggplot(gridA, aes(cell_line, assay, fill = state)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  scale_fill_manual(values = PAL_PRESENT, guide = "none") +
  labs(title = "A · Keystone GSE230649 (hg19, 28 samples)",
       subtitle = "Grey = not deposited") +
  theme_matrix() +
  theme(axis.text.x = element_blank())

pA_strip <- ggplot(stripA, aes(cell_line, y = 1, fill = paralog_status)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  scale_fill_manual(values = PAL_AMP, name = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_matrix() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        legend.position = "bottom")

panelA <- pA_mat / pA_strip + plot_layout(heights = c(5, 1.1))

# =============================================================================
# Panel B — genome build by dataset
# =============================================================================
fetchable <- man[!man$status %in% c("api", "manual_required", "pending_triage"), ]
bld <- as.data.frame(table(dataset = fetchable$dataset_id,
                           build   = fetchable$genome_build))
bld <- bld[bld$Freq > 0, ]
# Collapse GRCh38 label variants for display; they are the same assembly.
bld$build_grp <- ifelse(grepl("^hg19$", bld$build), "hg19 (project build)",
                 ifelse(grepl("hg38|GRCh38", bld$build), "hg38 / GRCh38 — needs liftOver",
                        "reference resource"))
bld$dataset <- reorder(bld$dataset, bld$Freq, FUN = sum)

n_lift <- sum(fetchable$liftover_required)

pB <- ggplot(bld, aes(Freq, dataset, fill = build_grp)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(
    "hg19 (project build)"           = "#2C5985",
    "hg38 / GRCh38 — needs liftOver" = "#C1662F",
    "reference resource"             = "#8C8C8C"), name = NULL) +
  labs(title = "B · Genome build is split across deposits",
       subtitle = paste0(n_lift, " of ", nrow(fetchable),
                         " files are not in the project build"),
       x = "files", y = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_project() +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank())

# =============================================================================
# Panel C — lineage-TF coverage against the keystone (risk R-14)
# =============================================================================
myc_lines <- sort(unique(ks$cell_line[ks$assay %in% c("MYC", "MYCN", "MYCL1")]))
lt <- design[design$role == "lineage_tf", ]
lt_lines <- sort(unique(lt$cell_line))
extra <- setdiff(lt_lines, myc_lines)          # lineage-TF lines outside the keystone

rowsC <- c(myc_lines, extra)
COLS  <- c("MYC-family ChIP", "ASCL1", "NEUROD1", "POU2F3")
gridC <- expand.grid(cell_line = rowsC, track = COLS, stringsAsFactors = FALSE)

gridC$state <- "absent"
gridC$state[gridC$track == "MYC-family ChIP" & gridC$cell_line %in% myc_lines] <- "present"
for (tf in c("ASCL1", "NEUROD1", "POU2F3")) {
  ln <- unique(lt$cell_line[lt$assay == tf])
  gridC$state[gridC$track == tf & gridC$cell_line %in% ln] <- "present"
}

# Order rows by paralog group, keystone first, then the non-keystone lines.
amp_lookup <- setNames(line_amp$paralog_status, line_amp$cell_line)
ord <- data.frame(cell_line = rowsC, stringsAsFactors = FALSE)
ord$keystone <- ord$cell_line %in% myc_lines
ord$amp <- factor(as.character(amp_lookup[ord$cell_line]), levels = AMP_ORDER)
ord <- ord[order(!ord$keystone, ord$amp, ord$cell_line), ]
gridC$cell_line <- factor(gridC$cell_line, levels = rev(ord$cell_line))
gridC$track     <- factor(gridC$track, levels = COLS)

# Label non-keystone lines so it is obvious they cannot support a within-line test.
lab <- setNames(ord$cell_line, ord$cell_line)
lab[!ord$keystone] <- paste0(ord$cell_line[!ord$keystone], " †")

pC <- ggplot(gridC, aes(track, cell_line, fill = state)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  scale_fill_manual(values = PAL_PRESENT, guide = "none") +
  scale_y_discrete(labels = lab) +
  scale_x_discrete(position = "top") +
  annotate("segment", x = 1.5, xend = 1.5,
           y = 0.5, yend = length(levels(gridC$cell_line)) + 0.5,
           colour = "#1A1A1A", linewidth = 0.4) +
  labs(title = "C · Lineage-TF controls barely overlap the keystone lines",
       subtitle = "Rows: lines with MYC-family ChIP.  † = lineage-TF line absent from the keystone.",
       caption = paste0(
         "The confounding test central to risk R-01 asks whether a MYC-bound region is also lineage-TF-bound IN THE SAME CELLS.\n",
         "Available for POU2F3 (H1048, H211, H526 — two paralog groups, already hg19); marginal for ASCL1 (SHP77 only, n=1, hg38);\n",
         "unavailable for NEUROD1, whose only line H446 is not in the keystone. Resolution is therefore heterogeneous and any\n",
         "claim of lineage-independence must be qualified per transcription factor (risk R-14).")) +
  theme_matrix() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# =============================================================================
# assemble
# =============================================================================
fig <- (panelA | pB) / pC +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title = "Data landscape for a paralog-aware analysis of MYC regulatory hubs in SCLC",
    subtitle = paste0(
      "Inventory of verified public data as acquired at M4 — metadata only, no analytical results.\n",
      "Assay and cell-line assignments are parsed from deposited file names; genome builds are read from each\n",
      "series' own data-processing declaration; amplification status is author-declared (Plotnik et al. 2024)."),
    caption = paste0("Generated ", format(Sys.Date(), "%Y-%m-%d"),
                     " by scripts/07_figures/fig_s01_data_landscape.R"),
    theme = theme_project() +
      theme(plot.title    = element_text(size = FIG_BASE + 3, face = "bold"),
            plot.subtitle = element_text(size = FIG_BASE - 0.5, lineheight = 1.25),
            plot.margin   = margin(6, 8, 4, 6))
  )

save_fig(fig, "fig_s01_data_landscape.png", width = 200, height = 190)

# --- caption file, so the figure is interpretable without the script ----------
writeLines(c(
  "# Figure S1 — Data landscape",
  "",
  "**This figure reports metadata, not results.** No analytical output exists at M4.",
  "",
  paste0("**A.** Composition of the keystone deposit GSE230649 (hg19, 28 samples): ",
         "MYC ChIP in 5 lines, MYCN in 2, MYCL1 in 2, H3K27ac in 10, ATAC in 9. ",
         "Two gaps are carried forward rather than smoothed over — H211 has MYC ChIP ",
         "but no ATAC, and H196 has accessibility and H3K27ac but no MYC-family ChIP. ",
         "The strip beneath encodes author-declared paralog amplification status."),
  "",
  paste0("**B.** Files by genome build. The project build is hg19, fixed by the ",
         "keystone, but all three supporting ATAC deposits and the ASCL1 ChIP are ",
         "hg38 — so the consensus accessible-region universe, which requires support ",
         "from at least two independent ATAC datasets, cannot be built without ",
         "crossing builds. Continuous signal is never lifted; intervals are called ",
         "in the native build and only intervals are lifted, with loss rate reported."),
  "",
  paste0("**C.** Cell-line overlap between the keystone and the lineage-TF controls. ",
         "This is the constraint with the clearest scientific consequence: the ",
         "confounding test central to risk R-01 requires lineage-TF occupancy in the ",
         "same lines where MYC-family occupancy was measured. POU2F3 provides three ",
         "such lines across two paralog groups and is already in hg19; ASCL1 provides ",
         "one; NEUROD1 provides none. Conclusions about lineage-independence are ",
         "therefore qualified per transcription factor rather than stated uniformly."),
  "",
  "Sources: GEO series and sample records for GSE230649, GSE269424, GSE256345,",
  "GSE281523, GSE281524, GSE210113, GSE249362 (verified 2026-07-26); amplification",
  "status from Plotnik et al. 2024, Mol Cancer Res (PMID 38747975)."
), "figures/fig_s01_data_landscape_caption.md")
cat("wrote figures/fig_s01_data_landscape_caption.md\n")

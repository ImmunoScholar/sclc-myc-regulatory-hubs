# -----------------------------------------------------------------------------
# theme_project.R — one figure style for the whole project.
#
#   source("R/theme_project.R")
#
# Every figure parameter comes from config/params.yml, so changing figure sizing
# or the paralog palette is a config edit, not a hunt through plotting code.
# Established at M4 rather than M10 so the final figures are not a retrofit.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(yaml)
})

.fig_cfg <- yaml::read_yaml("config/params.yml")$figures

# Millimetres, as configured: 90 mm single-column, 180 mm double.
FIG_W_SINGLE <- .fig_cfg$width_single
FIG_W_DOUBLE <- .fig_cfg$width_double
FIG_DPI      <- .fig_cfg$dpi
FIG_BASE     <- .fig_cfg$base_font_size

# Paralog colours are fixed project-wide: the same gene must never change colour
# between figures, or readers will compare panels wrongly.
PAL_PARALOG <- unlist(.fig_cfg$palette_paralog)

# Amplification-status palette, derived from the paralog colours so the two
# encodings agree. Grey for non-amplified: it is a reference group, not a fourth
# category competing for attention.
PAL_AMP <- c(
  "MYC-amp"       = PAL_PARALOG[["MYC"]],
  "MYCN-amp"      = PAL_PARALOG[["MYCN"]],
  "MYCL-amp"      = PAL_PARALOG[["MYCL1"]],
  "non-amplified" = "#8C8C8C"
)

# Presence/absence fill used across the inventory panels. Deliberately low
# chroma: these panels are about structure, and the paralog colours carry
# the meaning.
PAL_PRESENT <- c("present" = "#2C5985", "absent" = "#EDEDED")

theme_project <- function(base_size = FIG_BASE) {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(colour = "#1A1A1A"),
      plot.title       = element_text(size = base_size + 1.5, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size - 0.5, colour = "#4D4D4D", hjust = 0),
      plot.caption     = element_text(size = base_size - 1.5, colour = "#666666", hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size - 0.5),
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(3.2, "mm"),
      strip.text       = element_text(size = base_size, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "#E8E8E8"),
      plot.margin      = margin(4, 4, 4, 4)
    )
}

# Matrix/heatmap panels: gridlines are noise when every cell is a tile.
theme_matrix <- function(base_size = FIG_BASE) {
  theme_project(base_size) +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title  = element_blank()
    )
}

# Single save path so DPI, device and units cannot drift between figures.
save_fig <- function(plot, filename, width = FIG_W_DOUBLE, height = 100) {
  dir.create("figures", showWarnings = FALSE)
  path <- file.path("figures", filename)
  ggplot2::ggsave(path, plot, width = width, height = height, units = "mm",
                  dpi = FIG_DPI, device = ragg::agg_png, bg = "white")
  cat("wrote ", path, "  (", width, "x", height, " mm @ ", FIG_DPI, " dpi)\n", sep = "")
  invisible(path)
}

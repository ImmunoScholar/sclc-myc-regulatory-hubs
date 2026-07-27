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

# Diverging scale for signed continuous values (correlation heatmaps). Config-
# driven for the same reason as the paralog colours: a hex literal typed into a
# plotting script is invisible to the accessibility check and outlives the palette
# it was copied from.
PAL_DIVERGING <- unlist(.fig_cfg$palette_diverging)

scale_fill_diverging <- function(limits, name = NULL, ...) {
  ggplot2::scale_fill_gradient2(
    low = PAL_DIVERGING[["low"]], mid = PAL_DIVERGING[["mid"]],
    high = PAL_DIVERGING[["high"]], midpoint = 0, limits = limits, name = name, ...)
}

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

# --- palette accessibility, CHECKED not asserted -------------------------------
# config figures$colourblind_check requires deutan/protan/tritan simulation with a
# minimum perceptual distance. "Colourblind safe" is a testable claim and an
# untested one is worth nothing, so this runs rather than being declared.
check_palette <- function(cols, min_dist = .fig_cfg$colourblind_check$min_perceptual_distance,
                          label = "palette") {
  if (!requireNamespace("colorspace", quietly = TRUE)) {
    warning("colorspace unavailable — accessibility NOT checked"); return(invisible(NA))
  }
  sims <- list(normal = identity,
               deutan = colorspace::deutan,
               protan = colorspace::protan,
               tritan = colorspace::tritan)
  worst <- Inf; worst_where <- ""
  for (nm in names(sims)) {
    sc <- sims[[nm]](unname(cols))
    rgb <- t(grDevices::col2rgb(sc))
    lab <- grDevices::convertColor(rgb / 255, from = "sRGB", to = "Lab")
    d <- as.matrix(stats::dist(lab))
    diag(d) <- Inf
    if (min(d) < worst) { worst <- min(d); worst_where <- nm }
  }
  ok <- worst >= min_dist
  cat(sprintf("  palette check [%s]: min CIE dE = %.1f under %s (need >= %s) -> %s\n",
              label, worst, worst_where, min_dist, if (ok) "PASS" else "FAIL"))
  if (!ok) warning(sprintf("%s fails colourblind check (dE %.1f < %s)", label, worst, min_dist))
  invisible(ok)
}

# --- legibility, CHECKED at final print size -----------------------------------
# config figures$min_text_pt. A figure shipped illegible is a failure, so the
# smallest theme text is verified against the configured floor.
check_text_size <- function(base_size = FIG_BASE, min_pt = .fig_cfg$min_text_pt) {
  smallest <- base_size - 1.5      # the smallest size used in theme_project()
  ok <- smallest >= min_pt
  cat(sprintf("  text check: smallest element %.1f pt (floor %s) -> %s\n",
              smallest, min_pt, if (ok) "PASS" else "FAIL"))
  if (!ok) warning("figure text falls below the configured minimum")
  invisible(ok)
}

# --- single save path ----------------------------------------------------------
# Emits PNG *and* vector PDF (config figures$vector_output): a raster-only figure
# cannot be rescaled for print. Records every figure in a manifest so none is
# untraceable (config figures$figure_manifest).
save_fig <- function(plot, filename, width = FIG_W_DOUBLE, height = 100,
                     script = NA_character_, caption = NA_character_) {
  dir.create("figures", showWarnings = FALSE)
  stem <- sub("\\.(png|pdf)$", "", filename)
  png_path <- file.path("figures", paste0(stem, ".png"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, units = "mm",
                  dpi = FIG_DPI, device = ragg::agg_png, bg = "white")
  cat("wrote ", png_path, "  (", width, "x", height, " mm @ ", FIG_DPI, " dpi)\n", sep = "")

  pdf_path <- NA_character_
  if (isTRUE(.fig_cfg$vector_output)) {
    pdf_path <- file.path("figures", paste0(stem, ".pdf"))
    ok <- tryCatch({
      ggplot2::ggsave(pdf_path, plot, width = width, height = height, units = "mm",
                      device = grDevices::cairo_pdf, bg = "white"); TRUE
    }, error = function(e) { warning("PDF export failed: ", conditionMessage(e)); FALSE })
    if (ok) cat("wrote ", pdf_path, "  (vector)\n", sep = "")
  }

  if (isTRUE(.fig_cfg$figure_manifest)) {
    mf <- "figures/figure_manifest.csv"
    row <- data.frame(figure = stem, png = basename(png_path),
                      pdf = ifelse(is.na(pdf_path), NA, basename(pdf_path)),
                      width_mm = width, height_mm = height, dpi = FIG_DPI,
                      script = script, caption = caption,
                      generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                      stringsAsFactors = FALSE)
    if (file.exists(mf)) {
      old <- utils::read.csv(mf, stringsAsFactors = FALSE)
      old <- old[old$figure != stem, , drop = FALSE]
      row <- rbind(old, row)
    }
    utils::write.csv(row, mf, row.names = FALSE)
  }
  invisible(png_path)
}

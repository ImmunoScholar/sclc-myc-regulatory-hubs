# -----------------------------------------------------------------------------
# 00_select_palette.R — choose a paralog palette that passes the accessibility check.
#
#   Rscript scripts/07_figures/00_select_palette.R
#
# The palette frozen in config since M1 — MYC #B2182B, MYCN #2166AC,
# MYCL1 #1B7837 — FAILS the project's own colourblind requirement: min CIE dE
# 12.7 under deuteranopia against a floor of 15. Red vs green, the most common
# confusion, in a palette labelled "colourblind_safe: true". Figure S1 already
# shipped with it.
#
# Rather than guess a replacement, candidates are SCORED with the same function
# that caught the failure, and the winner is the one with the largest worst-case
# perceptual separation across normal, deutan, protan and tritan vision.
#
# Candidates are drawn from established colourblind-safe families (Okabe-Ito;
# ColorBrewer diverging endpoints) rather than invented.
#
# TWO criteria, not one. Maximising separation-between-colours alone selects pale
# tints, because CIE distance rewards lightness spread: the best scorer on that
# metric is a near-white gold. The paralogs are drawn as thin lines, 1.7 pt points
# and coloured text, all on a white panel, so a colour must ALSO stand off the
# background. A category the reader cannot see is not accessible whatever its
# pairwise dE says. Both thresholds are hard filters; dE only ranks the survivors.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(colorspace) })

to_lab <- function(cols) {
  grDevices::convertColor(t(grDevices::col2rgb(unname(cols))) / 255,
                          from = "sRGB", to = "Lab")
}

worst_dE <- function(cols) {
  sims <- list(normal = identity, deutan = colorspace::deutan,
               protan = colorspace::protan, tritan = colorspace::tritan)
  w <- Inf; where <- ""
  for (nm in names(sims)) {
    d <- as.matrix(stats::dist(to_lab(sims[[nm]](unname(cols))))); diag(d) <- Inf
    if (min(d) < w) { w <- min(d); where <- nm }
  }
  list(dE = w, where = where)
}

# Lightness contrast against the white panel. Luminance, not hue, is what carries
# a thin mark; L* is on 0-100, so 100 - L* is the headroom each colour has.
min_bg_contrast <- function(cols) min(100 - to_lab(cols)[, "L"])

candidates <- list(
  "current (config, M1)"      = c(MYC="#B2182B", MYCN="#2166AC", MYCL1="#1B7837"),
  "Okabe-Ito verm/blue/orange"= c(MYC="#D55E00", MYCN="#0072B2", MYCL1="#E69F00"),
  "Okabe-Ito verm/blue/purple"= c(MYC="#D55E00", MYCN="#0072B2", MYCL1="#CC79A7"),
  "Okabe-Ito blue/orange/green"=c(MYC="#0072B2", MYCN="#E69F00", MYCL1="#009E73"),
  "red/blue/purple"           = c(MYC="#B2182B", MYCN="#2166AC", MYCL1="#762A83"),
  "dark red/mid blue/gold"    = c(MYC="#A50026", MYCN="#4575B4", MYCL1="#FDAE61"),
  "verm/navy/light-orange"    = c(MYC="#D55E00", MYCN="#08519C", MYCL1="#FDBF6F"),
  "purple/orange/teal"        = c(MYC="#762A83", MYCN="#E08214", MYCL1="#1B7837"),
  # conventional assignment (red = MYC) given a fair hearing, with the third
  # colour moved off green so it does not collide with red under deuteranopia
  "dark red/blue/dark orange" = c(MYC="#A50026", MYCN="#2166AC", MYCL1="#B35806"),
  "dark red/blue/teal"        = c(MYC="#A50026", MYCN="#2166AC", MYCL1="#01665E"),
  "crimson/teal/purple"       = c(MYC="#B2182B", MYCN="#01665E", MYCL1="#762A83")
)

MIN    <- 15   # config figures$colourblind_check$min_perceptual_distance
MIN_BG <- 30   # L* headroom vs the white panel, so thin marks stay visible

cat("=========== candidate paralog palettes ===========\n")
cat("required: min CIE dE >= ", MIN, " across normal/deutan/protan/tritan\n", sep = "")
cat("required: every colour >= ", MIN_BG, " L* from white\n\n", sep = "")
res <- data.frame()
for (nm in names(candidates)) {
  w  <- worst_dE(candidates[[nm]])
  bg <- min_bg_contrast(candidates[[nm]])
  pass_cb <- w$dE >= MIN; pass_bg <- bg >= MIN_BG
  res <- rbind(res, data.frame(palette = nm, worst_dE = round(w$dE, 1),
                               limiting = w$where, bg_contrast = round(bg, 1),
                               pass_colourblind = pass_cb, pass_background = pass_bg,
                               pass = pass_cb && pass_bg,
                               cols = paste(candidates[[nm]], collapse = " "),
                               stringsAsFactors = FALSE))
  cat(sprintf("  %-28s dE %5.1f (%-6s) %s   bg %5.1f %s\n", nm, w$dE, w$where,
              if (pass_cb) "PASS" else "FAIL", bg, if (pass_bg) "PASS" else "FAIL"))
}

ok <- res[res$pass, ]
cat("\n=========== selection ===========\n")
if (!nrow(ok)) {
  cat("NO candidate passes both criteria. Widen the set before shipping any figure.\n")
  quit(status = 1)
}
cat(nrow(ok), " of ", nrow(res), " candidates pass both; ranking survivors by dE.\n", sep = "")
best <- ok[which.max(ok$worst_dE), ]
cat("selected: ", best$palette, "  (worst dE ", best$worst_dE,
    " under ", best$limiting, "; background contrast ", best$bg_contrast, ")\n", sep = "")
cat("colours : ", best$cols, "\n\n", sep = "")

sel <- candidates[[best$palette]]
cat("config/params.yml figures$palette_paralog should read:\n")
for (p in names(sel)) cat(sprintf("    %-6s \"%s\"\n", paste0(p, ":"), sel[[p]]))

write.csv(res, "data/metadata/palette_selection.csv", row.names = FALSE)
cat("\nwrote data/metadata/palette_selection.csv\n")

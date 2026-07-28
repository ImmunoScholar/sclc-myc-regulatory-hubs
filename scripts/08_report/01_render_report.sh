#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 01_render_report.sh — render report.qmd to HTML.
#
#   bash scripts/08_report/01_render_report.sh
#
# Prefers the Quarto CLI. Falls back to knitr + pandoc, which are already pinned
# in renv.lock and present on this machine, so the report renders WITHOUT
# requiring a Quarto install. The fallback is not a silent substitution: it says
# which path it took, and the two differ in presentation only (Quarto adds the
# folded-code UI and its own theme), not in content or numbers.
#
# Every figure and every value in the report is read from files the analysis
# wrote, so a stale analysis cannot produce a fresh-looking report.
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

SRC=report.qmd
OUT=report.html
[ -f "$SRC" ] || { echo "missing $SRC"; exit 1; }

# The report reads from data/metadata and figures/. If those are stale or absent
# the render must fail rather than emit a plausible-looking document.
need=(data/metadata/m5_gate_results.csv data/metadata/moes_ranking.csv
      data/metadata/m9_panel_coverage.csv figures/figure_manifest.csv)
missing=0
for f in "${need[@]}"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; missing=1; }
done
[ $missing -eq 0 ] || { echo "refusing to render against missing analysis outputs"; exit 1; }

if command -v quarto >/dev/null 2>&1; then
  echo "=== rendering with Quarto $(quarto --version) ==="
  quarto render "$SRC" --to html
  rc=$?
else
  echo "=== Quarto CLI not found — rendering with knitr + pandoc ==="
  echo "    (content identical; Quarto adds folded-code UI and its own theme)"
  Rscript -e '
    suppressPackageStartupMessages(library(knitr))
    opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
    knitr::knit("report.qmd", output = ".report_knit.md", quiet = TRUE)
  ' || { echo "knit failed"; exit 1; }
  pandoc .report_knit.md -o "$OUT" \
      --standalone --toc --toc-depth=3 --number-sections \
      --metadata title="Paralog-resolved MYC regulatory programmes in small cell lung cancer" \
      --highlight-style=tango
  rc=$?
  rm -f .report_knit.md
fi

if [ $rc -eq 0 ] && [ -s "$OUT" ]; then
  echo
  echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
else
  echo "RENDER FAILED (exit $rc)"; exit 1
fi

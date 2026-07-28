#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run_all.sh — run the pipeline in order.
#
#   bash scripts/run_all.sh              # analysis stages (assumes data present)
#   bash scripts/run_all.sh --with-data  # including the ~6 h download
#   bash scripts/run_all.sh --figures    # figures and report only
#
# Run from the repository root.
#
# STAGES ARE ORDERED BY DEPENDENCY, not by preference. M5 -> M6 -> M8 is the
# critical path; M8 must not run before M6, because the lineage result changes
# what "a MYC-associated hub" means and integrating evidence first would bake a
# confounder into the framework.
#
# This stops at the first failing stage. A pipeline that continues past a failed
# stage produces outputs that look complete and are not.
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

MODE="${1:---analysis}"
started=$(date +%s)
run() {
  local label="$1"; shift
  printf '\n=============================================================\n'
  printf '>>> %s\n' "$label"
  printf '=============================================================\n'
  if ! "$@"; then
    printf '\nSTAGE FAILED: %s\n' "$label"
    printf 'Stopping. Later stages depend on this one.\n'
    exit 1
  fi
}

if [ "$MODE" = "--with-data" ]; then
  echo "NOTE: the download stage fetches ~12 GB and takes roughly 6 hours."
  run "M4 · build dataset manifest"  Rscript scripts/01_data/01_build_manifest.R
  run "M4 · download"                bash    scripts/01_data/02_download.sh
  run "M4 · verify checksums"        Rscript scripts/01_data/03_verify.R
  run "M4 · QC report"               Rscript scripts/01_data/04_qc_report.R
fi

if [ "$MODE" != "--figures" ]; then
  # M5 regulatory layer. The numbered scripts in 02_regulatory run in sequence;
  # each writes its own QC and several are gated.
  for f in scripts/02_regulatory/*.R; do
    run "M5 · $(basename "$f")" Rscript "$f"
  done

  for f in scripts/03_tumour/*.R;    do run "M6 · $(basename "$f")" Rscript "$f"; done
  for f in scripts/04_functional/*.R; do
    case "$(basename "$f")" in 01_inspect_depmap.sh) continue;; esac
    run "M7 · $(basename "$f")" Rscript "$f"
  done
  for f in scripts/05_integration/*.R; do run "M8 · $(basename "$f")" Rscript "$f"; done
  for f in scripts/06_spatial/*.R;     do run "M9 · $(basename "$f")" Rscript "$f"; done
fi

run "M10 · figures"        bash scripts/07_figures/99_build_all.sh
run "M10 · output lint"    Rscript scripts/08_report/02_lint_outputs.R
run "M10 · report"         bash scripts/08_report/01_render_report.sh
run "audit · decisions"    bash scripts/00_setup/audit_decisions.sh

elapsed=$(( $(date +%s) - started ))
printf '\n=============================================================\n'
printf 'PIPELINE COMPLETE in %dh %dm\n' $((elapsed/3600)) $(((elapsed%3600)/60))
printf '=============================================================\n'

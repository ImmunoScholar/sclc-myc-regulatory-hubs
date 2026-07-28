#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 99_build_all.sh — regenerate every figure, in order, and verify the result.
#
#   bash scripts/07_figures/99_build_all.sh
#
# Run from the repository root. This is the entry point the M11 clean-clone check
# uses: a figure set that cannot be rebuilt from a fresh clone is not reproducible
# regardless of what the manifest says.
#
# It FAILS LOUDLY. Any script that errors, any accessibility check that reports
# FAIL, and any figure listed in the manifest without a file on disk stops the
# run with a non-zero exit. A build script that reports success while a check
# failed is worse than no build script.
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

FIGS=(
  00_select_palette.R
  fig_s01_data_landscape.R
  fig01_lineage_dominance.R
  fig02_m5_gate.R
  fig03_functional_domain.R
  fig_s02_regulatory_pipeline.R
  fig04_moes_null.R
  fig_s03_lineage_confound.R
)

rc=0
fail_checks=0
log=$(mktemp)

for f in "${FIGS[@]}"; do
  printf '=== %s ===\n' "$f"
  if ! Rscript "scripts/07_figures/$f" 2>&1 | grep -v out-of-sync | tee -a "$log"; then
    echo "  ERROR: $f exited non-zero"
    rc=1
  fi
  echo
done

# The accessibility checks print PASS/FAIL rather than stopping, so that a run
# reports every problem instead of only the first. This turns them into a gate.
if grep -q -- '-> FAIL' "$log"; then
  echo "=============================================================="
  echo "ACCESSIBILITY CHECKS FAILED:"
  grep -- '-> FAIL' "$log"
  echo "=============================================================="
  fail_checks=1
fi

echo "=== manifest vs disk ==="
missing=0
while IFS=, read -r fig png pdf rest; do
  fig=${fig//\"/}; png=${png//\"/}; pdf=${pdf//\"/}
  [ "$fig" = "figure" ] && continue
  for want in "figures/$png" "figures/$pdf"; do
    if [ ! -s "$want" ]; then echo "  MISSING: $want"; missing=1; fi
  done
  printf '  ok  %-30s %s + %s\n' "$fig" "$png" "$pdf"
done < figures/figure_manifest.csv

n_fig=$(( $(wc -l < figures/figure_manifest.csv) - 1 ))
echo
echo "=============================================================="
echo "figures in manifest : $n_fig"
echo "script errors       : $rc"
echo "failed checks       : $fail_checks"
echo "missing files       : $missing"
if [ $rc -eq 0 ] && [ $fail_checks -eq 0 ] && [ $missing -eq 0 ]; then
  echo "RESULT: all figures rebuilt and verified"
else
  echo "RESULT: FAILED"
fi
echo "=============================================================="
rm -f "$log"
[ $rc -eq 0 ] && [ $fail_checks -eq 0 ] && [ $missing -eq 0 ]

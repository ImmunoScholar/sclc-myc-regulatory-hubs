#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_builds2.sh — genome build + assay per series, using REAL sample IDs.
#
# Supersedes the representative-sample list in probe_builds.sh, where three IDs
# were guessed rather than read from the series records and consequently
# returned unrelated samples (mouse liver, a kidney biopsy, MCF7). Sample IDs
# are now extracted from each series' own !Series_sample_id lines. Never guess
# an accession.
#
# Read-only. Usage: bash scripts/01_data/probe_builds2.sh
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

OUT=data/metadata/geo_soft
mkdir -p "$OUT"

gsm_dir () { local g="$1"; echo "${g:0:$(( ${#g} - 3 ))}nnn"; }

for acc in GSE249362 GSE261348 GSE261345 GSE230649; do
  sf="$OUT/${acc}.soft.txt"
  [ -s "$sf" ] || { echo "missing $sf"; continue; }

  echo "============================================================="
  echo "[$acc]  total samples: $(grep -c '^!Series_sample_id' "$sf")"

  # Take the first three real sample IDs from the series record itself.
  ids=$(grep '^!Series_sample_id' "$sf" | sed 's/.*= *//' | tr -d '\r' | head -3)
  for gsm in $ids; do
    f="$OUT/${gsm}.soft.txt"
    if [ ! -s "$f" ]; then
      curl -4 -sS --max-time 90 --retry 3 --retry-delay 8 \
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm}&targ=self&form=text&view=quick" \
        -o "$f"
      sleep 3
    fi
    echo "  -- $gsm"
    grep -m1 '^!Sample_title'            "$f" | cut -c1-110
    grep -m1 '^!Sample_library_strategy' "$f"
    grep -iE '^!Sample_characteristics_ch1.*(antibody|cell line|segment|roi)' "$f" | head -2 | cut -c1-110
    grep -iE '^!Sample_data_processing.*(Assembly|Genome_build)' "$f" | head -2 | cut -c1-90
    grep -m1 '^!Sample_supplementary_file' "$f" | cut -c1-140
  done
done

echo
echo "============================================================="
echo "=== GSE249362: how many samples are actually POU2F3/ncBAF ChIP? ==="
sf="$OUT/GSE249362.soft.txt"
grep '^!Series_sample_id' "$sf" | sed 's/.*= *//' | tr -d '\r' | wc -l
echo "(full per-sample triage deferred: 125 samples, needs the series matrix)"

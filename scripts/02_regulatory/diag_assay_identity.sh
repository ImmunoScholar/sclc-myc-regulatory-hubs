#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# diag_assay_identity.sh — WHAT ARE THE NINE UNLABELLED GSE230649 TRACKS?
#
# They were identified as ATAC by ABSENCE of a target token in the file name.
# That was never verified against the GEO records. Given GSE269424 was labelled
# "ATAC-seq" and turned out to be a TF-overexpression experiment, an unverified
# assay assignment is not acceptable as the foundation of Aim 1.
#
# Read-only. Usage: bash scripts/02_regulatory/diag_assay_identity.sh
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

OUT=data/metadata/geo_soft
mkdir -p "$OUT"

# The nine tracks with no target token in the file name.
GSMS="GSM7230512 GSM7230513 GSM7230514 GSM7230515 GSM7230516 GSM7230517 GSM7230518 GSM7230519 GSM7230520"

for gsm in $GSMS; do
  f="$OUT/${gsm}.soft.txt"
  if [ ! -s "$f" ]; then
    curl -4 -sS --max-time 90 --retry 3 --retry-delay 6 \
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm}&targ=self&form=text&view=quick" \
      -o "$f"
    sleep 2
  fi
done

echo "=============================================================="
echo "THE NINE UNLABELLED TRACKS"
echo "=============================================================="
for gsm in $GSMS; do
  f="$OUT/${gsm}.soft.txt"
  [ -s "$f" ] || { echo "[$gsm] FETCH FAILED"; continue; }
  printf '\n[%s]\n' "$gsm"
  grep -m1 '^!Sample_title'            "$f" | cut -c1-100
  grep -m1 '^!Sample_library_strategy' "$f"
  grep -m1 '^!Sample_library_selection' "$f"
  grep -iE '^!Sample_characteristics_ch1' "$f" | head -4 | cut -c1-100
  grep -m1 '^!Sample_data_processing.*[Aa]ssembly' "$f" | cut -c1-70
done

echo
echo "=============================================================="
echo "CONTRAST: a known ChIP track from the same series"
echo "=============================================================="
f="$OUT/GSM7230493.soft.txt"
if [ -s "$f" ]; then
  grep -m1 '^!Sample_title'            "$f"
  grep -m1 '^!Sample_library_strategy' "$f"
  grep -m1 '^!Sample_library_selection' "$f"
  grep -iE '^!Sample_characteristics_ch1' "$f" | head -3
fi

echo
echo "=============================================================="
echo "SERIES-LEVEL: how many samples of each library strategy?"
echo "=============================================================="
grep -h '^!Sample_library_strategy' "$OUT"/GSM7230*.soft.txt 2>/dev/null \
  | sed 's/.*= *//' | tr -d '\r' | sort | uniq -c

echo
echo "=============================================================="
echo "DECISIVE: does any unlabelled track declare a ChIP antibody?"
echo "(an input/control track usually says 'input', 'IgG' or 'none')"
echo "=============================================================="
for gsm in $GSMS; do
  f="$OUT/${gsm}.soft.txt"
  ab=$(grep -i 'antibody' "$f" 2>/dev/null | head -1 | sed 's/.*= *//' | tr -d '\r')
  printf '  %-12s antibody: %s\n' "$gsm" "${ab:-<none declared>}"
done

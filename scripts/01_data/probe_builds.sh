#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_builds.sh — determine the GENOME BUILD of every series from GSM-level
# data-processing metadata, which is the authoritative declaration.
#
# Why this matters more than anything else in M4: the project build is hg19
# (mandated by GSE230649). Risk R-05 is that an hg38 file gets mixed in silently
# and produces plausible-looking but entirely wrong peak-to-gene links. At least
# one series (GSE281524) has "hg38" in its file names. Guessing is not
# acceptable here.
#
# Also resolves: GSE269424 is titled "ATAC-seq" but its files are named
# "*_ASCL1_*". One of those is misleading and the manifest must record which.
#
# Read-only. Usage: bash scripts/01_data/probe_builds.sh
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

OUT=data/metadata/geo_soft
mkdir -p "$OUT"

# One representative GSM per series (first sample listed in each series record).
declare -A REP=(
  [GSE230649]=GSM7230493
  [GSE269424]=GSM8315374
  [GSE256345]=GSM8094116
  [GSE281523]=GSM8622646
  [GSE281524]=GSM8622658
  [GSE210113]=GSM6421552
  [GSE249362]=GSM7942899
  [GSE60052]=GSM1464282
  [GSE261348]=GSM8140436
  [GSE261345]=GSM8140261
)

for acc in "${!REP[@]}"; do
  gsm="${REP[$acc]}"
  f="$OUT/${gsm}.soft.txt"
  if [ ! -s "$f" ]; then
    curl -4 -sS --max-time 90 --retry 3 --retry-delay 8 \
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm}&targ=self&form=text&view=quick" \
      -o "$f"
    sleep 3
  fi
done

echo "=== genome build / assay per series (from GSM data-processing fields) ==="
for acc in GSE230649 GSE269424 GSE256345 GSE281523 GSE281524 GSE210113 GSE249362 GSE60052 GSE261348 GSE261345; do
  gsm="${REP[$acc]}"
  f="$OUT/${gsm}.soft.txt"
  echo "-------------------------------------------------------------"
  printf '[%s]  rep sample %s\n' "$acc" "$gsm"
  if [ ! -s "$f" ]; then echo "  FETCH FAILED"; continue; fi
  grep -m1 '^!Sample_title' "$f" | cut -c1-120
  grep -m1 '^!Sample_library_strategy' "$f"
  grep -m1 '^!Sample_source_name_ch1' "$f" | cut -c1-120
  # antibody tells us ChIP target; absent for ATAC
  grep -i '^!Sample_characteristics_ch1.*\(antibody\|chip\|cell line\)' "$f" | head -3 | cut -c1-120
  echo "  --- build evidence ---"
  grep -iE '^!Sample_data_processing.*(hg19|hg38|GRCh37|GRCh38|genome build|assembly)' "$f" | head -4 | cut -c1-200
  grep -m1 '^!Sample_supplementary_file' "$f" | cut -c1-160
done

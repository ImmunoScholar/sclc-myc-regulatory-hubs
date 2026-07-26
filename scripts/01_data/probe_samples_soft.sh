#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_samples_soft.sh — full per-sample metadata via the series family.soft
# file, for series that publish no series_matrix (ChIP-only series often don't).
#
# One gzipped SOFT file carries every sample's title, library strategy and
# antibody, so this replaces fetching 125 GSM pages one at a time.
#
# Read-only. Usage: bash scripts/01_data/probe_samples_soft.sh
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

OUT=data/metadata/geo_soft
mkdir -p "$OUT"

fetch_family () {
  local acc="$1"
  local num="${acc#GSE}"
  local d="GSE${num:0:$(( ${#num} - 3 ))}nnn"
  local f="$OUT/${acc}_family.soft"
  if [ ! -s "$f" ]; then
    curl -4 -sS --max-time 300 --retry 3 --retry-delay 8 \
      "https://ftp.ncbi.nlm.nih.gov/geo/series/${d}/${acc}/soft/${acc}_family.soft.gz" \
      | gunzip -c 2>/dev/null \
      | grep -E '^\^SAMPLE|^!Sample_title|^!Sample_library_strategy|^!Sample_characteristics_ch1|^!Sample_supplementary_file' \
      > "$f"
    sleep 2
  fi
  echo "$f"
}

summarise () {
  local acc="$1" label="$2"
  echo "#################################################################"
  echo "### $label"
  local f; f=$(fetch_family "$acc")
  if [ ! -s "$f" ]; then echo "  (family.soft unavailable)"; return; fi
  echo "  samples: $(grep -c '^\^SAMPLE' "$f")"
  echo "  --- library strategies ---"
  grep '^!Sample_library_strategy' "$f" | sed 's/.*= *//' | sort | uniq -c
  echo "  --- ChIP antibodies ---"
  grep -i '^!Sample_characteristics_ch1 *= *chip antibody' "$f" \
    | sed 's/.*antibody: *//' | sed 's/ *(.*//' | sort | uniq -c | head -15
  echo "  --- sample titles ---"
  grep '^!Sample_title' "$f" | sed 's/.*= *//' | head -30
}

summarise GSE281524 "GSE281524 — ASCL1 ChIP (hg38), 10 samples"
summarise GSE210113 "GSE210113 — NEUROD1 / H3K27ac (hg19), 14 samples"

echo
echo "#################################################################"
echo "### GSE249362 — POU2F3/ncBAF (hg19), 125 samples: TRIAGE"
f=$(fetch_family GSE249362)
if [ -s "$f" ]; then
  echo "  samples: $(grep -c '^\^SAMPLE' "$f")"
  echo "  --- library strategies ---"
  grep '^!Sample_library_strategy' "$f" | sed 's/.*= *//' | sort | uniq -c
  echo "  --- ChIP antibodies (all) ---"
  grep -i '^!Sample_characteristics_ch1 *= *chip antibody' "$f" \
    | sed 's/.*antibody: *//' | sed 's/ *(.*//' | sort | uniq -c
  echo "  --- cell lines ---"
  grep -i '^!Sample_characteristics_ch1 *= *cell line' "$f" \
    | sed 's/.*line: *//' | sort | uniq -c | head -12
  echo "  --- titles containing POU2F3 ---"
  grep '^!Sample_title' "$f" | sed 's/.*= *//' | grep -i pou2f3 | head -20
else
  echo "  (family.soft unavailable)"
fi

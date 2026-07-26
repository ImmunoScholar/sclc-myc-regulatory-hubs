#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_lineoverlap.sh — do the lineage-TF control datasets share CELL LINES
# with the keystone?
#
# Why this is not a detail: the lineage-TF confounding analysis is a PRIMARY
# analysis (risk R-01). MYC amplification correlates with subtype, so
# "paralog-specific" hubs may be lineage-TF targets in disguise. Testing that at
# the binding level requires ASCL1/NEUROD1/POU2F3 occupancy in the SAME lines
# where MYC/MYCN/MYCL1 occupancy was measured. If the line sets barely overlap,
# the comparison has to be made differently, and it is better to know that now
# than after 6 GB of downloads.
#
# Read-only. Usage: bash scripts/01_data/probe_lineoverlap.sh
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

OUT=data/metadata/geo_soft
mkdir -p "$OUT"

get_matrix_head () {
  local acc="$1"
  local num="${acc#GSE}"
  local d="GSE${num:0:$(( ${#num} - 3 ))}nnn"
  local f="$OUT/${acc}_matrix_head.txt"
  if [ ! -s "$f" ]; then
    curl -4 -sS --max-time 150 --retry 3 --retry-delay 6 \
      "https://ftp.ncbi.nlm.nih.gov/geo/series/${d}/${acc}/matrix/${acc}_series_matrix.txt.gz" \
      | gunzip -c 2>/dev/null | head -70 > "$f"
    sleep 2
  fi
  echo "$f"
}

show_titles () {
  local acc="$1" label="$2"
  echo "############ $label ############"
  local f; f=$(get_matrix_head "$acc")
  if [ ! -s "$f" ]; then echo "  (matrix unavailable)"; echo; return; fi
  grep -m1 '^!Sample_title' "$f" | tr '\t' '\n' | tail -n +2 | tr -d '"'
  echo
}

show_titles GSE230649 "GSE230649 keystone (hg19) — MYC/MYCN/MYCL1/H3K27ac/ATAC"
show_titles GSE281524 "GSE281524 ASCL1 ChIP (hg38)"
show_titles GSE269424 "GSE269424 ATAC (hg38)"
show_titles GSE210113 "GSE210113 NEUROD1/H3K27ac (hg19)"

echo "############ GSE249362 POU2F3 triage (hg19) ############"
f=$(get_matrix_head GSE249362)
if [ -s "$f" ]; then
  echo "-- titles mentioning POU2F3 --"
  grep -m1 '^!Sample_title' "$f" | tr '\t' '\n' | tail -n +2 | tr -d '"' \
    | grep -i 'pou2f3' | head -20 || echo "  (none in first 70 lines)"
  echo "-- antibody characteristics present --"
  grep '^!Sample_characteristics_ch1' "$f" | tr '\t' '\n' | tr -d '"' \
    | grep -i 'antibody' | sort -u | head -15
else
  echo "  (matrix unavailable)"
fi

#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_last_two.sh — resolve GSE281524 (ASCL1) and GSE249362 (POU2F3).
#
# Neither publishes a series_matrix or a family.soft that this network can pull.
# GSE281524 has only 10 samples so per-GSM fetching is cheap. GSE249362 has 125,
# but its supplementary FILE NAMES encode cell line and antibody
# (e.g. ..._ChIP_NCIH1048_DMSO_RPB1_...), so filelist.txt is enough to triage
# without 125 separate requests.
#
# Read-only. Usage: bash scripts/01_data/probe_last_two.sh
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

OUT=data/metadata/geo_soft
mkdir -p "$OUT"
gsm_dir () { local g="$1"; echo "${g:0:$(( ${#g} - 3 ))}nnn"; }

echo "############ GSE281524 — ASCL1 ChIP (hg38): all 10 samples ############"
ids=$(grep '^!Series_sample_id' "$OUT/GSE281524.soft.txt" | sed 's/.*= *//' | tr -d '\r')
for gsm in $ids; do
  f="$OUT/${gsm}.soft.txt"
  if [ ! -s "$f" ]; then
    curl -4 -sS --max-time 90 --retry 3 --retry-delay 8 \
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm}&targ=self&form=text&view=quick" \
      -o "$f"
    sleep 2
  fi
  t=$(grep -m1 '^!Sample_title' "$f" | sed 's/.*= *//')
  ab=$(grep -im1 'chip antibody' "$f" | sed 's/.*antibody: *//' | sed 's/ *(.*//')
  printf '  %-12s %-34s ab=%s\n' "$gsm" "$t" "${ab:-none}"
done

echo
echo "############ GSE249362 — POU2F3/ncBAF (hg19): filelist triage ############"
fl="$OUT/GSE249362_filelist.txt"
if [ ! -s "$fl" ]; then
  curl -4 -sS --max-time 120 --retry 4 --retry-delay 15 \
    -A "Mozilla/5.0 (X11; Linux x86_64) research-data-fetch" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249362/suppl/filelist.txt" \
    -o "$fl" || true
fi
if [ -s "$fl" ] && head -1 "$fl" | grep -q 'Archive'; then
  echo "  total entries: $(grep -c '^File' "$fl")"
  echo "  --- distinct file types ---"
  awk -F'\t' 'NR>1 {print $5}' "$fl" | sort | uniq -c
  echo "  --- antibody token frequency (parsed from file names) ---"
  awk -F'\t' 'NR>1 {print $2}' "$fl" \
    | grep -oiE '_(POU2F3|BRD9|H3K27ac|RPB1|SMARCA4|SMARCA2|IgG|input)_' \
    | tr 'A-Z' 'a-z' | sort | uniq -c | sort -rn
  echo "  --- cell-line token frequency ---"
  awk -F'\t' 'NR>1 {print $2}' "$fl" \
    | grep -oiE '(NCIH1048|H1048|NCIH211|H211|H526|H69|H847|SHP77|COLO668|H889|H524|H196|OS3|H1836|H446)' \
    | tr 'a-z' 'A-Z' | sort | uniq -c | sort -rn
  echo "  --- POU2F3 files (candidate downloads) ---"
  awk -F'\t' 'NR>1 {print $2"\t"$4}' "$fl" | grep -i pou2f3 | head -20
else
  echo "  filelist STILL unavailable (NCBI returned an error page or 403)."
  echo "  GSE249362 stays 'pending_triage' in the manifest; nothing is guessed."
fi

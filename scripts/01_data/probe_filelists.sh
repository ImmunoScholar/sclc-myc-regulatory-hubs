#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_filelists.sh — inspect GEO filelist.txt for each series.
#
# Purpose: find out whether GEO publishes MD5 checksums and per-file sizes. If
# it does, the manifest can carry REAL checksums rather than trust-on-first-use
# values we computed ourselves (which prove only that a file did not change
# since we fetched it, not that we fetched the right file).
#
# Read-only. Usage: bash scripts/01_data/probe_filelists.sh
# -----------------------------------------------------------------------------
set -uo pipefail

series_dir () {  # GSE230649 -> GSE230nnn
  local acc="$1" num prefix
  num="${acc#GSE}"
  if [ "${#num}" -le 3 ]; then prefix="GSEnnn"; else prefix="GSE${num:0:$(( ${#num} - 3 ))}nnn"; fi
  echo "$prefix"
}

echo "=== GSE230649 filelist.txt (keystone) ==="
curl -4 -sS --max-time 45 \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230649/suppl/filelist.txt" | head -35

echo
echo "=== supplementary directory listing per series (name / size) ==="
for acc in GSE269424 GSE256345 GSE281523 GSE281524 GSE210113 GSE249362 GSE60052 GSE261348 GSE261345; do
  d=$(series_dir "$acc")
  echo "--- $acc ---"
  curl -4 -sS --max-time 45 "https://ftp.ncbi.nlm.nih.gov/geo/series/${d}/${acc}/suppl/" \
    | sed -e 's/<[^>]*>//g' \
    | grep -vE '^\s*$|Index of|Parent Directory|Last modified|HHS Vulnerability|^Name' \
    | head -12
done

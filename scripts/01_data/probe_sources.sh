#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_sources.sh — verify every data host is reachable BEFORE writing download
# code that assumes it is. bioconductor.org turned out to be unroutable from
# this network (D-006); assume nothing.
#
# Read-only. Downloads nothing except tiny index pages.
# Usage: bash scripts/01_data/probe_sources.sh
# -----------------------------------------------------------------------------
set -uo pipefail

probe () {
  local label="$1" url="$2"
  local code
  code=$(curl -4 -sS -o /dev/null -w '%{http_code}' --max-time 30 -L "$url" 2>&1)
  printf '  %-28s %-6s %s\n' "$label" "$code" "$url"
}

echo "=== hosts ==="
probe "NCBI (GEO web)"       "https://www.ncbi.nlm.nih.gov/geo/"
probe "NCBI FTP (https)"     "https://ftp.ncbi.nlm.nih.gov/geo/"
probe "cBioPortal API"       "https://www.cbioportal.org/api/studies/sclc_ucologne_2015"
probe "DepMap"               "https://depmap.org/portal/"
probe "figshare (DepMap DL)" "https://ndownloader.figshare.com/"
probe "SCLC-CellMiner"       "https://discover.nci.nih.gov/SclcCellMinerCDB/"
probe "MSigDB"               "https://www.gsea-msigdb.org/gsea/msigdb/"

echo
echo "=== GEO supplementary file listings (sizes, no download) ==="
for acc in GSE230649 GSE269424 GSE256345 GSE281523 GSE281524 GSE210113 GSE249362 GSE60052 GSE261348 GSE261345; do
  url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=${acc}"
  code=$(curl -4 -s -o /dev/null -w '%{http_code}' --max-time 30 "$url")
  printf '  %-12s HTTP %s\n' "$acc" "$code"
done

echo
echo "=== FTP-style supplementary directory index for the keystone dataset ==="
# GEO exposes a browsable suppl directory; this is where file names and sizes live.
curl -4 -sS --max-time 45 \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230649/suppl/" \
  | sed -e 's/<[^>]*>//g' | grep -vE '^\s*$' | head -40

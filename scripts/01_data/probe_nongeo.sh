#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_nongeo.sh — verify the non-GEO resources and the liftOver chain file.
#
# The chain file has become load-bearing: four of the seven coordinate-level
# datasets turned out to be hg38 while the project build is hg19. Confirm the
# chain is fetchable before designing around it.
#
# Also checks DepMap. The Bioconductor `depmap` package routes through
# ExperimentHub -> bioconductor.org, which is unroutable from this network
# (D-006), so a direct download path is required.
#
# Read-only. Usage: bash scripts/01_data/probe_nongeo.sh
# -----------------------------------------------------------------------------
set -uo pipefail

head_info () {
  local label="$1" url="$2"
  local out
  out=$(curl -4 -sSI --max-time 45 -L "$url" 2>&1 \
        | awk 'BEGIN{IGNORECASE=1} /^HTTP/{c=$2} /^Content-Length:/{l=$2} /^Content-Type:/{t=$2} END{gsub(/\r/,"",l);gsub(/\r/,"",t);print c" "(l==""?"-":l)" "(t==""?"-":t)}')
  printf '  %-34s %s\n' "$label" "$out"
}

echo "=== UCSC liftOver chain files (hg38 -> hg19) ==="
head_info "hg38ToHg19.over.chain.gz" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
head_info "hg19ToHg38.over.chain.gz" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
echo "  (mirror check)"
head_info "hgdownload2 mirror" \
  "https://hgdownload2.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"

echo
echo "=== UCSC hg19 chrom sizes (for build assertions) ==="
head_info "hg19.chrom.sizes" \
  "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes"

echo
echo "=== cBioPortal: what does sclc_ucologne_2015 actually contain? ==="
curl -4 -sS --max-time 45 "https://www.cbioportal.org/api/studies/sclc_ucologne_2015/molecular-profiles" \
  | tr ',' '\n' | grep -iE 'molecularProfileId|molecularAlterationType|datatype' | head -24

echo
echo "=== DepMap: locate a downloadable release ==="
head_info "depmap downloads page" "https://depmap.org/portal/api/download/files"
echo "  first 600 chars of the download index:"
curl -4 -sS --max-time 60 "https://depmap.org/portal/api/download/files" | cut -c1-600

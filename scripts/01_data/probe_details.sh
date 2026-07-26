#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_details.sh — resolve three questions before the manifest is written:
#   1. Can individual GSM supplementary files be fetched directly? If so we
#      avoid the monolithic 8.6 GB tar entirely (28 resumable downloads instead
#      of one, and no need for ~17 GB of transient disk to untar).
#   2. Do the GEO endpoints support HTTP range requests (i.e. is `curl -C -`
#      resume real, or would a dropped connection restart from zero)?
#   3. What are the FULL supplementary listings for the two GeoMx series? The
#      first probe was truncated and showed only large PNGs; the count matrices
#      must be identified so ~1.5 GB of images can be skipped.
#
# Read-only: HEAD requests and index pages only.
# -----------------------------------------------------------------------------
set -uo pipefail

GSM_BASE="https://ftp.ncbi.nlm.nih.gov/geo/samples"

gsm_dir () { local g="$1"; echo "${g:0:$(( ${#g} - 3 ))}nnn"; }

echo "=== 1. individual GSM supplementary file access + range support ==="
for spec in \
  "GSM7230503 GSM7230503_MYC_H1048_R1_treat_pileup.bedgraph.gz" \
  "GSM7230510 GSM7230510_MYCN_H526_R1_treat_pileup.bedgraph.gz" \
  "GSM7230512 GSM7230512_COLO668_treat_pileup.bedgraph.gz" ; do
  set -- $spec
  gsm="$1"; fn="$2"
  url="${GSM_BASE}/$(gsm_dir "$gsm")/${gsm}/suppl/${fn}"
  read -r code len ranges < <(curl -4 -sS -I --max-time 45 -L "$url" \
    | awk 'BEGIN{IGNORECASE=1} /^HTTP/{c=$2} /^Content-Length:/{l=$2} /^Accept-Ranges:/{r=$2} END{gsub(/\r/,"",l); gsub(/\r/,"",r); print c, (l==""?"-":l), (r==""?"none":r)}')
  printf '  %-12s HTTP %-4s bytes=%-12s accept-ranges=%s\n' "$gsm" "$code" "$len" "$ranges"
done

echo
echo "  explicit range request test (first 100 bytes):"
u="${GSM_BASE}/GSM7230nnn/GSM7230503/suppl/GSM7230503_MYC_H1048_R1_treat_pileup.bedgraph.gz"
got=$(curl -4 -sS --max-time 45 -r 0-99 -o /dev/null -w '%{http_code} %{size_download}' -L "$u")
echo "    HTTP+bytes: $got   (expect '206 100' if resume is genuinely supported)"

echo
echo "=== 2. GeoMx FULL supplementary listings ==="
for acc in GSE261348 GSE261345; do
  num="${acc#GSE}"; d="GSE${num:0:$(( ${#num} - 3 ))}nnn"
  echo "--- $acc (non-PNG entries only) ---"
  curl -4 -sS --max-time 60 "https://ftp.ncbi.nlm.nih.gov/geo/series/${d}/${acc}/suppl/" \
    | sed -e 's/<[^>]*>//g' \
    | grep -vE '^\s*$|Index of|Parent Directory|Last modified|HHS Vulnerability|^Name' \
    | grep -viE '\.png' || echo "  (none)"
done

echo
echo "=== 3. exact byte sizes from filelist.txt for the remaining series ==="
for acc in GSE269424 GSE256345 GSE281523 GSE281524 GSE210113 GSE249362; do
  num="${acc#GSE}"; d="GSE${num:0:$(( ${#num} - 3 ))}nnn"
  echo "--- $acc ---"
  curl -4 -sS --max-time 60 "https://ftp.ncbi.nlm.nih.gov/geo/series/${d}/${acc}/suppl/filelist.txt" \
    | head -6
done

echo
echo "=== 4. GSE60052 processed matrix header (what form are the data in?) ==="
u60="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60052/suppl/GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv.gz"
curl -4 -sS --max-time 60 -r 0-3000 "$u60" | gunzip -c 2>/dev/null | cut -c1-300 | head -3

#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# probe_geo_metadata.sh — fetch authoritative GEO series metadata.
#
# Filenames are not evidence. Two of the series in the frozen dataset inventory
# look, from their file names alone, like they may be a different assay or a
# different genome build than recorded. This pulls the actual GEO series records
# so the manifest is built from the source of truth.
#
# Saves brief SOFT records under data/metadata/geo_soft/ (a few KB each, tracked
# as provenance) and prints the fields that matter.
#
# Read-only w.r.t. the analysis. Usage: bash scripts/01_data/probe_geo_metadata.sh
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

OUT=data/metadata/geo_soft
mkdir -p "$OUT"

ACCS="GSE230649 GSE269424 GSE256345 GSE281523 GSE281524 GSE210113 GSE249362 GSE60052 GSE261348 GSE261345"

for acc in $ACCS; do
  f="$OUT/${acc}.soft.txt"
  if [ ! -s "$f" ]; then
    curl -4 -sS --max-time 90 --retry 3 --retry-delay 5 \
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${acc}&targ=self&form=text&view=brief" \
      -o "$f"
    sleep 2   # be polite to NCBI; unthrottled bursts get 403s
  fi
done

echo "=== series summary ==="
for acc in $ACCS; do
  f="$OUT/${acc}.soft.txt"
  echo "-------------------------------------------------------------"
  echo "[$acc]"
  if [ ! -s "$f" ]; then echo "  FETCH FAILED"; continue; fi
  grep -m1 '^!Series_title'            "$f" | cut -c1-160
  grep -m1 '^!Series_type'             "$f"
  grep    '^!Series_type'              "$f" | tail -n +2
  grep -m1 '^!Series_platform_organism' "$f"
  grep -m1 '^!Series_pubmed_id'        "$f" || echo "!Series_pubmed_id = (none listed)"
  # Genome build is usually declared in the data-processing lines.
  grep -i 'genome build\|hg19\|hg38\|GRCh37\|GRCh38' "$f" | head -3
  grep -m1 '^!Series_sample_id' "$f" > /dev/null && \
    echo "  n samples: $(grep -c '^!Series_sample_id' "$f")"
done

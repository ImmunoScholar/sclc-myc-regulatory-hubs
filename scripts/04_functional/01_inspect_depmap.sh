#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 01_inspect_depmap.sh — what is actually in the DepMap downloads?
#
#   bash scripts/04_functional/01_inspect_depmap.sh
#
# Read-only. Establishes structure BEFORE any analysis is written: identifier
# format, dimensions, whether cell-line metadata came embedded, and whether the
# keystone lines are present. Filenames from Custom Downloads differ from the
# All Data names the manifest expects, so nothing is assumed.
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
D=data/raw/depmap

for f in depmap_crispr_chronos_26Q1.csv depmap_expression_26Q1.csv depmap_cn_log2_26Q1.csv; do
  p="$D/$f"
  [ -f "$p" ] || { echo "MISSING $f"; continue; }
  echo "=============================================================="
  echo "$f"
  echo "  size : $(du -h "$p" | cut -f1)"
  echo "  rows : $(( $(wc -l < "$p") - 1 ))  (excluding header)"
  ncol=$(head -1 "$p" | awk -F',' '{print NF}')
  echo "  cols : $ncol"
  echo
  echo "  --- first 8 column names ---"
  head -1 "$p" | tr ',' '\n' | head -8 | sed 's/^/    /'
  echo
  echo "  --- first data row, first 4 fields ---"
  sed -n '2p' "$p" | cut -d',' -f1-4 | sed 's/^/    /'
  echo
  echo "  --- metadata columns embedded? ---"
  head -1 "$p" | tr ',' '\n' | grep -inE 'lineage|primary.?disease|subtype|oncotree|cell.?line.?name|stripped|depmap.?id|model.?id|sex|age' \
    | head -12 | sed 's/^/    /' || echo "    none found"
  echo
done

echo "=============================================================="
echo "KEYSTONE LINE LOOKUP"
echo "=============================================================="
echo "Our lines: COLO668 H1048 H196 H211 H524 H526 H69 H847 H889 SHP77"
echo "DepMap uses ACH-###### ModelIDs and/or stripped names (NCIH1048, SHP77)."
echo
p="$D/depmap_crispr_chronos_26Q1.csv"
if [ -f "$p" ]; then
  echo "  --- first column values, first 5 ---"
  cut -d',' -f1 "$p" | sed -n '2,6p' | sed 's/^/    /'
  echo
  echo "  --- rows matching keystone name patterns ---"
  cut -d',' -f1-3 "$p" | grep -iE 'COLO ?668|NCIH1048|NCIH196|NCIH211|NCIH524|NCIH526|NCIH69|NCIH847|NCIH889|SHP ?77' \
    | head -15 | sed 's/^/    /' || echo "    no name matches in the first columns (likely ACH IDs only)"
fi

echo
echo "=============================================================="
echo "SCLC LINE COUNT (needs a lineage column, or Model.csv)"
echo "=============================================================="
if [ -f "$p" ]; then
  if head -1 "$p" | grep -qiE 'lineage|oncotree|primary.?disease'; then
    echo "  a lineage-like column exists — SCLC lines countable from this file"
  else
    echo "  NO lineage column in the CRISPR file."
    echo "  Model.csv (or the metadata checkbox) is REQUIRED: selective dependency"
    echo "  is defined by comparing SCLC against other lineages, so lineage"
    echo "  assignment is not optional for the primary M7 analysis."
  fi
fi

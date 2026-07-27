#!/usr/bin/env bash
# Verify whether recent edits actually landed on the WSL filesystem.
cd "$(dirname "$0")/../.." || exit 1
F=scripts/03_tumour/05_occupancy_confounding.R
echo "file  : $F"
echo "size  : $(stat -c%s "$F") bytes"
echo "mtime : $(stat -c%y "$F")"
echo "lines : $(wc -l < "$F")"
echo
echo "marker 'peak-set QC'      : $(grep -c 'peak-set QC' "$F" || true)"
echo "marker 'POU_SUBTYPE_LINE' : $(grep -c 'POU_SUBTYPE_LINE' "$F" || true)"
echo "marker 'DECISIVE CONTRAST': $(grep -c 'DECISIVE CONTRAST' "$F" || true)"
echo "old text 'reaches the CHROMATIN' : $(grep -c 'reaches the CHROMATIN' "$F" || true)"
echo
echo "--- interpretation block as it exists on disk ---"
sed -n '/interpretation ===/,/ASCL1 NOT TESTED/p' "$F" | head -25

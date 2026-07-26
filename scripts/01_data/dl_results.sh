#!/usr/bin/env bash
# Inspect per-file worker results so far. Read-only.
cd "$(dirname "$0")/../.." || exit 1
D=logs/dl_results
[ -d "$D" ] || { echo "no results yet"; exit 0; }
n=$(find "$D" -type f | wc -l)
echo "worker results written: $n"
echo
echo "=== by status ==="
cat "$D"/* 2>/dev/null | cut -f1 | sort | uniq -c
echo
echo "=== any failures ==="
cat "$D"/* 2>/dev/null | awk -F'\t' '$1=="FAIL" {printf "  %s/%s :: %s\n", $2, $3, $4}'
echo "  (nothing above = no failures yet)"
echo
echo "=== .part files with sizes (in-flight or awaiting resume) ==="
find data/raw -name '*.part' -printf '%10s  %p\n' | sort -k2

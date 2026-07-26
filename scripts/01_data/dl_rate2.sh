#!/usr/bin/env bash
# Live throughput of the parallel downloader. Read-only.
cd "$(dirname "$0")/../.." || exit 1

echo "active curl streams : $(pgrep -c curl 2>/dev/null || echo 0)"
b1=$(du -sb data/raw 2>/dev/null | cut -f1)
sleep 20
b2=$(du -sb data/raw 2>/dev/null | cut -f1)
rate=$(( (b2 - b1) / 20 ))
total=$(awk -F'\t' 'NR>1 {s+=$4} END {print s}' data/metadata/download_list.tsv)

echo "on disk             : $(numfmt --to=iec "$b2")"
echo "target              : $(numfmt --to=iec "$total")"
echo "progress            : $(awk -v b="$b2" -v t="$total" 'BEGIN{printf "%.1f%%", 100*b/t}')"
echo "current rate        : $(numfmt --to=iec "$rate")/s"
if [ "$rate" -gt 0 ]; then
  eta=$(( (total - b2) / rate ))
  echo "projected remaining : $(( eta / 3600 ))h $(( (eta % 3600) / 60 ))m"
fi
echo "files complete      : $(find data/raw -type f ! -name '*.part' ! -name '.gitkeep' | wc -l) / 66"
echo "in flight (.part)   : $(find data/raw -name '*.part' | wc -l)"
if pgrep -f 02_download.sh > /dev/null; then echo "process             : running"; else echo "process             : NOT running"; fi

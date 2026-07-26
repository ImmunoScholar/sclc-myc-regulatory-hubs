#!/usr/bin/env bash
# Measure actual download throughput and project a completion time. Read-only.
cd "$(dirname "$0")/../.." || exit 1
LOG=logs/02_download.log

start_line=$(grep -m1 'download run:' "$LOG" | sed 's/.*download run: //')
start_epoch=$(date -d "$start_line" +%s 2>/dev/null)
now_epoch=$(date +%s)
elapsed=$(( now_epoch - start_epoch ))

bytes=$(du -sb data/raw 2>/dev/null | cut -f1)
total=$(awk -F'\t' 'NR>1 {s+=$4} END {print s}' data/metadata/download_list.tsv)

echo "started        : $start_line"
echo "now            : $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
echo "elapsed        : $(( elapsed / 60 )) min"
echo
echo "on disk        : $(numfmt --to=iec "$bytes")"
echo "target         : $(numfmt --to=iec "$total")"
echo "progress       : $(awk -v b="$bytes" -v t="$total" 'BEGIN{printf "%.1f%%", 100*b/t}')"
echo

if [ "$elapsed" -gt 0 ]; then
  rate=$(( bytes / elapsed ))
  echo "mean rate      : $(numfmt --to=iec "$rate")/s"
  remain=$(( total - bytes ))
  if [ "$rate" -gt 0 ]; then
    eta=$(( remain / rate ))
    echo "projected ETA  : $(( eta / 3600 ))h $(( (eta % 3600) / 60 ))m remaining"
  fi
fi

echo
echo "instantaneous check (10 s sample):"
b1=$(du -sb data/raw 2>/dev/null | cut -f1)
sleep 10
b2=$(du -sb data/raw 2>/dev/null | cut -f1)
inst=$(( (b2 - b1) / 10 ))
echo "  current rate : $(numfmt --to=iec "$inst")/s"

echo
echo "files done     : $(grep -c '^  ok ' "$LOG")"
echo "curl retries   : $(grep -ci 'transient\|retry\|Connection' "$LOG" 2>/dev/null || echo 0)"
if pgrep -f 02_download.sh > /dev/null; then echo "process        : running"; else echo "process        : NOT running"; fi

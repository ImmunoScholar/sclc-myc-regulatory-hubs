#!/usr/bin/env bash
# Compact download progress. Usage: bash scripts/01_data/dl_status.sh
cd "$(dirname "$0")/../.." || exit 1
LOG=logs/02_download.log
[ -f "$LOG" ] || { echo "no log yet"; exit 0; }
echo "completed files : $(grep -c '^  ok ' "$LOG")"
echo "failures        : $(grep -c '^  FAIL ' "$LOG")"
echo "on disk         : $(du -sh data/raw 2>/dev/null | cut -f1)"
echo "expected total  : 12.08 GB across 66 files"
echo
echo "currently:"
tail -n 2 "$LOG"
echo
if pgrep -f 02_download.sh > /dev/null; then echo "STATUS: running"; else echo "STATUS: not running"; fi

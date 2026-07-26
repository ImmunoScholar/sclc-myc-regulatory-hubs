#!/usr/bin/env bash
# What outputs actually exist right now? Read-only.
cd "$(dirname "$0")/../.." || exit 1
echo "=== figures/ ==="
ls -la figures/
echo
echo "=== results/ (all files) ==="
find results -type f | sort
echo
echo "=== image files anywhere (excluding .git and renv) ==="
find . -path ./.git -prune -o -path ./renv -prune -o -type f \
     -name '*.png' -print -o -type f -name '*.pdf' -print \
     -o -type f -name '*.svg' -print | sort
echo "(none above means none exist)"
echo
echo "=== data downloaded so far ==="
du -sh data/raw 2>/dev/null
find data/raw -type f -name '*.part' | wc -l | xargs echo "in-progress .part files:"
find data/raw -type f ! -name '*.part' ! -name '.gitkeep' | wc -l | xargs echo "completed files:"

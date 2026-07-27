#!/usr/bin/env bash
# Objective snapshot of project state. Read-only.
cd "$(dirname "$0")/../.." || exit 1
echo "=== commits ==="; git rev-list --count HEAD
echo "=== tracked files ==="; git ls-files | wc -l
echo
echo "=== scripts by stage ==="
for d in scripts/*/; do
  n=$(find "$d" -type f \( -name '*.R' -o -name '*.sh' \) | wc -l)
  printf '  %-28s %s\n' "$(basename "$d")" "$n"
done
echo
echo "=== analysis outputs so far ==="
printf '  figures      : %s\n' "$(find figures -type f -name '*.png' -o -name '*.pdf' 2>/dev/null | wc -l)"
printf '  result tables: %s\n' "$(find results/tables -type f 2>/dev/null | wc -l)"
printf '  result objects: %s\n' "$(find results/objects -type f ! -name '.gitkeep' 2>/dev/null | wc -l)"
echo
echo "=== data acquired ==="
du -sh data/raw 2>/dev/null | cut -f1 | xargs printf '  raw    : %s\n'
printf '  files  : %s\n' "$(find data/raw -type f ! -name '.gitkeep' | wc -l)"
printf '  processed: %s\n' "$(find data/processed -type f ! -name '.gitkeep' 2>/dev/null | wc -l)"
echo
echo "=== decisions / risks logged ==="
printf '  decisions: %s\n' "$(grep -c '^### D-' docs/decision_log.md)"
printf '  risks    : %s\n' "$(grep -c '^### R-' docs/06_risk_log.md)"
printf '  open HIGH: %s\n' "$(grep -c '^### R-.*HIGH.*OPEN' docs/06_risk_log.md)"
echo
echo "=== config: key frozen parameters ==="
grep -E 'signal_quantile|min_datasets_supporting|stitch_distance|max_distance|permutations' config/params.yml | sed 's/^/  /'

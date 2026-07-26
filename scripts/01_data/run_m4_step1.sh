#!/usr/bin/env bash
# Rebuild the manifest, then dry-run the downloader. Usage:
#   bash scripts/01_data/run_m4_step1.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
mkdir -p logs

echo "=== 01_build_manifest.R ==="
Rscript scripts/01_data/01_build_manifest.R > logs/01_manifest.log 2>&1
rc=$?
echo "exit code: $rc"
tail -n 8 logs/01_manifest.log
[ $rc -ne 0 ] && exit $rc

echo
echo "=== 02_download.sh --dry-run ==="
bash scripts/01_data/02_download.sh --dry-run > logs/02_dryrun.log 2>&1
echo "exit code: $?"
grep -E '^considered|^skipped|^fetched|^failed|^RESULT|WOULD FETCH' logs/02_dryrun.log | tail -20

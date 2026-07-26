#!/usr/bin/env bash
# Wrapper: run the package install with full logging.
# Usage: bash scripts/00_setup/run_install.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
mkdir -p logs
LOG="logs/03_install_packages.log"
echo "log: $(pwd)/$LOG"
Rscript scripts/00_setup/03_install_packages.R > "$LOG" 2>&1
rc=$?
echo "exit code: $rc"
tail -20 "$LOG"
exit $rc

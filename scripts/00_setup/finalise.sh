#!/usr/bin/env bash
# Re-snapshot the environment, verify the lockfile, and run the test suite.
# Usage: bash scripts/00_setup/finalise.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

echo "=== renv snapshot ==="
Rscript scripts/00_setup/04_snapshot.R > logs/04_snapshot.log 2>&1
rc=$?
echo "snapshot exit code: $rc"
grep -E 'packages locked|R version in lock|RESULT' logs/04_snapshot.log

echo
echo "=== lockfile reconciliation ==="
Rscript scripts/00_setup/diag_lockfile.R 2>&1 | grep -E 'locked in|present in|locked but not|RESULT'

echo
echo "=== GSVA loads and is usable ==="
Rscript -e 'suppressPackageStartupMessages(library(GSVA)); cat("GSVA", as.character(packageVersion("GSVA")), "OK\n"); cat("gsvaParam available:", exists("gsvaParam"), "\n")' 2>&1 | tail -3

echo
echo "=== environment test suite ==="
Rscript tests/test_environment.R 2>&1 | grep -E 'Test passed|Failure|Error|All environment'

#!/usr/bin/env bash
# Wrapper so the long download can run detached with a single log.
# Usage: bash scripts/01_data/run_download.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
mkdir -p logs
: > logs/02_download.log
bash scripts/01_data/02_download.sh
rc=$?
echo "download exit code: $rc"
exit $rc

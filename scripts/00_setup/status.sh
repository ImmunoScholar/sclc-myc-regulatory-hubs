#!/usr/bin/env bash
# Compact status of the package install.
# NOTE: renv writes progress with carriage returns and no newline, which can
# prepend to our own status lines. Translate \r to \n before grepping, or
# genuine [ok]/[FAIL] lines are silently missed.
cd "$(dirname "$0")/../.." || exit 1
LOG=logs/03_install_packages.log
[ -f "$LOG" ] || { echo "no log yet"; exit 0; }
echo "=== groups / results so far ==="
tr '\r' '\n' < "$LOG" | grep -E '^####|\[ok\]|\[have\]|\[FAIL\]|^RESULT|^FAILED'
echo
echo "=== currently ==="
tail -c 400 "$LOG" | tr '\r' '\n' | tail -3
echo
if pgrep -f 03_install_packages.R > /dev/null; then echo "STATUS: running"; else echo "STATUS: not running"; fi

#!/usr/bin/env bash
cd "$(dirname "$0")/../.." || exit 1
echo "=== all risks, with status ==="
grep -E '^### R-' docs/06_risk_log.md | sed 's/\*\*//g' | sed 's/^### //'
echo
echo "=== open, HIGH severity ==="
grep -E '^### R-' docs/06_risk_log.md | sed 's/\*\*//g' | grep -i 'HIGH' | grep -i 'OPEN'
echo
echo "=== closed ==="
grep -E '^### R-' docs/06_risk_log.md | sed 's/\*\*//g' | grep -i 'CLOSED'

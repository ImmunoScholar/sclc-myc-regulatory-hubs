#!/usr/bin/env bash
# Stage, commit from .git/COMMIT_MSG_TMP, push, and report.
# Usage: bash scripts/00_setup/commit_push.sh
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
git add -A
if git diff --cached --quiet; then
  echo "nothing to commit"
else
  git commit -q -F .git/COMMIT_MSG_TMP || exit 1
  rm -f .git/COMMIT_MSG_TMP
fi
git push -q origin main || exit 1
echo "=== log ==="
git --no-pager log --oneline
echo
echo "=== working tree ==="
if [ -z "$(git status --porcelain)" ]; then echo "clean"; else git status --short; fi
echo
echo "=== remote sync ==="
git --no-pager log --oneline -1 origin/main

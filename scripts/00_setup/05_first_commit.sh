#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 05_first_commit.sh — stage, verify, and make the initial commit.
#
# Run from the project root:  bash scripts/00_setup/05_first_commit.sh
#
# Verifies BEFORE committing that no data file has crept into the index. The
# pre-commit hook is a second line of defence; this is the first.
# -----------------------------------------------------------------------------
set -euo pipefail

cd "$(dirname "$0")/../.." || exit 1
[ -f .gitignore ] || { echo "FATAL: cannot locate the project root" >&2; exit 1; }

git add -A

echo "=== staged files ==="
git diff --cached --name-only | sort

echo
echo "=== pre-commit verification ==="

fail=0

# 1. Nothing from the data payload directories except .gitkeep placeholders.
leaked=$(git diff --cached --name-only | grep -E '^data/(raw|processed)/' | grep -v '\.gitkeep$' || true)
if [ -n "$leaked" ]; then
  echo "FAIL: data payload staged:"; echo "$leaked"; fail=1
else
  echo "  ok: no data payload staged"
fi

# 2. No large-data file extensions anywhere in the index.
badext=$(git diff --cached --name-only | grep -Ei '\.(bedGraph|bw|bam|fastq\.gz|rds|RData|h5|h5ad|tar|zip|gz)$' || true)
if [ -n "$badext" ]; then
  echo "FAIL: data-format file staged:"; echo "$badext"; fail=1
else
  echo "  ok: no data-format extensions staged"
fi

# 3. Nothing over 5 MB.
big=0
while IFS= read -r f; do
  [ -f "$f" ] || continue
  sz=$(stat -c%s "$f")
  if [ "$sz" -gt $((5*1024*1024)) ]; then echo "FAIL: $f is $sz bytes"; big=1; fi
done < <(git diff --cached --name-only)
[ "$big" -eq 0 ] && echo "  ok: no file over 5 MB" || fail=1

# 4. No credential files.
cred=$(git diff --cached --name-only | grep -E '(^|/)(\.Renviron|\.env|.*\.pem|.*\.key)$' || true)
if [ -n "$cred" ]; then echo "FAIL: credential file staged:"; echo "$cred"; fail=1
else echo "  ok: no credential files staged"; fi

# 5. renv library must not be tracked (renv.lock must be).
lib=$(git diff --cached --name-only | grep -E '^renv/library/' || true)
if [ -n "$lib" ]; then echo "FAIL: renv library staged"; fail=1
else echo "  ok: renv library not staged"; fi
if git diff --cached --name-only | grep -qx 'renv.lock'; then
  echo "  ok: renv.lock IS staged"
else
  echo "  WARNING: renv.lock is not staged — has the snapshot been run?"
fi

echo
if [ "$fail" -ne 0 ]; then
  echo "ABORTING: verification failed. Nothing committed." >&2
  exit 1
fi

total=$(git diff --cached --name-only | wc -l)
bytes=$(git diff --cached --name-only | xargs -r -I{} stat -c%s {} 2>/dev/null | awk '{s+=$1} END {print s+0}')
echo "verification PASSED — ${total} files, ${bytes} bytes total"
echo

git commit -q -F - <<'MSG'
Initial commit: project scaffold, environment, and specification

Repository structure, locked R environment, and the frozen project
specification for an analysis of MYC-paralog regulatory hubs in SCLC.

Environment:
- R 4.6.1 / Bioconductor 3.23, pinned via renv.lock
- CRAN from Posit P3M noble binaries; Bioconductor from the GWDG mirror
  (bioconductor.org is unroutable from this network)
- Pure R; no Python dependency

Data hygiene:
- .gitignore written before any download and verified adversarially:
  11 decoy data files confirmed ignored, 9 deliverable paths confirmed
  still trackable
- pre-commit hook rejects any staged file over 5 MB

Documentation:
- docs/project_contract.md  scope freeze, aims, pre-registered controls
- docs/decision_log.md      numbered decisions with rationale
- docs/research_journal.md  dated working record
- docs/01-07                gap statement, dataset inventory, architecture,
                            dependencies, risk log, milestone roadmap

No data files are tracked. All datasets are fetched by scripts and
remain subject to their original repository terms.
MSG

echo "=== commit created ==="
git --no-pager log --stat --oneline -1 | head -40
echo
echo "Branch: $(git rev-parse --abbrev-ref HEAD)"
echo "Remote: $(git remote get-url origin)"
echo
echo "NOT PUSHED. Review, then push with:  git push -u origin main"

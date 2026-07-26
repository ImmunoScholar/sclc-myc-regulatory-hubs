#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 01_init_git.sh — initialise Git, install a large-file pre-commit guard, and
# ADVERSARIALLY TEST that .gitignore blocks data before any download happens.
#
# Run from the project root:  bash scripts/00_setup/01_init_git.sh
# Idempotent: safe to re-run.
# -----------------------------------------------------------------------------
set -euo pipefail

REMOTE_URL="git@github.com:ImmunoScholar/sclc-myc-regulatory-hubs.git"
MAX_MB=5

if [ ! -f .gitignore ]; then
  echo "FATAL: .gitignore not found. It must exist before git init." >&2
  exit 1
fi

# --- directory placeholders so a fresh clone has the expected tree ------------
for d in data/raw data/processed results/objects logs figures notebooks tests; do
  mkdir -p "$d"
  [ -f "$d/.gitkeep" ] || touch "$d/.gitkeep"
done

# --- git init -----------------------------------------------------------------
if [ -d .git ]; then
  echo "git repository already initialised"
else
  git init -q
  echo "git initialised"
fi
git symbolic-ref HEAD refs/heads/main
git config core.autocrlf false
git config core.filemode false

if git remote get-url origin >/dev/null 2>&1; then
  git remote set-url origin "$REMOTE_URL"
else
  git remote add origin "$REMOTE_URL"
fi
echo "origin -> $(git remote get-url origin)"

# --- pre-commit guard: block large files even if .gitignore is edited later ---
cat > .git/hooks/pre-commit <<'HOOK'
#!/usr/bin/env bash
# Refuse any staged file larger than MAX_BYTES. Defence in depth: .gitignore is
# a policy, this is enforcement. Bypass deliberately with --no-verify only if
# you are certain.
set -euo pipefail
MAX_BYTES=$((5 * 1024 * 1024))
fail=0
while IFS= read -r f; do
  [ -f "$f" ] || continue
  size=$(stat -c%s "$f")
  if [ "$size" -gt "$MAX_BYTES" ]; then
    printf 'BLOCKED: %s is %s bytes (limit %s)\n' "$f" "$size" "$MAX_BYTES" >&2
    fail=1
  fi
done < <(git diff --cached --name-only --diff-filter=ACM)
if [ "$fail" -ne 0 ]; then
  echo "pre-commit: large file(s) staged. Data belong in data/ (git-ignored)," >&2
  echo "reproduced by download scripts, not committed." >&2
  exit 1
fi
HOOK
chmod +x .git/hooks/pre-commit
echo "pre-commit large-file guard installed (limit ${MAX_MB} MB)"

# --- ADVERSARIAL TEST ---------------------------------------------------------
# Plant decoy files that mimic real downloads and assert git ignores every one.
echo
echo "=== .gitignore adversarial test ==="
DECOYS=(
  "data/raw/GSE230649_H1048_MYC.bedGraph"
  "data/raw/GSE230649_RAW.tar"
  "data/raw/IMfirst_DSP_rawcounts.xlsx"
  "data/processed/consensus_regions_chr1.rds"
  "data/processed/signal_matrix.h5"
  "results/objects/genie3_network.rds"
  "logs/download.log"
  "stray_download.bedGraph"
  "renv/library/placeholder.txt"
  ".Renviron"
  "data/metadata/big_object.rds"
)
mkdir -p data/raw data/processed results/objects logs renv/library data/metadata
fail=0
for f in "${DECOYS[@]}"; do
  printf 'decoy' > "$f"
  if git check-ignore -q "$f"; then
    printf '  IGNORED  %s\n' "$f"
  else
    printf '  LEAKED   %s   <-- .gitignore does NOT cover this\n' "$f"
    fail=1
  fi
  rm -f "$f"
done

# Files that MUST remain trackable — an over-broad .gitignore is also a failure.
echo
echo "=== must-be-tracked test ==="
KEEPERS=(
  "README.md"
  "renv.lock"
  "docs/decision_log.md"
  "config/params.yml"
  "scripts/01_data/download.R"
  "results/tables/moes_ranking.csv"
  "figures/fig01_regions.png"
  "data/metadata/manifest.csv"
  "data/raw/.gitkeep"
)
mkdir -p docs config scripts/01_data results/tables figures data/metadata
for f in "${KEEPERS[@]}"; do
  created=0
  if [ ! -f "$f" ]; then printf 'placeholder' > "$f"; created=1; fi
  if git check-ignore -q "$f"; then
    printf '  BLOCKED  %s   <-- .gitignore is TOO BROAD\n' "$f"
    fail=1
  else
    printf '  trackable %s\n' "$f"
  fi
  [ "$created" -eq 1 ] && rm -f "$f"
done
rmdir renv/library renv 2>/dev/null || true

echo
if [ "$fail" -ne 0 ]; then
  echo "RESULT: FAIL — fix .gitignore before proceeding. No downloads until this passes." >&2
  exit 1
fi
echo "RESULT: PASS — .gitignore blocks all data decoys and blocks no deliverables."

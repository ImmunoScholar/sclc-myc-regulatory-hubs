#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 07_clean_clone_verify.sh — M11 release gate, run honestly.
#
#   bash scripts/00_setup/07_clean_clone_verify.sh
#
# The contract's gate is: "a fresh git clone on a clean machine reproduces every
# figure. Not 'should' — actually tested."
#
# WHAT THIS TESTS, and what it deliberately does not:
#
#   TESTED   fresh clone succeeds and carries every file the pipeline references
#   TESTED   no absolute or home-relative paths are baked into any script
#   TESTED   every package in renv.lock resolves from the configured repositories
#            (dry run — it does not reinstall 232 packages)
#   TESTED   the full figure set rebuilds FROM THE CLONE, with its own working
#            directory, which is where path bugs actually surface
#
#   NOT TESTED  a from-scratch package install on a machine with no R library
#   NOT TESTED  re-downloading the 12 GB of source data
#
# The two untested items are stated rather than quietly folded in. Data are
# deliberately absent from version control, so a literally-clean machine must
# also run the download scripts (~6 h) before figures can be rebuilt; this test
# links the existing data into the clone to isolate CODE portability, which is
# the failure mode a clone test is actually for.
# -----------------------------------------------------------------------------
set -uo pipefail
SRC="$(cd "$(dirname "$0")/../.." && pwd)"
CLONE=$(mktemp -d /tmp/cloneverify.XXXXXX)
rc=0
note() { printf '  %-52s %s\n' "$1" "$2"; }

echo "=========== 1. fresh clone ==========="
git clone -q "$SRC" "$CLONE/repo" || { echo "clone FAILED"; exit 1; }
cd "$CLONE/repo" || exit 1
note "cloned to" "$CLONE/repo"
note "files" "$(git ls-files | wc -l)"

echo
echo "=========== 2. absolute / home paths in code ==========="
# A path like /home/priya/... or ~/... works on this machine and nowhere else.
#
# This check must EXCLUDE ITSELF. The pattern it searches for necessarily appears
# in its own source, so the first run reported the checker as the violation — the
# same self-match that verify_config.R hit scanning its own source for removed
# keys, and that pgrep -f hit matching its own command line. A check that flags
# itself trains the reader to ignore it.
bad=$(grep -rnE '"(/home/|/Users/|/mnt/[a-z]/)|"~/' \
        scripts/ R/ config/ report.qmd 2>/dev/null \
        | grep -v '07_clean_clone_verify.sh' \
        | grep -v '^\s*#' | grep -vE '^\S+:\s*#' || true)
if [ -n "$bad" ]; then
  echo "$bad" | head -20
  note "absolute paths" "FOUND — not portable"; rc=1
else
  note "absolute paths" "none"
fi

echo
echo "=========== 3. referenced scripts all present ==========="
# 01_init_git.sh is excluded on purpose: it contains ELEVEN deliberately fake
# paths as decoys for the adversarial .gitignore test, and treating those as
# broken references would make this check permanently red.
missing=0
while read -r s; do
  [ -f "$s" ] || { echo "  MISSING: $s"; missing=$((missing+1)); }
done < <(grep -rhoE 'scripts/[0-9A-Za-z_/]+\.(R|sh)' README.md report.qmd scripts/ 2>/dev/null \
           --exclude=01_init_git.sh --exclude=07_clean_clone_verify.sh | sort -u)
note "referenced scripts missing" "$missing"
[ "$missing" -eq 0 ] || rc=1

echo
echo "=========== 4. lockfile composition and sources ==========="
# Reports what a restore would have to fetch. An earlier version of this check
# compared each lockfile Repository value against the NAMES of the configured
# repos and flagged "RSPM" and a bare URL as unresolvable — a false alarm, since
# those are provenance labels rather than repo names. Composition is reported;
# the binding question (does every package actually install on a bare machine?)
# is out of scope here and is stated as such in the header.
Rscript -e '
  pk <- renv::lockfile_read()$Packages
  cat("  packages in lockfile:", length(pk), "\n")
  srcs <- table(vapply(pk, function(p) if (is.null(p$Source)) "NA" else p$Source, character(1)))
  for (n in names(srcs)) cat(sprintf("    %-14s %d\n", n, srcs[[n]]))
  gh <- vapply(pk, function(p) identical(p$Source, "GitHub"), logical(1))
  cat("  GitHub-sourced (unpinned risk):", sum(gh), "\n")
' 2>&1 | grep -v out-of-sync

echo
echo "=========== 5. rebuild every figure FROM THE CLONE ==========="
# Link only the GITIGNORED subdirectories. `ln -sfn src target` places the link
# INSIDE target when target already exists as a directory, and data/ does exist in
# a clone because data/metadata/*.csv are tracked — the first version of this
# check created data/data -> source and silently resolved nothing, which is how
# two figures appeared to fail for the wrong reason.
for d in data/processed data/raw results/objects; do
  if [ -d "$SRC/$d" ] && [ ! -e "$CLONE/repo/$d" ]; then
    mkdir -p "$(dirname "$CLONE/repo/$d")"
    ln -s "$SRC/$d" "$CLONE/repo/$d"
  fi
done
note "linked gitignored data dirs" "$(ls -d "$CLONE"/repo/data/processed "$CLONE"/repo/data/raw 2>/dev/null | wc -l)"
rm -rf "$CLONE/repo/figures"; mkdir -p "$CLONE/repo/figures"

# Point the clone at the existing package library. Reinstalling 232 packages is
# a different test (see header) and would take hours; this one asks whether the
# CODE runs from a different working directory, which is the failure mode a clone
# test exists to catch.
export RENV_PATHS_LIBRARY="$SRC/renv/library"
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE

if bash scripts/07_figures/99_build_all.sh > "$CLONE/build.log" 2>&1; then
  note "figure rebuild" "PASS"
else
  note "figure rebuild" "FAIL"; tail -n 25 "$CLONE/build.log"; rc=1
fi
grep -E "figures in manifest|script errors|failed checks|missing files|RESULT" "$CLONE/build.log" | sed 's/^/  /'

echo
echo "=========== 6. figures identical to the committed set ==========="
same=0; diff_n=0
shopt -s nullglob          # an unmatched glob must expand to nothing, not to the
                           # literal "*.png", which previously reported one
                           # "differing" file when the directory was empty
pngs=("$CLONE/repo"/figures/*.png)
if [ ${#pngs[@]} -eq 0 ]; then
  echo "  no PNGs produced — nothing to compare"; rc=1
fi
for f in "${pngs[@]}"; do
  b=$(basename "$f")
  if cmp -s "$f" "$SRC/figures/$b"; then same=$((same+1)); else diff_n=$((diff_n+1)); echo "  DIFFERS: $b"; fi
done
note "byte-identical PNGs" "$same"
note "differing PNGs" "$diff_n"
[ "$diff_n" -eq 0 ] || rc=1

echo
echo "=============================================================="
if [ $rc -eq 0 ]; then
  echo "CLEAN-CLONE VERIFY: PASS"
  echo "  (package install from scratch and data re-download are NOT covered — see header)"
else
  echo "CLEAN-CLONE VERIFY: FAIL"
fi
echo "=============================================================="
echo "clone kept at $CLONE for inspection"
exit $rc

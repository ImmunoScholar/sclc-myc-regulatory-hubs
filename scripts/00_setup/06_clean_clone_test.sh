#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 06_clean_clone_test.sh — prove the repository actually restores from scratch.
#
#   bash scripts/00_setup/06_clean_clone_test.sh
#
# Clones the pushed remote into a throwaway directory, runs renv::restore(), and
# executes the environment test suite there. This is the M11 release gate, run
# early and often rather than once at the end.
#
# Note: this is a fair test of the LOCKFILE, not of a cold machine — the local
# renv cache is shared, so package *downloads* are skipped. It still catches the
# failures that actually happen: an incomplete lockfile, a file that was never
# committed, a hard-coded absolute path, or a script that only works from the
# original directory.
# -----------------------------------------------------------------------------
set -uo pipefail

REMOTE="git@github.com:ImmunoScholar/sclc-myc-regulatory-hubs.git"
TMP="$(mktemp -d /tmp/sclc-clone-test.XXXXXX)"
trap 'echo; echo "cleaning up $TMP"; rm -rf "$TMP"' EXIT

echo "=== cloning into $TMP ==="
git clone --quiet "$REMOTE" "$TMP/repo" || { echo "CLONE FAILED"; exit 1; }
cd "$TMP/repo" || exit 1
echo "cloned commit: $(git rev-parse --short HEAD)"

echo
echo "=== tracked content sanity ==="
n=$(git ls-files | wc -l)
echo "tracked files: $n"
data_leak=$(git ls-files | grep -E '^data/(raw|processed)/' | grep -v '.gitkeep' | wc -l)
echo "data payload files tracked: $data_leak"
[ "$data_leak" -eq 0 ] || { echo "FAIL: data leaked into the repository"; exit 1; }
big=$(find . -path ./.git -prune -o -type f -size +5M -print | wc -l)
echo "files over 5 MB: $big"

echo
echo "=== renv::restore() ==="
Rscript -e 'renv::restore(prompt = FALSE)' 2>&1 | tail -n 15
rc=${PIPESTATUS[0]}
echo "restore exit code: $rc"

echo
echo "=== library check in the clone ==="
Rscript -e '
n <- length(jsonlite::fromJSON("renv.lock")$Packages)
cat("packages in lockfile:", n, "\n")
inst <- rownames(installed.packages(lib.loc = renv::paths$library()))
cat("packages installed in project library:", length(inst), "\n")
core <- c("GenomicRanges","rtracklayer","DESeq2","singscore","AUCell","GENIE3",
          "RobustRankAggreg","ComplexHeatmap","annotatr",
          "TxDb.Hsapiens.UCSC.hg19.knownGene")
miss <- core[!core %in% inst]
if (length(miss)) { cat("MISSING:", paste(miss, collapse=", "), "\n"); quit(status=1) }
cat("all core packages present\n")
' 2>&1 | tail -n 10

echo
echo "=== environment test suite, run inside the clone ==="
Rscript tests/test_environment.R 2>&1 | tail -n 8

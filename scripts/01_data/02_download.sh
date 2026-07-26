#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 02_download.sh — resumable, size-verified acquisition of every manifest file.
#
#   bash scripts/01_data/02_download.sh [--dataset ID] [--jobs N] [--dry-run] [--force]
#
# Guarantees:
#   * RESUMABLE      — curl -C - continues partial files. Every endpoint confirmed
#                      to return Accept-Ranges: bytes / HTTP 206.
#   * IDEMPOTENT     — a file already at its exact expected size is skipped.
#   * FAILS LOUDLY   — wrong final size is an error, never a warning. Partial
#                      files are kept as .part so the next run resumes rather
#                      than restarting 8.6 GB.
#   * INTEGRITY      — every .gz is gzip -t tested; catches truncation that a
#                      correct byte count can mask on a resumed transfer.
#
# CONCURRENCY (default 2) — arrived at the hard way; see decision D-018.
#
#   serial          : ~171 KB/s mean, 0 failures, ~19 h projected for 12 GB
#   4 concurrent    : 403 rate-limit responses within ~2 minutes, on 4 of the
#                     first 8 files, three of which left 980-byte HTML error
#                     pages sitting in .part files
#
# A bandwidth probe had suggested an 8.75x aggregate gain from 4 streams, but it
# measured byte-RANGE requests against a single file, which does not trip NCBI's
# per-IP limit the way four requests for four distinct files does. The probe was
# correct about available bandwidth and wrong about what was constraining us.
#
# So: 2 jobs, staggered starts, and a worker that treats 403 as transient with
# exponential backoff while discarding any .part that is really an error page.
# Raising this is not a free speed-up — it trades throughput for 403s.
#
# Reads:  data/metadata/download_list.tsv
# Writes: data/raw/<dataset_id>/<file>   (git-ignored)
#         logs/02_download.log
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

LIST=data/metadata/download_list.tsv
LOG=logs/02_download.log
RESULT_DIR=logs/dl_results
ONLY=""; DRY=0; FORCE=0; JOBS=2

while [ $# -gt 0 ]; do
  case "$1" in
    --dataset) ONLY="$2"; shift 2 ;;
    --jobs)    JOBS="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    --force)   FORCE=1; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

[ -f "$LIST" ] || { echo "FATAL: $LIST not found. Run 01_build_manifest.R first." >&2; exit 1; }
mkdir -p logs "$RESULT_DIR"
rm -f "$RESULT_DIR"/*

log () { printf '%s\n' "$*" | tee -a "$LOG"; }

log "==============================================================="
log "download run: $(date -u '+%Y-%m-%dT%H:%M:%SZ')  jobs=$JOBS"
[ -n "$ONLY" ] && log "restricted to dataset: $ONLY"
[ "$DRY" -eq 1 ] && log "DRY RUN — nothing will be fetched"
log "==============================================================="

# --- build the worklist -------------------------------------------------------
WORK=$(mktemp); trap 'rm -f "$WORK"' EXIT
n_total=0; n_skip=0; bytes_todo=0

while IFS=$'\t' read -r ds fname url size dest; do
  [ "$ds" = "dataset_id" ] && continue
  [ -z "${ds:-}" ] && continue
  [ -n "$ONLY" ] && [ "$ds" != "$ONLY" ] && continue
  n_total=$((n_total + 1))

  if [ -f "$dest" ] && [ "$FORCE" -eq 0 ]; then
    have=$(stat -c%s "$dest")
    if [ "$have" = "$size" ]; then n_skip=$((n_skip + 1)); continue; fi
  fi
  printf '%s\t%s\t%s\t%s\t%s\n' "$ds" "$fname" "$url" "$size" "$dest" >> "$WORK"
  bytes_todo=$((bytes_todo + size))
done < "$LIST"

n_queued=$(wc -l < "$WORK" | tr -d ' ')
log "considered : $n_total"
log "complete   : $n_skip"
log "queued     : $n_queued  ($(numfmt --to=iec "$bytes_todo" 2>/dev/null || echo "$bytes_todo"))"
log ""

if [ "$DRY" -eq 1 ]; then
  while IFS=$'\t' read -r ds fname url size dest; do
    log "  WOULD FETCH  $ds/$fname  ($size bytes)"
  done < "$WORK"
  log "RESULT: dry run complete."
  exit 0
fi

if [ "$n_queued" -eq 0 ]; then
  log "RESULT: PASS — nothing to do; all requested files already complete."
  exit 0
fi

# --- run workers --------------------------------------------------------------
export DL_RESULT_DIR="$RESULT_DIR"
start=$(date +%s)

# xargs, not GNU parallel: xargs is in coreutils and needs no extra dependency.
# -P runs JOBS workers; each worker is its own process writing its own result
# file, so concurrent progress cannot corrupt a shared log.
tr '\t' '\n' < "$WORK" | xargs -d '\n' -n 5 -P "$JOBS" \
  bash scripts/01_data/_dl_one.sh 2>>"$LOG"

elapsed=$(( $(date +%s) - start ))

# --- collect results ----------------------------------------------------------
n_ok=0; n_fail=0; n_sk=0; bytes_got=0
FAILED=()
for f in "$RESULT_DIR"/*; do
  [ -f "$f" ] || continue
  IFS=$'\t' read -r st ds fname detail < "$f"
  case "$st" in
    OK)   n_ok=$((n_ok + 1)); bytes_got=$((bytes_got + ${detail:-0})); log "  ok    $ds/$fname" ;;
    SKIP) n_sk=$((n_sk + 1)) ;;
    FAIL) n_fail=$((n_fail + 1)); FAILED+=("$ds/$fname :: $detail"); log "  FAIL  $ds/$fname — $detail" ;;
  esac
done

log ""
log "==============================================================="
log "fetched    : $n_ok"
log "skipped    : $((n_skip + n_sk))"
log "failed     : $n_fail"
log "downloaded : $(numfmt --to=iec "$bytes_got" 2>/dev/null || echo "$bytes_got")"
log "elapsed    : $((elapsed / 60))m $((elapsed % 60))s"
if [ "$elapsed" -gt 0 ] && [ "$bytes_got" -gt 0 ]; then
  log "mean rate  : $(numfmt --to=iec $(( bytes_got / elapsed )) 2>/dev/null)/s (aggregate, $JOBS jobs)"
fi
log "==============================================================="

if [ "$n_fail" -gt 0 ]; then
  log ""
  log "FAILURES:"
  for f in "${FAILED[@]}"; do log "  - $f"; done
  log ""
  log "RESULT: FAIL — re-run to resume. Downstream analysis must not proceed."
  exit 1
fi

log "RESULT: PASS — all requested files present at their expected size."
log "Next: Rscript scripts/01_data/03_verify.R"

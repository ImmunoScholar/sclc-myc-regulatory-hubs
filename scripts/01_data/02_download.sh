#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 02_download.sh — resumable, size-verified acquisition of every manifest file.
#
#   bash scripts/01_data/02_download.sh [--dataset ID] [--dry-run] [--force]
#
# Guarantees:
#   * RESUMABLE      — curl -C - continues a partial file. Every endpoint was
#                      confirmed to return Accept-Ranges: bytes / HTTP 206.
#   * IDEMPOTENT     — a file already present at the exact expected byte size is
#                      skipped, so re-running costs nothing.
#   * FAILS LOUDLY   — a wrong final size is an error, never a warning. Partial
#                      files are left in place (named .part) so the next run
#                      resumes rather than restarting 8.6 GB.
#   * INTEGRITY      — every .gz is tested with `gzip -t` after download, which
#                      catches truncation that a size check alone can miss.
#
# Deliberately NOT parallel. NCBI rate-limits bursts and starts returning 403
# HTML error pages that a naive downloader will happily save as "data".
#
# Reads:  data/metadata/download_list.tsv   (from 01_build_manifest.R)
# Writes: data/raw/<dataset_id>/<file>      (all git-ignored)
#         logs/02_download.log
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

LIST=data/metadata/download_list.tsv
LOG=logs/02_download.log
ONLY=""; DRY=0; FORCE=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dataset) ONLY="$2"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    --force)   FORCE=1; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

[ -f "$LIST" ] || { echo "FATAL: $LIST not found. Run 01_build_manifest.R first." >&2; exit 1; }
mkdir -p logs

log () { printf '%s\n' "$*" | tee -a "$LOG"; }

log "==============================================================="
log "download run: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
[ -n "$ONLY" ] && log "restricted to dataset: $ONLY"
[ "$DRY" -eq 1 ] && log "DRY RUN — nothing will be fetched"
log "==============================================================="

n_total=0; n_skip=0; n_get=0; n_fail=0; bytes_got=0
FAILED_FILES=()

# Skip the header; fields: dataset_id, file_name, source_url, size_expected, dest
while IFS=$'\t' read -r ds fname url size dest; do
  [ "$ds" = "dataset_id" ] && continue
  [ -z "${ds:-}" ] && continue
  if [ -n "$ONLY" ] && [ "$ds" != "$ONLY" ]; then continue; fi
  n_total=$((n_total + 1))

  mkdir -p "$(dirname "$dest")"

  # --- already complete? ------------------------------------------------------
  if [ -f "$dest" ] && [ "$FORCE" -eq 0 ]; then
    have=$(stat -c%s "$dest")
    if [ "$have" = "$size" ]; then
      n_skip=$((n_skip + 1))
      continue
    fi
    log "  size differs for $fname (have $have, want $size) — resuming"
  fi

  if [ "$DRY" -eq 1 ]; then
    log "  WOULD FETCH  $ds/$fname  ($size bytes)"
    n_get=$((n_get + 1))
    continue
  fi

  # --- fetch (resumable) ------------------------------------------------------
  log "  fetching $ds/$fname ($(numfmt --to=iec "$size" 2>/dev/null || echo "$size") )"
  part="${dest}.part"
  # If a completed-but-wrong-size file exists, move it aside to resume against.
  [ -f "$dest" ] && [ ! -f "$part" ] && mv "$dest" "$part"

  http=$(curl -4 -sS -L \
      --retry 5 --retry-delay 10 --retry-connrefused \
      --connect-timeout 30 --max-time 7200 \
      --speed-limit 1024 --speed-time 120 \
      -C - -o "$part" \
      -w '%{http_code}' \
      -A "sclc-myc-regulatory-hubs/1.0 (academic research; contact via repository)" \
      "$url" 2>>"$LOG") || true

  # 416 = range beyond end, i.e. the .part is already complete. Not an error.
  if [ "$http" != "200" ] && [ "$http" != "206" ] && [ "$http" != "416" ]; then
    log "  FAIL  $fname — HTTP $http"
    FAILED_FILES+=("$ds/$fname (HTTP $http)")
    n_fail=$((n_fail + 1))
    continue
  fi

  # --- verify size ------------------------------------------------------------
  got=$(stat -c%s "$part" 2>/dev/null || echo 0)
  if [ "$got" != "$size" ]; then
    log "  FAIL  $fname — size mismatch: got $got, expected $size (kept as .part for resume)"
    FAILED_FILES+=("$ds/$fname (size $got != $size)")
    n_fail=$((n_fail + 1))
    continue
  fi

  # --- verify gzip integrity --------------------------------------------------
  case "$fname" in
    *.gz)
      if ! gzip -t "$part" 2>>"$LOG"; then
        log "  FAIL  $fname — gzip integrity test failed (corrupt despite correct size)"
        FAILED_FILES+=("$ds/$fname (gzip -t failed)")
        n_fail=$((n_fail + 1))
        rm -f "$part"       # a size-correct but corrupt file must not be resumed
        continue
      fi
      ;;
  esac

  mv "$part" "$dest"
  n_get=$((n_get + 1))
  bytes_got=$((bytes_got + got))
  log "  ok    $fname"
done < "$LIST"

log ""
log "==============================================================="
log "considered : $n_total"
log "skipped    : $n_skip (already complete)"
log "fetched    : $n_get"
log "failed     : $n_fail"
log "downloaded : $(numfmt --to=iec "$bytes_got" 2>/dev/null || echo "$bytes_got")"
log "==============================================================="

if [ "$n_fail" -gt 0 ]; then
  log ""
  log "FAILURES:"
  for f in "${FAILED_FILES[@]}"; do log "  - $f"; done
  log ""
  log "RESULT: FAIL — re-run to resume. Downstream analysis must not proceed."
  exit 1
fi

if [ "$DRY" -eq 1 ]; then
  log "RESULT: dry run complete."
  exit 0
fi

log "RESULT: PASS — all requested files present at their expected size."
log "Next: Rscript scripts/01_data/03_verify.R"

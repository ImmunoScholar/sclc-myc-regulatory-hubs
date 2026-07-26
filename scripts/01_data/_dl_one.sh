#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# _dl_one.sh — fetch and validate ONE file. Worker for 02_download.sh.
#
#   bash scripts/01_data/_dl_one.sh <dataset> <file_name> <url> <expected_size> <dest>
#
# Not called directly. Separate script because parallel workers must be separate
# processes, and because each writes its own result file — interleaved concurrent
# writes to a shared log produce corrupted lines that look like data errors.
#
# Guarantees:
#   resumable (curl -C -) · size-verified · gzip-tested · atomic rename on success
#
# HANDLES NCBI's 403 RATE LIMITING. This is the failure that actually happens
# under concurrency, and it is nasty for two reasons:
#   1. curl does NOT retry 403 by default (it is not classed as transient), so a
#      single rate-limit response killed the file outright.
#   2. the 403 HTML BODY gets written into the .part file. A later `curl -C -`
#      resumes from the end of that markup, yielding a file that is a few hundred
#      bytes of HTML followed by real data. Observed: 980-byte .part files whose
#      first two bytes were 3c3f ("<?") instead of the gzip magic 1f8b.
# So: any .part that is not valid gzip, or begins with markup, is DISCARDED
# before retrying, and 403 is retried with exponential backoff.
# -----------------------------------------------------------------------------
set -uo pipefail

ds="$1"; fname="$2"; url="$3"; size="$4"; dest="$5"
RESULT_DIR="${DL_RESULT_DIR:-logs/dl_results}"
MAX_TRIES="${DL_MAX_TRIES:-5}"
mkdir -p "$RESULT_DIR" "$(dirname "$dest")"

key=$(printf '%s' "$dest" | md5sum | cut -d' ' -f1)
res="$RESULT_DIR/$key"
emit () { printf '%s\t%s\t%s\t%s\n' "$1" "$ds" "$fname" "${2:-}" > "$res"; }

part="${dest}.part"

# Discard a .part that cannot be a genuine partial transfer.
scrub_part () {
  [ -f "$part" ] || return 0
  local sz magic head_txt
  sz=$(stat -c%s "$part")
  magic=$(head -c 2 "$part" | od -An -tx1 | tr -d ' \n')
  head_txt=$(head -c 128 "$part" 2>/dev/null | tr -d '\0')
  case "$fname" in
    *.gz)
      if [ "$magic" != "1f8b" ]; then rm -f "$part"; return 1; fi ;;
  esac
  case "$head_txt" in
    '<!DOCTYPE'*|'<html'*|'<HTML'*|'<?xml'*) rm -f "$part"; return 1 ;;
  esac
  if [ "$sz" -lt 1024 ]; then rm -f "$part"; return 1; fi
  return 0
}

# --- already complete? --------------------------------------------------------
if [ -f "$dest" ]; then
  have=$(stat -c%s "$dest")
  if [ "$have" = "$size" ]; then emit SKIP "already complete"; exit 0; fi
  [ -f "$part" ] || mv "$dest" "$part"
fi
scrub_part || true

# --- attempt loop -------------------------------------------------------------
attempt=1
backoff=20
last=""
while [ "$attempt" -le "$MAX_TRIES" ]; do

  http=$(curl -4 -sS -L \
      --retry 3 --retry-delay 8 --retry-connrefused \
      --connect-timeout 30 --max-time 10800 \
      --speed-limit 512 --speed-time 300 \
      -C - -o "$part" \
      -w '%{http_code}' \
      -A "sclc-myc-regulatory-hubs/1.0 (academic research; contact via repository)" \
      "$url" 2>/dev/null) || true

  case "$http" in
    200|206|416)
      got=$(stat -c%s "$part" 2>/dev/null || echo 0)
      if [ "$got" = "$size" ]; then
        case "$fname" in
          *.gz)
            if ! gzip -t "$part" 2>/dev/null; then
              rm -f "$part"      # size-correct but corrupt: never resume this
              last="gzip -t failed (corrupt; .part discarded)"
              attempt=$((attempt + 1)); sleep "$backoff"; backoff=$((backoff * 2))
              continue
            fi ;;
        esac
        mv "$part" "$dest"
        emit OK "$got"
        exit 0
      fi
      last="size $got != expected $size"
      ;;
    403|429|503)
      # Rate limiting. The response body is now sitting in .part; remove it.
      scrub_part || true
      last="HTTP $http (rate limited)"
      ;;
    *)
      scrub_part || true
      last="HTTP $http"
      ;;
  esac

  attempt=$((attempt + 1))
  if [ "$attempt" -le "$MAX_TRIES" ]; then
    sleep "$backoff"
    backoff=$((backoff * 2))
  fi
done

emit FAIL "$last after $MAX_TRIES attempts"
exit 1

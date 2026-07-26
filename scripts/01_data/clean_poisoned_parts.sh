#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# clean_poisoned_parts.sh — remove .part files that are HTTP error pages.
#
# A 403 response body gets written to the .part file. curl -C - would then resume
# FROM THE END OF THAT HTML, producing a file that is a few hundred bytes of
# markup followed by real data. The size check would eventually reject it, but
# only after re-downloading hundreds of megabytes, and a partial resume chain
# could in principle land on the right size with wrong leading bytes.
#
# Any .part that begins with markup, or is implausibly small for a genomics
# file, is deleted so the next attempt starts clean.
#
# Usage: bash scripts/01_data/clean_poisoned_parts.sh [--apply]
# Without --apply it only reports.
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

echo "=== scanning .part files ==="
n_bad=0; n_ok=0
while IFS= read -r p; do
  sz=$(stat -c%s "$p")
  head_txt=$(head -c 256 "$p" 2>/dev/null | tr -d '\0')
  is_html=0
  case "$head_txt" in
    '<!DOCTYPE'*|'<html'*|'<HTML'*|'<?xml'*|*'Access forbidden'*|*'<title>'*) is_html=1 ;;
  esac
  # gzip members begin with 0x1f 0x8b; a real .gz partial must start that way.
  gz_magic=$(head -c 2 "$p" | od -An -tx1 | tr -d ' \n')
  expect_gz=0
  case "$p" in *.gz.part) expect_gz=1 ;; esac
  bad_magic=0
  if [ "$expect_gz" -eq 1 ] && [ "$gz_magic" != "1f8b" ]; then bad_magic=1; fi

  if [ "$is_html" -eq 1 ] || [ "$bad_magic" -eq 1 ] || [ "$sz" -lt 4096 ]; then
    reason=""
    [ "$is_html" -eq 1 ]   && reason="HTML error page"
    [ "$bad_magic" -eq 1 ] && reason="${reason:+$reason; }not gzip (magic=$gz_magic)"
    [ "$sz" -lt 4096 ]     && reason="${reason:+$reason; }implausibly small ($sz bytes)"
    echo "  POISONED  $p"
    echo "            reason: $reason"
    n_bad=$((n_bad + 1))
    [ "$APPLY" -eq 1 ] && rm -f "$p" && echo "            deleted"
  else
    echo "  ok        $(numfmt --to=iec "$sz")  $(basename "$p")"
    n_ok=$((n_ok + 1))
  fi
done < <(find data/raw -name '*.part')

echo
echo "resumable partials : $n_ok"
echo "poisoned partials  : $n_bad"
[ "$APPLY" -eq 0 ] && [ "$n_bad" -gt 0 ] && echo "(re-run with --apply to delete)"

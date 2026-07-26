#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# diag_bandwidth.sh — is the download slow because of PER-CONNECTION throttling
# or because the link itself is saturated?
#
# This distinction decides whether parallelism is worth doing:
#   * if 4 streams give ~4x the aggregate of 1 stream -> per-connection limit,
#     parallelism fixes it
#   * if aggregate stays flat -> the link (or the route) is the ceiling and
#     parallelism just splits the same bandwidth into smaller pieces
#
# Uses byte-range requests against a real target file, time-boxed, so it
# measures without downloading anything we keep.
#
# Read-only. Usage: bash scripts/01_data/diag_bandwidth.sh
# -----------------------------------------------------------------------------
set -uo pipefail

URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7230nnn/GSM7230496/suppl/GSM7230496_H3K27Ac_H211_R1_treat_pileup.bedgraph.gz"
SECS=15
CHUNK=$((60 * 1024 * 1024))   # 60 MB ceiling per stream; time is the real limit

probe_stream () {   # $1 = byte offset, $2 = output marker
  curl -4 -sS -o /dev/null \
    --max-time "$SECS" \
    -r "${1}-$(( $1 + CHUNK ))" \
    -w '%{size_download}\n' \
    "$URL" 2>/dev/null
}

echo "NOTE: the main download is still running, so these figures are measured"
echo "alongside it and understate what a free link would give."
echo

echo "=== 1 stream, ${SECS}s ==="
s1=$(probe_stream 0)
r1=$(( ${s1:-0} / SECS ))
echo "  bytes: ${s1:-0}   rate: $(numfmt --to=iec "$r1")/s"

echo
echo "=== 4 concurrent streams, ${SECS}s (distinct offsets) ==="
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
off=0
for i in 1 2 3 4; do
  probe_stream "$off" > "$tmp/s$i" &
  off=$(( off + 70 * 1024 * 1024 ))
done
wait
agg=0
for i in 1 2 3 4; do
  v=$(cat "$tmp/s$i" 2>/dev/null | tr -d '[:space:]')
  v=${v:-0}
  echo "  stream $i: $(numfmt --to=iec $(( v / SECS )))/s"
  agg=$(( agg + v ))
done
r4=$(( agg / SECS ))
echo "  AGGREGATE : $(numfmt --to=iec "$r4")/s"

echo
echo "=== verdict ==="
if [ "$r1" -gt 0 ]; then
  ratio=$(awk -v a="$r4" -v b="$r1" 'BEGIN{printf "%.2f", a/b}')
  echo "  aggregate / single = ${ratio}x"
  awk -v x="$ratio" 'BEGIN{
    if (x > 2.2)      print "  -> PER-CONNECTION throttling dominates. Parallelism will help substantially.";
    else if (x > 1.4) print "  -> partial gain available. Modest parallelism (3-4) is worthwhile.";
    else              print "  -> link or route is the ceiling. Parallelism will NOT help; do not bother.";
  }'
fi

echo
echo "=== route sanity: a non-NCBI host for comparison ==="
for u in "https://packagemanager.posit.co/cran/__linux__/noble/latest/src/contrib/PACKAGES.gz" \
         "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"; do
  out=$(curl -4 -sS -o /dev/null --max-time 25 -w '%{speed_download} %{size_download}' "$u" 2>/dev/null)
  sp=$(echo "$out" | awk '{printf "%d", $1}')
  echo "  $(numfmt --to=iec "${sp:-0}")/s  <- $(echo "$u" | cut -c1-58)"
done

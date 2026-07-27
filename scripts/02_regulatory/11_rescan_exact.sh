#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 11_rescan_exact.sh — exact per-line region calls at four thresholds.
#
#   bash scripts/02_regulatory/11_rescan_exact.sh
#
# Replaces the approximation in 09_. There, one scan at x2 was re-thresholded by
# each run's MEAN signal to derive the stricter sets. That keeps whole runs
# including their low-signal flanks, so x2.5/x3/x4 came out systematically wider
# than true thresholding — which invalidated the threshold sensitivity analysis,
# the very thing meant to make the conclusion robust.
#
# This does the real thing: a single pass per file maintaining FOUR independent
# run accumulators, one per threshold, writing each to its own file. Nine passes
# rather than thirty-six.
#
# Thresholds are per-line multiples of that line's mean signal over covered bases
# (read from the histograms), so differing depth cannot let deep tracks dominate.
#
# Output: data/processed/regions/exact/<line>_x<mult>.bed
# -----------------------------------------------------------------------------
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1

OUT=data/processed/regions/exact
mkdir -p "$OUT"

MULTS=(2.5 3 4 5)
CHROMS='($1=="1"||$1=="2"||$1=="3"||$1=="4"||$1=="5"||$1=="6"||$1=="7"||$1=="8"||$1=="9"||$1=="10"||$1=="11"||$1=="12"||$1=="13"||$1=="14"||$1=="15"||$1=="16"||$1=="17"||$1=="18"||$1=="19"||$1=="20"||$1=="21"||$1=="22"||$1=="X")'

# mean signal over covered bases, from the cached histogram
line_mean () {
  awk -F'\t' '!/^#/ { s=$1/2; if (s>0) { num+=s*$2; den+=$2 } } END { printf "%.6f", num/den }' \
    "data/processed/atac_profile/$1.hist.tsv"
}

for f in data/raw/GSE230649/GSM*_treat_pileup.bedgraph.gz; do
  base=$(basename "$f")
  # ATAC tracks only: no target token after the GSM id
  case "$base" in
    GSM*_H3K27Ac_*|GSM*_MYC_*|GSM*_MYCN_*|GSM*_MYCL1_*) continue ;;
  esac
  line=$(echo "$base" | sed -E 's/^GSM[0-9]+_([^_]+)_treat_pileup.*/\1/')

  # skip if every output already exists
  done_all=1
  for m in "${MULTS[@]}"; do
    [ -s "$OUT/${line}_x${m}.bed" ] || done_all=0
  done
  if [ "$done_all" -eq 1 ]; then echo "  [cached] $line"; continue; fi

  gm=$(line_mean "$line")
  t1=$(awk -v g="$gm" 'BEGIN{printf "%.6f", g*2.5}')
  t2=$(awk -v g="$gm" 'BEGIN{printf "%.6f", g*3}')
  t3=$(awk -v g="$gm" 'BEGIN{printf "%.6f", g*4}')
  t4=$(awk -v g="$gm" 'BEGIN{printf "%.6f", g*5}')

  start=$(date +%s)
  zcat "$f" | awk -F'\t' \
    -v T1="$t1" -v T2="$t2" -v T3="$t3" -v T4="$t4" \
    -v O1="$OUT/${line}_x2.5.bed" -v O2="$OUT/${line}_x3.bed" \
    -v O3="$OUT/${line}_x4.bed"   -v O4="$OUT/${line}_x5.bed" \
    'BEGIN { t[1]=T1; t[2]=T2; t[3]=T3; t[4]=T4;
             o[1]=O1; o[2]=O2; o[3]=O3; o[4]=O4 }
     '"$CHROMS"' {
       v=$4+0
       for (j=1;j<=4;j++) {
         if (v >= t[j]) {
           if (open[j] && cc[j]==$1 && ce[j]==$2) { ce[j]=$3 }
           else {
             if (open[j]) print cc[j]"\t"cs[j]"\t"ce[j] > o[j]
             cc[j]=$1; cs[j]=$2; ce[j]=$3; open[j]=1
           }
         } else if (open[j]) {
           print cc[j]"\t"cs[j]"\t"ce[j] > o[j]; open[j]=0
         }
       }
     }
     END { for (j=1;j<=4;j++) if (open[j]) print cc[j]"\t"cs[j]"\t"ce[j] > o[j] }'

  el=$(( $(date +%s) - start ))
  printf '  %-9s mean=%s  ' "$line" "$gm"
  for m in "${MULTS[@]}"; do
    printf 'x%s:%s  ' "$m" "$(wc -l < "$OUT/${line}_x${m}.bed")"
  done
  printf ' (%dm%02ds)\n' $((el/60)) $((el%60))
done

echo
echo "=== output ==="
ls -1 "$OUT" | wc -l | xargs echo "files:"
du -sh "$OUT" | cut -f1 | xargs echo "size :"
echo
echo "RESULT: exact scans complete."
echo "Next: Rscript scripts/02_regulatory/12_finalise_universe.R"

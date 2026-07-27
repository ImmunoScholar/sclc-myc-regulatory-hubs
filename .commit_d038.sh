#!/usr/bin/env bash
set -e
cd /home/priya/projects/sclc-myc-regulatory-hubs

rm -f .run_s01.sh .run_figs.sh

echo "=== clock check (host vs WSL) ==="
date

echo "=== status ==="
git status --short

echo "=== figure manifest ==="
cat figures/figure_manifest.csv

cat > .git/COMMIT_MSG_TMP <<'EOF'
M10: replace paralog palette after it failed its own accessibility check (D-038)

The colourblind check written at M4 to enforce figures$colourblind_check was
never actually run until the first results figure. It failed immediately on the
project's own palette: min CIE dE 12.7 under deuteranopia against a floor of 15,
MYC #B2182B against MYCL1 #1B7837. Red vs green, carried since M1, already
shipped in the M4 inventory figure.

Replacement selected by score rather than taste, in a committed re-runnable
script, on two criteria: the configured colourblind floor AND >= 30 L*
separation from the white panel. The second was added after the first ranking
returned a pale gold on a dE of 44.3 -- maximising pairwise distance alone
rewards near-white tints, and the paralogs are drawn as thin lines, small points
and coloured text. Three of eleven candidates pass the colourblind test and fail
the background test, including the top scorer.

Selected MYC #762A83, MYCN #E08214, MYCL1 #1B7837: dE 30.5, background 36.7.
No scientific result changes; colour carries no quantitative claim.

Also fixed, all exposed by this:
- the diverging heatmap scale was a hardcoded hex in the plotting script, where
  the accessibility check could not see it; now config-driven
- panel A borrowed a paralog colour for a variance-source encoding
- fig_s01 never ran the accessibility gate at all; it does now, and PAL_AMP
  passes at 29.1 having never previously been tested

Four rendering defects passed all automated checks and were caught only by
inspecting the PNG: clipped main title, colliding panel titles, panel C labels
drawn white-on-white outside downward bars, and panel A labels collapsing to the
group centre because subsetting the layer left position_dodge one fill level.
The checks verify the palette and the theme, not the finished plot.
EOF

git add -A
git commit -F .git/COMMIT_MSG_TMP
rm -f .git/COMMIT_MSG_TMP
echo "=== committed ==="
git log --oneline -3

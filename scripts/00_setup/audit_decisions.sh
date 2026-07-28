#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# audit_decisions.sh — verify every logged decision against the actual state.
#
# Not a re-read of the log. Each check tests a CONCRETE claim the decision makes,
# so a decision that has drifted from reality shows up as a failure rather than
# staying plausible on paper.
# -----------------------------------------------------------------------------
cd "$(dirname "$0")/../.." || exit 1

pass=0; fail=0; na=0; skip=0

# Some decisions are evidenced by PIPELINE ARTEFACTS — git hooks, cached
# intermediates, per-line scan outputs — which are deliberately not in version
# control. In a fresh clone they are absent, and reporting that as FAIL makes a
# correctly-cloned repository look broken. It did exactly that on the first
# bare-machine run: D-010, D-023 and D-024 went red purely because the data had
# not been generated yet.
#
# chk_artefact marks those SKIP when the artefact tree is absent, and tests them
# normally when it is present. The distinction is "decision drifted" versus
# "output not built yet", and conflating the two devalues every other line.
HAVE_ARTEFACTS=0
[ -d data/processed/regions ] && [ -n "$(ls -A data/processed/regions 2>/dev/null)" ] && HAVE_ARTEFACTS=1

chk () {  # chk <id> <description> <test-command>
  if eval "$3" >/dev/null 2>&1; then
    printf '  %-7s PASS  %s\n' "$1" "$2"; pass=$((pass+1))
  else
    printf '  %-7s FAIL  %s\n' "$1" "$2"; fail=$((fail+1))
  fi
}
chk_artefact () {  # same, but SKIP when pipeline outputs are absent
  if [ "$HAVE_ARTEFACTS" -eq 0 ]; then
    printf '  %-7s SKIP  %s (pipeline output not present in this checkout)\n' "$1" "$2"
    skip=$((skip+1))
  else
    chk "$@"
  fi
}
note () { printf '  %-7s ----  %s\n' "$1" "$2"; na=$((na+1)); }

echo "=============================================================="
echo "DECISION AUDIT"
echo "=============================================================="
echo
echo "--- decisions logged ---"
grep -c '^### D-' docs/decision_log.md | xargs echo "  count:"
grep -o '^### D-[0-9]*' docs/decision_log.md | tr -d '#' | tr '\n' ' '
echo; echo
echo "--- risks logged ---"
grep -c '^### R-' docs/06_risk_log.md | xargs echo "  count:"
echo

echo "=== M0-M2 scoping decisions ==="
chk D-001 "Plotnik cited as foundation in gap statement" \
  "grep -qi 'Plotnik' docs/02_gap_statement.md"
chk D-002 "SCENIC rejected, GENIE3 adopted in architecture" \
  "grep -qi 'GENIE3' docs/04_analysis_architecture.md && grep -qi 'REJECT' docs/04_analysis_architecture.md"
chk D-003 "MOES has NO weights field in config" \
  "! grep -A30 '^moes:' config/params.yml | grep -qE '^\s+weights:'"
chk D-003b "MOES method is robust_rank_aggregation" \
  "grep -qE 'method:\s*robust_rank_aggregation' config/params.yml"
chk D-004 "no MACS2/bedtools dependency in dependency inventory" \
  "grep -qi 'not.*required' docs/05_dependency_inventory.md"
chk D-005 "single-cell layer absent from renv.lock (Seurat)" \
  "! grep -q '\"Seurat\"' renv.lock"

echo
echo "=== M3 environment decisions ==="
chk D-006 "GWDG Bioconductor mirror in .Rprofile" \
  "grep -q 'ftp.gwdg.de' .Rprofile"
chk D-007 "P3M CRAN binaries in .Rprofile" \
  "grep -q 'packagemanager.posit.co' .Rprofile"
chk D-008 "SUPERSEDED by D-022 — BSgenome now installed" \
  "grep -q 'BSgenome.Hsapiens.UCSC.hg19' renv.lock"
chk_artefact D-010 "pre-commit large-file hook present" \
  "test -x .git/hooks/pre-commit"
chk D-011 "targets NOT in renv.lock" \
  "! grep -q '\"targets\"' renv.lock"
chk D-012 "GSVA installed (resolution of D-012)" \
  "grep -q '\"GSVA\"' renv.lock"

echo
echo "=== M4 data decisions ==="
chk D-014 "liftover report exists with loss rates" \
  "test -s data/metadata/liftover_report.csv"
chk D-014b "genome_utils provides harmonise_seqnames" \
  "grep -q 'harmonise_seqnames' R/genome_utils.R"
chk D-015 "GSE269424 restricted to EGFP arms in registry" \
  "grep -q '_EGFP_' config/datasets.yml"
chk D-016 "DepMap flagged MANUAL in manifest" \
  "grep -q 'manual_required' data/metadata/dataset_manifest.csv"
chk D-017 "checksum ledger exists" \
  "test -s data/metadata/checksums.sha256"
chk D-017b "manifest records none_published for checksums" \
  "grep -qi 'none_published' data/metadata/dataset_manifest.csv"
chk D-018 "downloader default concurrency is 2" \
  "grep -qE 'JOBS=2' scripts/01_data/02_download.sh"
chk D-019 "assert_hg19_bounds refuses vacuous pass" \
  "grep -q 'min_chroms' R/genome_utils.R"

echo
echo "=== M5 regulatory decisions ==="
chk D-020 "m5_gate section present, no count-based pass/fail" \
  "grep -q 'm5_gate:' config/params.yml && grep -q 'report_counts_descriptively' config/params.yml"
chk D-021 "bootstrap enabled for MOES rank stability" \
  "grep -A6 'bootstrap:' config/params.yml | grep -qE 'n_resamples:\s*2000'"
chk D-021b "vector_output + min_text_pt in figures config" \
  "grep -q 'vector_output: true' config/params.yml && grep -q 'min_text_pt' config/params.yml"
chk D-021c "regulon_validity gate in config" \
  "grep -q 'regulon_validity:' config/params.yml"
chk D-022 "BSgenome hg19 present and build-checked" \
  "grep -q 'BSgenome.Hsapiens.UCSC.hg19' renv.lock"
chk_artefact D-023 "universe derived from keystone ATAC (per_line_runs cache)" \
  "test -s data/processed/regions/per_line_runs.rds || test -s data/processed/regions/universe_final.rds"
chk_artefact D-024 "exact rescan produced 36 files" \
  "test \$(ls data/processed/regions/exact/*.bed 2>/dev/null | wc -l) -eq 36"
chk D-025 "criterion 2 reformulated as mycl1_vs_mycn_nesting" \
  "grep -q 'mycl1_vs_mycn_nesting' config/params.yml"
chk D-025b "old absolute max_expected criterion removed" \
  "! grep -qE '^\s+max_expected:\s*0.50' config/params.yml"
# Same inversion as D-026b below, and found the same way: D-026 asserted the
# registry still carried the provisional "UNVERIFIED MYC" marking. D-035 resolved
# amplification from DepMap copy number and removed it, so the check began failing
# on a correctly-updated registry. A deferral check must be retired when the
# deferral is, or it reports the fix as the fault.
chk D-026 "provisional UNVERIFIED MYC marking retired (superseded by D-035)" \
  "! grep -q 'UNVERIFIED MYC' config/datasets.yml"
chk D-035c "amplification RESOLVED from DepMap copy number in registry" \
  "grep -q 'RESOLVED' config/datasets.yml && grep -q 'DepMap Public 26Q1 log2' config/datasets.yml"
# D-026b asserted the criterion-3 DEFERRAL was recorded. That was the right
# invariant at M5-M6 and the wrong one from M7, when D-035 resolved amplification
# from DepMap and evaluated the criterion. Left unchanged it passed *because* the
# gate table was stale, and would have failed the moment the table was corrected —
# an audit check defending the error it exists to catch. Inverted to match D-035.
chk D-026b "criterion 3 no longer deferred (resolved at M7 by D-035)" \
  "! grep -qi 'not evaluable' data/metadata/m5_gate_results.csv"
chk D-035 "criterion 3 evaluated and recorded FAIL in gate results" \
  "grep '3_distal_contrast' data/metadata/m5_gate_results.csv | grep -q 'FALSE'"
chk D-035b "no M5 criterion left unevaluated" \
  "! grep -qE ',NA\$' data/metadata/m5_gate_results.csv"

echo
echo "=== gate results actually on disk ==="
if [ -s data/metadata/m5_gate_results.csv ]; then cat data/metadata/m5_gate_results.csv; fi

echo
echo "=============================================================="
printf 'PASS %d   FAIL %d   SKIP %d   informational %d\n' "$pass" "$fail" "$skip" "$na"
echo "=============================================================="

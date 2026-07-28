# MYC-paralog regulatory hubs in small cell lung cancer

Paralog-resolved MYC, MYCN and MYCL1 regulatory programmes are built from SCLC
cell-line chromatin, scored in patient tumours, and combined across independent
evidence layers by a weight-free rank-aggregation framework.

Plotnik et al. 2024 (PMID 38747975) established that the three paralogs occupy
distinct enhancer repertoires in SCLC cell lines. That is prior work. This
repository asks whether those programmes retain paralog identity outside the
cell lines they were derived from, and whether their downstream genes converge
across evidence layers.

## Results

**Paralog identity does not survive in patient tumours.** Neuroendocrine lineage
state accounts for the regulon scores; the paralog's own expression does not.
Across two independent cohorts (GSE60052 n=79, George 2015 n=81), lineage
explains 0.244–0.428 of the unique variance in regulon score against 0.000–0.068
for the paralog. No paralog association survives lineage adjustment in the
replication cohort (0 of 3). The same pattern recurs in GeoMx spatial tissue,
where lineage explains 0.474 and 0.508 against 0.019 and 0.010 for MYC.

**The integration framework returns no ranked hub list.** Two of four specified
evidence domains were admitted; the transcriptional domain was dropped as
lineage-confounded and the network domain on both of its pre-conditions. The two
surviving domains overlap at the top of their rankings by up to 2.0–2.2× chance
(global p = 0.019 MYC, 0.008 MYCN, 0.065 MYCL1), but across 10,387 genes no gene
reaches FDR < 0.05 for any paralog; the lowest FDR obtained is 0.358. The
convergence is detectable in aggregate and not attributable to individual genes.

**MYC and MYCN regulons differ in programme.** The MYC regulon is
neurogenesis-weighted (OR 1.85, p = 7.0e-6), reproducing the enhancer-to-neurogenesis
link of the source study through an independent pipeline. The MYCN regulon is not
(OR 1.09, p = 0.33) and is Hallmark/housekeeping-weighted instead (OR 2.99,
p = 3.6e-4). The source study grouped MYC and MYCN together.

**Exploratory.** MYCL expression tracks BET-inhibitor AUC in the Sanger screens
(GDSC2 birabresib rho +0.621, lineage-adjusted +0.631; GDSC2 molibresib +0.617,
adjusted +0.587; GDSC1 PFI-1 +0.426, adjusted +0.429). Positive rho means
MYCL-high lines are *less* sensitive. It does not appear in PRISM (rho +0.046,
p = 0.84, n = 21), and GDSC1 and GDSC2 share a platform and cell lines, so they
are not two independent replications. Provisional.

Both negative outcomes were registered as reportable findings in the project
contract before any result existed.

## Aims and outcomes

| Aim | Scope | Outcome |
|---|---|---|
| 1 · Regulatory layer | Consensus ATAC region universe, bedGraph signal, active regions, super-enhancers, peak-to-gene linking, paralog regulons | 3 of 4 pre-registered gate criteria pass |
| 2 · Patient translation | Regulon scoring in two bulk cohorts; lineage-TF confounding as a primary analysis | Lineage dominates; paralog identity lost |
| 3 · Functional evidence | DepMap selective dependency; BET-inhibitor drug response; GENIE3 network importance | Dependency enriched pooled, not per paralog; drug association provisional; network excluded |
| 4 · Evidence integration | Two-stage Robust Rank Aggregation across evidence domains | Runs on 2 of 4 domains; no ranked hub list |
| 5 · Spatial coherence | GeoMx within-tumour coherence, no prognostic claims | Only the MYC regulon clears panel coverage |

The M5 gate criterion that fails is the distal-fraction contrast: the direction
reproduces (MYC-amplified 47.5% vs MYC-expressing 44.2%) but the magnitude is
+3.3 points against a published +27. The shared ATAC-defined universe is ~84%
distal by construction, which compresses the achievable range, where the source
study called peaks de novo per line. The failure is recorded in the gate table
rather than reformulated.

## Method notes

**Regulons are anchored on measured binding, not motifs.** All three paralogs
bind the same canonical E-box, so motif-based inference cannot separate them and
would collapse the independent variable into one MYC regulon. SCENIC and its
variants are excluded for that reason; GENIE3 was retained as an expression-only
secondary view and then excluded on its own pre-conditions.

**Rank aggregation, not a weighted score.** ATAC, H3K27ac, ChIP and
super-enhancer calls are correlated by construction, so a weighted sum would
count chromatin evidence several times over. Rank aggregation needs no weights
and handles missing layers natively.

**Peak-to-gene linking is not nearest-gene.** 48,756 links are retained from
127,019 candidate pairs; 32.6% agree with the nearest gene, and the link set is
enriched 6.06× over a distance-matched null.

## Interpretation and limits

MOES (Multi-Omics Evidence Score) is a heuristic prioritisation, not a predictive
or clinically validated model. It orders candidates by weight of independent
evidence. It does not establish causation and makes no claim about clinical
utility. Treatment-outcome, survival and immune-exclusion analyses are out of
scope.

- **No single-cell layer.** The suitable patient atlas is access-restricted; the
  accessible alternative is CDX-derived. It was dropped rather than substituted,
  so patient translation rests on two bulk cohorts.
- **Chromatin evidence is cell-line derived**, two lines per paralog for MYCN and
  MYCL1. The chromatin layer is descriptive; inference is reserved for the
  tumour cohorts.
- **Region boundaries are re-derived.** The keystone deposit ships no peak files,
  so regions come from a shared ATAC-defined universe. This is a declared
  deviation and is the likely source of two discordances with the source study.
- **Lineage-TF resolution is uneven.** POU2F3 occupancy is available in three
  keystone lines, ASCL1 in one, NEUROD1 in none. No blanket claim of lineage
  control is made.
- **The spatial panel is targeted** (~1,800 genes). Only the MYC regulon clears
  the coverage gate, at 56 of 500 members (11.2%); MYCN reaches 8.4% and MYCL1
  9.6%. The spatial result is a 56-gene proxy.
- **Pharmacogenomic coverage is thin** and the drug association appears on one
  platform only.
- **Package installation from source is untested.** The lockfile resolves all 232
  packages, but every clean-clone run so far reused a warm renv cache. Compiling
  from scratch on a machine that has never built these packages remains unproven.

## Data

All data are public and none are redistributed here. Everything is fetched by
scripts in `scripts/01_data/` against a checksum manifest, with download dates
recorded in `data/metadata/`.

| Layer | Source | Notes |
|---|---|---|
| MYC/MYCN/MYCL1 + H3K27ac ChIP-seq, ATAC-seq | GSE230649 | Keystone, 10 SCLC lines, hg19, bedGraph only (8.6 GB) |
| ATAC-seq support | GSE269424, GSE256345, GSE281523 | hg38, lifted; GSE269424 restricted to EGFP control arms |
| Lineage-TF ChIP-seq | GSE281524 (ASCL1), GSE210113 (NEUROD1), GSE249362 (POU2F3) | Confounder controls |
| Bulk tumour RNA-seq | GSE60052 (79 tumours + 7 normal); George et al. 2015 via cBioPortal (81) | Discovery and replication |
| CRISPR dependency, expression, copy number | DepMap Public 26Q1 | 25 SCLC vs 1,183 other lines |
| Pharmacogenomics | PRISM Repurposing Secondary, Sanger GDSC1, GDSC2 | BET inhibitors, named not pattern-matched |
| Spatial | GSE261348, GSE261345 (GeoMx DSP) | Targeted ~1,800-gene panel |

Genome build is hg19 throughout, fixed by the keystone series. The build is
asserted at every coordinate-handling step; lifts are explicit, logged, and
reported with loss rates.

## Reproducing

R ≥ 4.4 on Linux or macOS. No Python, no cloud services, no paid tools.

System libraries, required before `renv::restore()` on a machine that has never
built these packages:

```bash
sudo apt-get update && sudo apt-get install -y build-essential cmake gsfonts libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libfreetype-dev libharfbuzz-dev libfribidi-dev libpng-dev libtiff-dev libjpeg-dev libcairo2-dev libmagick++-dev pandoc
```

Clone outside any existing checkout, restore the environment, and reinstall the
large-file hook (git hooks are not cloned):

```bash
git clone git@github.com:ImmunoScholar/sclc-myc-regulatory-hubs.git
cd sclc-myc-regulatory-hubs
Rscript -e 'renv::restore()'
bash scripts/00_setup/01_init_git.sh
```

Run the whole pipeline in dependency order, stopping at the first failing stage:

```bash
bash scripts/run_all.sh
```

`--with-data` includes the ~6 h download; `--figures` rebuilds figures and the
report only. Individual stages can be run in numeric order instead; each writes
outputs to `results/` and QC to `logs/`.

Rebuild and check the figure set, then render the report:

```bash
bash scripts/07_figures/99_build_all.sh
bash scripts/08_report/01_render_report.sh
```

Every value in the report is read at render time from files the analysis wrote,
and the render stops if those files are absent. Rendering prefers the Quarto CLI
and falls back to knitr with pandoc, which differ in presentation only.

### Checks

```bash
bash scripts/00_setup/audit_decisions.sh        # decisions still match the code
Rscript scripts/08_report/02_lint_outputs.R     # contract vocabulary and limits
bash scripts/00_setup/07_clean_clone_verify.sh  # fresh clone rebuilds every figure
```

The clean-clone check clones the repository, rebuilds all eight figures from the
clone, and compares them byte-for-byte against the committed set. It does not
cover package installation from source or data re-download.

### Requirements

About 40 GB of free disk. Memory is the binding constraint: bedGraph signal is
processed chromosome by chromosome and the pipeline runs in roughly 10 GB of RAM.
Loading signal genome-wide will not work on a normal machine.

## Repository layout

```
config/           analysis parameters (YAML); thresholds live here, not in code
data/raw/         downloaded data                      (git-ignored)
data/processed/   intermediates                        (git-ignored)
data/metadata/    manifests, checksums, figure inputs   (tracked)
docs/             contract, gap statement, decision log, risk log, journal
scripts/00_setup/       environment, git, audit and verification gates
scripts/01_data/        manifest, download, verify, QC
scripts/02_regulatory/  Aim 1 — regions, signal, super-enhancers, regulons
scripts/03_tumour/      Aim 2 — regulon scoring, lineage confounding
scripts/04_functional/  Aim 3 — dependency, amplification, drug response
scripts/05_integration/ Aim 4 — evidence build, MOES, controls, benchmark
scripts/06_spatial/     Aim 5 — spatial coherence
scripts/07_figures/     figure generation and palette selection
scripts/08_report/      report render, output lint, publication tables
results/tables/   stage tables and six publication tables  (tracked)
results/objects/  serialised intermediates              (git-ignored)
figures/          eight figures, captions, manifest      (tracked)
logs/             run and QC logs                       (git-ignored)
tests/            assertions on intermediate outputs
```

Do not clone this repository inside its own working tree. A nested clone is a
separate git repository, and `git add -A` records it as a gitlink — a pointer to
a commit nobody else can fetch. `.gitignore` now excludes the path.

`docs/project_contract.md` holds the scope and the pre-registered controls;
`docs/decision_log.md` records every methodological decision and deviation.

## Outputs

Figures 1–4 carry the main results: lineage dominance over paralog identity, the
M5 gate, functional evidence with domain accounting, and the MOES null with its
power analysis. Figures S1–S4 cover the data landscape, region construction, the
lineage confound at three levels, and spatial coherence. Each has a caption file
beside it.

Six publication tables in `results/tables/publication/` correspond to the M5
gate, lineage dominance, dependency, MOES, spatial coherence and drug response.

## Prior work

- Plotnik JP et al. (2024) *Mol Cancer Res* — paralog-resolved MYC/MYCN/MYCL1
  binding in SCLC; the foundation for Aim 1. PMID 38747975
- Dammert MA et al. (2019) *Nat Commun* — MYC paralog-specific apoptotic priming.
- Ireland AS et al. (2020) *Cancer Cell* — MYC-driven subtype plasticity.
- Jung LA et al. (2017) *Cancer Res* — 18-gene MYC-activity signature, used as a
  paralog-blind benchmark.
- Zhang et al. (2025) *Clin Cancer Res* — MYC-paralog amplification and the
  spatial immune microenvironment; the reason the spatial aim is restricted.

## License

Code is released under the MIT License. The underlying datasets remain subject to
the terms of their original repositories and should be cited directly.

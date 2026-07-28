# MYC-paralog regulatory hubs in small cell lung cancer

A multi-omics analysis investigating whether MYC, MYCN, and MYCL1 maintain distinct regulatory programs in patient SCLC tumours, despite occupying distinct enhancer repertoires in cell lines.

## Main findings

The analysis builds paralog-resolved regulons from cell-line chromatin data and tests their activity across three independent patient cohorts (tumour RNA-seq, spatial transcriptomics) plus functional genomics.

**Primary result:** Paralog-specific regulatory programs lose paralog identity in patient tumours. Neuroendocrine lineage state—not paralog expression—explains regulon activity. Lineage explains 24–43% of regulon-score variance; paralog expression explains 0–7%. This pattern holds across two bulk tumour cohorts (n=79, n=81) and is replicated in spatial tissue, with consistent effect sizes.

**Evidence integration:** A two-stage rank aggregation framework identifies weak aggregate convergence between cis-regulatory and functional evidence (~2× enrichment at the top of both rankings), but no individual gene survives multiple-testing correction across 10,387 genes (best FDR 0.358). The aggregated signal is detectable; the per-gene signal is not.

**Exploratory finding:** MYCL expression tracks BET-inhibitor resistance in two independent pharmacogenomic screens (Sanger GDSC1, GDSC2; rho +0.59–+0.62 after lineage adjustment), but does not replicate in PRISM. Effect is provisional and platform-dependent.

## Scientific context

Plotnik et al. (2024) mapped distinct enhancer repertoires for the three MYC paralogs across ten SCLC cell lines using ChIP-seq and ATAC-seq. That foundational work is cited but not repeated here.

This project extends beyond *in vitro* by:
- Constructing paralog-resolved regulatory programs (super-enhancer-aware, peak-to-gene linked)
- Scoring them in independent human SCLC tumour cohorts
- Testing whether paralog signals remain separable from lineage-transcription-factor programs
- Integrating independent evidence layers (chromatin, CRISPR dependency, drug response, spatial)

## Methodology

| Component | Details |
|---|---|
| **Chromatin layer** | GSE230649 (10 SCLC lines): MYC/MYCN/MYCL1 ChIP-seq, H3K27ac, ATAC-seq, hg19. Consensus accessible-region universe from three independent ATAC datasets. Signal quantified, active regions called, super-enhancers stitched, peak-to-gene links inferred. |
| **Patient translation** | GSE60052 (79 samples) and George et al. 2015 via cBioPortal (81 samples). Regulon scoring with explicit confounding test: paralog signals vs. lineage-TF programs (ASCL1, NEUROD1, POU2F3). |
| **Functional domains** | CRISPR dependency (DepMap): selective essentiality in SCLC vs. non-SCLC lines. Drug response: BET inhibitors across three pharmacogenomic screens (PRISM, GDSC1, GDSC2), neuroendocrine-adjusted. Network inference: excluded (59 lines underpowered for tree-ensemble GRN; MYC expression confounded by NE state). |
| **Evidence integration** | Two-stage Robust Rank Aggregation: within-domain ranking → between-domain aggregation → permutation-based FDR control. Designed without arbitrary weights to handle correlated layers (chromatin/ATAC/ChIP/super-enhancer). |
| **Spatial validation** | GeoMx DSP (GSE261348, GSE261345): targeted ~1,800-gene panel. Coverage gate: minimum 20 measured genes AND ≥10% of regulon members. Only MYC regulon qualifies (56/500 members, 11.2%). |

## Key analyses and their scope

### Regulatory layer (Aim 1)
Three paralog-resolved regulons built from measured ChIP binding and peak-to-gene linking (not motif-based inference). Validation:
- MYCN occupancy shows expected ~84% nesting within MYC regions (Plotnik's reported 84% overlap confirmed)
- MYCL1 shows independent occupancy pattern
- MYC regulon enriched for neurogenesis (p=7e-6, independent replication)
- MYCN regulon enriched for Hallmark housekeeping (p=3.6e-4; published sources grouped MYCN with MYC)
- MYCL1 regulon does not reach significance

### Patient translation (Aim 2)
Regulon scores do not track paralog expression after adjustment for neuroendocrine lineage state:
- Lineage explains 24–43% of regulon-score variance
- Paralog expression explains 0–7% (overlapping confidence intervals)
- Effect replicates identically in both cohorts
- Result holds in spatial tissue (third independent modality)

### Functional evidence (Aim 3)
Selective CRISPR dependency present at aggregate level (pooled regulons, OR 2.71, p=0.0048) but not separable per paralog. Drug response: MYCL-BET-inhibitor association in GDSC1/GDSC2 (not PRISM). MYC trends opposite direction but does not survive confounding adjustment.

### Evidence integration (Aim 4)
Four evidence domains specified: cis-regulatory (admitted), functional (admitted), transcriptional (dropped; lineage-confounded at aggregate and per-gene level), network (dropped; NE-confounded and underpowered). Two-domain MOES degenerate (no hierarchy). Weak inter-domain concordance. No hub list reported.

### Spatial coherence (Aim 5)
MYC regulon shows within- and between-tumour variance; effect sizes match bulk cohorts. Limited to MYC due to targeted panel coverage. No prognostic, immune-exclusion, or outcome claims.

## Interpretation

The evidence-integration framework is a **heuristic prioritisation tool**, not a predictive or clinically validated model. No individual gene identified through this analysis should be treated as a validated therapeutic target without experimental corroboration.

The primary finding—paralog-program collapse under lineage dominance—is a robust null result. It arose as a pre-registered reportable outcome and is reported with power analysis rather than as a study failure. The phenomenon appears at three biological levels (chromatin occupancy, gene expression, spatial transcripts) and across two independent tumour cohorts.

### Known limitations

- **No single-cell layer.** The suitable patient-derived atlas (HTAN/Chan 2021) is access-restricted; the public alternative is CDX-derived and cannot support a patient-translation claim. Dropped rather than substituted.
- **Chromatin from cell lines, n=2 per paralog.** Descriptive layer only; no per-line statistical claims. Inference reserved for tumour cohorts where sample size is adequate.
- **Consensus region universe.** Boundaries are re-derived from ATAC-defined consensus, not Plotnik's de-novo peak calls. Introduces structure bias (distal-fraction criterion fails) and is reported as such.
- **Lineage-TF coverage uneven.** POU2F3 occupancy available in 3 keystone lines; ASCL1 in 1; NEUROD1 in 0. Confounding-resolution heterogeneous.
- **Spatial panel targeted.** ~1,800 genes measurable; most regulon members absent. MYC regulon barely clears the coverage gate (by 1.2 percentage points).
- **Pharmacogenomics sparse and platform-dependent.** MYCL-BET-inhibitor association appears in Sanger data only (GDSC1, GDSC2; same platform, overlapping lines) and not in PRISM. Pre-registration did not cover this specific comparison.

## Data sources

All data are public. No datasets are redistributed; all are fetched by scripts with checksum verification and download dates recorded.

| Source | Accession | Role | Notes |
|---|---|---|---|
| GSE230649 | Plotnik et al. ChIP-seq/ATAC | Keystone: MYC/MYCN/MYCL1, H3K27ac, ATAC (10 lines) | hg19 bedGraph only; 8.6 GB |
| GSE269424, GSE256345, GSE281523 | ATAC (support) | Consensus accessible-region universe | hg38; liftOver to hg19 |
| GSE281524, GSE210113, GSE249362 | Lineage-TF ChIP | ASCL1, NEUROD1, POU2F3 controls | Confounder identification |
| GSE60052 | Bulk SCLC RNA-seq | Primary tumour cohort | n=79 + 7 normal |
| George et al. 2015 | Bulk SCLC RNA-seq | Replication cohort | n=81; via cBioPortal |
| DepMap Public 26Q1 | CRISPR dependency | Selective essentiality | Expression, copy number, dependency |
| GDSC1, GDSC2, PRISM | Pharmacogenomics | Drug-response association | Sanger and Broad platforms; BET inhibitors |
| GSE261348, GSE261345 | GeoMx DSP | Spatial transcriptomics | Targeted ~1,800-gene panel |

**Genome build:** hg19 throughout (mandated by GSE230649). Build asserted at every coordinate-handling step; any liftOver logged and loss rate QC'd.

## Reproducing the analysis

**Requirements:** R ≥ 4.4 on Linux/macOS. All dependencies pinned in `renv.lock`. No Python, no cloud services, no paid tools.

**System libraries:**
```bash
sudo apt-get update && sudo apt-get install -y \
  build-essential cmake gsfonts \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libfreetype-dev libharfbuzz-dev libfribidi-dev \
  libpng-dev libtiff-dev libjpeg-dev libcairo2-dev libmagick++-dev \
  pandoc
```

**Setup:**
```bash
git clone git@github.com:ImmunoScholar/sclc-myc-regulatory-hubs.git
cd sclc-myc-regulatory-hubs
Rscript -e 'renv::restore()'
bash scripts/00_setup/01_init_git.sh
```

**Run pipeline:**
```bash
bash scripts/run_all.sh [--with-data] [--figures]
```

The pipeline processes data → regulons → patient translation → functional evidence → integration in dependency order, writing outputs to `results/` and logs to `logs/`. Stops at first error rather than producing incomplete outputs.

**Rebuild figures:**
```bash
bash scripts/07_figures/99_build_all.sh
```

**Render report:**
```bash
bash scripts/08_report/01_render_report.sh
```

(Prefers Quarto CLI; falls back to knitr+pandoc if unavailable. Quarto changes presentation only, not content.)

**System resources:** ~40 GB disk free. RAM is the binding constraint (~10 GB); pipeline processes bedGraph signal chromosome-by-chromosome rather than genome-wide to fit.

## Repository structure

```
scripts/           Numbered analysis stages (00–08)
  01_data/         Download, verify, QC all source datasets
  02_m5/           Regulatory layer (Aim 1)
  03_m6/           Patient translation (Aim 2)
  04_m7/           Functional evidence (Aim 3)
  05_m8/           Evidence integration (Aim 4)
  06_m9/           Spatial validation (Aim 5)
  07_figures/      Generate publication figures
  08_report/       Render Quarto report

config/            YAML analysis parameters (thresholds, methods)
data/
  raw/             Downloaded public datasets (git-ignored)
  processed/       Intermediates (git-ignored)
  metadata/        Manifests, checksums, provenance (tracked)

results/
  tables/          Output tables, metric summaries (tracked)
  objects/         Serialized R objects (git-ignored)

figures/           Publication figures and captions (tracked)
logs/              Execution logs and QC reports (git-ignored)
tests/             Unit tests on intermediate outputs
docs/
  project_contract.md          Frozen scope and pre-registered controls
  02_gap_statement.md          Novelty assessment vs. Plotnik et al.
  03_dataset_inventory.md      Every dataset, role, access status
  04_analysis_architecture.md  Pipeline logic and method choices
  decision_log.md              Detailed record of every deviation

R/                 Shared utility functions
renv/              Package version management
```

## Figures

- **Figure 1:** Lineage dominance over paralog identity (variance partition, expression correlation, MYC/NE antagonism)
- **Figure 2:** M5 regulatory-layer gate (four pre-registered criteria: MYCN nesting, differential occupancy, distal-fraction, E-box motif specificity)
- **Figure 3:** Functional evidence and evidence-domain accounting (CRISPR dependency, regulon enrichment, MOES domain status)
- **Figure 4:** Convergence weak and unlocalisable to genes (top-K overlap, permutation FDR, power analysis, paralog-specificity bias)
- **Figure S1:** Data landscape (composition, build mix, cell-line overlap with controls)
- **Figure S2:** Regulatory pipeline (ATAC threshold validation, stitching robustness, peak-to-gene enrichment)
- **Figure S3:** Lineage confound at every level (chromatin occupancy, expression, genome-wide association)
- **Figure S4:** Spatial coherence (coverage gate, variance partitioning, regulon vs. lineage effect)

## License

Code released under MIT License (see [`LICENSE`](LICENSE)). Underlying datasets subject to their original repositories' terms; cite directly.

## References

**Foundational work:**
- Plotnik JP et al. (2024) Mol Cancer Res. Paralog-resolved MYC/MYCN/MYCL1 binding in SCLC. PMID 38747975

**Positive/negative controls and benchmarks:**
- Dammert MA et al. (2019) Nat Commun. MYC paralog-specific apoptotic priming.
- Tlemsani C et al. (2020) Cell Rep. SCLC-CellMiner pharmacogenomics.
- Ireland AS et al. (2020) Cancer Cell. MYC-driven lineage plasticity.
- Zhang et al. (2025) Clin Cancer Res. MYC-paralog amplification and spatial immune microenvironment.
- Jung LA et al. (2017) Cancer Res 77(4). MYC-activity signature for benchmark.

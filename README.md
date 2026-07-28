# MYC-paralog regulatory hubs in small cell lung cancer

A reproducible, multi-omics analysis that translates paralog-resolved MYC
regulatory programs from SCLC cell lines into patient tumours, and prioritises
candidate regulatory hubs using a transparent, weight-free evidence-integration
framework.

> **Status: analysis complete, release in preparation.** All five aims have been
> run to a conclusion. Two of those conclusions are **negative results**, and they
> are the main findings rather than a shortfall — see
> [`report.html`](report.html) for the full write-up and
> [`docs/decision_log.md`](docs/decision_log.md) for why each decision was made.

## What this project found

1. **Paralog-resolved regulons are internally coherent in cell-line chromatin**
   (3 of 4 pre-registered gate criteria pass) but **lose paralog identity in
   patient tumours.** Neuroendocrine lineage state explains the regulon scores;
   the paralog's own expression does not. Replicated across two cohorts
   (160 tumours) with a power analysis, and reproduced again in spatial tissue.
2. **The evidence-integration framework returns no prioritised hub list.** Three
   of four planned evidence domains were excluded on their own evidence. The two
   that remain converge weakly in aggregate (~2× top-K overlap) but no individual
   gene survives testing across 10,387 — aggregate detectable, per-gene
   unattributable.
3. Along the way: the MYC regulon independently reproduces the published
   MYC-enhancer→neurogenesis link, while the **MYCN regulon does not** — a
   discordance with the source study, which grouped them together.

Both negative results were **pre-registered as reportable findings** before any
result existed. Nothing here was adjusted until it passed.

---

## The question

MYC, MYCN and MYCL1 are amplified in largely mutually exclusive subsets of small
cell lung cancer, and they are not interchangeable — they occupy **distinct
enhancer repertoires**. Plotnik et al. (2024) established this directly, with
ChIP-seq of all three paralogs plus H3K27ac and ATAC-seq across ten SCLC cell
lines.

That work is the foundation of this project and is **not repeated here.**

What it did not do was leave the cell line. Every result was *in vitro*; there
was no super-enhancer calling, no peak-to-gene linking beyond nearest-gene, no
network inference, no patient tumour data, and no prioritisation of the
downstream targets it identified.

**This project asks the next question:**

> Do paralog-specific MYC regulatory programs operate in patient tumours — and
> which of their downstream hubs are robust enough, across independent evidence
> layers, to warrant functional investigation?

## Approach

Five aims, one coherent argument rather than five separate analyses:

| Aim | What it does |
|---|---|
| **1 · Regulatory layer** | Build a consensus accessible-region universe from ATAC-seq, quantify MYC/MYCN/MYCL1 and H3K27ac bedGraph signal over it, call active regions and super-enhancers, and link regions to genes — producing **paralog-resolved regulons** |
| **2 · Patient translation** | Score those regulons in two independent human SCLC tumour cohorts, and test explicitly whether "paralog-specific" signal is separable from lineage-TF (ASCL1 / NEUROD1 / POU2F3) programs |
| **3 · Functional evidence** | Selective CRISPR dependency (DepMap). *GENIE3 network importance was **excluded** on both its pre-conditions — 59 lines is underpowered and MYC expression tracks lineage across them. Drug-response association was **not carried out**.* |
| **4 · Evidence integration** | Two-stage Robust Rank Aggregation across evidence domains. *Ran on **two** of four domains; returns **no ranked hub list** — see findings above.* |
| **5 · Spatial coherence** | A restricted check of within-tumour coherence. No prognostic claims. *Only the MYC regulon clears the panel-coverage gate.* |

Aims are shown as specified, with their actual outcome in italics. An aim that
quietly changes shape between specification and write-up is the easiest kind of
overstatement to commit and the hardest for a reader to detect.

### Two design choices worth explaining

**Why not SCENIC?** All three MYC paralogs bind the same canonical E-box.
Motif-based network inference is therefore structurally unable to tell them
apart, and would collapse the project's independent variable into a single "MYC
regulon". Regulons here are anchored on **measured paralog binding** instead,
with GENIE3 supplying an independent, expression-only view.

**Why rank aggregation instead of a weighted score?** Any weights would be
arbitrary. Worse, the evidence layers are not independent — ATAC, H3K27ac, ChIP
and super-enhancer calls are correlated by construction, so a weighted sum would
silently count chromatin evidence four times. Rank aggregation is weight-free and
handles missing layers natively.

## Data

All data are public. Nothing in this repository is redistributed — everything is
fetched by scripts in `scripts/01_data/`, with a checksum manifest and recorded
download dates in `data/metadata/`.

| Layer | Source | Notes |
|---|---|---|
| MYC/MYCN/MYCL1 + H3K27ac ChIP-seq, ATAC-seq | GEO **GSE230649** | 10 SCLC lines. **hg19**, bedGraph only (8.6 GB) |
| ATAC-seq (support) | GSE269424, GSE256345, GSE281523 | |
| Lineage-TF ChIP-seq (confounder controls) | GSE281524 (ASCL1), GSE210113 (NEUROD1), GSE249362 (POU2F3) | |
| Bulk tumour RNA-seq | GSE60052 (79 tumours + 7 normal); George et al. 2015 via cBioPortal | |
| CRISPR dependency & expression | DepMap (release pinned in the manifest) | |
| ~~Pharmacogenomics~~ | ~~SCLC-CellMiner CDB~~ | **Not used.** The drug-response association under Aim 3 was never carried out; listed here because it was specified |
| Spatial | GSE261348, GSE261345 (GeoMx DSP) | Targeted ~1,800-gene panel |

**Genome build is hg19 throughout**, mandated by the source alignment of
GSE230649. Build is asserted at every coordinate-handling step; any lift is
explicit, logged, and QC'd for loss rate.

## Reproducing

Requires R ≥ 4.4 on Linux/macOS. No paid tools, no cloud services, no Python.

```bash
git clone git@github.com:ImmunoScholar/sclc-myc-regulatory-hubs.git
cd sclc-myc-regulatory-hubs
```

Restore the exact package versions:

```bash
Rscript -e 'renv::restore()'
```

Reinstall the large-file guard (git hooks are not cloned):

```bash
bash scripts/00_setup/01_init_git.sh
```

Then run the numbered scripts in `scripts/` in order. Each stage writes its
outputs to `results/` and its QC to `logs/`.

Rebuild every figure and verify the set:

```bash
bash scripts/07_figures/99_build_all.sh
```

Render the report:

```bash
bash scripts/08_report/01_render_report.sh
```

The report source is [`report.qmd`](report.qmd) and the rendered output is
[`report.html`](report.html). Rendering prefers the Quarto CLI; if Quarto is not
installed it falls back to `knitr` + `pandoc`, which are already pinned in
`renv.lock`. The two paths differ in presentation only — the fallback loses
Quarto's folded-code UI and theme, not any content or number. To use the Quarto
path, install the CLI from <https://quarto.org/docs/get-started/>.

Every value in the report is read at render time from files the analysis wrote,
so a stale analysis cannot produce a fresh-looking report, and the render refuses
to run if those files are missing.

Or run the whole pipeline in dependency order:

```bash
bash scripts/run_all.sh
```

It stops at the first failing stage — a pipeline that continues past a failure
produces outputs that look complete and are not. Pass `--with-data` to include
the ~6 h download, or `--figures` for figures and report only.

### System requirements

Roughly 40 GB of free disk. **Memory is the real constraint:** the pipeline is
written to process bedGraph signal chromosome by chromosome and runs in about
10 GB of RAM. Loading the signal genome-wide will not work on a normal machine.

## Repository layout

```
config/          analysis parameters (YAML) — thresholds live here, not in code
data/raw/        downloaded data          (git-ignored)
data/processed/  intermediates            (git-ignored)
data/metadata/   manifests, checksums, provenance   (tracked)
docs/            project contract, gap statement, decision log, journal
scripts/         numbered analysis stages
results/tables/  output tables            (tracked)
results/objects/ serialised intermediates (git-ignored)
figures/         publication figures      (tracked)
logs/            run and QC logs          (git-ignored)
tests/           assertions on intermediate outputs
```

Start with [`docs/project_contract.md`](docs/project_contract.md) for scope and
pre-registered controls, and [`docs/decision_log.md`](docs/decision_log.md) for
why the project is built the way it is.

## Interpretation and limits

**MOES is a heuristic prioritisation, not a predictive or clinically validated
model.** It orders candidates by weight of independent
evidence so that finite experimental effort can be aimed sensibly. It does not
establish causation, and it is not a claim about clinical utility.

Known limitations, stated up front rather than buried:

- **No single-cell layer.** The only suitable patient atlas is access-restricted,
  and the accessible alternative is cell-line/xenograft-derived. It was dropped
  rather than substituted, so patient translation rests on two bulk cohorts
  (n=79 and n≈81). Effective sample size is modest and is reported as such.
- **Chromatin evidence comes from cell lines**, n=2 per paralog. No per-line
  statistical claims are made from that; the chromatin layer is descriptive and
  inference is reserved for the tumour cohorts.
- **Region boundaries are re-derived**, not the original peak calls, because the
  source deposit contains no peak files. This is a declared methodological
  deviation, checked against the published region counts before proceeding.
- **Paralog amplification correlates with SCLC subtype**, so paralog-specific
  signal may be partly lineage-driven. This is tested directly as a primary
  analysis rather than assumed away — and if the programs prove inseparable from
  lineage identity, that negative result is the reported finding.
- **The spatial panel is targeted (~1,800 genes)**, so most regulon members are
  not measurable spatially. Coverage fraction is reported on every spatial figure.

## Prior work this builds on

- Plotnik JP et al. (2024) *Mol Cancer Res* — paralog-resolved MYC/MYCN/MYCL1
  binding in SCLC. The foundation for Aim 1. [PMID 38747975](https://pubmed.ncbi.nlm.nih.gov/38747975/)
- Dammert MA et al. (2019) *Nat Commun* — MYC paralog-specific apoptotic priming.
  Used as an orthogonal positive control.
- Tlemsani C et al. (2020) *Cell Rep* — SCLC-CellMiner.
- Ireland AS et al. (2020) *Cancer Cell* — MYC-driven subtype plasticity.
- Zhang et al. (2025) *Clin Cancer Res* — MYC-paralog amplification and the
  spatial immune microenvironment. The reason the spatial aim here is restricted.

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). The underlying
datasets remain subject to the terms of their original repositories and should be
cited directly.

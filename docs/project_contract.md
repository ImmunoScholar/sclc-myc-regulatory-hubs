# Project Contract

**Status:** ACTIVE · **Frozen:** 2026-07-26 · **Applies from:** M3 onward

This document is the binding scope agreement for the project. M0 (scoping) and
M1 (specification) are **frozen**. The biological question, research gap,
dataset selection, GRN strategy, and evidence-integration framework are settled
and are not to be redesigned. Any deviation requires a numbered entry in
[`decision_log.md`](decision_log.md) stating the new evidence and the rationale.

---

## 1. The question

> How do MYC-associated transcriptional and epigenetic regulatory programs shape
> SCLC state and regulatory vulnerability — and which downstream regulatory hubs
> are robust enough, across independent evidence layers, to warrant functional
> investigation?

## 2. The gap (frozen — see `02_gap_statement.md` for the full argument)

Plotnik et al. 2024 (Mol Cancer Res, PMID 38747975) established that MYC, MYCN
and MYCL1 occupy **distinct enhancer repertoires** in SCLC cell lines. That
finding is **prior work and is never re-claimed here.**

What remains unaddressed, and what this project does:

1. Convert paralog-resolved MYC binding into **enhancer-anchored,
   super-enhancer-aware regulons** (Plotnik did a 1 kb promoter/distal split only).
2. Test those regulons in **independent human SCLC tumour cohorts** (Plotnik was
   entirely *in vitro*).
3. Prioritise hubs with a **transparent, weight-free evidence-integration
   framework** (no prioritisation of any kind exists in the prior work).

**Novelty grade: MODERATE-HIGH.** The regulatory-characterisation step is prior
work; the integration, translation and prioritisation are not.

## 3. Aims

| Aim | Content | Milestone |
|---|---|---|
| **Aim 1** | Paralog-resolved regulatory layer: consensus regions → bedGraph signal → active regions → super-enhancers → peak-to-gene linking → regulons | M5 |
| **Aim 2** | Patient-tumour translation: regulon scoring in GSE60052, replication in George 2015, **lineage-TF confounding analysis** | M6 |
| **Aim 3** | Functional & network evidence: DepMap selective dependency, drug-response association, GENIE3 network importance | M7 |
| **Aim 4** | MOES integration: two-stage Robust Rank Aggregation across four evidence domains | M8 |
| **Aim 5** | Spatial coherence, **restricted** — no prognostic or immune-exclusion claims | M9 |

## 4. Frozen methodological decisions

These were settled in M1 with reasoning recorded. Do not re-litigate.

| # | Decision | Core reason |
|---|---|---|
| 1 | **SCENIC / pySCENIC / SCENIC+ rejected** as the GRN method | All three MYC paralogs bind the same canonical E-box. Motif-based inference is structurally paralog-blind and would collapse the project's independent variable into one "MYC regulon". SCENIC+ additionally impossible — no public SCLC scATAC-seq exists. |
| 2 | **Primary GRN = ChIP-anchored regulon construction**; secondary = GENIE3 | Uses measured paralog binding, not inferred motifs. Only approach that preserves paralog identity. |
| 3 | **MOES = two-stage Robust Rank Aggregation, not a weighted sum** | Evidence layers are non-independent — ATAC, H3K27ac, ChIP and SE evidence are correlated by construction. A weighted sum would quadruple-count chromatin. Four domains → RRA within, then across. Precedent: NetICS (Dimitrakopoulos 2018). |
| 4 | **No de novo peak calling** | GSE230649 ships bedGraph only (8.6 GB, no peak files); re-aligning 28 SRA samples is not feasible on this hardware. Signal is quantified over an ATAC-defined consensus region universe in R. Removes any need for MACS2 or bedtools. |
| 5 | **Single-cell layer excluded** | The only suitable patient atlas (Chan 2021) is HTAN/Synapse-gated and unobtainable; the GEO alternative is CDX-derived and cannot support a patient-translation claim. Dropping beats substituting. |
| 6 | **Spatial restricted to coherence/heterogeneity** | Zhang et al. 2025 already published the MYC-paralog ↔ spatial immune-exclusion link. |
| 7 | **Methylation, proteomics, miRNA excluded** | Available but off-hypothesis. Inclusion for completeness would violate biological coherence. |

## 5. Pre-registered controls

Declared **before** any result exists, so they cannot be tuned after the fact.

- **Negative control (hard gate at M8):** Plotnik's paralog-*shared* housekeeping
  promoters (e.g. `RPS26`) must **not** rank as paralog-specific hubs. Failure
  means the framework lacks specificity and is reported as failed — not quietly
  adjusted until it passes.
- **Positive control:** top MYC-specific hubs should be consistent with the
  MYC→BCL2/MCL1 apoptotic-priming axis of Dammert et al. 2019.
- **Benchmark comparators:** MSigDB `HALLMARK_MYC_TARGETS_V1/V2` and the Jung
  et al. 2017 18-gene MYC-activity signature — both paralog-blind, which is the
  point of comparison.
- **Sanity gate at M5 (hard):** re-derived active-region counts must be within a
  defensible range of Plotnik's reported 18,823 MYC / 4,017 MYCN / 5,688 MYCL1,
  with MYCN/MYC overlap near their reported ~84%. **Substantial discordance
  halts the project for method review; it is not explained away.**

## 6. Falsifiability

The project has a defined failure mode and a defined response to it.

If paralog-specific regulons show no differential activity in the tumour
cohorts, or prove inseparable from lineage-TF programs (risk R-01), the honest
conclusion becomes *"MYC-paralog regulatory programs in SCLC are largely
inseparable from lineage identity"* — reported as a negative result with a power
analysis. The framework and the audit trail remain the contribution. **The
project does not collapse; the conclusion changes.**

## 7. Interpretive limits — binding on all outputs

MOES is a **heuristic prioritisation, not a predictive or clinically validated
model.** This statement appears in the README, in every output table header, and
in every relevant figure legend.

Prohibited vocabulary in all reports, figures and code comments:
`predicts` · `driver of` · `therapeutic target` · `causes` · `validated`
(unless describing someone else's experimentally validated result, with citation).

Provisional conclusions are labelled provisional.

## 8. Reproducibility contract

- Primary language **R**; Python only if genuinely required (currently: none).
- `renv.lock` pins every package version. `renv::restore()` from a clean clone
  must work.
- **No data files in version control.** Data are reproduced by download scripts
  plus a checksum manifest. Enforced by `.gitignore` *and* a pre-commit
  large-file hook (5 MB limit).
- Genome build is **asserted at every coordinate-handling step**. Project build
  is hg19 (mandated by GSE230649); any lift is explicit, logged, and QC'd for
  loss rate.
- Every external resource records its **download date and version/release ID**
  in `data/metadata/`. DepMap is pinned to a named release.
- Release gate (M11): a fresh `git clone` on a clean machine reproduces every
  figure. Actually tested, not assumed.

## 9. Out of scope

Treatment-outcome and survival analysis · immune-exclusion claims · any *in vivo*
or wet-lab interpretation beyond citing published work · clinical
recommendations of any kind · methylation, proteomic and miRNA layers.

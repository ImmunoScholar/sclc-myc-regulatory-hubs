# Project Brief

**Project:** SCLC MYC-Paralog Regulatory Hubs
**Owner:** Priya
**Started:** 2026-07-26
**Status:** Phase 2 (specification) — coding not yet authorised

---

## Working title

**Evidence-integrated prioritisation of MYC-paralog-associated regulatory hubs in small cell lung cancer: from cell-line enhancer maps to patient-tumour validation**

---

## Central hypothesis

MYC, MYCN, and MYCL1 drive **non-interchangeable** enhancer-anchored regulatory programs in SCLC. A subset of the genes downstream of these programs constitutes **regulatory hubs** whose importance is corroborated independently by chromatin accessibility, paralog binding, enhancer/super-enhancer activity, network centrality, transcriptional behaviour in patient tumours, cellular dependency, and drug-response association. These convergently supported hubs are stronger candidates for functional investigation than genes nominated by any single evidence layer — including the paralog-blind MYC signatures currently in general use.

**Null hypothesis (must be reportable):** paralog-specific regulons show no differential activity across patient tumours stratified by MYC-paralog status, and multi-layer convergence identifies no hubs beyond those recoverable from expression alone.

---

## Specific aims

**Aim 1 — Build paralog-resolved, enhancer-anchored regulons.**
Re-derive MYC/MYCN/MYCL1 binding and H3K27ac landscapes from GSE230649; call super-enhancers (absent from the source publication); link regulatory regions to target genes using distance-weighted and expression-correlated methods rather than nearest-gene; replicate accessibility support in two independent ATAC datasets.
*Deliverable:* three paralog regulons + one shared/core regulon, with per-gene evidence provenance.

**Aim 2 — Test whether paralog programs operate in patient tumours.**
Score regulon activity in two independent bulk SCLC tumour cohorts. Test differential activity by MYC-paralog status and by A/N/P/I subtype. Explicitly test whether paralog regulons are separable from lineage-TF programs, or merely a proxy for them.
*Deliverable:* cross-cohort regulon activity, with an explicit lineage-TF confounding analysis.

**Aim 3 — Add functional and pharmacogenomic evidence.**
Test whether regulon members show selective CRISPR dependency in MYC-paralog-amplified vs non-amplified SCLC lines (DepMap), and whether their expression associates with drug response (SCLC-CellMiner, PRISM).
*Deliverable:* dependency and drug-association statistics per candidate hub.

**Aim 4 — Integrate evidence into a transparent prioritisation framework (MOES).**
Two-stage rank aggregation across four conceptually independent evidence domains, with leave-one-domain-out stability analysis, permutation-based FDR, and built-in positive and negative controls.
*Deliverable:* ranked hub table with full per-layer audit trail, presented as a **heuristic prioritisation framework, not a predictive or clinically validated model.**

**Aim 5 — Orthogonal spatial coherence check (restricted).**
In two independent GeoMx cohorts, ask only whether prioritised hubs show coherent, regionally heterogeneous expression in patient tissue. **No prognostic or predictive claims** — that space is occupied (Zhang et al. 2025).

---

## What this project is not

- Not a drug-discovery project. No therapeutic claims.
- Not a claim to have discovered paralog-specific enhancer preferences — that is Plotnik et al. 2024, cited as foundation.
- Not a predictive/clinical model. MOES is a prioritisation heuristic.
- Not a maximal-technology showcase. Layers that do not serve the hypothesis are excluded by design (see gap statement §4).

---

## Success criteria

1. Every figure and table regenerable from a clean checkout with one command.
2. Every prioritised hub traceable to its supporting evidence per layer.
3. Discovery / validation / replication roles declared *before* analysis and never mixed.
4. At least one pre-registered negative control passes (paralog-shared housekeeping promoters must NOT surface as paralog-specific hubs).
5. Honest reporting if the central hypothesis fails.

---

## Persistent artefacts

| Document | Purpose |
|---|---|
| `01_project_brief.md` | This file |
| `02_gap_statement.md` | Gap + novelty assessment (FINALISED) |
| `03_dataset_inventory.md` | Every dataset, role, access status, traps |
| `04_analysis_architecture.md` | Pipeline, GRN method choice, MOES design |
| `05_dependency_inventory.md` | R/Python/system/GitHub requirements |
| `06_risk_log.md` | Open risks and mitigations |
| `07_milestone_roadmap.md` | Phased plan with gates |
| `08_qc_checklist.md` | To be written in Phase 3 |
| `09_figure_list.md` | To be written in Phase 3 |

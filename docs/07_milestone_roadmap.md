# Milestone Roadmap

Each milestone has a **gate** — an explicit, checkable condition. Do not pass a gate on optimism.

---

## ✅ M0 — Scoping (COMPLETE, 2026-07-26)
Phase 1 candidate gaps → Gap A selected → Phase 1b novelty validation.
**Outcome:** original framing partially superseded by Plotnik et al. 2024; gap re-scoped and finalised.
**Gate passed:** a defensible, non-overlapping gap exists with verified public data.

---

## ✅ M1 — Specification (COMPLETE, 2026-07-26)
Project brief, gap statement, dataset inventory, analysis architecture, dependency inventory, risk log, roadmap.
**Gate passed:** every omics layer justified or excluded; GRN method chosen with reasoning; MOES designed without arbitrary weights.

---

## ✅ M2 — Decisions before any code (COMPLETE, 2026-07-26)
1. **Single-cell layer: DROPPED.** Synapse cannot be authorised; the pre-agreed decision rule was invoked rather than substituting CDX data (R-03 closed). Must appear in README limitations.
2. **Hardware: measured, not assumed.** 948 GB disk free, 10 GiB RAM, 6 cores, R 4.6.1 — WSL2 Ubuntu-24.04. RAM is the binding constraint; chromosome-wise bedGraph processing becomes mandatory (R-13 opened).
3. **Repo: `sclc-myc-regulatory-hubs`**, public, SSH. Keys already work (`ImmunoScholar`); R-12 closed. Build root `/home/priya/projects/sclc-myc-regulatory-hubs` on the native Linux FS.
**Gate passed:** all three resolved; two risks closed, one opened.

---

## M3 — Repository & environment (Phase 3)
Scaffold, `renv` init + snapshot, `.gitignore`, README skeleton, apt dependencies installed, first commit and push.
**Gate:** `renv::restore()` works from a clean clone; `sessionInfo()` captured; no data files tracked by git.

---

## M4 — Data acquisition & QC (Phase 4)
Download scripts with checksums, manifest (incl. pinned DepMap release), file-existence / dimension / sample-count / genome-build assertions, a rendered QC report.
**Gate:** every dataset in the inventory either downloaded and QC-passed, or formally dropped with a logged reason. Genome build asserted for every coordinate file.

---

## M5 — Regulatory layer & regulons (Aim 1)
Consensus regions → signal quantification → active regions → super-enhancers → peak-to-gene linking → paralog regulons.
**Gate (hard):** our active-region counts are within a defensible range of Plotnik's reported 18,823 / 4,017 / 5,688. **Substantial discordance halts the project for method review — it is not explained away.** MYCN/MYC overlap should approximate their reported ~84%.

---

## M6 — Patient-tumour translation (Aim 2)
Regulon scoring in GSE60052; replication in George 2015; **lineage-TF confounding analysis**.
**Gate:** the lineage-TF analysis is complete and reported *whatever it shows*. If paralog regulons are inseparable from lineage programs, the project's conclusion changes here (R-01) and the remaining aims are reframed accordingly rather than abandoned.

---

## M7 — Functional & network evidence (Aim 3 + GRN)
DepMap selective dependency, drug-response association, GENIE3 network importance.
**Gate:** each evidence layer produces a per-gene score with documented missingness; no layer silently drops genes.

---

## M8 — MOES integration (Aim 4)
Two-stage RRA, layer-correlation matrix, leave-one-domain-out stability, permutation FDR, positive/negative controls, Hallmark + Jung 2017 benchmark.
**Gate (hard):** the **pre-registered negative control passes** — paralog-shared housekeeping promoters (e.g. `RPS26`) must not rank as paralog-specific hubs. Failure means the framework lacks specificity and must be fixed or reported as failed, not quietly adjusted until it passes.

---

## M9 — Spatial coherence (Aim 5, restricted)
GeoMx coherence and regional heterogeneity, with panel-coverage fraction reported.
**Gate:** no prognostic, predictive, or immune-exclusion claims appear anywhere in the output.

---

## M10 — Figures, tables, report (Phase 6)
Publication-grade figures, publication-ready tables, Quarto manuscript-style report.
**Gate:** every figure interpretable without external explanation; every table has notes where interpretation is non-obvious; the "heuristic, not predictive" statement present throughout.

---

## M11 — Release readiness (Phase 7)
Clean checkout test, path verification, dependency docs, reproducibility notes, commit history tidy, release-ready README, optional GitHub Actions.
**Gate (hard):** a fresh `git clone` on a clean machine reproduces every figure. Not "should" — actually tested.

---

## Ordering rule

M5 → M6 → M8 is the critical path. M7 and M9 can proceed in parallel with M6. **M8 must not begin before M6 completes** — the lineage-TF result may change what "a MYC-associated hub" means, and integrating evidence before that question is settled would bake a confounder into the framework.

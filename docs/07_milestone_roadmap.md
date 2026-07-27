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

## ✅ M3 — Repository & environment (Phase 3) — COMPLETE, 2026-07-26
Scaffold, `.gitignore` (written first, verified adversarially), git init + remote, `renv` init + snapshot, package install, project contract / decision log / research journal, environment tests, commits pushed.
**Done:** **35/35 packages** installed, **231 locked** in `renv.lock`; `data/metadata/session_info.txt` and `environment.yml` written; `tests/test_environment.R` passing (64 assertions, including config invariants: hg19, chromosome-wise processing, no MOES weights, GSVA as the secondary scorer); pushed to `github.com/ImmunoScholar/sclc-myc-regulatory-hubs`; no data tracked.
**System deps:** `libmagick++-dev` + `pandoc` installed — the latter satisfies the M10 reporting dependency ahead of schedule.
**Gate PASSED:** clean-clone test run and passing — fresh `git clone` → `renv::restore()` → full test suite green, executed inside the clone (`scripts/00_setup/06_clean_clone_test.sh`). Caveat: the renv cache is shared, so this validates the lockfile, not a cold-machine source build.

---

## ✅ M4 — Data acquisition & QC (COMPLETE, 2026-07-26)
Manifest generated from live GEO filelists, resumable checksum-aware downloads, verify gate, QC pipeline.
**Delivered:** 66 automated files / **12 GB** acquired in 6 h 22 m (2 concurrent, zero failures); `03_verify.R` 66/66 pass with the SHA256 ledger established; `04_qc_report.R` **167/167 checks pass** including full-stream hg19 build assertions on all 28 bedGraphs.
**Gate PASSED.** Every dataset downloaded and QC-passed, or formally dropped with a logged reason. Genome build asserted for every hg19 coordinate file — and the assertion is now capable of failing (see R-16).

**Five findings that changed the analysis, none of which raised an error on their own:**
1. **Ensembl seqnames in the keystone** (`20`, not `chr20`) against UCSC annotation everywhere else. Mixing them returns zero overlaps, not an error. `R/genome_utils.R` now mandatory for every coordinate join (D-019, R-16).
2. **The build assertion was passing on nothing** — it compared `20` to `chr20`, matched zero rows, found zero violations, reported success. It guarded R-05, the highest-consequence technical risk. Now refuses to pass on fewer than 20 compared chromosomes.
3. **GSE269424 is a TF-overexpression experiment**, not native ATAC. Only the four EGFP control arms are used; the ASCL1/NEUROD1 arms have deliberately remodelled chromatin (D-015).
4. **Lineage-TF controls barely overlap the keystone** — POU2F3 3 lines, ASCL1 1, NEUROD1 none. R-01's confounding test has heterogeneous resolution (R-14).
5. **GSE261348 has a deposit gap** — 9 high-depth unflagged segments annotated but not deposited, while a 40×-under-sequenced flagged one was. Characterised and pinned; slide `IMF-001/002` excluded at M9 (R-17).

**Outstanding, non-blocking:** DepMap manual download + release pin, needed only at M7 (D-016).

---

## M5 — Regulatory layer & regulons (Aim 1)
liftOver of hg38 intervals → consensus regions → chunked signal quantification → active regions → super-enhancers → peak-to-gene linking → paralog regulons.

**Gate (hard) — REBUILT 2026-07-26, see D-020.** The original gate (counts near Plotnik's 18,823 / 4,017 / 5,688) was **circular**: our count is `|universe| × P(signal>quantile) × P(H3K27ac⁺)`, and a 0.75 quantile passes 25% of regions by construction, so any universe near 125k reproduces those numbers whether or not the pipeline works. It could not fail.

Replaced with **relative, threshold-invariant** criteria:
1. **MYCN ⊂ MYC at 0.84 ± 0.15**, held across the sensitivity range.
2. **MYCL1 not a MYC subset** (overlap < 0.50).
3. **Distal-fraction contrast**, MYC-amplified > MYC-expressing, difference ≥ 0.10 (published 0.39 vs 0.12).
4. **Motif validation**: paralog E-box central dinucleotides enriched over shuffled background — the only criterion independent of every parameter we chose.

Counts are reported across `sensitivity_quantiles`, descriptively, never as pass/fail.
**Substantial discordance on 1–4 halts the project for method review — it is not explained away.**

**Additional gate before M6 — regulon internal validity.** Reformulated (D-028): the original spec was uncomputable (needed unverified amplification and absent expression) *and* circular (a regulon built from a paralog's own lines scores high there by construction).

**M5 STATUS 2026-07-26 — substantially complete.**

| step | result |
|---|---|
| liftOver | 0.91% loss, 687,962 intervals |
| universe | 102,334 → **97,106** after blacklist; TSS 8.55×, H3K27ac 2.83× |
| signal | 28 tracks quantified; fold-over-background (D-027 rationale) |
| **gate criterion 1** | MYCN⊂MYC **0.912**, spread 0.037 across 50× set-size change — PASS |
| **gate criterion 2** | MYCN 5.50× vs MYCL1 4.63×, OR 3.14, p=1.1e-60 — PASS (magnitude modest) |
| gate criterion 3 | **BLOCKED** → M7 (needs DepMap copy number, D-026) |
| **gate criterion 4** | paralog E-box specificity **3/3**, χ² p=1.6e-101 — PASS |
| super-enhancers | 170–699 per line (1.9–6.6% of stitched); amplicon-resident negligible |
| peak-to-gene | 48,756 links; null enrichment **6.06×**; only **32.6%** are nearest-gene |
| **regulons** | LOO 8/9 · Jaccard 0.059–0.130 · **2 of 3 programmes validated** |

**Regulon programmes:** MYC → **neurogenesis** (p=7.0e-6, reproduces Plotnik from an independent pipeline) · MYCN → **Hallmark/housekeeping** (p=3.6e-4) · MYCL1 → **none** (best p=0.094, FAIL).

**Two findings to carry forward:** (i) Plotnik grouped MYC *and* MYCN enhancer targets under neurogenesis; we reproduce it for MYC only — MYCN is housekeeping-weighted (D-029). (ii) MYCL1 proceeds flagged, not validated, and is re-gated at M7 with real expression (D-030).

---

## ✅ M6 — Patient-tumour translation (Aim 2) — COMPLETE 2026-07-27
Regulon scoring in GSE60052; replication in George 2015; lineage-TF confounding at both expression and chromatin level.

**Gate PASSED as written:** the lineage-TF analysis was completed and reported *whatever it showed*. It showed R-01 materialising, and the project's conclusion changed here rather than the analysis being abandoned.

| result | GSE60052 (n=79) | George 2015 (n=81) |
|---|---|---|
| lineage > paralog unique variance | **3/3** | **3/3** |
| unique lineage R² | 0.355–0.428 | 0.244–0.426 |
| unique paralog R² | 0.001–0.068 | 0.000–0.019 |
| paralog associations surviving adjustment | 1/3 (not robust) | **0/3** |
| MYC vs NE score (Ireland 2020) | **−0.590** | **−0.466** |

**Conclusion — replicated across 160 patients:** paralog-resolved regulatory programmes do **not** retain paralog identity in patient tumours independently of neuroendocrine lineage state. The mechanism (MYC/NE antagonism) reproduces independently in both cohorts.

**Per-gene rescue also fails** (D-033): regulon genes survive lineage adjustment at background rate (enrichment 0/3; MYCN zero of 33,683 genes). The transcriptional domain is dropped from MOES.

**Chromatin-level test (R-14), suggestive not established:** paralog-*specific* regions are ~3× depleted for lineage-TF binding relative to shared regions (OR 0.29–0.33) — which would locate the failure in the expression readout rather than in shared occupancy. But the POU2F3 peak sets fail a biology sanity check (the POU2F3-driven line H1048 has the fewest peaks), so this needs CCLE corroboration at M7.

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

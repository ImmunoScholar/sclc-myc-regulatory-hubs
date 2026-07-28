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

## ✅ M5 — Regulatory layer & regulons (Aim 1) — COMPLETE 2026-07-28
liftOver of hg38 intervals → consensus regions → chunked signal quantification → active regions → super-enhancers → peak-to-gene linking → paralog regulons.

**Gate (hard) — REBUILT 2026-07-26, see D-020.** The original gate (counts near Plotnik's 18,823 / 4,017 / 5,688) was **circular**: our count is `|universe| × P(signal>quantile) × P(H3K27ac⁺)`, and a 0.75 quantile passes 25% of regions by construction, so any universe near 125k reproduces those numbers whether or not the pipeline works. It could not fail.

Replaced with **relative, threshold-invariant** criteria:
1. **MYCN ⊂ MYC at 0.84 ± 0.15**, held across the sensitivity range.
2. **MYCL1 not a MYC subset** (overlap < 0.50).
3. **Distal-fraction contrast**, MYC-amplified > MYC-expressing, difference ≥ 0.10 (published 0.39 vs 0.12).
4. **Motif validation**: paralog E-box central dinucleotides enriched over shuffled background — the only criterion independent of every parameter we chose.

Counts are reported descriptively across the fold grid, never as pass/fail.

**FINAL GATE RESULT (closed 2026-07-28): 3 PASS / 1 FAIL.**

| criterion | observed | result |
|---|---|---|
| 1 · MYCN⊂MYC ≈ 0.84 | **0.912**, spread 0.037 across a 50× set-size change | **PASS** |
| 2 · MYCN vs MYCL1 differential nesting | 5.50× vs 4.63×, OR 3.14, p 1.1e-60 | **PASS** (magnitude modest) |
| 3 · distal-fraction contrast | +3.3 pts vs published +27; within-group spread (12.5) exceeds it | **FAIL** |
| 4 · paralog E-box specificity | 3/3 own-motif top share, χ² p 1.6e-101 | **PASS** |

Criterion 3 was unblocked at M7 by DepMap copy number, which established **exactly two** MYC-amplified lines (H524 log2 6.73, H211 2.35) against MYC-expressing H1048 and SHP77 — confirming Plotnik's stated design and correcting the project registry, which had listed all five MYC ChIP lines as amplified (D-035).

**The failure is reported, not reformulated.** Unlike criteria 2 and 4 — where the *test* was demonstrably mis-specified and was rebuilt — criterion 3's test is sound and the data simply do not support the published magnitude. Most likely structural: our shared ATAC universe is ~84% distal by construction, compressing the achievable range, where Plotnik called peaks de novo per line. This is the second discordance tracing to the shared-grid design (D-023) after MYCL1's inflated occupancy overlap.

**Additional gate before M6 — regulon internal validity.** Reformulated (D-028): the original spec was uncomputable (needed unverified amplification and absent expression) *and* circular (a regulon built from a paralog's own lines scores high there by construction).

**M5 STATUS — COMPLETE 2026-07-28** (criterion 3 closed at M7, D-035).

| step | result |
|---|---|
| liftOver | 0.91% loss, 687,962 intervals |
| universe | 102,334 → **97,106** after blacklist; TSS 8.55×, H3K27ac 2.83× |
| signal | 28 tracks quantified; fold-over-background (D-027 rationale) |
| **gate criterion 1** | MYCN⊂MYC **0.912**, spread 0.037 across 50× set-size change — PASS |
| **gate criterion 2** | MYCN 5.50× vs MYCL1 4.63×, OR 3.14, p=1.1e-60 — PASS (magnitude modest) |
| gate criterion 3 | **FAIL** — evaluated at M7 once DepMap copy number resolved amplification (D-035); +3.3 pts vs published +27 |
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

## ✅ M7 — Functional & network evidence (Aim 3 + GRN) — COMPLETE 2026-07-28
DepMap selective dependency, drug-response association, GENIE3 network importance.
**Gate PASSED:** every admitted layer produces a per-gene score with documented missingness.

| layer | outcome |
|---|---|
| CRISPR selective dependency | Enriched **pooled** across the three regulons: OR 2.71 (95% CI 1.29–5.15), p = 0.0048, on 11 selective genes in 984 tested. **Not separable per paralog** — ORs 2.45 / 2.21 / 3.09 on 4, 4 and 5 genes, intervals almost fully overlapping |
| Paralogs' own dependency | None is SCLC-selective. MYC is essential in both groups but *less* so in SCLC (+0.503, FDR 0.18) |
| Drug response (D-047) | MYCL expression tracks BET-inhibitor AUC in **Sanger only** — GDSC2 birabresib rho +0.621 (NE-adj +0.631), GDSC2 molibresib +0.617 (+0.587), GDSC1 PFI-1 +0.426 (+0.429). **Does not replicate in PRISM** (+0.046, p = 0.84). GDSC1/GDSC2 share a platform, so they are not two independent replications. Recorded exploratory and provisional |
| GENIE3 network (D-037) | **EXCLUDED on both pre-conditions** — 59 SCLC lines is below what tree-ensemble inference needs, and MYC expression tracks NE state across them (rho −0.441), the same confound that barred a tumour-based network |
| Amplification (D-035) | DepMap copy number established **exactly two** MYC-amplified lines (H524 log2 6.73, H211 2.35), correcting the project registry and unblocking M5 criterion 3 |

---

## ✅ M8 — MOES integration (Aim 4) — COMPLETE 2026-07-28
Two-stage RRA, layer-correlation matrix, leave-one-domain-out stability, permutation FDR, positive/negative controls, Hallmark + Jung 2017 benchmark.

**Result: no prioritised hub list.** MOES ran on **2 of 4** domains — transcriptional dropped as lineage-confounded at aggregate *and* per-gene level (D-033), network excluded (D-037). With two domains the two-stage hierarchy is degenerate and is reported as a single aggregation. Across **10,387** genes, **no gene reaches FDR < 0.05 for any paralog**; the lowest FDR obtained is **0.358**. Top-of-list overlap between domains is real but weak — max 2.03× (MYC), 2.19× (MYCN), 1.72× (MYCL1), global p 0.019 / 0.008 / 0.065 — while genome-wide Spearman between domains is ≈ 0 (−0.023, +0.012, −0.008). Aggregate detectable, per-gene unattributable.

**Gate (hard) — NOT MET, and recorded as uninformative rather than passed.**

| control | outcome |
|---|---|
| Negative (hard gate) | **VACUOUS.** `RPS26` is not in the MOES universe, and MOES produces no hubs, so nothing *can* fail the control. A control that cannot fail is no evidence of specificity. Calling it "passed" would repeat the vacuous-pass error caught at R-16, where a build assertion compared `20` to `chr20`, matched zero rows and reported PASS |
| Positive (MYC→BCL2/MCL1) | **NOT MET.** BCL2 sits in the *MYCN* regulon; best MYC rank 128 at FDR 0.93. MCL1 is in no regulon. Consistent with the null rather than evidence against it |
| HALLMARK_MYC_TARGETS | MET |
| Jung 2017 signature (D-046) | MET — curated from the source table and run. Regulon overlap OR 4.49 / 6.62 / 7.33; the signature agrees with the regulon scores (strongest for MYCN, +0.438) and does **not** track MYC expression in tumours (rho −0.185, n.s.) |

---

## ✅ M9 — Spatial coherence (Aim 5, restricted) — COMPLETE 2026-07-28
GeoMx coherence and regional heterogeneity, with panel-coverage fraction reported.
**Gate PASSED:** no prognostic, predictive or immune-exclusion claim appears in any output (enforced by `02_lint_outputs.R`).

Only the **MYC** regulon clears the coverage gate, at 56 of 500 members (11.2%); MYCN reaches 8.4% and MYCL1 9.6%, and neither is scored. MYC clears by 1.2 points and MYCL1 misses by 0.4 — knife-edge verdicts against a threshold raised mid-project, reported as such. Everything downstream is a **56-gene proxy**, not the regulon. Variance splits roughly half between slides and half within (0.554 / 0.446); restricted to the three unambiguously single-patient slides it is 0.701 / 0.299, so about a third of the score is regional variation inside one tumour. Lineage dominance replicates in this third modality: lineage explains 0.474 and 0.508 against 0.019 and 0.010 for MYC.

---

## ✅ M10 — Figures, tables, report (Phase 6) — COMPLETE 2026-07-28
**Gate PASSED.** 8 figures (4 main, 4 supplementary), each with a caption file; 6 publication tables in `results/tables/publication/`; Quarto report reading every value from files the analysis wrote, refusing to render if they are absent. The "heuristic, not predictive" statement is present in the report, README and MOES figure legend, and is machine-checked by `scripts/08_report/02_lint_outputs.R` along with the contract's prohibited vocabulary.

The paralog palette was replaced after failing its own accessibility check (D-041): the M1 palette reached only ΔE 12.7 under deuteranopia against a floor of 15.

---

## ✅ M11 — Release readiness (Phase 7) — COMPLETE 2026-07-28
**Gate PASSED for what it covers.** `scripts/00_setup/07_clean_clone_verify.sh` clones the repository, rebuilds all 8 figures **from the clone with no data present**, and compares them byte-for-byte against the committed set: 8 byte-identical PNGs, 0 differing. Decision audit 34 PASS / 0 FAIL / 0 SKIP; output lint 0 failures.

**What the gate does not cover, stated rather than implied.** It does not install 232 packages on a bare machine and does not re-download the 12 GB of source data. Every clean-clone run so far was made inside a WSL instance with a warm renv cache — `renv::restore()` reported 215 packages "linked from cache" in 0.7 s. So the lockfile is proven complete and internally consistent; **compiling these packages from source on a machine that has never built them remains unproven** (D-048). The run did surface two previously undocumented system packages, `cmake` (for `fs`) and `gsfonts` (for `magick`), now in the README install line. Closing this properly needs a container or CI runner with no R library and no renv cache.

Three defects found by running the gate rather than trusting it (D-048): a hardcoded `/mnt/c/Users/Priya/Downloads` path in `06_drug_response.R`; an audit that reported 3 FAIL on a correct clone because it tested deliberately gitignored artefacts; and a gitlink to a nested verification clone committed by `git add -A`. All three fixed.

---

## Ordering rule

M5 → M6 → M8 is the critical path. M7 and M9 can proceed in parallel with M6. **M8 must not begin before M6 completes** — the lineage-TF result may change what "a MYC-associated hub" means, and integrating evidence before that question is settled would bake a confounder into the framework.

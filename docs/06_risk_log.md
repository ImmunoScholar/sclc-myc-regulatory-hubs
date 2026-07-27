# Risk Log

Opened 2026-07-26. Update in place — do not delete entries; mark them CLOSED with a date and outcome.

Severity: **HIGH** = threatens the central claim · **MED** = threatens an aim · **LOW** = nuisance

---

### R-01 · Scientific · HIGH · OPEN
**"Paralog-specific" hubs turn out to be lineage-TF targets in disguise.**
MYC-amplification correlates with SCLC subtype (Plotnik: MYC-amp correlates with NEUROD1 expression; MYCL with ASCL1). A hub attributed to MYC may simply be an ASCL1/NEUROD1 target in a co-varying background.
**Mitigation:** the lineage-TF confounding analysis is a *primary* planned analysis, not a robustness check — using ASCL1 (GSE281524), NEUROD1 (GSE210113) and POU2F3 (GSE249362) ChIP references. Report the overlap fraction openly.
**If it fails:** the honest finding becomes "MYC-paralog regulatory programs in SCLC are largely inseparable from lineage identity" — a legitimate, publishable negative result. The project does not collapse; the conclusion changes.

---

### R-02 · Technical · HIGH · OPEN
**GSE230649 ships bedGraph only (8.6 GB), no peak files.**
Re-alignment of 28 samples from SRA is not laptop-feasible.
**Mitigation:** quantify bedGraph signal over a consensus ATAC-defined region universe instead of de novo peak calling (architecture Part B/E).
**Residual risk:** our region boundaries will not match Plotnik's MACS2 calls. Must be declared as a methodological deviation, with a region-count comparison against their reported 18,823 / 4,017 / 5,688 as a sanity check. If our counts are wildly discordant, the quantification approach needs revisiting before proceeding.

---

### R-03 · Data access · HIGH · **CLOSED 2026-07-26 — single-cell layer DROPPED**
**No good patient-tumour single-cell dataset is accessible.**
Chan et al. 2021 atlas (the right dataset: 155k cells, 21 patient specimens) lists no GEO accession — HTAN/Synapse only. The Synapse connector cannot be authorised in a non-interactive session. The GEO fallback GSE138474 is CDX/cell-line-derived on pre-10x chemistry.
**Outcome:** the pre-agreed decision rule was invoked. The single-cell layer is **formally excluded** from the project rather than filled with a CDX substitute — a model-derived layer does not support a patient-translation argument, and including it for completeness would violate the project's biological-coherence principle.
**Consequences:** Aim covering cell-state specificity/plasticity is removed; `pyscenic`, `synapseclient`, `Seurat`, `SingleCellExperiment`, `scater`, `zellkonverter` are dropped from the dependency inventory. Patient translation now rests on the two bulk tumour cohorts (R-09 applies with more force as a result).
**Reversibility:** if Synapse/HTAN access is obtained later, the layer can be reinstated as an additive aim without disturbing M5–M8. This exclusion must be stated plainly in the README limitations section — not omitted.

---

### R-04 · Data · MED · OPEN
**George 2015 has no copy-number data on cBioPortal.**
Verified via API: only mutations, mRNA expression, structural variants. MYC-paralog amplification status is therefore unavailable from that route for the replication cohort.
**Mitigation:** source amplification calls from the paper's supplementary tables; if unavailable, define an expression-based "paralog-high" call and label it as such everywhere. **Never silently substitute expression for amplification** — the two are not equivalent and Plotnik's own data show MYC-amplified behaves differently from MYC-expressing.

---

### R-05 · Technical · **HIGH** (escalated 2026-07-26 from MED) · OPEN
**Genome build mismatch (hg19 vs hg38) — now confirmed, not hypothetical.**
Builds were verified at M4 from each series' own `!Sample_data_processing` declaration. The keystone (GSE230649), NEUROD1 (GSE210113) and POU2F3 (GSE249362) are hg19. **All three supporting ATAC datasets (GSE269424, GSE256345, GSE281523) and the ASCL1 ChIP (GSE281524) are hg38.** **23 of 66 manifest files require lifting.**

Escalated because the consensus region universe requires ≥2 independent ATAC datasets and only one of them is in the project build — so the central Aim 1 object cannot be built without crossing builds. The lineage-TF controls are also split across builds, and that analysis is *primary* (R-01), not a robustness check.

**Mitigation (policy D-014):** never lift continuous signal; call intervals in the native build and lift only intervals; report per-dataset lift loss rate; prefer sources publishing intervals over signal; assert build at load — `04_qc_report.R` fails any hg19-declared file containing an interval past its chromosome's length in `hg19.chrom.sizes`.
**Residual risk:** liftOver loses intervals in regions restructured between builds, non-uniformly across the genome. Loss is reported per dataset, and a lifted region universe is never presented as equivalent to a natively-hg19 one.

---

### R-14 · Scientific · **HIGH** · OPEN *(opened 2026-07-26)*
**The lineage-TF controls barely overlap the keystone cell lines, so the confounding analysis cannot be fully within-line.**

Verified cell-line composition:

| Control | Lines | Overlap with keystone MYC-family ChIP lines | Build |
|---|---|---|---|
| POU2F3 (GSE249362) | NCIH1048, NCIH211, NCIH526 | **3 lines, 2 paralogs** (H1048/H211 MYC-amp, H526 MYCN-amp) | hg19 ✓ |
| ASCL1 (GSE281524) | H1836, SHP77 | **1 line** (SHP77, MYC-amp) | hg38 |
| NEUROD1 (GSE210113) | H446 only | **none** | hg19 |

Keystone MYC ChIP lines are H1048, H211, H524, H847, SHP77; MYCL1 COLO668, H889; MYCN H526, H69.

**Why it matters.** R-01 is that "paralog-specific" hubs may be lineage-TF targets in disguise. The cleanest test is within-line co-occupancy: is a MYC-bound region also lineage-TF-bound *in the same cells*? That is available for POU2F3 (3 lines) and marginally for ASCL1 (1 line), and **not at all for NEUROD1**.

**Mitigation.** The confounding analysis is stratified by the evidence actually available, and its resolution is stated per TF rather than implied uniform:
- **POU2F3** — within-line co-occupancy in 3 keystone lines. The strongest arm, and in the project's own build.
- **ASCL1** — within-line in SHP-77 only; n=1, descriptive, no inference from it.
- **NEUROD1** — **subtype-level only.** H446 is SCLC-N and the MYC-amplified keystone lines are NEUROD1-high, so the comparison is "does the MYC-specific region set overlap the NEUROD1 repertoire measured in a different SCLC-N line", which cannot separate line-specific from subtype-specific effects.

**This does not change the frozen design** — the lineage-TF confounding analysis remains a primary analysis using ASCL1/NEUROD1/POU2F3 references (D-002 era decision). What changes is that its **resolution is heterogeneous**, and any claim of lineage-independence must be qualified per TF. A uniform "we controlled for lineage TFs" statement would be false.
**If it bites:** if paralog-specific hubs cannot be distinguished from POU2F3/ASCL1 targets in the lines where within-line data exist, that is the R-01 negative result and is reported as such.

---

### R-17 · Data · MED · OPEN *(opened 2026-07-26)*
**GSE261348 appears to have a deposit error: nine high-quality segments are annotated but their counts were never deposited.**

Found by QC comparing the workbook's `SegmentProperties` sheet (184 annotated segments) against `TargetCountMatrix` (175 deposited ROI columns). The nine were characterised rather than written off as a QC exclusion, and the pattern is **inverted from what a QC exclusion looks like**:

| | 9 absent segments | median deposited segment |
|---|---|---|
| Raw reads | **28.9 M** | 7.7 M |
| Aligned reads | 27.7 M | 7.4 M |
| Deduplicated | 10.4 M | 2.9 M |
| Sequencing saturation | 61.4 % | 60.9 % |
| Nuclei count | 2,706 | 2,503 |
| QC flags | **none** | — |

All nine are ROIs 002–010 of slide `IMF-001/002`. The **one** segment from that slide that *was* deposited (ROI 001) has 694,776 raw reads — roughly 40× fewer than its slide-mates — carries the `Low Negative Probe Count` flag, and has 241 nuclei.

So the deeply-sequenced, unflagged segments were dropped and the shallow flagged one was kept. That is not a filtering decision; the most plausible explanation is an upstream deposition error. Direction of the mismatch was checked both ways — 0 matrix columns are absent from `SegmentProperties` — so it is genuine absence, not a naming mismatch.

**Consequence.** Slide `IMF-001/002` (from the naming, likely patients IMF-001 and IMF-002) effectively has no usable spatial data. The cohort's usable ROI count is 175, of which 1 is unreliable.

**Mitigation.** Declared in `config/params.yml` under `spatial.known_segment_gaps`, including the M9 action: **exclude slide `IMF-001/002` entirely** rather than retain its single under-sequenced ROI. QC now verifies the gap still matches the characterised shape — an undeclared gap, or a declared one whose shape changed on re-download, fails. The discrepancy is not tolerated, it is pinned.

**Impact on conclusions.** Low. Spatial is already restricted to a coherence and heterogeneity check with no prognostic claims (D-001). But the reduced patient count must appear in the spatial figure legends alongside the panel-coverage fraction.

---

### R-16 · Technical · **HIGH** · OPEN *(opened 2026-07-26)*
**Sequence-name conventions differ between the keystone data and every annotation resource — and the QC check meant to catch build problems silently passed on nothing.**

Discovered when QC ran against real data for the first time. Two distinct problems:

**(a) The naming mismatch.** GSE230649 bedGraphs use **Ensembl** sequence names — bare `1`, `20`, `X` (verified: first line is `20\t0\t60000\t0.00000`; full set is 1–22, X, Y with no scaffolds). Every annotation resource in the project uses **UCSC** names: `TxDb.Hsapiens.UCSC.hg19.knownGene`, `hg19.chrom.sizes`, and the UCSC liftOver chains. The coordinates are identical hg19; only the convention differs.

This does not raise an error. `GenomicRanges` operations between objects with disjoint seqlevels return **zero overlaps**, which reads as a biological result — "MYC binds none of these regions" — rather than as a bug.

**(b) The vacuous assertion — the more serious of the two.** `04_qc_report.R` asserted hg19 bounds by merging bedGraph seqnames against `hg19.chrom.sizes`. Because `20` never matches `chr20`, the merge returned zero rows, found zero violations, and **recorded a PASS**. The check guarding R-05, the project's highest-consequence technical risk, could not fail.

A second flaw compounded it: assertions ran on a 200,000-line head sample, but these bedGraphs are written in **BAM-header order** (`20 21 22 1 3 2 5 4 …`), not sorted, so a head sample contains exactly one chromosome. Each file is 45.5M lines.

**Mitigation.**
- `R/genome_utils.R` provides `harmonise_seqnames()`, and **every coordinate join must pass through it**. Established at M4 so M5 inherits it rather than rediscovering the problem.
- `assert_hg19_bounds()` **refuses to return a pass** when fewer than 20 chromosomes were actually compared. A check that had nothing to compare now reports failure, not success.
- Build assertions stream each file in full by default. `--fast` skips them and records that they were skipped as a *failure*, so a fast run can never be mistaken for a clean one.
- QC additionally asserts that all 23 analysis chromosomes are present, which catches truncation.

**Lesson recorded.** Two QC defects in this project have now been of the same kind: a check that returns a pass because it silently examined nothing (this one, and the `!isTRUE(vector)` resumability check at M4). **A QC check must be tested against data known to violate it**, or it is decoration. Any future assertion added to this project should be demonstrated failing before it is trusted.

---

### R-15 · Data · LOW · **CLOSED 2026-07-26**
**GSE60052 deposits no raw counts.**
Verified by reading the file header: the only matrix is `normalized.log2`, gene-symbol rows, 86 columns.
**Outcome:** DESeq2 cannot be used on this cohort (it requires raw integer counts). Differential analysis uses **limma** on the log2 matrix. Both packages were already in the dependency inventory, so this selects between existing options rather than changing the design. Two parsing traps recorded and handled at load: sample headers carry leading whitespace, and normals are encoded by a `.normal` suffix (7 normal / 79 tumour, asserted in QC).

---

### R-06 · Scope · MED · OPEN
**Spatial layer partially scooped.**
Zhang et al. 2025 (Clin Cancer Res) already links ecDNA MYC-paralog amplification to spatial immune exclusion.
**Mitigation:** spatial re-scoped to a restricted coherence/heterogeneity check with no prognostic or immune-exclusion claims (gap statement §4). Cite Zhang explicitly as prior art.

---

### R-07 · Data · MED · OPEN
**GeoMx panel covers only ~1,800 genes.**
Most regulon members will not be measurable spatially.
**Mitigation:** report the measurable intersection fraction on every spatial figure; treat spatial as supporting evidence with declared partial coverage; never let a gene's absence from the panel count against it in MOES (RRA handles missing layers natively — this is one reason it was chosen).

---

### R-08 · Reproducibility · MED · OPEN
**DepMap re-releases quarterly; gene-effect values change between releases.**
**Mitigation:** pin and record the exact release string in the download manifest. Any figure citing DepMap names the release. Never write "DepMap" without a release ID.

---

### R-09 · Statistical · MED · OPEN
**Small effective sample sizes.** Two lines per paralog in the ChIP data; 79 and ~81 tumours in the bulk cohorts.
**Mitigation:** no per-line statistical claims from n=2; treat paralog groups as descriptive at the chromatin level and reserve inference for the tumour cohorts. Report confidence intervals, not just p-values. Power considerations stated up front in the report.

---

### R-10 · Interpretive · MED · OPEN
**MOES could be over-read as a predictive model.**
**Mitigation:** the "heuristic prioritisation, not a predictive or clinically validated model" statement appears in the README, every output table header, and every relevant figure legend. Prohibited vocabulary list enforced at report review (`predicts`, `driver of`, `therapeutic target`).

---

### R-11 · Technical · LOW · OPEN
**SCLC-CellMiner has no documented API or R package.**
`rcellminer` covers NCI-60, not the SCLC CDB.
**Mitigation:** manual/scripted UI download with a manually recorded provenance note (date, file, portal version). Acceptable — this is a one-off download.

---

### R-12 · Infrastructure · LOW · **CLOSED 2026-07-26**
**GitHub push blocked by missing SSH credentials.**
**Outcome:** verified in WSL2 Ubuntu-24.04 — `ssh -T git@github.com` authenticates as GitHub user `ImmunoScholar` using an existing `~/.ssh/id_ed25519`. `git` 2.43.0 present; `user.name`/`user.email` already set globally. No key generation needed.
**Residual:** the `.gitignore`-before-first-commit requirement still stands and is enforced at M3.

---

### R-13 · Technical · MED · OPEN *(opened 2026-07-26)*
**RAM is 10 GB, not 16 GB.**
Measured on the target machine: 10 GiB RAM, 8 GiB swap, 6 cores, 948 GB free on the WSL2 Linux filesystem. Disk is a non-issue; memory is the binding constraint.
**Consequence:** chromosome-by-chromosome bedGraph processing is **mandatory, not optional**. Whole-genome `rtracklayer::import()` of an 8.6 GB bedGraph set will exhaust memory.
**Mitigation:** every coordinate-handling script iterates over seqnames with explicit `which=` range restriction, writes per-chromosome intermediates to `data/processed/`, and frees objects between iterations. Peak RSS asserted in QC. Swap exists but relying on it means unacceptable runtimes — treat swap use as a failure signal, not a safety net.

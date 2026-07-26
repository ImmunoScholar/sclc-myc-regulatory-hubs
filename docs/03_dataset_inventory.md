# Dataset Inventory

All accessions verified against GEO / cBioPortal / source portals on **2026-07-26**.
Legend: ✅ verified this session · ⚠️ verified with a caveat · ❓ not yet verified

---

## Role definitions (declared before analysis — do not change post hoc)

- **DISCOVERY** — used to construct regulons and nominate hubs.
- **VALIDATION** — independent data used to test discovery-derived claims.
- **EXTERNAL REPLICATION** — different platform/population; last line of evidence.
- **REFERENCE** — annotation or prior signatures; not a test set.

---

## D1 · Regulatory / chromatin layer — DISCOVERY

| Accession | Content | Role | Status |
|---|---|---|---|
| **GSE230649** | ChIP-seq: **MYC, MYCN, MYCL1, H3K27ac** + ATAC-seq, 10 SCLC lines, 28 samples. Paralog amplification declared by the authors: MYC-amp — H1048, H211, H524, H847, SHP77; MYCN-amp — H526, H69; MYCL-amp — COLO668, H889 | DISCOVERY (keystone) | ⚠️ |
| **GSE269424** | ATAC-seq, SCLC cell lines, 8 samples | VALIDATION (accessibility) | ✅ |
| **GSE256345** | ATAC-seq, SMARCA4 / state plasticity, 13 samples | VALIDATION (accessibility) | ✅ |
| **GSE281523** | FFPE-ATAC + PDX-ATAC, SCLC, 12 samples | EXTERNAL REPLICATION (tissue-derived accessibility) | ✅ |
| **GSE281524** | ASCL1 ChIP-seq, SCLC lines, 10 samples | REFERENCE (lineage-TF confounder control) | ✅ |
| **GSE210113** | NEUROD1 / BET ChIP-seq, 14 samples | REFERENCE (lineage-TF confounder control) | ✅ |
| **GSE249362** | POU2F3 / ncBAF, 125 samples | REFERENCE (lineage-TF confounder control) | ✅ |

> ⚠️ **GSE230649 trap — the single biggest technical constraint in this project.**
> Processed data are distributed as **bedGraph only (8.6 GB tar), with no peak files.** Raw reads are in SRA (hg19-aligned in the source paper).
> **Decision:** do NOT re-align from SRA. Instead quantify bedGraph signal over a **consensus ATAC-defined region set**, then apply signal-threshold-based region calling. This avoids ~28 realignments, is laptop-feasible, and is arguably cleaner (regions defined by accessibility, not by per-antibody peak-caller behaviour). Consequence: our region boundaries will differ from Plotnik's MACS2 calls — this must be stated as a methodological deviation and the region-count comparison reported.
> Genome build: source is **hg19**. All coordinates must be handled in hg19 or explicitly lifted to hg38. **Mixing builds is the most likely silent failure in this project.**

---

## D2 · Bulk tumour transcriptome — VALIDATION / EXTERNAL REPLICATION

| Source | Content | Role | Status |
|---|---|---|---|
| **GSE60052** | 79 primary SCLC tumours + 7 normal lung. Processed matrix: `GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv.gz` (9.4 MB). Chinese patient cohort, Illumina HiSeq 2000 | VALIDATION (primary tumour) | ✅ |
| **cBioPortal `sclc_ucologne_2015`** (George et al. 2015) | Mutations (MAF), mRNA expression (RNA-seq, continuous + z-score), structural variants | EXTERNAL REPLICATION (different population/platform) | ⚠️ |
| **CCLE / DepMap expression** | SCLC cell-line expression, for regulon construction and correlation | DISCOVERY (support) | ❓ |

> ⚠️ **George 2015 trap:** cBioPortal exposes **no copy-number profile** for this study (verified via API — only MUTATION_EXTENDED, MRNA_EXPRESSION ×2, STRUCTURAL_VARIANT). MYC-paralog amplification status for this cohort therefore **cannot** come from cBioPortal CNA. Options: (a) source from the paper's supplementary tables, (b) derive a expression-based paralog-high call and state it as such. Must be resolved before Aim 2. **Do not silently substitute expression for amplification.**

---

## D3 · Functional / pharmacogenomic layer — VALIDATION

| Source | Content | Role | Status |
|---|---|---|---|
| **DepMap** | CRISPR gene effect (Chronos), expression, copy number, model metadata, PRISM drug repurposing | VALIDATION (dependency) | ⚠️ |
| **SCLC-CellMiner CDB** (`discover.nci.nih.gov/SclcCellMinerCDB/`) | Drug activity, DNA copy number, mutation, methylation, RNA (microarray + RNA-seq log2 FPKM+1), protein (RPPA, SWATH-MS), microRNA, histone, CRISPR, metabolomics. Cell-line sets: NCI-DTP, UTSW, CCLE-Broad-MIT, GDSC-MGH-Sanger, CTRP-Broad-MIT, MD Anderson, Achilles, CSHL | VALIDATION (drug response) | ⚠️ |

> ⚠️ **DepMap:** public and free, no login for downloads, but the portal is bot-protected so the release name could not be auto-read. **Pin the exact release string at download time and record it in the manifest** — DepMap re-releases quarterly and gene-effect values change between releases. Never cite "DepMap" without a release ID.
> ⚠️ **SCLC-CellMiner:** confirmed to expose tab-delimited downloads via the web UI. **No documented REST API or R package was found on the portal.** Assume manual/scripted UI download; record file provenance manually. `rcellminer` on Bioconductor covers NCI-60, **not** the SCLC CDB — do not assume it works here without checking.
> **From CellMiner we will use ONLY expression, copy number, and drug activity.** Methylation, proteomics, microRNA and metabolomics are explicitly excluded — they do not serve the hypothesis.

---

## D4 · ~~Single-cell layer~~ — **EXCLUDED 2026-07-26**

**Decision: the single-cell layer is dropped from the project.** Recorded here rather than deleted, because the exclusion and its reasoning are part of the scientific record and must be reproduced in the README limitations section.

| Source | Content | Why not used |
|---|---|---|
| **Chan et al. 2021 SCLC atlas** | 155,098 transcriptomes, 21 SCLC specimens / 19 patients; SCLC-A ×14, SCLC-N ×6, SCLC-P ×1 | The scientifically correct dataset, but **inaccessible.** Data availability names no GEO accession — only HTAN and `github.com/dpeerlab`, requiring HTAN/Synapse authorisation that is not available here. |
| **GSE138474** | SuperSeries (GSE138267 + GSE138418), 43 samples | Accessible but **unfit for purpose.** Predominantly CDX and cell-line-derived models, not primary patient tumours, on Illumina HiSeq 2000 (pre-10x chemistry, low cell numbers). |

> **Rationale for exclusion over substitution:** the project's whole argument is *translation into patient tumours*. A CDX/cell-line-derived single-cell layer cannot support that argument, and adding it would be inclusion-for-completeness — precisely the failure mode the project charter forbids. Dropping a layer honestly is stronger than carrying a weak one.
> **Consequence:** patient-tumour translation now rests entirely on the two bulk cohorts (D3). Risk R-09 (small effective sample size) applies with more force; state this in the limitations.
> **Reversible:** if HTAN/Synapse access is obtained later, this becomes an additive aim without disturbing M5–M8. See risk R-03 (CLOSED).

---

## D5 · Spatial layer — RESTRICTED ORTHOGONAL VALIDATION

| Accession | Content | Role | Status |
|---|---|---|---|
| **GSE261348** | GeoMx DSP, **175 ROIs / 32 ES-SCLC patients**, IMfirst phase IIIb trial (EudraCT 2019-002784-10). Files: `IMfirst_DSP_rawcounts.xlsx`, `IMfirst_DSP_normalizedcounts.xlsx` + per-patient PNGs | EXTERNAL REPLICATION | ✅ |
| **GSE261345** | GeoMx DSP, 121 ROIs, CANTABRICO cohort | EXTERNAL REPLICATION | ✅ |
| GSE263196 | Spatial transcriptomics, human SCLC, 5 samples | Optional | ✅ |

> ⚠️ **Panel constraint:** GeoMx here is a **targeted ~1,800-gene cancer/immune panel** plus custom SCLC subtype probes (NEUROD1, POU2F3, YAP1). Most regulon members will be **absent from the panel.** Any spatial analysis must be restricted to the measurable intersection, and the coverage fraction reported explicitly on every spatial figure.
> ⚠️ **Scope constraint:** treatment-outcome analysis in these cohorts is **out of scope** (see gap statement §4 — Zhang et al. 2025). Use for spatial coherence and regional heterogeneity only.

---

## D6 · Reference sets

| Source | Use |
|---|---|
| MSigDB `HALLMARK_MYC_TARGETS_V1` / `V2` | Benchmark comparator — the paralog-blind status quo our regulons are tested against |
| Jung et al. 2017 18-gene Myc activity signature | Benchmark comparator |
| Plotnik et al. 2024 reported region counts (18,823 / 4,017 / 5,688) | Sanity check on our re-derivation |
| Plotnik "MYC shared" housekeeping promoters (e.g. `RPS26`) | **Pre-registered negative control** — must NOT rank as paralog-specific hubs |
| Dammert et al. 2019 MYC→BCL2/MCL1 axis | Orthogonal positive control for MYC-specific hubs |

---

## Access actions — **ALL RESOLVED 2026-07-26**

1. ~~Connect the Synapse connector for the Chan atlas~~ → **not obtained; single-cell layer dropped** (D4, R-03).
2. ~~Confirm disk and RAM~~ → **measured:** 948 GB free, 10 GiB RAM, 6 cores (dependency inventory). RAM forces chromosome-wise bedGraph processing (R-13).
3. No credentials are needed for GEO, cBioPortal, DepMap, or SCLC-CellMiner. ✅

**Remaining manual step, at M3:** create the empty public GitHub repo `sclc-myc-regulatory-hubs` (no auto-README). SSH already authenticates as `ImmunoScholar`.

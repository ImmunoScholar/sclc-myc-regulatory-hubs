# Analysis Architecture

---

## Part A — Per-layer justification

The instruction was explicit: justify each layer or recommend excluding it. Layers included for completeness alone are removed.

| Layer | Why included | Question it answers | Contribution to hypothesis | Verdict |
|---|---|---|---|---|
| MYC/MYCN/MYCL1 ChIP-seq | Only direct evidence of paralog identity at a locus | *Which paralog binds here?* | Defines regulon membership. Irreplaceable — no motif method can substitute (see Part B) | **INCLUDE — core** |
| H3K27ac ChIP-seq | Distinguishes bound-and-active from bound-and-inert | *Is this a functioning enhancer?* | Filters binding to transcriptionally productive regions; enables SE calling | **INCLUDE — core** |
| Super-enhancer calling | Absent from Plotnik et al.; SEs mark lineage/identity genes | *Which targets sit under disproportionate regulatory investment?* | Primary hub-nomination signal | **INCLUDE — novel** |
| ATAC-seq | Plotnik showed paralog enhancer preference is *determined by* accessibility | *Is the region open, and in which lines?* | Defines the consensus region universe; independent replication across 3 datasets | **INCLUDE — core** |
| Bulk tumour RNA-seq | The translation step the field is missing | *Do these programs operate in patients?* | Tests Aim 2. Without it the project stays in vitro like the prior work | **INCLUDE — core** |
| GRN inference | Independent, expression-only view of network importance | *Is this gene central independent of chromatin evidence?* | Supplies an evidence domain not derived from chromatin | **INCLUDE — supporting** |
| DepMap CRISPR | Converts correlation into functional consequence | *Does removing it matter, selectively, in paralog-amplified lines?* | Strongest single non-transcriptional evidence | **INCLUDE — core** |
| Drug response (CellMiner/PRISM) | Pharmacological corroboration | *Does expression track drug sensitivity?* | Weakest layer (noisy, confounded by lineage) — down-weighted by design | **INCLUDE — supporting, low confidence** |
| Spatial (GeoMx) | Patient-tissue orthogonality | *Is hub expression spatially coherent and regionally variable?* | Restricted role only; panel covers ~1,800 genes | **INCLUDE — restricted** |
| scRNA-seq | Cell-state resolution of MYC programs | *Do hubs mark a discrete tumour-cell state?* | Valuable only if patient-derived — and no patient-derived atlas is accessible | **EXCLUDED 2026-07-26** — Chan atlas unobtainable (HTAN/Synapse only); GSE138474 is CDX-derived and cannot support a patient-translation claim. See R-03, dataset inventory D4 |
| Methylation | — | — | Does not test the hypothesis | **EXCLUDE** |
| Proteomics (RPPA/SWATH) | — | — | Coverage too sparse for regulon-scale scoring | **EXCLUDE** |
| microRNA / metabolomics | — | — | Off-hypothesis | **EXCLUDE** |

---

## Part B — Gene regulatory network method: decision and justification

### The decisive technical constraint

**MYC, MYCN and MYCL1 all bind the same canonical E-box (`CACGTG`).** Plotnik et al. found only subtle central-dinucleotide preferences at enhancers (`CAGATG` / `CACATG` / `CACCTG`).

Every motif-driven GRN method — SCENIC, pySCENIC, SCENIC+ — resolves TF→target links via motif enrichment against cisTarget databases. **Such methods are structurally incapable of distinguishing MYC paralogs.** They would collapse all three into a single "MYC regulon," destroying the variable the entire project is about.

This is the single most important methodological finding of Phase 1b and it dictates the design.

### Method assessment

| Method | Verdict | Reasoning |
|---|---|---|
| **SCENIC+** | **REJECT** | Requires paired scATAC + scRNA (multiome). A GEO survey found **no SCLC scATAC-seq dataset** — all public SCLC ATAC is bulk. The method cannot be run on our data. Also motif-based → paralog-blind. |
| **SCENIC / pySCENIC** | **REJECT outright** | Motif-based → paralog-blind (fatal for the primary aim). Its one conditional use — defining MYC-program-high cell states — required the scRNA layer, which is now excluded (R-03). No remaining role in this project. |
| **ARACNe / ARACNe-AP** | **REJECT** | Mutual-information based; needs large sample counts for stable MI estimates (our tumour cohorts are n≈79–81). Java toolchain conflicts with the R-first requirement. No accuracy advantage over tree-based methods in benchmarks. |
| **GRNBoost2** | **SECONDARY (optional)** | Same algorithmic family as GENIE3, substantially faster. Python (`arboreto`). Reserve as a speed fallback only if GENIE3 runtime becomes limiting. |
| **GENIE3** | **ADOPT as the co-expression method** | Bioconductor R package — fits the R-first mandate. Random-forest based, handles continuous expression, no pseudotime requirement. Per BEELINE (Pratapa et al. 2019, Nat Methods, 692 citations), methods **not requiring pseudotime-ordered cells are generally more accurate.** |
| **ChIP-anchored regulon construction** (custom) | **ADOPT as primary** | Uses measured paralog binding rather than inferred motifs. The only approach that preserves paralog identity. |

### Final recommendation — a justified combination

1. **Primary: ChIP-anchored, enhancer-linked regulon construction.**
   Consensus ATAC regions → bedGraph signal quantification for MYC/MYCN/MYCL1 and H3K27ac → active-region calling → SE calling (stitching + signal ranking) → **peak-to-gene linking by distance-weighted assignment cross-checked against expression correlation across SCLC lines**, replacing Plotnik's nearest-gene assignment.

2. **Secondary: GENIE3** on bulk expression (CCLE SCLC lines + tumour cohorts), with MYC/MYCN/MYCL1 as candidate regulators, producing a **network-importance score** as an independent evidence domain.

3. ~~Optional tertiary: pySCENIC~~ — **removed 2026-07-26** with the single-cell layer. The GRN component is therefore entirely R-based: ChIP-anchored construction plus GENIE3.

**Framing caveat, per BEELINE:** GRN inference accuracy is *moderate* even for the best methods. GRN evidence therefore enters MOES as **one domain among four**, never as the arbiter.

---

## Part C — The Multi-Omics Evidence Score (MOES)

### Why not a weighted sum

Three reasons a weighted linear score is the wrong instrument here:

1. **Weights would be arbitrary.** The instruction was explicit that they must not be.
2. **The layers are not independent.** ATAC accessibility, H3K27ac activity, paralog binding and SE status are correlated *by construction* — a region is called active partly because it is open. A weighted sum silently multiplies this shared signal, inflating chromatin evidence roughly fourfold relative to dependency evidence.
3. **Coverage is ragged.** GeoMx measures ~1,800 genes; DepMap covers only expressed/targeted genes. A sum penalises genes for missing layers rather than treating absence as missing data.

### Recommended strategy — two-stage rank aggregation

Grounded in the cancer-gene-prioritisation literature, where **rank aggregation is the established weight-free integration strategy** (NetICS, Dimitrakopoulos et al. 2018, *Bioinformatics* — graph diffusion plus robust rank aggregation; related approaches MinNetRank and iRank use minimum-strategy and PageRank integration respectively).

**Stage 1 — collapse correlated layers into four conceptually independent evidence domains.**

| Domain | Constituent layers |
|---|---|
| **D1 Cis-regulatory** | Paralog ChIP signal at linked regions · H3K27ac activity · SE association · ATAC accessibility |
| **D2 Transcriptional** | DE in MYC-high vs MYC-low tumours (GSE60052) · replication in George 2015 · correlation with paralog expression across SCLC lines. *(The planned single-cell state-enrichment component is removed — see R-03.)* |
| **D3 Network** | GENIE3 importance · regulon membership · network centrality |
| **D4 Functional** | DepMap selective CRISPR dependency (paralog-amp vs non-amp) · drug-response association |

Within each domain: convert every layer to a within-layer normalised rank, aggregate by **Robust Rank Aggregation** (Kolde et al. 2012; R package `RobustRankAggreg`).

**Stage 2 — aggregate the four domain ranks by RRA again**, yielding a per-gene p-value under the null that a gene's ranks are drawn uniformly at random. Benjamini–Hochberg adjust across genes.

### Why RRA specifically
- **No weights required** — satisfies the stated constraint.
- Returns a **calibrated significance statement**, not an uninterpretable composite number.
- **Handles missing layers natively** — critical given GeoMx and DepMap coverage gaps.
- Robust to a single layer behaving badly.
- Two-stage grouping directly addresses the non-independence problem that defeats weighted sums.

### Mandatory rigour components

| Component | Purpose |
|---|---|
| **Layer–layer Spearman correlation matrix** (published as a figure) | Transparency about non-independence; justifies the domain grouping empirically rather than by assertion |
| **Leave-one-domain-out stability** | Recompute MOES dropping each domain; report rank correlation and top-N persistence. A hub that vanishes when one domain is removed is not robust |
| **Permutation null** (≥1,000 gene-label permutations within layers) | Empirical FDR |
| **Positive control** | Established MYC targets and the Dammert et al. 2019 MYC→BCL2/MCL1 axis should score highly for MYC specifically |
| **Negative control (pre-registered)** | Plotnik's paralog-**shared** housekeeping promoters (ribosome biogenesis, e.g. `RPS26`) must **NOT** surface as paralog-specific hubs. If they do, the specificity of the framework has failed and must be reported as such |
| **Benchmark comparison** | MOES-derived hubs vs `HALLMARK_MYC_TARGETS_V1/V2` and Jung et al. 2017's 18-gene signature — does SCLC-native, paralog-resolved evidence add anything over the paralog-blind status quo? |
| **Lineage-TF confounding analysis** | Using ASCL1/NEUROD1/POU2F3 ChIP references, test whether "paralog-specific" hubs are actually lineage-TF targets. **This is the most likely way the project's central claim fails, so it is tested directly, not defensively.** |
| **Sensitivity analysis** | Equal-weight sum reported alongside RRA — if conclusions diverge, both are shown |

### Reporting stance
MOES is a **heuristic prioritisation framework**. It is not predictive, not validated, and not clinical. Every output table and figure legend carries this statement. Language such as "predicts", "identifies drivers of", or "therapeutic target" is prohibited in the manuscript-style report.

---

## Part D — Discovery / validation / replication map

| Stage | Data | May be used for |
|---|---|---|
| **Discovery** | GSE230649 (ChIP + ATAC), CCLE/DepMap expression | Regulon construction, hub nomination |
| **Validation 1** | GSE60052 (79 tumours) | Testing regulon activity in patients |
| **Validation 2** | DepMap CRISPR, SCLC-CellMiner drug activity | Functional corroboration |
| **Replication** | George 2015 (cBioPortal), GSE269424 + GSE256345 + GSE281523 (ATAC) | Independent confirmation, different platform/population |
| **Orthogonal** | GSE261348 + GSE261345 (GeoMx) | Spatial coherence only |

**Rule:** no dataset changes role during the project. Any change is a protocol deviation and must be logged in the risk log with a date and reason.

---

## Part E — Pipeline stages (build order)

```
00_setup        environment, renv snapshot, directory scaffold
01_download     retrieval + checksums + manifest (records DepMap release ID)
02_qc_raw       file existence, dimensions, sample counts, genome-build assertions
03_regions      consensus ATAC region universe (hg19)
04_signal       bedGraph quantification: MYC/MYCN/MYCL1/H3K27ac over regions
05_active       active-region calling; comparison to Plotnik's reported counts
06_se           super-enhancer calling (stitch + rank)
07_link         peak-to-gene linking (distance-weighted + expression-correlated)
08_regulons     paralog regulons + shared/core regulon
09_bulk         tumour cohort scoring, DE, lineage-TF confounding analysis
10_grn          GENIE3 network importance
11_depmap       selective dependency analysis
12_drug         CellMiner / PRISM association
13_moes         two-stage RRA, stability, permutation, controls
14_spatial      restricted GeoMx coherence check
15_figures      publication figures
16_report       Quarto manuscript-style report
```

Each stage: one script, explicit inputs/outputs, hard-fail assertions at entry and exit, logged session info.

# Gap Statement & Novelty Assessment

**Status:** FINALISED — Phase 1b complete, 2026-07-26
**Verdict:** Original Gap A framing is **partially superseded** and has been **re-scoped**. The re-scoped gap survives.

---

## 1. The paper that changes the framing

**Plotnik JP et al. (2024)** *MYC Family Amplification Dictates Sensitivity to BET Bromodomain Protein Inhibitor Mivebresib (ABBV-075) in Small-Cell Lung Cancer.* Mol Cancer Res.
DOI: [10.1158/1541-7786.MCR-23-0599](https://doi.org/10.1158/1541-7786.MCR-23-0599) | PMID 38747975 | PMC11294817
This is the paper behind **GSE230649**, the dataset identified in Phase 1 as our keystone.

Retrieved and read in full via PubMed Central. This was the single most important finding of Phase 1b: **the first analytical step of the original Gap A had already been performed by the authors of the dataset we intended to use.**

### What Plotnik et al. already did (NO LONGER NOVEL — do not re-claim)

| Analysis | Their result |
|---|---|
| ChIP-seq of MYC, MYCN, MYCL1 in paralog-amplified SCLC lines (2 lines each) | First direct genome-wide comparison of all three paralogs in one tissue type |
| Intersection with H3K27ac to define "transcriptionally active" bound regions | 18,823 MYC / 4,017 MYCN / 5,688 MYCL1 active regions |
| Promoter vs. distal partition (1 kb cutoff) | MYC-amplified lines: 39% distal/enhancer vs. 12% in MYC-expressing lines |
| Motif analysis of enhancer-bound regions | Paralog-specific E-box central dinucleotides: MYC `CAGATG`, MYCN `CACATG`, MYCL `CACCTG` |
| Peak overlap between paralogs | 84% of MYCN regions fall inside a MYC peak; MYCL1 regions largely non-overlapping |
| ATAC-seq, 9 lines, BETi-sensitive vs. resistant | Differential accessibility concentrated in introns/enhancers (only 2.3% in promoters) |
| PCA of H3K27ac and ATAC landscapes | Lines separate along PC1 by MYC amplification status |
| Nearest-gene + g:Profiler ontology | MYC/MYCN enhancers → neurogenesis; MYCL1 + shared promoters → housekeeping/ribosome biogenesis |

**Consequence:** the claim *"MYC paralogs have distinct enhancer-binding preferences in SCLC"* is published. Our project must **build on it as established prior work, not rediscover it.**

### What Plotnik et al. explicitly did NOT do (the surviving gap)

Verified by full-text reading — none of the following appear anywhere in the paper:

1. **No super-enhancer calling.** Only a 1 kb promoter/distal binary split. No ROSE-style stitching, no SE ranking, no SE-associated gene sets.
2. **No peak-to-gene linking beyond "nearest gene."** No distance-weighted or expression-correlated enhancer→target assignment.
3. **No gene regulatory network inference of any kind.** No SCENIC, GENIE3, ARACNe, or network centrality.
4. **No patient tumour data whatsoever.** Every result is cell-line-derived. No bulk tumour cohort was scored.
5. **No single-cell analysis.**
6. **No spatial analysis.**
7. **No DepMap / CRISPR dependency integration.** CCLE is used only to look up copy-number for subgrouping.
8. **No external pharmacogenomic integration.** Drug response is their own in-house mivebresib/JQ1 IC50 assays only — no PRISM, no CTRP, no SCLC-CellMiner.
9. **No prioritisation or evidence-integration framework.** No ranking of candidate genes or hubs.
10. Their own stated limitation: analysis is entirely *in vitro*; they defer *in vivo* and downstream-pathway work to "future work."

---

## 2. Other near-neighbour literature checked

| Study | Overlap | Effect on our scope |
|---|---|---|
| **Dammert et al. 2019**, Nat Commun — MYC paralog-dependent apoptotic priming | Establishes paralog-specific *vulnerabilities* (MYC→BCL2 repression via MIZ1/DNMT3a; MCL1 dependency) via CRISPR-activation | Confirms paralog-specific functional divergence is real. Different axis (apoptotic machinery, not enhancer networks). **Use as orthogonal validation** — our top MYC-specific hubs should be consistent with their apoptotic-priming findings. |
| **Tlemsani et al. 2020**, Cell Rep — SCLC-CellMiner | States: "analyses reveal transcription networks linking SCLC subtypes with MYC and its paralogs" | Closest prior work on the *transcriptional network* side, but expression-correlation-based only — no chromatin, no binding, no enhancers. **Our regulatory anchoring is the differentiator.** Also our primary pharmacogenomic resource. |
| **Zhang et al. 2025**, Clin Cancer Res — ecDNA MYC paralog amplification shapes immunosuppressive TME | Uses IMC/mIHC spatial profiling; ecMYC+ tumours show reduced T-cell infiltration and spatial immune exclusion | **Materially reduces the novelty of the spatial-immune angle proposed in Phase 1 Gap B.** Spatial is therefore demoted from a discovery aim to a restricted orthogonal validation layer (see §4). |
| **Ireland et al. 2020**, Cancer Cell — MYC drives temporal evolution of SCLC subtypes | MYC drives ASCL1→NEUROD1→YAP1 plasticity | Foundational context. Motivates why paralog identity, not just "MYC-high," matters. Not overlapping. |
| **Jung et al. 2017**, Cancer Res — 18-gene Myc activity signature | Pan-cancer MYC activity signature with prognostic value | Precedent that MYC-activity signatures exist — but literature-curated, pan-cancer, **paralog-blind**, and not SCLC-native. Directly motivates our alternative. **Use as a benchmark comparator.** |
| **Miyakawa et al. 2022**, Cancer Sci — ASCL1 regulates SE-associated miRNAs | SE landscape in SCLC exists in the literature for the *lineage TF* axis | Confirms the SE-analysis gap is specific to the MYC-paralog axis, not SCLC generally. |

---

## 3. The re-scoped gap (FINAL)

> Plotnik et al. established that MYC, MYCN, and MYCL1 occupy **distinct enhancer repertoires** in SCLC cell lines. What remains entirely unaddressed is whether those paralog-specific regulatory programs **operate in patient tumours**, and **which of their downstream target hubs are robust enough to warrant functional investigation** when judged against independent chromatin, transcriptional, network, dependency, and pharmacogenomic evidence.
>
> No study has (i) converted paralog-resolved MYC binding into enhancer-anchored, super-enhancer-aware regulons, (ii) tested those regulons in independent human SCLC tumour cohorts, or (iii) applied a transparent, weight-free evidence-integration framework to prioritise MYC-associated regulatory hubs.

**Why this is defensible:**
- It does not re-claim anything Plotnik published; it explicitly cites them as the foundation.
- The bottleneck it addresses is real: cell-line enhancer maps are abundant, patient-tumour translation is absent.
- It is falsifiable. If paralog-specific regulons show no differential activity in tumour cohorts, that is a publishable negative result and we report it as such.
- The output (a ranked, evidence-audited hub list with per-layer provenance) is a genuine research product, not a gene list.

**Novelty grade: MODERATE-HIGH.** Downgraded from Phase 1's provisional "HIGH" because the regulatory-characterisation step is now prior work. The integration, translation, and prioritisation are novel.

---

## 4. Scope decisions forced by this review

| Decision | Rationale |
|---|---|
| **Cite Plotnik as foundation; re-derive their peaks rather than re-claim them** | Honest positioning. We use GSE230649 as input, not as a discovery. |
| **Add super-enhancer calling** | Genuine white space — they did a 1 kb binary split only. |
| **Add proper peak-to-gene linking** | "Nearest gene" is a known weak assignment; improving it is a real methodological contribution. |
| **Spatial demoted to restricted orthogonal validation** | Zhang 2025 covers MYC-paralog↔immune-exclusion spatially. We use spatial only to ask whether hub expression is spatially coherent and regionally heterogeneous — no prognostic claims. |
| **scRNA-seq layer is CONDITIONAL** | See risk log R-03. Best available patient-tumour atlas (Chan 2021) is not in GEO; the GEO alternative (GSE138474) is CDX/cell-line-derived on an old platform. |
| **Methylation, proteomics, miRNA EXCLUDED** | Available in SCLC-CellMiner but do not serve the hypothesis. Included only for completeness would violate the biological-coherence principle. |

---

## 5. If the gap collapses later

Trigger conditions and pre-agreed responses:

- **A paper appears doing regulon-based MYC-paralog scoring in SCLC tumours** → pivot to the methodological contribution (the MOES framework itself, benchmarked against Jung 2017's signature) and reframe as a reproducibility/benchmarking study.
- **Paralog-specific regulons show no separation in tumour cohorts** → report as a negative result with power analysis; the framework and the audit trail remain the contribution.
- **GSE230649 bedGraph proves unusable** → fall back to the independent ATAC sets (GSE269424, GSE256345) plus lineage-TF ChIP, and restrict the binding layer to what can be recovered.

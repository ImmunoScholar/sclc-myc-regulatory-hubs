# Figure S1 — Data landscape

**This figure reports metadata, not results.** No analytical output exists at M4.

**A.** Composition of the keystone deposit GSE230649 (hg19, 28 samples): MYC ChIP in 5 lines, MYCN in 2, MYCL1 in 2, H3K27ac in 10, ATAC in 9. Two gaps are carried forward rather than smoothed over — H211 has MYC ChIP but no ATAC, and H196 has accessibility and H3K27ac but no MYC-family ChIP. The strip beneath encodes author-declared paralog amplification status.

**B.** Files by genome build. The project build is hg19, fixed by the keystone, but all three supporting ATAC deposits and the ASCL1 ChIP are hg38 — so the consensus accessible-region universe, which requires support from at least two independent ATAC datasets, cannot be built without crossing builds. Continuous signal is never lifted; intervals are called in the native build and only intervals are lifted, with loss rate reported.

**C.** Cell-line overlap between the keystone and the lineage-TF controls. This is the constraint with the clearest scientific consequence: the confounding test central to risk R-01 requires lineage-TF occupancy in the same lines where MYC-family occupancy was measured. POU2F3 provides three such lines across two paralog groups and is already in hg19; ASCL1 provides one; NEUROD1 provides none. Conclusions about lineage-independence are therefore qualified per transcription factor rather than stated uniformly.

Sources: GEO series and sample records for GSE230649, GSE269424, GSE256345,
GSE281523, GSE281524, GSE210113, GSE249362 (verified 2026-07-26); amplification
status from Plotnik et al. 2024, Mol Cancer Res (PMID 38747975).

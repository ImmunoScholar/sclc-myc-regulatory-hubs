# Figure S3 — The lineage confound at every level tested

**A.** Chromatin-level co-occupancy. In the three keystone lines with POU2F3
ChIP, regions active for MYC or MYCN are POU2F3-bound at 2.5–3.1x the rate of
the ATAC universe they were drawn from (all p < 1e-15). The entanglement
between MYC-family occupancy and lineage-factor occupancy is therefore present
in the very data the regulons were built from, before any tumour is scored.
Resolution is not uniform across lineage TFs and was never pooled: POU2F3 has
three keystone lines, ASCL1 one, NEUROD1 none (risk R-14).

**B.** Cell-line level. Across 59 CCLE SCLC lines, MYC expression correlates
negatively with neuroendocrine score (rho -0.441, FDR 0.008) and with ASCL1
(rho -0.409, FDR 0.010). MYCN and MYCL show no comparable structure. This is
the independent reproduction of Ireland et al. 2020 in a third dataset, and it
is also one of the two grounds on which the network domain was excluded: a
cell-line network would inherit exactly this confound (D-037).

**C.** Gene level, genome-wide. Of 33,682 genes tested per paralog, the number
still associated with the paralog after adjusting for NE score and all four
lineage TFs is 206 for MYC, 62 for MYCL1, and zero for MYCN. Crucially, the
survivors are not concentrated in the regulons: enrichment inside versus
outside is OR 1.21 (p = 0.57) for MYC and OR 1.57 (p = 0.48) for MYCL1 —
background rate. Survivors exist, but they are not the regulon genes, so they
do not rescue the paralog-specific claim.

Taken together these three levels rule out the most obvious objection to
Figure 1 — that the null is an artefact of tumour-level regulon scoring.

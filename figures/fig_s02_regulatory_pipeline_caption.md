# Figure S2 — Constructing the regulatory layer

**A.** Intrinsic validation of the ATAC threshold. TSS enrichment and H3K27ac
fold both rise monotonically with the threshold, which is the evidence that the
retained regions are regulatory signal rather than background. Because the
relationship holds across the whole range rather than at one point, the chosen
3x cutoff is not load-bearing (D-024): conclusions do not depend on it.

**B.** The MYCN-in-MYC nesting fraction behind gate criterion 1, across four
thresholds and three replicate conventions. Within the pre-specified
≥2-replicate rule the value varies by only 0.037 across a 2.7-fold change in
threshold. The union and all-replicate conventions sit at materially different
values, so the convention is a real choice and is shown rather than hidden.
The circled point is the value reported in the gate.

**C.** Super-enhancer calling. The first stitching rule — first slope ≥ 1 —
returned 8,238 of 8,238 stitched enhancers as super-enhancers for H1048, i.e.
all of them, which is not a super-enhancer call at all. It was replaced by a
knee-point cutoff with a hard guard at 30% of stitched enhancers. Every line
now lands between 1.9% and 6.6%, far below the guard.

**D.** Peak-to-gene linking. Distance-weighted, activity-correlated linking
retains 48,756 of 127,019 candidate pairs. Only 32.6% of retained links agree with the
nearest gene, so this is not nearest-gene assignment under another name, and the
link set is enriched 6.1x over a distance-matched null.
An earlier version of this step assigned distal regions only, which left
promoter-bound genes out of every regulon and understated the enrichment
(D-029).

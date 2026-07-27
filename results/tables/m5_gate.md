# M5 Gate

Generated: 2026-07-27 11:25:14 UTC

No criterion is a region count: counts are determined by our own threshold
and universe size, so matching a published count would validate nothing (D-020).

| criterion | observed | expected | result |
|---|---|---|---|
| 1_mycn_in_myc | 0.886 (spread 0.012) | 0.84 +/- 0.15, stable | PASS |
| 2_differential_nesting | MYCN 5.12x vs MYCL1 4.48x, OR 2.27, p=3e-41 | enrichment(MYCN) > enrichment(MYCL1), p < 0.05 | PASS |
| 3_distal_contrast | not evaluable | MYC-amplified > MYC-expressing by >=0.10 | **BLOCKED** |
| 4_motif_specificity | 3/3 own-motif top share (p=0.037); chi-sq p=7.7e-100 | 3/3 own-motif highest share, motif x set dependence p<0.05 | PASS |

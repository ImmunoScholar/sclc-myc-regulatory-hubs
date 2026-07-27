# M5 Gate

Generated: 2026-07-27 13:30:01 UTC

No criterion is a region count: counts are determined by our own threshold
and universe size, so matching a published count would validate nothing (D-020).

| criterion | observed | expected | result |
|---|---|---|---|
| 1_mycn_in_myc | 0.912 (spread 0.037) | 0.84 +/- 0.15, stable | PASS |
| 2_differential_nesting | MYCN 5.50x vs MYCL1 4.63x, OR 3.14, p=1.1e-60 | enrichment(MYCN) > enrichment(MYCL1), p < 0.05 | PASS |
| 3_distal_contrast | not evaluable | MYC-amplified > MYC-expressing by >=0.10 | **BLOCKED** |
| 4_motif_specificity | 3/3 own-motif top share (p=0.037); chi-sq p=1.6e-101 | 3/3 own-motif highest share, motif x set dependence p<0.05 | PASS |

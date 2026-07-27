# M5 — accessible-region universe

Generated: 2026-07-27 08:43:20 UTC

Derived from the nine keystone ATAC tracks (GSE230649, hg19, cells matched to the ChIP). Primary threshold **x2.5** each line's mean signal over covered bases (D-024).

Four universes are retained so every M5 gate metric can be reported across
thresholds. The conclusion should not depend on the choice; if it does,
that fragility is the finding.

| x mean | regions | % genome | distal % | multi-line % | ext. corrob % |
|---|---|---|---|---|---|
| 2 | 897,764 | 12.619 | 96.6 | 34.8 | 28.2 |
| 2.5 **(primary)** | 478,619 | 6.995 | 94.7 | 33.8 | 42.4 |
| 3 | 331,613 | 5.014 | 93.3 | 35 | 50.8 |
| 4 | 174,430 | 3.13 | 90.3 | 36 | 62.7 |

Validated intrinsically: H524/chr1 TSS enrichment 8.9x and matched-cell
H3K27ac fold 2.83x at the primary threshold (D-024). Agreement with an
external laboratory's ATAC peak calls was abandoned as a target — it
measured inter-laboratory reproducibility, not fitness for purpose.

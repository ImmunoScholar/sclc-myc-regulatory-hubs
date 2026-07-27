# M5 — MYC-family amplification inferred from coverage

Generated: 2026-07-27 11:43:31 UTC

Plotnik profiled "two cell lines harboring the alteration for each
amplification type" but GSE230649 has five MYC ChIP samples, and neither
the paper nor the GEO records name which are amplified. Measured here
instead, as amplicon coverage relative to each track's genome-wide median.

| line | MYC | MYCN | MYCL1 | strongest | ratio |
|---|---|---|---|---|---|
| COLO668 | 1.31 | 0.61 | 12.21 | MYCL1 | 12.21 |
| H1048 | 1.12 | 0.78 | 1.11 | MYC | 1.12 |
| H196 | 1.66 | 0.71 | 1.42 | MYC | 1.66 |
| H211 | 4.24 | 0.71 | 1.31 | MYC | 4.24 |
| H524 | 1.21 | 0.68 | 1.47 | MYCL1 | 1.47 |
| H526 | 0.94 | 6.04 | 1.81 | MYCN | 6.04 |
| H69 | 0.86 | 0.52 | 0.71 | MYC | 0.86 |
| H847 | 0.8 | 1.04 | 1.48 | MYCL1 | 1.48 |
| H889 | 0.56 | 0.65 | 3.73 | MYCL1 | 3.73 |
| SHP77 | 1.81 | 1.1 | 1.01 | MYC | 1.81 |

Coverage reflects copy number and accessibility together, so a call is
treated as solid only where ATAC and H3K27ac agree and the elevation holds
across the broad (+/-1 Mb) window.

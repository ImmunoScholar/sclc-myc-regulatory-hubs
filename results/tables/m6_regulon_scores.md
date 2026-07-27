# M6 — regulon scores in patient tumours (GSE60052, n=79)

Generated: 2026-07-27 13:54:59 UTC

Primary scorer singscore (rank-based, cohort-independent); GSVA as the
methodologically distinct sensitivity check.

| paralog | genes used | rho with own expression | 95% CI | FDR | programme |
|---|---|---|---|---|---|
| MYC | 368/500 | -0.202 | -0.405 to 0.02 | 0.111 | validated |
| MYCN | 398/500 | 0.135 | -0.089 to 0.346 | 0.236 | validated |
| MYCL1 | 373/500 | 0.268 | 0.05 to 0.462 | 0.0507 | **unvalidated** |

## Specificity matrix (Spearman rho)

| regulon | MYC | MYCN | MYCL1 |
|---|---|---|---|
| MYC | -0.202 | 0.175 | 0.242 |
| MYCN | -0.416 | 0.135 | 0.118 |
| MYCL1 | -0.468 | 0.075 | 0.268 |

Diagonal dominance: 2/3 regulons correlate most strongly with their own paralog (P = 0.259 by chance).

A null result here is a reportable finding, not a failure: it would mean
paralog-resolved programmes do not retain paralog identity in patient
tumours (gap statement section 5).

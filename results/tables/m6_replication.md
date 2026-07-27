# M6 — replication in an independent cohort

Generated: 2026-07-27 15:28:56 UTC

GSE60052 (n=79) vs George et al. 2015 (n=81, cBioPortal). Independent
patients, sequencing and processing.

**Scoring**: mean per-gene z-score in both cohorts. singscore was not used
for the comparison because its ranking background differed between cohorts
(33,683 vs 1,004 genes) and would have compared scoring artefacts.

**Gene naming**: GSE60052 carries `MYCL1`, George carries `MYCL`. Resolved
explicitly.

| cohort | paralog | genes | rho raw | unique paralog R2 | unique lineage R2 | rho partial | FDR |
|---|---|---|---|---|---|---|---|
| GSE60052 | MYC | 368 | 0.07 | 0.001 | **0.428** | 0.022 | 0.846 |
| GSE60052 | MYCN | 398 | 0.297 | 0.016 | **0.377** | 0.204 | 0.107 |
| GSE60052 | MYCL1 | 373 | 0.39 | 0.068 | **0.355** | 0.317 | 0.0136 |
| George2015 | MYC | 379 | -0.1 | 0.019 | **0.269** | 0.038 | 0.967 |
| George2015 | MYCN | 418 | 0.099 | 0.001 | **0.244** | 0.005 | 0.967 |
| George2015 | MYCL1 | 382 | 0.186 | 0 | **0.426** | 0.114 | 0.936 |

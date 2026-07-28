# M9 — spatial coherence and orthogonal validation

Generated: 2026-07-28 08:49:48 UTC

## Scope and what limits it

Restricted orthogonal validation only (datasets.yml D5). Two GeoMx DSP cohorts on a targeted CTA panel: GSE261348 (IMfirst) and GSE261345 (CANTABRICO).

**Only the MYC regulon can be scored.** Panel coverage against the gate
(>= 20 genes AND >= 10% of members):

| paralog | measured / size | coverage | scoreable |
|---|---|---|---|
| MYC | 56 / 500 | 11.2% | **yes** |
| MYCN | 42 / 500 | 8.4% | no |
| MYCL1 | 48 / 500 | 9.6% | no |

MYCL1 fails by 0.4 percentage points and MYC clears by 1.2. These are knife-edge verdicts against a threshold that was itself raised during the project, and should be read as such rather than as clean separations. The MYC result below is a **56-gene proxy**, not the 500-gene regulon.

Two further limits apply to every number here:

1. **Slide is not tumour.** Most slides carry two patient identifiers with no per-ROI patient label, so within-slide variance is an upper bound on within-tumour variance. Only 3 IMfirst slides are unambiguously one patient; CANTABRICO has none.
2. **The NE score rests on 3 of 10 markers** (NCAM1, ASCL1, DLL3); the rest are off-panel.

## Coherence — is the score a property of the tumour or the region?

| cohort | subset | ROIs | slides | between-slide R² | within-slide R² |
|---|---|---|---|---|---|
| GSE261348 | all slides | 174 | 18 | 0.554 | 0.446 |
| GSE261348 | single-patient slides | 21 | 3 | 0.701 | 0.299 |
| GSE261345 | all slides | 121 | 14 | 0.528 | 0.471 |
| GSE261345 | single-patient slides | 0 | 0 | — | — |

## Replication — lineage versus the paralog's own expression

| cohort | ROIs | unique R² MYC | unique R² lineage | full R² | lineage dominates |
|---|---|---|---|---|---|
| GSE261348 | 174 | 0.0189 | 0.4744 | 0.4744 | **yes** |
| GSE261345 | 121 | 0.0097 | 0.5076 | 0.5704 | **yes** |

## Associations

| cohort | variable | Spearman rho | p |
|---|---|---|---|
| GSE261348 | MYC | +0.113 | 0.137 |
| GSE261348 | NE | +0.169 | 0.0259 |
| GSE261348 | ASCL1 | +0.084 | 0.268 |
| GSE261348 | NEUROD1 | +0.520 | 1.99e-13 |
| GSE261348 | POU2F3 | +0.274 | 0.000256 |
| GSE261348 | YAP1 | -0.175 | 0.0209 |
| GSE261348 | MYC (lineage-adjusted) | -0.022 | 0.771 |
| GSE261345 | MYC | +0.142 | 0.12 |
| GSE261345 | NE | -0.312 | 0.000492 |
| GSE261345 | ASCL1 | -0.435 | 6e-07 |
| GSE261345 | NEUROD1 | +0.647 | 1.1e-15 |
| GSE261345 | POU2F3 | +0.327 | 0.000254 |
| GSE261345 | YAP1 | +0.195 | 0.0324 |
| GSE261345 | MYC (lineage-adjusted) | -0.172 | 0.0594 |

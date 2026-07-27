# M7 — network domain feasibility

Generated: 2026-07-27 18:45:01 UTC
DepMap Public 26Q1 expression, 59 SCLC cell lines.

Tested BEFORE building, so a confounded domain cannot be admitted while
appearing independent.

- **Power**: n = 59 SCLC lines. GENIE3 conventionally needs hundreds. **Inadequate.**
- **Confounding**: MYC vs NE score rho = -0.441 (tumours gave -0.590 and -0.466). **Confounded.**

**Verdict: EXCLUDE** the network domain.

| paralog | target | rho | FDR |
|---|---|---|---|
| MYC | NE_score | -0.441 | 0.00813 |
| MYC | ASCL1 | -0.409 | 0.00982 |
| MYC | NEUROD1 | 0.096 | 0.584 |
| MYC | POU2F3 | 0.078 | 0.595 |
| MYC | YAP1 | 0.001 | 0.996 |
| MYCN | NE_score | -0.184 | 0.306 |
| MYCN | ASCL1 | -0.199 | 0.306 |
| MYCN | NEUROD1 | 0.187 | 0.306 |
| MYCN | POU2F3 | 0.261 | 0.136 |
| MYCN | YAP1 | 0.083 | 0.595 |
| MYCL | NE_score | 0.278 | 0.136 |
| MYCL | ASCL1 | 0.266 | 0.136 |
| MYCL | NEUROD1 | 0.12 | 0.549 |
| MYCL | POU2F3 | -0.102 | 0.584 |
| MYCL | YAP1 | -0.159 | 0.378 |

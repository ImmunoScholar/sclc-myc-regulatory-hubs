# Aim 3 — drug-response association (BET inhibitors)

Generated: 2026-07-28 11:09:56 UTC

The last open commitment from the project contract. Deliberately narrow and
on-thesis: the source study reports that MYC-family amplification dictates
sensitivity to BET bromodomain inhibition in SCLC, so the question asked of
public pharmacogenomic data is whether MYC-family status tracks BET-inhibitor
sensitivity, and whether it survives neuroendocrine lineage adjustment.

Three independent screens, so no single-screen artefact can carry the result.

## Coverage

| screen | lines | SCLC (exact) | SCLC (grepl) | compounds | BET |
|---|---|---|---|---|---|
| PRISM | 727 | 22 | 29 | 1482 | 5 |
| GDSC1 | 948 | 51 | 63 | 316 | 2 |
| GDSC2 | 947 | 49 | 61 | 286 | 3 |

The grepl column is shown to make the trap visible: "Non-Small Cell Lung Cancer" contains "Small Cell Lung Cancer", so pattern-matching the lineage silently adds NSCLC lines to an SCLC set. Matching is exact.

## Associations

AUC is a sensitivity metric on which **lower means more sensitive**, so a
positive rho means higher expression tracks *less* sensitivity.

| screen | drug | gene | n | rho | p | rho (NE-adj) | p (NE-adj) |
|---|---|---|---|---|---|---|---|
| PRISM | MOLIBRESIB | MYC | 21 | -0.055 | 0.811 | +0.029 | 0.902 |
| PRISM | MOLIBRESIB | MYCN | 21 | +0.306 | 0.178 | +0.266 | 0.243 |
| PRISM | MOLIBRESIB | MYCL | 21 | +0.046 | 0.843 | -0.029 | 0.902 |
| GDSC1 | MOLIBRESIB (GDSC1:275) | MYC | 31 | +0.119 | 0.522 | +0.169 | 0.364 |
| GDSC1 | MOLIBRESIB (GDSC1:275) | MYCN | 31 | -0.351 | 0.0529 | -0.307 | 0.0927 |
| GDSC1 | MOLIBRESIB (GDSC1:275) | MYCL | 31 | -0.007 | 0.969 | -0.075 | 0.69 |
| GDSC1 | PFI-1 (GDSC1:1219) | MYC | 29 | -0.279 | 0.142 | -0.183 | 0.341 |
| GDSC1 | PFI-1 (GDSC1:1219) | MYCN | 29 | +0.119 | 0.538 | +0.091 | 0.638 |
| GDSC1 | PFI-1 (GDSC1:1219) | MYCL | 29 | +0.426 | 0.0213 | +0.429 | 0.0204 |
| GDSC2 | BIRABRESIB (GDSC2:1626) | MYC | 30 | -0.390 | 0.033 | -0.314 | 0.0912 |
| GDSC2 | BIRABRESIB (GDSC2:1626) | MYCN | 30 | -0.189 | 0.318 | -0.160 | 0.397 |
| GDSC2 | BIRABRESIB (GDSC2:1626) | MYCL | 30 | +0.621 | 0.000253 | +0.631 | 0.000187 |
| GDSC2 | MOLIBRESIB (GDSC2:1624) | MYC | 30 | -0.326 | 0.0784 | -0.231 | 0.219 |
| GDSC2 | MOLIBRESIB (GDSC2:1624) | MYCN | 30 | -0.115 | 0.545 | -0.111 | 0.561 |
| GDSC2 | MOLIBRESIB (GDSC2:1624) | MYCL | 30 | +0.617 | 0.000281 | +0.587 | 0.000648 |
| GDSC2 | PFI-1 (GDSC2:2173) | MYC | 34 | -0.026 | 0.882 | +0.004 | 0.982 |
| GDSC2 | PFI-1 (GDSC2:2173) | MYCN | 34 | -0.211 | 0.231 | -0.130 | 0.465 |
| GDSC2 | PFI-1 (GDSC2:2173) | MYCL | 34 | +0.221 | 0.21 | +0.178 | 0.313 |

Benjamini-Hochberg across all 18 tests: 2 survive at FDR < 0.05 (minimum FDR 0.003).

## Cross-screen replication — read this before the table above

The MYCL association is the **only** signal anywhere in this project that survives neuroendocrine lineage adjustment. That makes it the one most in need of scepticism, not the least.

It appears in **2 of 3 screens** (GDSC1, GDSC2) and **not in PRISM** (rho +0.046, p = 0.84, n = 21).

**GDSC1 and GDSC2 are not independent.** They are the same Sanger platform over heavily overlapping cell lines, so agreement between them is much weaker evidence than agreement between GDSC and PRISM would be. Counting them as two replications would overstate the result.

The direction is also worth stating plainly: positive rho means MYCL-high lines are **less** sensitive to BET inhibition. The source study's claim concerns MYC-family *amplification* driving *sensitivity*; MYC here trends the opposite way to MYCL (negative rho, i.e. more sensitive) but does **not** survive lineage adjustment. Nothing here reproduces or refutes that study, which tested amplification rather than expression and a different compound.

Reported as **exploratory and provisional**.

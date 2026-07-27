# M6 — lineage confounding (PRIMARY analysis, risk R-01)

Generated: 2026-07-27 13:58:04 UTC
GSE60052, n = 79 tumours.

## Subtype distribution

`SCLC-A=26  SCLC-N=18  SCLC-P=13  SCLC-Y=22`

## Variance partitioning

| regulon | R2 paralog only | R2 lineage only | R2 both | unique paralog | unique lineage |
|---|---|---|---|---|---|
| MYC | 0.048 | 0.386 | 0.386 | 0 | **0.338** |
| MYCN | 0.012 | 0.424 | 0.426 | 0.001 | **0.414** |
| MYCL1 | 0.098 | 0.49 | 0.506 | 0.016 | **0.408** |

## Paralog association before and after adjusting for lineage

| regulon | rho raw | rho partial | FDR |
|---|---|---|---|
| MYC | -0.202 | 0.029 | 0.796 |
| MYCN | 0.135 | 0.035 | 0.796 |
| MYCL1 | 0.268 | 0.138 | 0.675 |

## Per-TF resolution (R-14)

Occupancy-level confounding cannot be tested uniformly and no blanket
"we controlled for lineage TFs" claim is permitted:

- **POU2F3** — within-line in 3 keystone lines, hg19. Strongest arm.
- **ASCL1** — within-line in SHP-77 only (n=1), hg38. Descriptive.
- **NEUROD1** — H446 only, zero keystone overlap. Subtype-level only.
- **YAP1** — no ChIP data in any acquired dataset. Expression only.

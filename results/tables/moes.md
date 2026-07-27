# Reduced two-domain MOES

Generated: 2026-07-27 20:04:35 UTC

## What this is, and is not

MOES was specified at M1 as a **four-domain, two-stage, paralog-specific** prioritisation. It runs here on **two** domains, of which **one** attributes evidence to a paralog.

| domain | status | attributes to a paralog? |
|---|---|---|
| cis-regulatory | admitted | yes |
| functional | admitted (D-036) | no — gene level only |
| transcriptional | excluded (D-033) | — |
| network | excluded (D-037) | — |

**Two-stage RRA is degenerate here.** With two domains, stage 2 has nothing to aggregate that stage 1 did not; the reported score is a single aggregation, not a hierarchy. Leave-one-domain-out likewise reduces to single-domain rankings and is reported as such.

## Layer correlation

| paralog | Spearman(cis, functional) |
|---|---|
| MYC | -0.0229 |
| MYCN | +0.0119 |
| MYCL1 | -0.0075 |

Near zero, so the two domains are close to independent and RRA is not double-counting one signal.

## Leave-one-domain-out (degenerate)

| paralog | combined vs cis only | combined vs functional only |
|---|---|---|
| MYC | +0.621 | +0.654 |
| MYCN | +0.635 | +0.692 |
| MYCL1 | +0.699 | +0.572 |

## Between-paralog difference is chromatin alone

| pair | MOES rho | cis-only rho |
|---|---|---|
| MYC vs MYCN | +0.576 | +0.169 |
| MYC vs MYCL1 | +0.497 | +0.174 |
| MYCN vs MYCL1 | +0.590 | +0.330 |

The functional vector is identical for all three paralogs, so every difference between their MOES rankings originates in the cis layer. MOES cannot deliver multi-layer paralog-specific prioritisation on this evidence, and this table is why.

## Result

Universe: 10,387 genes with evidence in both domains. Permutations: 10,000. Bootstrap: 2,000 resamples of cis links.

| paralog | genes at FDR < 0.05 | best FDR achieved | rank-stable |
|---|---|---|---|
| MYC | 0 | 0.381 | 0 |
| MYCN | 0 | 0.371 | 0 |
| MYCL1 | 0 | 0.358 | 0 |

**No gene reaches FDR < 0.05 for any paralog.** The best FDR achieved is 0.358. This is not a marginal miss: the two admitted domains are close to independent (layer rho between -0.023 and +0.012), so genes ranking highly in one do not rank highly in the other more often than chance predicts. RRA has nothing to aggregate.

MOES therefore returns **no prioritised hub list**. Reporting a top-N table here would present the head of an unranked distribution as a result. The ranking is written to `data/metadata/moes_ranking.csv` with its FDR column so the absence is inspectable, not hidden.

See the sensitivity section below before reading this null as uninformative.

**Bootstrap caveat.** Intervals resample cis-regulatory links only; the functional domain is held fixed because its per-gene values are precomputed summaries over cell lines. Reported rank intervals therefore UNDERSTATE total uncertainty.

## Method sensitivity (positive control)

The result above is a null, so the method was tested for its ability to detect the alternative. A synthetic functional layer was blended with the real cis ranking at increasing strength and passed through the identical FDR machinery.

| injected signal | layer rho | min FDR | genes FDR < 0.05 | genes FDR < 0.20 |
|---|---|---|---|---|
| 0.00 | -0.000 | 0.878 | 0 | 0 |
| 0.05 | +0.028 | 0.721 | 0 | 0 |
| 0.10 | +0.101 | 0.084 | 0 | 8 |
| 0.20 | +0.209 | 0.018 | 1 | 16 |
| 0.40 | +0.508 | 0.004 | 6 | 32 |
| 1.00 | +1.000 | 0.000 | 288 | 1095 |

The response is graded and monotone: nothing on pure noise, discoveries appearing once the two layers correlate at rho ~ 0.21, and 288 genes when the layers are identical.

**The observed layer correlation is -0.023 and the best FDR achieved on real data is 0.381.** The convergence MOES was built to find is not merely non-significant here; it is absent at a level the method demonstrably detects. The null is a property of the evidence, not of the test.

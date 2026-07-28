# Reduced two-domain MOES

Generated: 2026-07-27 20:28:58 UTC

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
| MYC | 0 | 0.383 | 0 |
| MYCN | 0 | 0.375 | 0 |
| MYCL1 | 0 | 0.358 | 0 |

**No gene reaches FDR < 0.05 for any paralog.** The best FDR achieved is 0.358.

## Aggregate concordance — weak, real, and not localisable

The per-gene null above should NOT be read as 'the two domains share nothing'. Asked as an aggregate question — do the top K genes of each domain overlap more than the K²/N expected by chance? — the answer is yes for two of three paralogs:

| paralog | max overlap ratio | global p |
|---|---|---|
| MYC | 2.03 | 0.019 |
| MYCN | 2.19 | 0.008 |
| MYCL1 | 1.72 | 0.065 |

The global p compares the maximum ratio over K against the permutation distribution of that same maximum, so it controls family-wise error across the whole curve. Pointwise exceedances are recorded in `moes_concordance.csv` but are not a test: the K values are nested, so a single excursion produces a run of them. K is restricted to where at least 5 genes are expected to overlap by chance.

This sits alongside near-zero Spearman correlation between the domains (-0.023 to +0.012) without contradiction. Spearman measures monotone association across all 10,387 genes; concordance-at-K measures agreement specifically at the TOP of both lists. The evidence is that the two domains agree weakly about which genes are most important, and not at all about the ordering of the rest.

## What MOES can and cannot deliver

MOES returns **no prioritised hub list**, but the reason is more specific than an absence of signal. The convergence is real at roughly 2x and too DIFFUSE to attribute: spread thinly across the top of both rankings, no individual gene carries enough evidence to survive multiple testing over 10,387 genes. Aggregate detectable, per-gene unattributable.

Reporting a top-N table would therefore be wrong twice over: it would present genes whose individual FDR is ~0.36 as prioritised hits, and it would imply a paralog attribution that panel D of Figure 4 shows comes from chromatin alone. The full ranking with its FDR column is in `data/metadata/moes_ranking.csv` so this is inspectable, not hidden.

See the sensitivity section below before reading the per-gene null as uninformative.

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

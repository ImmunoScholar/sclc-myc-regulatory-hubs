# M5 — signal normalisation

Generated: 2026-07-27 10:04:10 UTC

**Adopted: fold-over-own-background.**

Paralog groups are perfectly confounded with cell lines (MYCN is only
H526+H69; MYCL1 only COLO668+H889), so the normalisation choice is
load-bearing for the headline comparison and is recorded explicitly.

| method | MYC regions | MYCN | MYCL1 | MYCN-in-MYC | MYCL1-in-MYC |
|---|---|---|---|---|---|
| raw | 4,977 | 4,575 | 7,956 | **0.51** | 0.395 |
| fold_over_background **(adopted)** | 4,977 | 4,575 | 7,956 | **0.51** | 0.395 |
| quantile_within_assay | 5,269 | 4,768 | 8,263 | **0.513** | 0.41 |
| median_scaled | 4,977 | 4,575 | 7,956 | **0.51** | 0.395 |

Quantile normalisation is shown for comparison and is **not** used.
It forces identical signal distributions, so if MYCN genuinely binds
fewer regions than MYC it erases that difference — it would reproduce the
published 0.84 whether or not the data support it. Agreement produced that
way is an artefact of the method, not evidence.

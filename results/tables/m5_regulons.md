# M5 — paralog regulons

Generated: 2026-07-27 13:49:02 UTC

Top-ranked genes by aggregate link score, capped at 500. Links are high+moderate tier only, from an H3K27ac activity proxy (D-027) — not expression.

| paralog | active regions | regulon genes |
|---|---|---|
| MYC | 16,115 | 500 |
| MYCN | 2,600 | 500 |
| MYCL1 | 5,536 | 500 |

## Validity gate

| check | detail | result |
|---|---|---|
| leave_one_line_out | 8/9 held-out lines pass | PASS |
| cross_paralog_distinctness | max Jaccard 0.130 (limit 0.60) | PASS |
| coherent_programme | MYC:NEUROGENESIS MYCN:HALLMARK MYCL1:NONE | **FAIL** |

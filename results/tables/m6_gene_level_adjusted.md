# M6 — per-gene paralog association, lineage-adjusted

Generated: 2026-07-27 14:02:22 UTC
GSE60052, n = 79 tumours. Adjusted for NE score + ASCL1, NEUROD1, POU2F3, YAP1 (df = 72).

The regulon-level test is null (R-01). MOES operates per gene, so this asks
whether individual genes carry lineage-independent paralog signal. The test
is ENRICHMENT relative to non-regulon genes, not a raw survivor count —
with ~33,700 genes some survive FDR by chance.

| paralog | regulon genes | survivors (FDR<0.05, +) | in regulon | rate in | rate out | OR | p |
|---|---|---|---|---|---|---|---|
| MYC | 367 | 76 | 1 | 0.27% | 0.23% | 1.21 | 0.566 |
| MYCN | 397 | 0 | 0 | 0% | 0% | NA | NA |
| MYCL1 | 373 | 58 | 1 | 0.27% | 0.17% | 1.57 | 0.476 |

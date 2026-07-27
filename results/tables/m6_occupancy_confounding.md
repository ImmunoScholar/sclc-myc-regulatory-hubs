# M6 — lineage confounding at the chromatin level (R-14)

Generated: 2026-07-27 15:19:44 UTC

Second, independent test of R-01. M6 found the confound in tumour
EXPRESSION; this asks whether paralog-bound regions are also
lineage-TF-bound IN THE SAME CELLS.

## POU2F3 within-line co-occupancy (the only strong arm)

| line | paralog | active regions | universe POU2F3 | active POU2F3 | enrichment | OR | p |
|---|---|---|---|---|---|---|---|
| H1048 | MYC | 12,403 | 9.3% | 23.8% | 2.57x | 4.06 | 0 |
| H526 | MYCN | 9,924 | 28.6% | 72.3% | 2.52x | 8.4 | 0 |
| H211 | MYC | 10,643 | 22.9% | 71% | 3.1x | 11.97 | 0 |

## Resolution per TF (R-14) — not uniform, never pooled

- **POU2F3** — 3 keystone lines (H1048, H211 = MYC; H526 = MYCN), hg19, peaks. Strong.
- **NEUROD1** — H446 only, zero keystone overlap. Subtype-level; cannot separate line from subtype.
- **ASCL1** — SHP-77 bigWig signal, hg38, n=1. NOT tested; would not support inference.

Drug-treated arms excluded: the FHD286 POU2F3 peak file is 4,606 bytes
against 181,493 for DMSO, so treated arms do not represent native binding.

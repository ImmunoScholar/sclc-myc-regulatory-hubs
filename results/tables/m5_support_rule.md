# M5 — accessible-region support rule

Generated: 2026-07-27 08:06:51 UTC

Each candidate rule scored against H524's independent MACS2 peak set on
chr1, restricted to regions our own H524 ATAC track supports.

| rule | regions | % genome | ext. corroborated | precision | recall | F1 | Jaccard |
|---|---|---|---|---|---|---|---|
| all (>=1 line) | 297,443 | 4.397 | 53.9% | 0.3705 | 0.3505 | **0.3603** | 0.2197 |
| >=2 lines | 102,334 | 2.474 | 80.5% | 0.3982 | 0.3489 | **0.3719** | 0.2285 |
| >=3 lines | 55,528 | 1.703 | 91% | 0.429 | 0.343 | **0.3812** | 0.2355 |
| >=2 lines OR ext-corroborated | 180,223 | 3.286 | 88.9% | 0.3922 | 0.3505 | **0.3702** | 0.2272 |
| >=1 line AND ext-corroborated | 160,303 | 2.933 | 100% | 0.4501 | 0.3505 | **0.3941** | 0.2454 |

Agreement with MACS2 is a calibration target, not ground truth (D-023).
Recall is bounded by the fact that the reference is an independent ATAC
experiment in transduced rather than parental H524; published replicate
Jaccards between independent ATAC experiments run roughly 0.3-0.6.

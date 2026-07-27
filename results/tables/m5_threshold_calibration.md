# M5 — ATAC signal threshold calibration

Generated: 2026-07-27 07:53:07 UTC

Calibrated on **H524 / chr1** against its independent MACS2 peak set (GSE269424, lifted to hg19).
Reference peaks: 9832 covering 4,705,669 bp.

Thresholds are expressed as multiples of the line's mean signal over
covered bases, so the same multiple transfers across lines of differing depth.

| x mean | threshold | regions | median width | Mb | precision | recall | F1 | Jaccard |
|---|---|---|---|---|---|---|---|---|
| 2 | 3.359 | 17,413 | 253 | 6.5 | 0.255 | 0.352 | **0.2957** | 0.1735 |
| 3 | 5.039 | 5,325 | 318 | 2.45 | 0.467 | 0.2434 | **0.32** | 0.1905 |
| 4 | 6.719 | 3,957 | 336 | 1.85 | 0.527 | 0.2075 | **0.2977** | 0.1749 |
| 5 | 8.398 | 2,631 | 344 | 1.2 | 0.5969 | 0.1522 | **0.2425** | 0.138 |
| 6 | 10.078 | 1,961 | 333 | 0.84 | 0.649 | 0.1161 | **0.197** | 0.1092 |
| 8 | 13.437 | 1,303 | 308 | 0.51 | 0.685 | 0.0743 | **0.1341** | 0.0718 |
| 10 | 16.797 | 874 | 291 | 0.32 | 0.7188 | 0.0494 | **0.0924** | 0.0485 |
| 12 | 20.156 | 548 | 268.5 | 0.18 | 0.7478 | 0.0291 | **0.056** | 0.0288 |
| 15 | 25.195 | 290 | 250.5 | 0.09 | 0.7674 | 0.0145 | **0.0285** | 0.0145 |
| 20 | 33.593 | 107 | 251 | 0.03 | 0.787 | 0.005 | **0.0099** | 0.005 |

**Chosen: x3 mean covered signal** (F1 = 0.32).

Agreement with MACS2 is a calibration target, not ground truth. MACS2 with
an input control uses information these deposits do not contain. A high F1
means our regions behave like real peak calls in matched cells; it does not
make them peak calls, and the residual weakness is recorded in D-023.

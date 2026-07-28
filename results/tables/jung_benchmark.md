# Benchmark — Jung et al. 2017 18-gene Myc activity signature

Generated: 2026-07-28 10:48:08 UTC

Source: Jung LA et al., *Cancer Res* 77(4):971-981, 2017. PMID 27923830,
doi:10.1158/0008-5472.CAN-15-2906, Table 1. Transcribed from the paper;
17 up / 1 down, matching the paper's own description.

This comparator is **paralog-blind**, which is the point of using it.

## Signature

`AIMP2`, `APEX1`, `BMP7`, `CCND2`, `DKC1`, `EIF2S1`, `ENO1`, `HMGA1`, `LMO3`, `NAP1L1`, `NME1`, `ODC1`, `RPL5`, `RPS5`, `SIX1`, `SRM`, `UXT`, `NDRG1`

Down-regulated: `NDRG1`.

## 1. Regulon overlap — underpowered, reported as such

| paralog | hits | expected | OR | p |
|---|---|---|---|---|
| MYC | 2 | 0.50 | 4.49 | 0.088 |
| MYCN | 3 | 0.56 | 6.62 | 0.016 |
| MYCL1 | 3 | 0.51 | 7.33 | 0.013 |

An 18-gene signature against 500-gene regulons in a 10,387-gene universe expects fewer than one overlapping gene by chance. A null here means the test had no power, not that there is no relationship.

## 2. The signature tracks neither MYC expression nor lineage

| comparison | Spearman rho | p |
|---|---|---|
| vs NE score | +0.075 | 0.513 |
| vs MYC expression | -0.185 | 0.103 |
| vs ASCL1 | +0.052 | 0.651 |
| vs POU2F3 | +0.233 | 0.0389 |

Unique R2: MYC expression **0.0006**, lineage **0.0206**.

**Both are near zero.** Lineage is nominally the larger of the two, but on values of 0.0206 against 0.0006 that ordering carries no weight and is not reported as a finding. The honest reading is that this signature is largely unrelated to *either* covariate in this cohort.

What it does show is that a published, paralog-blind MYC-activity signature does **not** track MYC expression in SCLC tumours (rho -0.185, n.s.). That is consistent with this project's broader observation that MYC mRNA is a poor proxy for MYC regulatory activity in SCLC, and it is an independent line of support for it — but it is not evidence of the lineage confounding reported in M6, and must not be read as such.

## 3. Agreement with this project's regulon scores

| regulon | Spearman rho | p |
|---|---|---|
| MYC | +0.229 | 0.042 |
| MYCN | +0.438 | 5.32e-05 |
| MYCL1 | +0.289 | 0.00985 |

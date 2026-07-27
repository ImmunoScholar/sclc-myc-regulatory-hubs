# M5 — candidate accessible-region universes

Generated: 2026-07-27 07:35:33 UTC

**No universe is committed to by this script.** These are options measured
against real data so the support rule can be chosen deliberately.

| universe | regions | median width | total Mb | % genome |
|---|---|---|---|---|
| `union_all` | 290,862 | 459 | 212.9 | 7.012 |
| `ds_support_ge2` | 84,020 | 889 | 105.4 | 3.471 |
| `line_support_ge2` | 148,245 | 738 | 156.6 | 5.158 |
| `line_support_ge3` | 58,831 | 1130 | 87.7 | 2.889 |

## Cell-line composition problem

External ATAC peak sources cover H524, Lu139, H146, H82.
**Only H524 is also a keystone MYC-family ChIP line.**
The keystone's own ATAC covers 9 lines but ships as bedGraph signal with
no peak calls, so it contributes no intervals until 03_ derives regions
from signal.

A `>=2 independent datasets` rule therefore requires accessibility in two
largely non-keystone cell lines, selecting for shared elements and against
line-specific ones — the opposite of what a paralog-specificity question needs.

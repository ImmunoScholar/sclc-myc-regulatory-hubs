# Figure 2 — The M5 gate on the paralog-resolved regulatory layer

**A.** The four pre-registered gate criteria and their verdicts. Three pass,
one fails. Criterion 3 was blocked at M5 because amplification status could not
be verified, and was evaluated at M7 once DepMap Public 26Q1 copy number
resolved it (D-035); the verdict is recorded in the gate table itself (D-039).

**B.** Criterion 4, compositional motif share. Each column of the underlying
matrix is one paralog's region set split across the three E-box variants, so
shares sum to 1 within a set. CACCTG is the commonest variant in every set,
which is why an enrichment-over-background formulation failed on GC content.
The test asked instead whether each motif reaches its highest share in its own
paralog's regions. It does, in 3/3 (chi-square p = 1.6e-101).

**C.** Criterion 3, the failure. The direction reproduces — MYC-amplified lines
carry a higher distal fraction than MYC-expressing lines — but the magnitude is
+3.3 points against a published +27. Contributing factors, none separable at
2 lines vs 2: the shared ATAC-defined universe is ~84% distal by construction,
which compresses the achievable range; Plotnik called peaks de novo per line
with no shared grid; and H847 is absent from DepMap so the expressing group is
n=2. No significance test is applied, because none would be honest at this n.

**D.** Leave-one-line-out AUC. Each regulon is rebuilt with one contributing
line withheld and scored on that line. 8 of 9 reach the 0.70 threshold; the
exception is MYC/H1048 at 0.69. This is internal validity only — it shows the
programmes are not driven by a single line, not that they generalise to tumours.
Figure 1 shows that they do not.

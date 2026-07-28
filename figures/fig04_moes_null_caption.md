# Figure 4 — Convergence is real, weak, and not localisable to genes

**A.** Overlap between the top K genes of each domain, divided by the K²/N
overlap expected by chance, with the 95% permutation envelope in grey. This is
Robust Rank Aggregation's own question: do genes rank highly in both lists more
often than chance? They do, by up to ~2x, and the excursion clears the envelope
for two of three paralogs — global p 0.019 (MYC), 0.008 (MYCN), 0.065 (MYCL1).

The global p compares the maximum ratio over K against the permutation
distribution of that maximum, controlling family-wise error across the curve.
Counting pointwise exceedances would badly overstate this: the K values are
nested, so one excursion produces a run of them.

This coexists with near-zero Spearman correlation between the domains (-0.023,
+0.012, -0.008) without contradiction. Spearman measures monotone association
across all 10,387 genes; concordance-at-K measures agreement at the TOP of both
lists. The domains agree weakly about which genes matter most and not at all
about the ordering of the rest.

Plotted rather than a rank-rank scatter because the cis scores are heavily
tied — a gene with a single promoter link scores exactly 1 — so a scatter bands
into vertical stripes showing tie structure rather than association. Only the K
range where at least 5 genes are expected to overlap by chance is drawn; below
that the ratio is integer noise.

Coverage of the cis layer is also uneven. It is
uneven — 10,176 genes for MYC but 2,626 for MYCN and 4,586 for MYCL1, of 10,387 in the
MOES universe. Genes with no peak-to-gene link share a single tied rank, and
including them would draw a dense band that looks like structure but is only
the tie block. This unevenness is a real limit on what MOES could have found
for the two smaller paralogs, and it follows directly from MYCN and MYCL1
having two contributing ChIP lines each against MYC's five.

**B.** Empirical FDR by MOES rank, from 10,000 permutations. No curve approaches
the 0.05 threshold; the best FDR achieved across all three paralogs is 0.358.
This is not a marginal miss that a softer threshold would rescue.

Panels A and B answer different questions and the pair is the finding: the
aggregate overlap in A is detectable, and the per-gene evidence in B is not.
A signal of ~2x spread across the top of both rankings leaves no individual
gene strong enough to survive multiple testing over 10,387. Aggregate
detectable, per-gene unattributable — which is why no hub list is reported
even though the domains are not independent.

**C.** The power statement, without which panel B would be an absence
presented as a finding. A synthetic functional layer was blended with the real
cis ranking at increasing strength and passed through the identical FDR
machinery. The response is monotone: nothing on pure noise, discoveries
appearing once the domains correlate at rho ≈ 0.2, and 288 genes when the two
layers are identical. The observed correlation sits at the far left of this
axis, well below where the method demonstrably detects convergence.

**D.** Why any apparent paralog-specificity in a MOES ranking would be an
artefact. The functional layer is identical for all three paralogs, so it pulls
every ranking toward the same genes: paralogs agree at rho +0.50 to +0.59 in
MOES but only +0.17 to +0.33 in the chromatin that actually distinguishes them.
A top-N table drawn from this would present chromatin evidence under a
multi-omics label.

No prioritised hub list is reported. Doing so would be wrong twice over: it
would present genes whose individual FDR is ~0.36 as prioritised hits, and it
would imply a paralog attribution that panel D shows comes from chromatin
alone. The full ranking with its FDR column is in
`data/metadata/moes_ranking.csv` so this is inspectable rather than hidden.

MOES is a heuristic prioritisation, not a predictive or clinically validated
model (project contract, section 7).

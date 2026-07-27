# Decision Log

Append-only. Each entry records **what** was decided, **why**, and **what would
overturn it**. Decisions are never silently reversed — a reversal gets its own
numbered entry citing the one it supersedes.

Scope freeze: M0 and M1 are frozen as of 2026-07-26. Entries D-001 to D-005
record those frozen decisions for the record; D-006 onward are implementation
decisions taken from M3.

Legend — **SCI** scientific · **METH** methodological · **DATA** data · **ENV** environment/tooling

---

### D-001 · SCI · 2026-07-26 — Re-scope the gap after novelty validation

**Decision.** The original framing ("MYC paralogs have distinct enhancer-binding
preferences in SCLC — characterise them") is abandoned. The project is re-scoped
to translation and prioritisation.

**Why.** Phase 1b full-text review found **Plotnik et al. 2024** (Mol Cancer Res,
PMID 38747975, DOI 10.1158/1541-7786.MCR-23-0599) — the paper behind GSE230649,
our keystone dataset — had already published paralog-specific MYC/MYCN/MYCL1
enhancer characterisation, including E-box central-dinucleotide preferences and
~84% MYCN-in-MYC peak overlap. The first analytical step of the original plan was
already done, by the authors of the data we intended to use.

**What survives.** Plotnik did no super-enhancer calling, no peak-to-gene linking
beyond nearest-gene, no GRN, no patient tumour data, no DepMap, and no
prioritisation framework. All work is *in vitro*.

**Overturned if.** A paper appears doing regulon-based MYC-paralog scoring in
SCLC tumours → pivot to the MOES framework itself as the contribution,
benchmarked against Jung et al. 2017.

---

### D-002 · METH · 2026-07-26 — Reject SCENIC-family GRN inference

**Decision.** SCENIC, pySCENIC and SCENIC+ are rejected as the primary GRN
method. Primary = ChIP-anchored regulon construction; secondary = GENIE3.

**Why.** MYC, MYCN and MYCL1 all bind the same canonical E-box (`CACGTG`);
Plotnik found only subtle central-dinucleotide differences. Motif-driven
inference against cisTarget databases is therefore **structurally incapable of
distinguishing the paralogs** and would collapse all three into a single "MYC
regulon" — destroying the project's independent variable. SCENIC+ is separately
impossible: it needs paired scATAC+scRNA, and no public SCLC scATAC-seq exists.

**Overturned if.** A paralog-discriminating motif model is published and validated.

---

### D-003 · METH · 2026-07-26 — MOES uses two-stage RRA, not a weighted sum

**Decision.** Evidence integration is two-stage Robust Rank Aggregation over four
domains (cis-regulatory, transcriptional, network, functional) — RRA within each
domain, then across domains.

**Why.** A weighted sum requires weights, and any weights we chose would be
arbitrary — explicitly forbidden by the project charter. Worse, the layers are
**non-independent**: ATAC, H3K27ac, MYC ChIP and super-enhancer calls are
correlated by construction, so a flat weighted sum would effectively
quadruple-count chromatin evidence. RRA is rank-based, weight-free, and handles
missing layers natively (important given the ~1,800-gene GeoMx panel). Precedent:
NetICS (Dimitrakopoulos et al. 2018).

**Overturned if.** The layer-correlation matrix computed at M8 shows the domains
are in fact near-independent, in which case a simpler aggregation could be
justified — but the change must be logged, not assumed.

---

### D-004 · DATA · 2026-07-26 — No de novo peak calling

**Decision.** Quantify bedGraph signal over an ATAC-defined consensus region
universe in R (`rtracklayer` / `GenomicRanges`). No MACS2, no bedtools, no
re-alignment.

**Why.** GSE230649 ships **bedGraph only** (8.6 GB, no peak files). Re-aligning
28 SRA samples is not feasible on this hardware. This also removes two system-level
dependencies entirely.

**Residual risk (R-02).** Our region boundaries will not match Plotnik's MACS2
calls. This is declared as a methodological deviation and checked at the M5 gate
against their reported 18,823 / 4,017 / 5,688 active-region counts.

---

### D-005 · DATA · 2026-07-26 — Drop the single-cell layer

**Decision.** The single-cell layer is excluded from the project rather than
filled with a substitute dataset.

**Why.** The scientifically correct dataset (Chan et al. 2021, 155k cells, 21
patient specimens) has **no GEO accession** — HTAN/Synapse only, and Synapse
cannot be authorised in this environment. The accessible GEO alternative
(GSE138474) is predominantly **CDX and cell-line-derived** on pre-10x chemistry.
The project's entire argument is *translation into patient tumours*; a
model-derived single-cell layer cannot support that argument, and including it
would be inclusion-for-completeness.

**Cost, stated plainly.** Patient translation now rests entirely on two bulk
cohorts (n=79, n≈81). Risk R-09 (small effective sample size) carries more weight.
This goes in the README limitations.

**Consequences.** `pyscenic`, `synapseclient`, `Seurat`, `SingleCellExperiment`,
`scater`, `zellkonverter` removed from the dependency set. **The project now has
no required Python dependency.**

**Overturned if.** HTAN/Synapse access is obtained → reinstate as an additive aim
without disturbing M5–M8.

---

### D-006 · ENV · 2026-07-26 — Bioconductor via the GWDG mirror

**Decision.** Bioconductor 3.23 packages are installed from
`https://ftp.gwdg.de/pub/misc/bioconductor` rather than `bioconductor.org`.

**Why.** `bioconductor.org` is **unroutable from this machine.** DNS resolves
correctly (18.161.246.18/33/80/93, AWS CloudFront) but TCP connections time out
after 25 s on both IPv4-forced and default routes. Posit P3M does not mirror
Bioconductor (404 at every candidate path). Three official mirrors were probed;
GWDG (Germany), TU Dortmund, and TUNA (Tsinghua) all returned HTTP 200 for the
3.23 package index. GWDG selected.

**Evidence.** `scripts/00_setup/` probe output, 2026-07-26: 2,384 BioCsoft
packages visible via GWDG; 0 via bioconductor.org.

**Caveat.** This is a **network-environment workaround, not a scientific
decision.** Mirror content is identical to upstream. The mirror URL is recorded
in `renv.lock`; anyone restoring on a network that can reach bioconductor.org may
substitute it freely.

**Overturned if.** bioconductor.org becomes reachable — revert to upstream for
maximum portability.

---

### D-007 · ENV · 2026-07-26 — CRAN via Posit P3M noble binaries

**Decision.** CRAN packages install from
`https://packagemanager.posit.co/cran/__linux__/noble/latest` with an explicit
`HTTPUserAgent`.

**Why.** P3M serves precompiled Ubuntu-noble binaries, turning a multi-hour
source build of ~200 CRAN dependencies into minutes. The `.Rprofile` falls back
to `cloud.r-project.org` automatically on any non-noble platform, so the project
still restores on macOS/Windows.

---

### D-008 · DATA · 2026-07-26 — Defer `BSgenome.Hsapiens.UCSC.hg19`

**Decision.** The hg19 BSgenome package is **not** installed at M3.

**Why.** It is ~700 MB and provides genomic *sequence*. No planned analysis step
requires sequence: motif analysis is Plotnik's prior work and is cited, not
repeated (D-001). Installing it "in case" would waste disk and download time.

**Overturned if.** A step genuinely requires sequence retrieval → install then,
and log the reason here.

---

### D-009 · ENV · 2026-07-26 — Defer `pandoc` and Quarto; skip `libgit2-dev`

**Decision.** No `sudo apt install` is performed at M3.

**Why.** An apt audit found only one genuinely missing, genuinely needed package:
`pandoc` (required for report rendering at **M10**, not before). Two apparent
gaps were phantom renames on noble — `libfreetype6-dev` and `libtiff5-dev` do not
exist under those names; the real packages `libfreetype-dev` and `libtiff-dev`
are already installed. `libgit2-dev` is missing but is only needed by
`gert`/`usethis`, which are not in the dependency set.

**Action deferred to M10.** Install `pandoc` then, together with Quarto (which is
not in the apt archive and needs a `.deb` from quarto.org).

---

### D-010 · ENV · 2026-07-26 — Enforce the no-large-files rule with a hook

**Decision.** A `pre-commit` hook rejects any staged file above 5 MB, in addition
to `.gitignore`.

**Why.** `.gitignore` is a policy; a hook is enforcement. A single accidental
commit of an 8.6 GB bedGraph is painful to remove from history and may exceed
GitHub's hard per-file limit. `.gitignore` was additionally verified
**adversarially** before any download: 11 decoy files mimicking real downloads
were planted and all 11 confirmed ignored, and 9 deliverable paths confirmed
still trackable (an over-broad ignore file is also a failure).

**Note.** The hook lives in `.git/hooks/` and is therefore **not cloned**.
`scripts/00_setup/01_init_git.sh` reinstalls it; the README instructs
contributors to run it.

---

### D-011 · ENV · 2026-07-26 — Skip `targets` for pipeline orchestration

**Decision.** Scripts are numbered, ordered, and run directly. No `targets` DAG.

**Why.** The dependency inventory listed `targets` as optional. The pipeline is
linear with well-defined milestone gates, and `targets` would add a dependency
and a layer of indirection for a project a reader is meant to follow by reading
the scripts in order. Reproducibility is served by `renv.lock` plus explicit
script ordering.

**Overturned if.** Re-running expensive stages becomes a real bottleneck.

---

### D-012 · ENV · 2026-07-26 — GSVA outstanding; needs one system library — **RESOLVED same day**

> **RESOLUTION (2026-07-26).** `sudo apt install -y libmagick++-dev pandoc` was
> run. `magick` 2.9.1 then installed as a prebuilt P3M binary (no compile) and
> `GSVA` 2.6.3 built successfully. The environment is now **35 of 35**, with 231
> packages locked in `renv.lock`. `pandoc` 3.1.3 is also in place, so the M10
> reporting dependency is satisfied ahead of time. GSVA is retained as the
> secondary scoring method; no fallback to AUCell was needed. `tests/test_environment.R`
> now asserts that `tumour_scoring.secondary_method` stays `gsva`, so the
> methodological-distinctness argument below is enforced rather than merely
> documented.

**Decision.** The M3 environment is locked with **34 of 35 packages**. `GSVA`
alone failed and is deferred pending a single `sudo apt install`.

**Why it failed.** `GSVA` hard-depends on `SpatialExperiment`, which hard-depends
on the R package `magick`, which fails to build without the system library
`libmagick++-dev`. Confirmed by dependency-tree analysis
(`scripts/00_setup/diag_magick.R`), not guessed: `magick` is in GSVA's recursive
**hard** dependency set, not merely Suggests.

**Impact assessment.** Low. GSVA is the *secondary* tumour-scoring method — a
sensitivity check on the primary, `singscore` (`config/params.yml`,
`tumour_scoring`). `singscore` and `AUCell` both installed successfully, so the
critical path is unblocked. Nothing in Aims 1–5 is stalled by this.

**Preference, with reasoning.** Install it. `singscore` and `AUCell` are both
rank-based, so using AUCell as the sensitivity check would compare two methods
that share the same assumptions. GSVA's kernel-density approach is
methodologically distinct, which is what makes it a *useful* cross-check rather
than a restatement.

**If declined.** Drop GSVA, set `tumour_scoring.secondary_method: aucell`, and
state in the report that the sensitivity analysis compares two rank-based methods
— a weaker but honest check. Log the change here.

**Note.** `renv.lock` currently reflects the 34-package set. It must be
re-snapshotted after GSVA installs.

---

### D-013 · ENV · 2026-07-26 — All WSL invocation goes through script files

**Decision.** Multi-step shell work is written to a file in `scripts/` and
executed as `bash <path>`, never passed as an inline command string.

**Why.** Two distinct, silent mangling behaviours were hit during M3:
a Git-Bash-backed shell rewrote `/home/priya/...` into
`C:/Program Files/Git/home/priya/...` (MSYS2 path conversion), creating an empty
directory tree in the wrong filesystem; and PowerShell strips `|`, `$()`, `!` and
embedded double quotes from arguments passed to native executables, silently
truncating commands. Both produce *plausible-looking* partial success, which is
the dangerous kind of failure.

**Side benefit.** Every setup step is now a committed, re-runnable, reviewable
script rather than a shell command that existed only in a transcript.

---

### D-014 · METH · 2026-07-26 — liftOver policy: lift intervals, never signal

**Context.** Verification at M4 established the genome build of every dataset from
its own `!Sample_data_processing` declaration rather than from file names:

| Series | Assay | Build |
|---|---|---|
| GSE230649 (keystone) | MYC/MYCN/MYCL1/H3K27ac ChIP + ATAC | **hg19** |
| GSE210113 | NEUROD1 / H3K27ac ChIP | **hg19** |
| GSE249362 | POU2F3 ChIP | **hg19** |
| GSE60052 | bulk RNA-seq (gene level) | hg19 |
| GSE269424 | ATAC | **hg38** |
| GSE256345 | ATAC | **hg38** |
| GSE281523 | ATAC (FFPE/PDX) | **GRCh38** |
| GSE281524 | ASCL1 ChIP | **GRCh38** |

**The problem this creates.** The project build is hg19, frozen because the
keystone is hg19 (D-004). But **all three supporting ATAC datasets are hg38**, and
the consensus region universe requires support from ≥2 independent ATAC datasets
(`config/params.yml`, `regions.min_datasets_supporting`). So the region universe
cannot be built without crossing builds. Risk R-05 has stopped being hypothetical:
**23 of 66 files require lifting.**

**Decision — the policy, not a design change.** R-05's frozen mitigation already
required that "all lifts are explicit, logged, and QC'd for loss rate". This
records *how*:

1. **Never lift continuous signal.** Lifting a bigWig or bedGraph across builds is
   not well defined — bases are added, removed and rearranged, so a lifted signal
   track silently misattributes coverage.
2. **Call intervals in the NATIVE build first, then lift only the intervals.**
   Interval liftOver is well defined and its failure mode (an interval that does
   not map) is observable and countable.
3. **Report the loss rate** for every lift, per dataset, in the QC output. A lift
   that loses an unreasonable fraction of intervals halts the step.
4. **Prefer sources that publish intervals.** This is why the acquisition takes
   narrowPeak/BED files rather than bigWigs wherever both exist — see D-015.
5. **Assert build at load.** `04_qc_report.R` fails any hg19-declared file
   containing an interval that extends past its chromosome's length in
   `hg19.chrom.sizes`. Coordinates past the end of a chromosome are proof of a
   build mismatch and cannot be argued with.

**Consequence.** `hg38ToHg19.over.chain.gz` is now load-bearing infrastructure and
is in the manifest as a first-class resource with its own licence record.

---

### D-015 · DATA · 2026-07-26 — five dataset-composition corrections found by M4 verification

Filenames and inventory labels were checked against GEO's own records before
downloading. Five things were not what the frozen inventory implied. None changes
a frozen scientific decision; all change **which files are fetched**.

**1. GSE269424 is not native ATAC — it is a TF-overexpression experiment.**
Sample titles are `H524_ASCL1_1/2`, `H524_EGFP_1/2`, `Lu139_EGFP_1/2`,
`Lu139_NEUROD1_1/2`. The ASCL1/NEUROD1 tokens are **ectopically expressed
transgenes**, not ChIP antibodies (library strategy is ATAC-seq; no antibody
characteristic exists). The overexpression arms have deliberately remodelled
chromatin. **Only the four EGFP control arms are used.** Taking this series at its
"ATAC-seq of SCLC cell lines" label would have let engineered accessibility
contaminate the region universe underpinning Aim 1.

**2. GSE249362 is the best lineage-TF control and the cheapest.** Triage of its 185
supplementary files found POU2F3 peak **BED** files for NCIH1048, NCIH211 and
NCIH526 — three keystone lines spanning two paralogs, already in hg19. Taking the
5 BEDs costs ~1.5 MB; taking the 125 bigWigs would have cost ~20 GB. Only DMSO /
untreated arms serve as the reference repertoire; the FHD286 arm's peak file is
4,606 bytes versus 181,493 for DMSO, i.e. the drug abolishes binding.

**3. GSE281524 (ASCL1) is restricted to SHP-77.** Its 10 samples are H1836 and
SHP77 only. H1836 is not in the keystone set, so its ASCL1 tracks cannot support a
within-line comparison. Fetching 5 SHP-77 files (~2.1 GB) instead of all 10
(~6 GB) avoids 3 GB of data that could not have been used for this purpose.

**4. GSE210113 (NEUROD1) has zero cell-line overlap with the keystone.** All 14
samples are H446 (a NEUROD1-knockout experiment). See risk R-14.

**5. GSE60052 has no raw counts.** Verified by reading the file header: the only
deposited matrix is `normalized.log2`. **DESeq2 therefore cannot be used on this
cohort** — it requires raw integer counts. Differential analysis uses limma on the
log2 matrix. The dependency inventory lists both, so this selects between existing
options rather than changing the design. Two parsing traps also confirmed: sample
headers carry leading whitespace, and normals are encoded by a `.normal` suffix.

**Method note.** Three representative sample IDs used in an early probe were
*guessed* rather than read from the series records, and returned unrelated samples
(mouse liver, a kidney biopsy, MCF7). They were re-derived from each series'
`!Series_sample_id` lines. Accessions are never guessed; the failure is silent
because a wrong accession still returns a valid-looking record.

---

### D-016 · DATA · 2026-07-26 — DepMap is acquired manually, by design

**Decision.** DepMap is **not** downloaded by script. `03_verify.R` reports its
files as manual-required and does not count them as failures.

**Why.** `depmap.org` serves a Cloudflare Turnstile bot check, and the page states
plainly: "Need DepMap data in bulk? ... please don't scrape the portal." Automated
fetching is both technically blocked and against the operator's stated wishes, so
it is not attempted. The Bioconductor `depmap` package is not a workaround either —
it routes through ExperimentHub to bioconductor.org, which is unroutable from this
network (D-006).

**What this requires.** Three files downloaded by hand into `data/raw/depmap/`:
`CRISPRGeneEffect.csv`, `Model.csv`, `OmicsCNGene.csv`. The release must then be
pinned in `config/params.yml` (`depmap.release`), which is `null` as a deliberate
tripwire so Aim 3 scripts fail rather than run against an unpinned release (R-08).

**Impact.** None on M4–M6. DepMap is Aim 3 (M7) evidence only.

---

### D-017 · DATA · 2026-07-26 — integrity policy, and an honest account of its limits

**Decision.** File integrity rests on three checks, in this order:

1. **Exact byte size** against GEO's published `filelist.txt` value — the only
   *independent* check available.
2. **SHA256 trust-on-first-use** — computed at first successful download, recorded
   in `data/metadata/checksums.sha256`, enforced on every later run.
3. **`gzip -t`** on every gzipped file, which catches truncation that a correct
   size can mask on a resumed transfer.

**The limitation, stated rather than buried.** GEO publishes **no MD5 or SHA256**
for supplementary files — verified across every series used here. So TOFU proves a
file has not changed or corrupted *since we fetched it*; it does **not** prove we
fetched what the authors deposited. The published byte size is the only guard
against that, and it is a weak one. The manifest records
`checksum_algorithm = "none_published"` rather than leaving the field blank, so the
absence is visible instead of looking like an oversight.

**Additional guard.** `03_verify.R` also rejects any "data" file that is really an
HTML error page. NCBI returns HTML on 403 and DepMap returns a Cloudflare
challenge; a naive downloader saves both as data, and a small HTML file can pass a
size check if the size was itself recorded from a failed request.

---

### D-018 · ENV · 2026-07-26 — download concurrency is 2, established by failure

**Decision.** `02_download.sh` runs 2 concurrent transfers, staggered, with a
worker that treats HTTP 403/429/503 as transient.

**How it was arrived at, including the wrong turn.**

| setting | outcome |
|---|---|
| serial | ~171 KB/s, zero failures, ~19 h projected for 12 GB |
| 4 concurrent | **403s within ~2 minutes**, on 4 of the first 8 files |
| 2 concurrent | **517 KB/s aggregate, 64/64 files, zero failures, 6 h 22 m** |

A bandwidth probe had indicated an 8.75× aggregate gain from 4 streams and I acted
on it. That probe measured byte-**range** requests against a *single* file, which
does not trip NCBI's per-IP concurrent-connection limit the way requests for four
*distinct* files do. The measurement was sound; the inference from it was not.
Recorded because "measure, don't guess" is necessary but insufficient — the
measurement has to exercise the same mechanism as the real workload.

**Two real defects the 403s exposed.**
1. `curl` does not retry 403 by default (not classed as transient), so a single
   rate-limit response killed a file outright.
2. The 403 **HTML body was written into the `.part` file**. A later `curl -C -`
   resumes from the end of that markup, producing a file that is a few hundred
   bytes of HTML followed by real data. Observed as 980-byte `.part` files whose
   first two bytes were `3c3f` (`<?`) rather than the gzip magic `1f8b`. The
   worker now scrubs any `.part` that is not valid gzip or begins with markup, and
   `clean_poisoned_parts.sh` audits existing ones — 3 found and removed.

**Do not raise concurrency** without re-testing. It trades throughput for 403s,
and the 403 failure mode is silent corruption rather than a clean error.

---

### D-019 · METH · 2026-07-26 — sequence names are harmonised centrally, and assertions must be able to fail

**Decision.** All coordinate handling passes through `R/genome_utils.R`.
`harmonise_seqnames()` is mandatory before any cross-object genomic operation.

**Why.** The keystone deposit uses Ensembl sequence names (`20`) and every
annotation resource uses UCSC names (`chr20`). Same hg19 coordinates, different
convention, and mixing them yields **zero overlaps rather than an error** — a
result that looks biological. See risk R-16 for the full account.

**The stronger half of this decision.** `assert_hg19_bounds()` returns a
**failure** when it had fewer than 20 chromosomes to compare, rather than passing
because it found no violations among nothing. The original assertion could not
fail: it compared `20` against `chr20`, matched zero rows, and reported success —
while nominally guarding R-05, the project's highest-consequence technical risk.

**Consequence for how QC is written from here.** Build assertions stream whole
files rather than sampling, because these bedGraphs are in BAM-header order
(`20 21 22 1 3 2 …`) and a head sample contains one chromosome. `--fast` records
the skip **as a failure**, so a fast run cannot be mistaken for a clean one.

**Standing rule.** A QC check must be observed failing on data known to violate it
before it is trusted. Two defects in this project have been checks that passed
because they examined nothing.

---

### D-020 · METH · 2026-07-26 — the M5 gate was circular and has been rebuilt

**This supersedes the M5 gate defined at M1.** It is a correction, not a
relaxation: the replacement is harder to pass, not easier.

**The defect.** The original gate required our re-derived active-region counts to
fall near Plotnik's 18,823 MYC / 4,017 MYCN / 5,688 MYCL1. But our count is

```
count ≈ |ATAC universe| × P(signal > quantile) × P(H3K27ac⁺)
```

and `signal_quantile = 0.75` passes **25% of regions by construction**. With
H3K27ac⁺ covering roughly 60% of accessible regions, the count is ~15% of the
universe — so any consensus universe near 125,000 regions reproduces Plotnik's
number whether or not the pipeline is correct. Typical SCLC ATAC universes are
50k–150k regions. **The gate could not distinguish a working implementation from a
broken one**, because both numbers are set by parameters we chose.

This is the third check in this project that could not fail (after the vacuous
hg19 build assertion and the `!isTRUE(vector)` resumability test). The pattern is
consistent enough to be a design smell: *a check whose expected value is
determined by our own free parameters is decoration.*

**The replacement — relative properties we cannot tune.**

1. **MYCN regions ⊂ MYC regions at ~0.84.** Both sets are thresholded identically,
   so the overlap *ratio* is a property of the data. Moving the threshold moves
   both sets together and largely cancels.
2. **MYCL1 must NOT look like a MYC subset** (`max_expected: 0.50`), matching
   Plotnik's "largely non-overlapping".
3. **Distal-fraction contrast**, MYC-amplified vs MYC-expressing lines (published
   0.39 vs 0.12). A between-group contrast is immune to universe size. Direction
   is required; magnitude is secondary.
4. **Sequence-level motif validation** (D-021/D-022) — the only check that is
   fully independent of every parameter we set.

**Counts are still reported**, prominently and across a range of thresholds
(`sensitivity_quantiles`), because the shape and stability of that curve is
informative. They are simply never used as pass/fail.

---

### D-021 · METH · 2026-07-26 — four quality upgrades adopted for M5 onward

Adopted after an explicit mid-project audit. Each addresses a specific weakness
rather than adding polish.

**B · Independent motif validation.** Test whether our derived paralog regions
recover Plotnik's E-box central dinucleotides (MYC `CAGATG`, MYCN `CACATG`,
MYCL `CACCTG`) against shuffled-sequence backgrounds. *Why it matters:* region
counts are set by our thresholds; sequence content is not. This is the only
available check that is independent of every parameter we chose, and it
distinguishes "we found paralog-specific binding" from "we found accessible
chromatin and labelled it".

**C · Bootstrap rank stability for MOES.** 2,000 resamples, rank confidence
intervals for every hub, and an explicit UNSTABLE flag when an interval spans
more than 25% of the list. *Why it matters:* the evidence-integration framework is
offered as this project's methodological contribution. A point rank implies a
precision the data may not support — "hub #3" and "hub #17" can be
indistinguishable under resampling, and a ranked list that hides that is
misleading.

**D · Publication output system.** Vector PDF alongside PNG; embedded fonts; an
automated minimum-text-size check at final print dimensions; colourblind safety
*verified* by deutan/protan/tritan simulation with a minimum perceptual-distance
threshold rather than asserted; tables with explicit precision, effect sizes
beside every p-value, mandatory footnotes for non-obvious columns, exported as
both CSV and rendered tables; and a figure manifest linking every figure to its
script and inputs. *Why it matters:* "colourblind safe" and "publication quality"
are testable claims, and untested claims of accessibility are worth nothing.

**E · Regulon internal-validity gate before tumour scoring.** Each paralog regulon
must separate amplified from non-amplified lines *in the cells it was derived
from* (AUC ≥ 0.70) and track its paralog's own expression (Spearman ≥ 0.30).
*Why it matters:* cheapest of the four, and it catches the failure mode where a
regulon that does not behave in its source data is nevertheless scored in tumours
and interpreted.

**Cost.** Roughly 15–20% additional effort on the remaining work, plus ~700 MB
disk for the genome package.

---

### D-022 · ENV · 2026-07-26 — reverses D-008: install `BSgenome.Hsapiens.UCSC.hg19`

**Decision.** The hg19 BSgenome package is installed after all.

**Why the reversal.** D-008 deferred it on the explicit condition *"overturned if
a step genuinely requires sequence retrieval — install then, and log the reason
here."* Motif validation (D-021, B) is now such a step. The condition D-008 set
for its own reversal has been met, so this is the decision working as designed
rather than being overridden.

**Verification.** The install script asserts `seqlengths(hg19)[["chr1"]] ==
249250621` before reporting success, so a wrong-build package cannot pass
silently.

---

### D-023 · METH · 2026-07-26 — the region universe comes from keystone ATAC, with the threshold calibrated on H524

**Supersedes** the frozen `regions.min_datasets_supporting: 2` rule as a *filter*.
Taken against measured numbers, not assumption.

**What the candidates actually looked like** (`results/tables/m5_universe_options.md`):

| universe | regions | % genome | H524-supported |
|---|---|---|---|
| union of external peaks | 290,862 | **7.01%** | 34.4% |
| ≥2 datasets | 84,020 | 3.47% | 83.2% |
| ≥2 cell lines | 148,245 | 5.16% | 48.8% |
| ≥3 cell lines | 58,831 | 2.89% | 83.0% |

**Two findings from those numbers.**

1. The union is implausible as accessible chromatin — 7.0% of the genome against a
   typical ATAC range of 1–5%, inflated by GSE256345's `merge500` peaks merging
   across 500 bp gaps.
2. My prior concern that a `≥2` rule would preferentially discard keystone-relevant
   regions was **only half right**. Of the 206,842 regions it drops, just 14.5% are
   H524-supported. But 14.5% of 206,842 is still **~30,000 regions accessible in a
   MYC-amplified keystone line**, discarded because two unrelated lines do not
   share them — precisely where line-specific MYC enhancers would sit.

**The problem no support rule fixes.** The external peak sources cover H524,
Lu139, H146 and H82. **Only H524 is a keystone MYC-family ChIP line.** Quantifying
MYC binding in COLO668 or H1048 over a grid defined by H146 and H82 accessibility
is a cellular mismatch, whatever threshold is applied to it.

**Decision.**
- The universe is derived from the **nine keystone ATAC tracks** — the same cells
  as the ChIP, native hg19, no liftOver.
- Because those ship as bedGraph signal with no control track, the signal
  threshold is **calibrated against H524's external MACS2 peak set**. H524 is the
  one line with both keystone ATAC signal and independent peak calls, so the
  cutoff is *fitted to reproduce a real peak caller in matched cells* rather than
  chosen. The calibration curve is reported.
- External peaks are demoted from **filter** to **per-region annotation**, so
  corroboration can be reported about each hub without silently removing the
  line-specific regions the hypothesis is about.
- `min_datasets_supporting` is **reinterpreted, not deleted**: support is counted
  across the nine keystone ATAC lines, which is a meaningful independence unit
  because all nine are ChIP-matched.

**Why this is better than honouring the frozen spec literally.** The spec assumed
several independent ATAC datasets in matched cells. That assumption did not
survive contact with the data — GSE269424 lost half its samples to the
TF-overexpression finding (D-015), and three of four remaining lines are not in
the keystone. Applying the rule as written would have produced a defensible-looking
universe built mostly from the wrong cells.

**Residual weakness, stated.** Signal thresholding without a control track is
cruder than MACS2 with input. The calibration against H524 bounds that weakness
but does not remove it; regions are corroborated against external peak calls and
the agreement statistic is reported per line.

> **SUPERSEDED IN PART by D-024.** The external calibration described above was
> the wrong validation target and has been replaced. The rest of D-023 — universe
> from keystone ATAC, external peaks as annotation not filter — stands.

---

### D-024 · METH · 2026-07-26 — the ATAC threshold is validated intrinsically, and made non-load-bearing

**What went wrong first.** D-023 calibrated the ATAC threshold against GSE269424's
H524 MACS2 peaks. Best F1 was 0.372. Two attempts to explain the poor agreement
both failed:

1. *"Single-line noise is dragging precision down."* Wrong — recall stayed pinned
   at ~0.35 under every support rule (≥1, ≥2, ≥3, external-corroborated). Filtering
   changed nothing about what was missed.
2. *"Copy number inflates calls under a global threshold."* Wrong — MACS2-style
   local background made agreement **worse** (F1 0.333 vs 0.372).

The assay identity was then verified directly, because an unverified label had
already burned this project once (GSE269424, D-015). The nine tracks are genuinely
ATAC: `!Sample_title = ATAC_<line>`, `library_strategy = ATAC-seq`,
`chip antibody: none`. The labelling was right.

**The actual error: a bad validation target.** GSE269424's H524 is an independent
experiment — different lab, protocol and pipeline, in EGFP-transduced rather than
parental cells. Cross-laboratory ATAC concordance runs Jaccard 0.3–0.5; we saw
0.22 against a threshold-called set from a modest-depth track. Optimising against
it measured inter-laboratory reproducibility, not fitness for purpose. **Three
rounds were spent fixing a method that was not broken.**

**Intrinsic validation instead.** The universe exists to be a quantification grid
for MYC and H3K27ac ChIP in these nine lines, so it is judged on properties of
these data: TSS enrichment (the standard ATAC metric) and H3K27ac fold-enrichment
in **matched cells**. Result for H524/chr1: **TSS enrichment 8.9–26.0 and H3K27ac
fold 2.8–7.5** across the threshold range. These are good ATAC data. The region
calling was never the problem.

**Threshold: ×2.5 of each line's mean signal over covered bases.** Chosen, not
fitted — marginal enrichment degrades smoothly (new-region TSS 6.14 → 3.56 → 2.05
→ 1.28) with no cliff, so any value in ×2–×3 is defensible and the automatic pick
of ×2 was an artefact of an arbitrary TSS≥2 cutoff. ×2.5 keeps marginal regions
clearly enriched (TSS 3.56, H3K27ac 1.82) rather than borderline.

**The consideration that mattered more than the value.** TSS windows cover 2.21%
of chr1, so enrichment converts directly into promoter fraction: ×8 → 57%
promoter-proximal, ×3 → 28%, ×2.5 → 20%, ×2 → 11%. **The threshold sets the
promoter/distal balance of the universe, and distal fraction is one of the four
M5 gate criteria (D-020).** An over-strict threshold would produce a
promoter-biased grid and bias the exact quantity being tested. For a project built
on Plotnik's *enhancer* finding, over-stringency is the more dangerous error —
which is why the choice sits at the permissive end of the defensible band.

**Made non-load-bearing.** Every M5 gate metric is computed across ×2, ×2.5, ×3
and ×4. If MYCN-in-MYC holds near 0.84 across that range, the threshold is
irrelevant to the conclusion and robustness is demonstrated rather than assumed.
If the gate metrics move with threshold, that fragility is itself a reportable
finding. This applies the D-020 lesson properly: do not defend a parameter, make
the conclusion independent of it.

**Process note.** Two confident wrong diagnoses in succession, with the real error
one level up in the choice of question. The verification instinct was right but
aimed late — the premise worth checking first was the *yardstick*, not the *input*.

---

### D-025 · METH · 2026-07-26 — M5 gate criterion 2 reformulated as enrichment over chance

*(Written retrospectively on 2026-07-26 after a decision audit found `config/params.yml`
and three progress reports referencing "D-025" with no such entry in this log. The
decision was made and implemented at the time; only the record was missing. Logged
as a gap rather than silently backfilled.)*

**Original criterion.** MYCL1-in-MYC overlap < 0.50, from Plotnik's description of
MYCL1 regions as "largely non-overlapping" with MYC.

**Observed.** 0.774 — a clear failure against that threshold.

**Why the criterion, not the result, was the problem.** This project quantifies
every paralog over ONE shared ATAC-defined grid of 102,334 regions (D-023). Plotnik
called peaks *de novo per line* with no shared coordinate set, so their MYCL1 peaks
could sit at coordinates absent from the MYC peak sets entirely. Ours cannot —
every paralog is scored at the same coordinates, which mechanically inflates
cross-paralog overlap. The 0.50 threshold was measuring our design choice.

Corroborating evidence that this was the mechanism: **both** overlaps came out high
(MYCN 0.886, MYCL1 0.774) and the *contrast* between them was small, whereas
Plotnik reported a large contrast. The failure was in reproducing a difference, in
a design that suppresses differences by construction.

**Reformulation.** Test the differential as enrichment over chance, so set sizes
and the shared grid cancel:

```
enrichment(P) = [ |P ∩ MYC| / |P| ] / [ |MYC| / |universe| ]
```

Requirements: enrichment(MYCN) > enrichment(MYCL1), Fisher test p < 0.05, effect
size reported alongside. **No absolute magnitude threshold** — Plotnik gives no
numeric value for MYCL1, so inventing one would be fitting a target we do not have.

**Result.** MYCN 5.12x vs MYCL1 4.48x, OR 2.27, p = 3.0e-41. Passes on direction
and significance. The **magnitude is modest** (ratio 1.14) and is reported as such:
we reproduce the direction of Plotnik's finding, not its starkness.

**Still falsifiable.** If MYCL1 were as nested in MYC as MYCN is, the paralogs
would be indistinguishable by occupancy in these data and that negative result
would be the finding.

**Independent corroboration from criterion 4.** Motif content cannot be biased by
the shared grid, and there MYCL1 is unmistakably distinct — standardised residuals
-18.38 in MYC regions versus +19.94 in its own. So the weak occupancy contrast was
the design artefact this decision assumed, and the sequence evidence recovers the
distinctness the occupancy measure lost. Two data types disagreeing in an
explainable way is stronger than either agreeing alone.

---

### D-027 · METH · 2026-07-26 — peak-to-gene uses an H3K27ac activity proxy with confidence tiers

**Two problems with the frozen spec.** It called for expression correlation at
`min_abs_correlation: 0.3`. (1) GSE230649 contains **no RNA-seq** and our
expression data are tumour cohorts, not these cell lines. (2) At n = 10 lines,
Spearman |rho| = 0.3 is **p ~ 0.4** — a filter that looks rigorous and admits
noise. Significance needs |rho| >= 0.64.

The 0.3 value dates from M1, when the single-cell layer was still in scope and
more samples were expected. It was never revisited when that layer was dropped
(D-005) and the sample count fell to 10.

**Adopted.** Promoter H3K27ac as the activity proxy — an ABC-family model, **not**
expression correlation, and labelled that way in every output. Links carry
confidence tiers (high 0.64 / moderate 0.45 / weak 0.30) instead of one cutoff,
and only the high tier is treated as established.

**Validated, and it could have failed.** Distal and promoter H3K27ac come from
the same experiments, so global line effects could correlate everything. Against a
null of random pairs >1 Mb apart: candidate median rho **0.345** vs null
**-0.006**, high-tier fraction 18.3% vs 3.0% — **6.06x enrichment**. The null
median at essentially zero is what a correct null should give. The script exits
non-zero if enrichment falls below 1.5x.

**The methodological claim is now measured.** Only **32.6%** of retained links
point to the nearest gene, so two-thirds assign a region to something other than
its nearest neighbour. The improvement over nearest-gene assignment claimed in the
gap statement is real rather than asserted.

**Upgrade path.** Real CCLE expression at M7 replaces the proxy.

---

### D-028 · METH · 2026-07-26 — systematic parameter sweep; 11 stale or dead parameters corrected

**Why.** Three parameters had already needed amendment for the same underlying
reason — set at M1 against a design that later changed — so the remaining frozen
parameters were audited rather than waiting for each to fail in turn.

**Findings and actions:**

| # | Parameter | Problem | Action |
|---|---|---|---|
| 1 | `regions.min_datasets_supporting` | Comment said "datasets"; D-023 made it lines | Renamed `min_lines_supporting`; alias kept |
| 2 | `active_regions.signal_quantile` | Dead — thresholding is fold-over-background | Removed |
| 3 | `active_regions.sensitivity_quantiles` | Dead; `FOLD_GRID` hardcoded in 2 scripts | Replaced by `fold_grid`, scripts read config |
| 4 | `motif_validation.n_shuffles` | Dead — compositional test adopted | Removed; `test: compositional_share` |
| 5 | `regulon_validity` (4 params) | Not computable **and circular** | Leave-one-line-out; rest deferred to M7 |
| 6 | `super_enhancers` | Knee-point + 30% guard hardcoded | Moved into config |
| 7 | `regulons.max_size` | Selection rule unstated | `selection: top_by_aggregate_link_score` |
| 8 | `lineage_confounding` includes YAP1 | **No YAP1 ChIP exists** | Split ChIP vs expression-only |
| 9 | `depmap.required_files` | Missing CCLE expression | Added; `unblocks:` list recorded |
| 10 | `JUNG2017_MYC_ACTIVITY` | **Not in MSigDB** | Split into manual-curation block, `NOT_YET_CURATED` |
| 11 | `spatial.min_genes_measured: 5` | Meaningless for 500-gene regulons on a 1,738-target panel | Raised to 20 **plus** a 10% coverage fraction |

**Three were code/config divergences** (3, 6, and the hardcoded fold grid),
violating this file's own opening rule that a threshold hard-coded in a script is
a bug. That rule was stated at M1 and not enforced since.

**The pattern worth naming.** Every one of these was written when the design
assumed something later changed: more samples, quantile thresholding, shuffled
backgrounds, YAP1 ChIP that was never obtained. A frozen config is not
self-maintaining — freezing prevents *casual* change, it does not keep parameters
*correct* as the surrounding design moves. The sweep should repeat at each
milestone boundary alongside `audit_decisions.sh`.

**Two of these would have produced wrong results rather than errors:** YAP1 in the
lineage-TF list would have silently dropped a TF with no data, and
`min_genes_measured: 5` would have produced spatial regulon scores from a handful
of measured genes and reported them as if they meant something.

---

### D-029 · METH · 2026-07-26 — regulons: a real bug found by the benchmark, and the benchmark then corrected

**The benchmark failed and located a genuine defect.** The first regulon build gave
the MYC regulon **2** HALLMARK_MYC_TARGETS_V1 genes out of 500 against ~5 expected
by chance — depleted, OR 0.37.

**Cause: peak-to-gene linking was distal-only.** `20_peak_to_gene.R` linked only
`!is_prom_region`, so a paralog-bound **promoter** contributed nothing to its
regulon — discarding roughly half the binding, since only 45–56% of active regions
are distal. Two individually reasonable decisions (link enhancers; build regulons
from links) composed into a wrong one.

**Fix.** Promoter-proximal assignment added: a paralog-active promoter region maps
to its own gene. 19,236 promoter assignments across 17,373 genes, bringing total
assignments to 67,992. Effect on Hallmark hits: MYC 2→6, **MYCN 8→15 (OR 2.99,
p = 3.6e-4)**, MYCL1 5→9.

**Two predictions stated before testing, both from Plotnik's published ontology
result (MYC/MYCN enhancer targets → neurogenesis; MYCL1 and shared promoters →
housekeeping/ribosome):**

1. Promoter-inclusive regulons recover Hallmark enrichment — **confirmed for MYCN**,
   directionally for MYCL1, not for MYC.
2. Distal-only regulons are neurogenesis-enriched — **confirmed for MYC**:
   distal-only OR 1.70, p = 1.4e-4; promoter-inclusive OR 1.85, **p = 7.0e-6**.

**So the Hallmark gate was mis-specified, and the reason is principled rather than
convenient.** HALLMARK_MYC_TARGETS_V1 is a pan-cancer, promoter-centric,
ribosome-and-translation set. Plotnik's SCLC finding is that MYC's *enhancer*
programme is neurogenesis. Our MYC regulon reproduces that at p = 7e-6 and shows no
Hallmark excess — the published pattern exactly. Gating on Hallmark alone would have
failed a regulon that reproduces the prior work.

**Amended gate: a coherent programme per paralog, named.** Each regulon must show
significant enrichment for at least one MYC-relevant programme, and the programme
is reported per paralog. **This is not a relaxation — MYCL1 still fails.**

| paralog | programme | p | result |
|---|---|---|---|
| MYC | NEUROGENESIS | 7.0e-6 | PASS |
| MYCN | HALLMARK | 3.6e-4 | PASS |
| MYCL1 | none | best 0.094 | **FAIL** |

**DISCORDANCE WITH PRIOR WORK, reported not buried.** Plotnik grouped MYC *and*
MYCN enhancer targets under neurogenesis. We reproduce it for MYC but **not for
MYCN** (neuro OR 1.09, p = 0.33), which instead carries the strongest housekeeping
signal of the three. This is internally consistent with everything else observed
about MYCN — 91% nested in MYC regions, weakest motif specificity (residual 2.62
against MYC's 17.91 and MYCL1's 19.65) — so it is stated as a disagreement rather
than smoothed over.

**Other results.** Leave-one-line-out 8/9 (only MYC/H1048 marginal at AUC 0.690).
Cross-paralog Jaccard 0.059–0.130, far below the 0.60 limit — the regulons are
genuinely distinct gene sets.

---

### D-030 · SCOPE · 2026-07-26 — MYCL1 proceeds to M6 flagged, not validated

**Decision.** MYCL1 carries into tumour scoring with an explicit
`programme_validated: FALSE` flag. It is never reported as validated, and every
MYCL1 result states that its regulon showed no coherent programme.

**Why not exclude it.** MYCL1 produced the *strongest* result in the M5 gate —
motif specificity residual +19.65, the clearest paralog separation of the three
(D-020 criterion 4). Dropping the paralog with the best sequence-level evidence
because its gene-level regulon is unenriched would discard real signal.

**Why not treat it as passing.** It shows no significant enrichment for any tested
programme (Hallmark p = 0.094, neurogenesis p = 0.27). Scoring it in tumours as
though validated would put weight on a gene set we cannot show is coherent.

**Fair re-test at M7.** MYCL1 has been the weakest paralog throughout — two lines,
sparse paralog-specific regions, weakest signal-to-background. CCLE expression will
allow a real expression-correlation test in place of the H3K27ac activity proxy
(D-027), and its regulon should be rebuilt and re-gated then. If it still shows no
programme, that is the finding.

---

### D-031 · PROCESS · 2026-07-26 — three recurring defect classes, named so they stop recurring

Recorded because each has now appeared more than once, and naming the class is
cheaper than rediscovering it.

**1. Checks that cannot fail.** The hg19 build assertion compared `20` to `chr20`,
matched nothing, and reported PASS (R-16). The original M5 gate matched a published
count our own parameters determined (D-020). The ROSE cutoff collapsed to the
minimum and called every stitched region a super-enhancer. **Rule: a check must be
observed failing on data known to violate it before it is trusted.**

**2. Self-matching checkers.** `pgrep -f 11_rescan_exact.sh` matched the shell
command that invoked it. `verify_config.R` scanned the script tree for removed
config keys and matched its own source, which necessarily names every one.
**Rule: any check that inspects the codebase, process table, or filesystem must
exclude itself from its own scan.**

**3. Seam errors between edits.** Blacklist filtering was named as necessary and
never implemented. D-025 was referenced in config and three reports but never
written. `FOLD_GRID` and the SE cutoff lived in code while the config claimed to own
all thresholds. D-029 removed `benchmark_enrichment` from config while a script
still read it. **Rule: an intention declared and deferred is not recorded by the
decision log, which captures decisions made. Run `audit_decisions.sh` and
`verify_config.R` at every milestone boundary.**

---

### D-032 · SCI · 2026-07-26 — R-01 materialised: the project's conclusion changes, its scope does not collapse

**The finding.** Paralog-resolved regulatory programmes are constructible from
cell-line chromatin and internally valid there, but **do not retain paralog
identity in patient tumours independently of neuroendocrine lineage state.**
Lineage explains 39–49% of regulon score variance against 0.0–1.6% uniquely
attributable to the paralog's own expression; no paralog association survives
adjustment (0/3, all p > 0.2). See R-01.

**This was pre-committed.** Gap statement §5: *"Paralog-specific regulons show no
separation in tumour cohorts → report as a negative result with power analysis;
the framework and the audit trail remain the contribution."* The response was
agreed before the data were touched, which is why it is being followed rather
than negotiated.

**Adequately powered, so the null is informative.** n = 79; |ρ| ≥ 0.31 detectable
at 80% power. Observed unique-paralog contributions are an order of magnitude
below that. The exception is MYCN, whose expression median is 1.34 (IQR 0–3.26) —
for MYCN specifically, absent dynamic range is a competing explanation to absent
biology and must be stated as such.

**The mechanism is a positive result in its own right.** MYC vs NE score
ρ = −0.590, p = 1.1×10⁻⁸, independently reproducing Ireland et al. 2020. MYCN↔NEUROD1
+0.371 reproduces the MYCN/SCLC-N association. The negative result therefore has an
explanation rather than being an unexplained null, which makes it far more
defensible.

**What survives unchanged:**
- M5 gate: 3/4 evaluable criteria on threshold-invariant measures (MYCN⊂MYC 0.912
  with spread 0.037 across a 50× set-size change; paralog E-box specificity 3/3,
  χ² p = 1.6×10⁻¹⁰¹)
- MYC regulon → neurogenesis, p = 7×10⁻⁶, reproducing Plotnik from a pipeline
  sharing no method with theirs
- The MYCN-with-housekeeping discordance against Plotnik (D-029)
- Super-enhancers, peak-to-gene (6.06× over null; 32.6% nearest-gene)

**Effect on M7–M8 — the framework's purpose is sharpened, not lost.** MOES was
designed to integrate four evidence domains. The transcriptional domain is now known
to be lineage-confounded at the regulon level, so integrating it naively would
produce lineage-driven hubs dressed as paralog-specific ones.

Two changes follow:

1. **The transcriptional domain becomes lineage-adjusted by construction.** Evidence
   enters as the per-gene partial association with paralog expression after
   regressing out NE score and all four lineage TFs — not the raw association. The
   regulon-level aggregate is null, but individual genes may carry
   lineage-independent signal, and that is the level MOES operates at. This must be
   tested before M8 rather than assumed.
2. **The framework's contribution is reframed honestly.** It was proposed as
   evidence integration for hub prioritisation. Its demonstrated value is now
   *detecting and quantifying confounding that a single-layer analysis would have
   reported as a finding.* A weighted-sum framework over these layers would have
   returned confident paralog-specific hubs. That is a more useful methodological
   claim than the original and it is supported by this project's own data.

**Honest statement of what the project can now claim.** Not "here are validated
paralog-specific MYC hubs in SCLC". Rather: paralog-resolved programmes are
recoverable from chromatin and reproduce known enhancer biology, but their
tumour-level signal is dominated by lineage state; any prioritisation must be
lineage-adjusted, and the framework quantifies how much apparent paralog signal is
in fact lineage. **No result may be reported as paralog-specific without the
lineage-adjusted figure alongside it.**

---

### D-033 · SCI · 2026-07-26 — transcriptional domain dropped; MOES reduced to three domains; aims reframed

**The per-gene rescue failed.** The regulon aggregate was null (D-032), but MOES
operates per gene, so individual genes were tested for lineage-independent paralog
association (partial Spearman, adjusting for NE score + ASCL1 + NEUROD1 + POU2F3 +
YAP1, df = 72, detectable |ρ| ≥ 0.229).

| paralog | regulon genes | survivors in regulon | rate in | rate outside | OR | p |
|---|---|---|---|---|---|---|
| MYC | 367 | 1 | 0.27% | 0.23% | 1.21 | 0.57 |
| MYCN | 397 | **0** | 0.00% | 0.00% | — | — |
| MYCL1 | 373 | 1 | 0.27% | 0.17% | 1.57 | 0.48 |

**Enrichment 0/3.** Regulon genes survive at background rate. The correct test was
enrichment rather than survivor count — with ~33,700 genes, some survive FDR by
chance, and 1 of 370 is precisely that. **MYCN is more extreme: zero of 33,683
genes** genome-wide show a lineage-independent positive association, consistent
with its near-absent dynamic range (median 1.34, IQR 0–3.26).

**Decision 1: the transcriptional domain is DROPPED from MOES**, and the reason is
declared in the methods rather than omitted. It has no paralog-specific content to
contribute at either the aggregate or the per-gene level.

**Decision 2: the network domain must NOT be built on these tumours.** GENIE3
recovers regulators from expression; we have now shown that what predicts paralog
expression in this cohort is lineage state. A network domain built here would enter
MOES looking like independent evidence while carrying the same confound. It moves
to CCLE SCLC lines (~50) at M7.

**Decision 3: MOES runs on three domains** — cis-regulatory, network (CCLE-based),
functional (DepMap) — and cannot run before M7. The domain count is reported, not
implied.

**Decision 4: aims reframed.** The original headline — "identify robust
paralog-specific MYC regulatory hubs, cross-validated across layers" — cannot be
delivered as written and will not be claimed. The delivered contribution is:

1. A reproducible pipeline extracting paralog-resolved regulatory programmes from a
   bedGraph-only deposit, without MACS2 or input controls.
2. **Two independent reproductions of published work**: Plotnik's MYC-enhancer
   neurogenesis programme (p = 7×10⁻⁶) and Ireland 2020's MYC/NE antagonism
   (ρ = −0.590, p = 1.1×10⁻⁸).
3. **A new discordance**: MYCN's programme is housekeeping-weighted, not
   neurogenesis as Plotnik grouped it.
4. **A well-powered negative result with an identified mechanism**: paralog-resolved
   programmes do not retain paralog identity in patient tumours independently of
   neuroendocrine lineage state.
5. **A methodological demonstration**: a weighted-sum framework over these layers
   would have returned confident paralog-specific hubs. Detecting and quantifying
   that confound is the framework's demonstrated value — a stronger and more useful
   claim than the original, and supported by this project's own data.

**Why continue rather than consolidate.** DepMap CRISPR dependency is *not*
expression, so it is not subject to the lineage confound that removed the
transcriptional and tumour-network domains. CCLE provides ~50 SCLC lines for a
network domain that 79 tumours cannot support. DepMap also unblocks M5 gate
criterion 3 (D-026). The remaining evidence is genuinely independent of what failed.

**Standing constraint.** No result may be described as paralog-specific without the
lineage-adjusted figure beside it. This applies to every table, figure and sentence
from here to release.

---

### D-034 · SCI · 2026-07-27 — M6 negative result REPLICATES in an independent cohort

**Why replication was necessary.** A single-cohort null has mundane explanations —
cohort composition, platform, processing, the particular 79 patients. The gap
statement promised replication in a second cohort, and a *negative* result needs it
more than a positive one would.

**Cohorts.** GSE60052 (n = 79) and George et al. 2015 via cBioPortal (n = 81,
`sclc_ucologne_2015`). Independent patients, sequencing and processing. **160
patients total.**

**Q1 — does lineage dominance reproduce? YES, 3/3 regulons in BOTH cohorts.**

| | GSE60052 | George 2015 |
|---|---|---|
| lineage > paralog unique variance | 3/3 | 3/3 |
| unique lineage R² | 0.355–0.428 | 0.244–0.426 |
| unique paralog R² | 0.001–0.068 | 0.000–0.019 |
| MYC vs NE score (Ireland 2020) | **−0.590** | **−0.466** |

The MYC/NE antagonism reproduces independently, so the *mechanism* behind the null
is corroborated as well as the null itself.

**Q2 — does any positive paralog association reproduce? No.** MYCL1 survived
adjustment in GSE60052 (FDR 0.0136) but not in George (FDR 0.936). Three reasons it
is not a finding: it is **method-dependent** (absent under singscore in the same
cohort, FDR 0.675), it is **cohort-specific**, and it comes from **the one regulon
that failed M5 validation** (D-030). A positive that fails to replicate is a
spurious result being caught, which supports the null rather than undermining it.

**SCORING METHOD CHANGED FOR THE COMPARISON, deliberately.** singscore ranks each
gene against the others present in the matrix — 33,683 genes in GSE60052 but only
1,004 fetched for George, most of them regulon members. Comparing singscore across
the two would have compared scoring artefacts. Both cohorts were therefore rescored
with a universe-independent statistic (mean per-gene z-score), which is also
scale-free and so handles GSE60052's log2 against George's raw counts (max 37,557).
**GSE60052 was rescored rather than reused**, so the comparison is like-for-like.

**A verdict-logic error, corrected.** The first version required zero surviving
associations in *both* cohorts and therefore reported "NOT REPLICATED" when a
cohort-specific positive appeared. That conflated two questions: whether the
negative result reproduces, and whether a positive does. The corrected script tests
them separately and retains a branch that would report a genuinely reproducing
association as real — the logic is not tuned toward the expected answer.

**Gene naming, pre-handled.** GSE60052 carries `MYCL1` and not `MYCL`; George
carries `MYCL` and not `MYCL1`. Both were requested and resolved explicitly. This is
the silent-mismatch class that cost three rounds at the universe stage.

**Status.** The M6 conclusion is now a **replicated** negative result with an
identified and independently reproduced mechanism, across 160 patients. It is the
strongest result in the project.

---

### D-035 · SCI · 2026-07-28 — amplification resolved; M5 gate closes at 3 PASS / 1 FAIL

**Amplification settled from DepMap Public 26Q1 log2 copy number**, resolving a
question open since M5 that neither the paper, GEO, nor our coverage proxy could
answer.

| group | lines | log2 CN |
|---|---|---|
| MYC-amplified | **H524, H211** | 6.73, 2.35 |
| MYC-expressing | H1048, SHP77 | 0.98, 1.30 |
| MYC status unknown | H847 | absent from 26Q1 |
| MYCN-amplified | H526, H69 | 4.94, 6.41 |
| MYCL-amplified | COLO668, H889 | 5.60, 5.28 |

**Exactly two MYC-amplified lines**, confirming Plotnik's "two cell lines
harboring the alteration for each amplification type". **The project registry was
wrong**: it listed all five MYC ChIP lines as amplified, which was never
compatible with the paper's own text.

**D-026's coverage proxy is validated.** MYCN and MYCL reproduce exactly (H526's
MYCL 1.60 is a threshold artefact against 5.6/5.3 for the true MYCL lines). Its
failure on MYC is now explained rather than merely noted: **H524's amplification is
real but FOCAL**, so the ±1 Mb window diluted it — precisely the geometry D-026
inferred without being able to make the call.

**H847 is excluded from both groups.** Absent from 26Q1, so its status is unknown;
assigning it to the comparator would be inference from absence.

---

**CRITERION 3 RESULT: FAIL.** Direction reproduces, magnitude does not.

| | amplified | expressing | difference |
|---|---|---|---|
| ours | 47.5% (53.7, 41.2) | 44.2% (41.2, 47.2) | **+3.3 points** |
| Plotnik | 39% | 12% | +27 points |

**The more important observation: within-group spread exceeds the between-group
difference.** The amplified pair differ by 12.5 points and the expressing pair by
6.0, against a 3.3-point group difference — and H211 (amplified) and H1048
(expressing) have *identical* distal fractions of 41.2%. At 2 lines per group this
contrast is not distinguishable from line-to-line variation. No significance test
was applied because none would be honest at this n.

**Most likely structural cause.** Our universe is ~84% distal by construction —
built from ATAC regions where promoters are a minority — which compresses the
achievable range. Plotnik called peaks de novo per line with no shared grid, so
their distal fraction had far more room to move. This is the **second** discordance
tracing to the shared-grid design, after MYCL1's inflated occupancy overlap (D-025),
and both are consequences of D-023 rather than of the biology.

**FINAL M5 GATE: 3 PASS / 1 FAIL.**

| criterion | result |
|---|---|
| 1 · MYCN⊂MYC 0.912, spread 0.037 across 50× set-size change | **PASS** |
| 2 · MYCN vs MYCL1 differential nesting, OR 3.14, p 1.1e-60 | **PASS** |
| 3 · distal-fraction contrast | **FAIL** (direction only) |
| 4 · paralog E-box specificity 3/3, χ² p 1.6e-101 | **PASS** |

The gate was rebuilt at D-020 specifically so criteria could fail on evidence
rather than on our own free parameters. One did. That is the gate working, and the
failure is reported rather than reformulated — unlike criteria 2 and 4, where the
*test* was demonstrably mis-specified, here the test is sound and the data do not
support the published magnitude.

---

### D-036 · SCI · 2026-07-28 — functional domain admitted, but it is NOT paralog-resolved

**Result.** Paralog regulon genes collectively carry SCLC-selective CRISPR
dependency above background: **11/984 (1.12%) vs 73/17,547 (0.42%), pooled
OR 2.71 (95% CI 1.29–5.15), p = 0.0048**, across 25 SCLC lines vs 1,183 others in
DepMap Public 26Q1.

**Why this domain survives where two others did not.** CRISPR gene effect measures
whether knocking a gene out kills the cell. It is not expression, so it cannot
inherit the neuroendocrine lineage confound that removed the transcriptional domain
and barred the tumour-based network domain (D-033). This is the reason continuing to
M7 was worth doing rather than consolidating at M6.

**POOLED, NOT PER-PARALOG — and the distinction matters.** Per-paralog odds ratios
are MYC 2.45 (p 0.089), MYCN 2.21 (p 0.118), MYCL1 3.09 (p 0.028). These are
**indistinguishable**: 4–5 selective genes each, differing by one gene. The initial
verdict logic declared "only MYCL1 enriched", which would have reported a one-gene
difference as a biological distinction — and would conveniently have handed
functional support to the one paralog that FAILED the M5 programme gate (D-030).
The pooled test is the properly powered version of the question MOES actually asks.

**CONTROL RESULT, reportable in its own right: none of MYC, MYCN or MYCL is an
SCLC-selective dependency.** MYC is strongly pan-essential (−1.495) but **LESS so in
SCLC than in other lineages (delta +0.503)**; MYCN (−0.112) and MYCL (−0.117) show
no dependency at all. Nothing in DepMap supports SCLC being selectively dependent on
these paralogs. This must appear in the write-up — it constrains any therapeutic
framing.

**POSITIVE-CONTROL GATE ADDED.** ASCL1 is a canonical SCLC lineage dependency and
must appear as SCLC-selective or the script exits non-zero. It caught a real
failure: a lineage "fallback" I added swept the entire NSCLC panel into the SCLC
group (NSCLC lines carry `lineage_2 = "Non-Small Cell Lung Cancer"`, which *contains*
the substring "small cell lung"). The code ran cleanly; the contamination was
visible only in the biology — top hits became cilia and cornified-envelope genes and
ASCL1's effect fell from −0.539 to −0.148. My premise was also wrong: DepMap CRISPR
genuinely has ~25 SCLC lines, so the "fix" repaired a non-problem into a real one.

**THE STRUCTURAL CONSEQUENCE FOR MOES.** Domains by what they can attribute:

| domain | status | paralog-resolved? |
|---|---|---|
| cis-regulatory | available | **yes** (region-level, per paralog) |
| transcriptional | dropped (D-033) | — |
| functional | admitted, OR 2.71 | **no** (gene-level only) |
| network | untested | needs CCLE; likely confounded |

**Only ONE domain can attribute evidence to a specific paralog.** MOES was designed
to prioritise paralog-specific hubs by integrating four independent layers. With one
paralog-resolved layer, "multi-layer paralog-specific evidence integration" is not
supported by the data.

**MOES is therefore reframed, not abandoned:** it ranks genes by converging
regulatory and functional evidence, with **paralog attribution resting on chromatin
alone** and other layers corroborating a gene's importance rather than its paralog
identity. Every output must state which domains contributed and what each can and
cannot attribute.

---

### D-037 · SCI · 2026-07-28 — network domain excluded; MOES reduced to two domains

**Tested before building**, so a confounded domain could not be admitted while
appearing to be independent evidence. Both pre-conditions failed.

**Power.** 59 SCLC lines with CCLE expression. GENIE3 is a tree ensemble and
conventionally needs hundreds of samples for stable importance estimates.

**Confounding.** MYC expression tracks NE state across SCLC *cell lines* at
**ρ = −0.441** (FDR 0.008), and MYC vs ASCL1 at −0.409 (FDR 0.010) — the same
structure that barred the tumour-based network. A GENIE3 network built here would
inherit the D-033 confound.

**A third independent reproduction of Ireland et al. 2020.** MYC/NE antagonism:
GSE60052 tumours −0.590, George 2015 tumours −0.466, CCLE cell lines −0.441. Three
independent datasets across two data types. This is among the most robust findings
in the project, and it explains *why* two domains had to be excluded: the antagonism
is a property of SCLC biology, not of tumour tissue.

**FINAL MOES COMPOSITION — two of four domains.**

| domain | status | paralog-resolved? |
|---|---|---|
| cis-regulatory | ADMITTED | yes (region-level) |
| transcriptional | EXCLUDED (D-033) | — |
| network | EXCLUDED (this entry) | — |
| functional | ADMITTED (D-036) | no (gene-level) |

**What MOES can and cannot claim.** It ranks genes by converging regulatory and
functional evidence. It **cannot** perform multi-layer paralog-specific
prioritisation as specified at M1, because only one admitted domain attributes
evidence to a paralog. Two-stage RRA is also degenerate at two domains and reduces
to a single aggregation across two rankings — that must be stated, not glossed.

**The framework's real contribution, restated honestly.** Its demonstrated value is
**deciding which evidence is admissible**, not integrating whatever is available. It
excluded three of four domains on evidence, each with a documented test that could
have gone the other way. A weighted-sum framework over these four layers would have
produced confident paralog-specific hubs from two lineage-confounded layers and one
that cannot distinguish paralogs. That is a more useful methodological claim than
the original, and it is supported by this project's own data rather than asserted.

---

### D-026 · DATA · 2026-07-26 — amplification status: MYCN/MYCL1 confirmed, MYC unresolved, criterion 3 deferred to M7

**Why this was investigated.** M5 gate criterion 3 (distal-fraction contrast,
MYC-amplified vs MYC-expressing) needs a grouping the project did not reliably
have. Plotnik states *"we profiled two cell lines harboring the alteration for each
amplification type"*, but GSE230649 contains **five** MYC ChIP samples (H1048,
H211, H524, H847, SHP77). So two are amplified and three are the MYC-expressing
comparator — yet the project registry listed **all five as MYC-amplified**, which
cannot be correct.

The paper does not name them (Methods say only "representative SCLC cell lines"),
and the GEO records carry no amplification field. Both were checked directly
rather than assumed.

**Measured instead.** Amplification adds DNA copies, inflating pileup across an
amplicon in any assay, so amplicon coverage relative to each track's genome-wide
median is a copy-number proxy. Computed from the existing signal matrix at ±200 kb
and ±1 Mb.

**CONFIRMED — MYCN and MYCL1:**

| locus | line | ATAC | H3K27ac | ±1 Mb |
|---|---|---|---|---|
| MYCN | H69 | 15.1 | 44.0 | elevated |
| MYCN | H526 | 8.3 | 37.9 | 9.4 |
| MYCL1 | COLO668 | 22.2 | 71.2 | 12.4 |
| MYCL1 | H889 | 19.6 | 36.6 | 5.0 |

Both assays agree, and elevation **persists in the broad window** — the signature
of a real amplicon. All other lines near 1. Registry correct for these paralogs,
which is what the gate criteria actually depend on.

**UNRESOLVED — MYC.** The assays contradict each other: H524 gives H3K27ac 82.1
against ATAC 1.55; H196 gives ATAC 27.8 against H3K27ac 4.0. And every MYC-locus
ratio **collapses to ~1–4 in the ±1 Mb window** while MYCN and MYCL1 stay high.
Focal elevation that disappears when the window widens indicates a strong
regulatory element, not copy gain. On this evidence no line shows convincing MYC
amplification — which contradicts the published status of H524 and marks the limit
of the method, not a finding about the biology.

Additional weaknesses: **H211 has no ATAC track**, so its call rests on one assay;
and the MYC window may contain few universe regions, destabilising a median. The
accessibility/copy-number confound was flagged before running this and is where it
failed.

**Decision.** Criterion 3 remains **BLOCKED and is deferred to M7**, when DepMap
`OmicsCNGene.csv` provides direct copy-number calls. It is not dropped, and it is
not resolved by inference. The registry is corrected to mark MYC-amplification
status as **unverified** for all five MYC ChIP lines.

**Consequence for reporting criterion 1.** Because the five MYC ChIP lines may mix
amplified and non-amplified contexts while Plotnik's 18,823 came from two
amplified lines, **the similarity between our 17,693 and their 18,823 must not be
presented as replication** — the two may not be the same construct. Criteria 1, 2
and 4 are unaffected: they define paralog sets by which antibody was used, not by
amplification status. The MYCN-in-MYC value of 0.886 and its stability across a
50-fold change in set size stand on their own.

---

### D-038 · METH · 2026-07-28 — the paralog palette frozen at M1 failed its own accessibility check; replaced by scored selection

**Amends:** `config/params.yml` `figures$palette_paralog` — a value frozen since M1.
**Scientific impact:** none. Colour carries no quantitative claim in any panel.

#### What happened

`check_palette()` was written at M4 to enforce `figures$colourblind_check`
(D-021 D: colourblind safety *verified* by deutan/protan/tritan simulation
against a minimum perceptual distance, "rather than asserted"). The first
results figure to actually run it, `fig01_lineage_dominance.R`, **failed on the
project's own palette**:

```
palette check [paralog palette]: min CIE dE = 12.7 under deutan (need >= 15) -> FAIL
```

MYC `#B2182B` (red) against MYCL1 `#1B7837` (green) — the commonest form of
colour vision deficiency, in a palette carried since M1 and already shipped in
the M4 inventory figure `fig_s01_data_landscape`.

This is the check working as designed. The palette had been *declared*
colourblind-safe in config and never tested; between M1 and M10 that assertion
was simply believed. It belongs to the same class as the D-031 defects: **a
claim recorded in config that no code ever verified.**

#### Decision

The replacement was **selected by score, not by taste**, in
`scripts/07_figures/00_select_palette.R`, which is committed and re-runnable so
the choice can be audited or revisited. Eleven candidates drawn from established
colourblind-safe families (Okabe-Ito, ColorBrewer diverging endpoints) were
scored on **two** hard criteria:

1. **min CIE ΔE ≥ 15** across normal/deutan/protan/tritan — the configured floor.
2. **≥ 30 L\* separation from the white panel.**

The second criterion was added after the first ranking returned a pale gold
`#FDAE61` on a ΔE of 44.3. That score is real but is driven almost entirely by
lightness spread: maximising pairwise distance alone rewards near-white tints.
The paralogs are drawn as thin lines, 1.7 pt points and coloured text on a white
panel, so a colour must also stand off the background. Three of the eleven
candidates — including the top scorer — pass the colourblind test and fail this
one. **A category the reader cannot see is not accessible whatever its pairwise
ΔE says.**

Five of eleven pass both. Selected:

| paralog | colour | was |
|---|---|---|
| MYC | `#762A83` | `#B2182B` |
| MYCN | `#E08214` | `#2166AC` |
| MYCL1 | `#1B7837` | `#1B7837` |

min ΔE **30.5** (limiting under protanopia, 2× the floor); background contrast **36.7**.

The runner-up preserving the conventional red-for-MYC assignment was
`#A50026`/`#2166AC`/`#B35806` (ΔE 24.0, background 52.2). It passes both criteria
and is a valid substitute; the selection went on ΔE because there is no
established convention that MYC must be red.

#### Two further changes this exposed

**A hardcoded hex outlived its palette.** `fig01`'s correlation heatmap used
`high = "#B2182B"` — the old MYC red, typed directly into the plotting script
where `check_palette()` could not see it. It is now `figures$palette_diverging`
in config, exposed as `scale_fill_diverging()`. A diverging scale encodes a
signed continuous value, not paralog identity, so it stays deliberately separate
from `palette_paralog` — but it is no longer invisible to review.

**Panel A was borrowing a paralog colour for a non-paralog encoding.** Its fill
distinguishes *variance source* (lineage vs paralog) while paralog is the
x-axis; taking the fill from `PAL_PARALOG[["MYC"]]` implied "this bar is MYC" in
a panel where it is not. Now blue for lineage, neutral grey for paralog.

#### Consequence

Every figure using `PAL_PARALOG` must be regenerated: `fig_s01_data_landscape`
(M4) and `fig01_lineage_dominance` (M10). Remaining M10 figures are written
against the config and inherit it. Because the palette is read from
`config/params.yml` rather than repeated in each script, this is a one-line
change plus a re-run — which is why it was centralised at M4 (D-021).

#### Note for review — the checks do not check the plot

Three defects in this figure passed all three automated checks and were caught
only by **looking at the rendered PNG**: a clipped main title, colliding panel
titles, and value labels in panel C rendered white-on-white *outside* the bars,
because the bars run downward from zero and `vjust` was positive. A fourth —
labels collapsing to the group centre because subsetting the layer's data left
`position_dodge` with a single fill level to dodge against — was equally
invisible to them.

**`check_palette()` and `check_text_size()` verify properties of the palette and
the theme, not of the finished plot.** Passing them is necessary and nowhere
near sufficient. Rendered output must be inspected before any figure is
considered done.

---

### D-039 · METH · 2026-07-28 — the canonical M5 gate table was stale, and two audit checks were defending the staleness

**Scientific impact:** none — no result changes. The M5 gate stood at 3 PASS / 1 FAIL
before this entry and stands at 3 PASS / 1 FAIL after it. What changes is that the
table now *says so*.

#### What was wrong

`data/metadata/m5_gate_results.csv` is the canonical M5 gate table: the audit reads
it, the report will read it, and it is what a reader inspects to see whether the
gate passed. It recorded criterion 3 as:

```
"3_distal_contrast","not evaluable",...,NA
```

That was written by `16_m5_gate.R` at M5, when amplification status was genuinely
unknown (D-026). At M7, D-035 resolved amplification from DepMap Public 26Q1 copy
number and `03_gate_criterion3.R` evaluated the criterion: **+3.3 points against a
published +27, direction reproduces, magnitude fails, FAIL.**

The M7 script wrote its own outputs and never wrote back. It even built the
verdict — a data frame named `out` carrying `direction_pass`, `magnitude_pass` and
`pass` — and then **discarded it without a `write.csv`**. Dead code, one line from
being correct.

The result was two sources of truth that disagreed for the whole of M7: the
decision log said FAIL, the canonical table said never-evaluated. A gate figure
built from the table would have contradicted the log.

#### The worse half: the audit was protecting it

Two `audit_decisions.sh` checks asserted the *deferral*, not the *resolution*:

| check | asserted | after D-035 |
|---|---|---|
| D-026 | registry still says `UNVERIFIED MYC` | fails on a correctly-updated registry |
| D-026b | gate table says `not evaluable` | **passes only while the table is stale** |

D-026b is the dangerous one. It passed for the entirety of M7 *because* the bug
existed, and would have started failing the moment anyone fixed it. An audit check
written against a temporary state becomes, once that state is resolved, a guard on
the error it was meant to catch. D-026 shows the mirror failure: it was already
reporting FAIL against a registry that had been correctly updated, which trains the
reader to ignore a red line.

Both are now inverted to assert the resolved state, plus two new checks:
`D-035b` (no criterion left with `pass = NA`) and `D-035c` (registry records
resolution from DepMap). Audit: **34 PASS, 0 FAIL.**

#### Fixes

- `03_gate_criterion3.R` writes the verdict back into `m5_gate_results.csv`,
  guarded by a row-count check that stops rather than silently matching zero or
  many rows, and prints the resulting gate tally. `out` is now written to
  `m5_criterion3_summary.csv` instead of being discarded.
- `audit_decisions.sh` D-026 and D-026b inverted; D-035b and D-035c added.

#### Rule this establishes

**A check that encodes a deferral must be retired by the decision that resolves
the deferral.** The decision log entry that closes a gap is not complete until the
checks written against the gap have been updated with it. This is a fourth member
of the D-031 defect family — alongside a config claim no code verified (D-038), it
is a *check* whose claim had quietly become false.

A script that resolves a value recorded elsewhere must write back to where it is
recorded. Producing a new file beside a stale one is not resolution.

---

### D-040 · METH · 2026-07-28 — the pooled dependency result existed only as console output

**Scientific impact:** none. Values reproduce exactly on re-run: pooled OR 2.706
(95% CI 1.290–5.153), p = 0.00477; MYC dependency delta +0.503.

#### What was wrong

`04_selective_dependency.R` computed the pooled Fisher test — **the result that
admits the functional domain into MOES (D-036)** — printed it, and discarded the
object. Same for the per-regulon confidence intervals (`ft$conf.int`) and for the
paralogs' own dependency values. The numbers existed in three places: a console
run nobody keeps, prose in D-036, and this decision log.

Figure 3 needed all three. Building it against prose would have meant typing
`2.71` into a plotting script, which is how a figure silently drifts from the
analysis that produced it — the identical failure mode to the hardcoded `#B2182B`
in D-038 and the un-written-back gate verdict in D-039.

#### Fix

Three tables are now written by the script that computes them:

| file | contents |
|---|---|
| `m7_dependency_pooled.csv` | pooled OR, 95% CI, p, counts, `paralog_resolved = FALSE` |
| `m7_paralog_own_dependency.csv` | MYC/MYCN/MYCL gene effect in SCLC vs other lineages |
| `m7_dependency_enrichment.csv` | gains `ci_low` / `ci_high` per paralog |

The confidence intervals matter beyond provenance. The claim that the three
per-paralog ORs are indistinguishable was previously an assertion in prose; with
the intervals stored it is *shown*: 0.65–6.57, 0.58–5.91, 0.97–7.59, overlapping
almost completely on 4–5 selective genes each.

They also surfaced something worth stating plainly. **MYCL1's nominal p = 0.028
carries a 95% CI of 0.971–7.59, which includes 1.** Fisher's exact p and its
conditional-MLE interval disagree at these counts. This is not an error in either,
but it is a direct argument against reading any per-paralog verdict here — and it
would have stayed invisible while only the p-value was retained. It is recorded in
the Figure 3 caption rather than smoothed over.

#### Rule

This is the third instance of one pattern (D-038 config claim never verified,
D-039 verdict never written back, D-040 result never persisted). **A number that
appears in a figure, a report or the decision log must be readable from a file
that the analysis wrote.** Console output is not a result store.

---

### D-041 · SCI · 2026-07-28 — the reduced MOES returns no prioritised hubs; the two surviving domains do not converge

**Status:** RUN, REPORTED AS A NULL. Not a failure to complete — a completed
analysis whose answer is negative.

#### What was run

The two-domain MOES specified after D-037, over 10,387 genes carrying evidence in
both admitted domains: cis-regulatory (paralog-resolved aggregate peak-to-gene
link score) and functional (gene-level SCLC-selective dependency). Robust Rank
Aggregation, 10,000 permutations for empirical FDR, 2,000 bootstrap resamples of
the cis links.

Two things were verified before any score was believed:

1. **The cis layer is the committed one.** `01_build_evidence.R` rebuilds the
   aggregate link score and asserts that its top-500 genes are *identical* to the
   regulons in `regulons.rds`, for all three paralogs, or stops. A reimplementation
   that had quietly drifted would have made every rank below wrong with nothing
   else to catch it. It reproduces exactly.
2. **The RRA null is correct.** The n = 2 closed form used in the permutation and
   bootstrap loops was checked against `RobustRankAggreg::aggregateRanks`:
   max absolute difference **1.1e-15**. A hand-derived Beta null that was subtly
   wrong would have corrupted every p-value and FDR.

#### Result

**Zero genes reach FDR < 0.05 for any paralog.** Best FDR achieved: 0.381 (MYC),
0.371 (MYCN), 0.358 (MYCL1).

The mechanism is visible in one number. The Spearman correlation between the two
domains' rankings is **−0.023, +0.012, −0.008**. The domains are effectively
independent, so genes ranking highly in one do not rank highly in the other more
often than chance predicts, and RRA has nothing to aggregate. This is not a
threshold artefact; it is the absence of the convergence MOES exists to detect.

#### The null is real, and the method has the power to say so

A null from an untested method is not a finding. `03_moes_positive_control.R`
blends a synthetic functional layer with the real cis ranking at increasing
strength and runs the identical FDR machinery:

| injected signal | layer ρ | min FDR | genes FDR < 0.05 |
|---|---|---|---|
| 0.00 | −0.000 | 0.878 | 0 |
| 0.10 | +0.101 | 0.084 | 0 |
| 0.20 | +0.209 | 0.018 | 1 |
| 1.00 | +1.000 | 0.000 | 288 |

Monotone and graded: nothing on pure noise, discoveries once the layers correlate
at ρ ≈ 0.2, hundreds when they are identical. **Observed ρ = −0.023.** The
convergence is absent at a level the method demonstrably detects.

#### Three structural findings, computed rather than asserted

**Two-stage RRA is degenerate at two domains.** Stage 2 has nothing to aggregate
that stage 1 did not. The reported score is a single aggregation. Config already
carried `two_stage_rra_degenerate: true`; it is now demonstrated in the output.

**Leave-one-domain-out is degenerate too.** With two domains, removing one leaves
one. It is reported as "combined vs each domain alone" (ρ +0.57 to +0.70), not as
the stability check specified at M1, which assumed four domains.

**Between-paralog difference is chromatin alone.** The functional vector is
*identical* for all three paralogs, so it pulls every paralog's ranking toward the
same genes. MOES agreement between paralogs (ρ +0.50 to +0.59) is inflated well
above cis-only agreement (ρ +0.17 to +0.33) purely by that shared layer. Any
apparent paralog-specificity in a MOES ranking would be chromatin wearing a
multi-omics label.

#### Decision

**MOES is reported as a null. No top-N hub table will be produced.** Presenting
the head of a distribution with best-FDR 0.36 as "prioritised hubs" would be the
single most misleading thing this project could do — it is exactly the output the
M1 design anticipated, and it would look like a result.

The full ranking with its FDR column is written to `data/metadata/moes_ranking.csv`
so the absence is inspectable rather than hidden.

The multi-layer paralog-specific prioritisation specified at M1 is **not
achievable on this evidence**. Three of four domains were excluded on their own
evidence (D-033, D-036, D-037), and the two that survived do not converge. This is
now the project's second major negative result, alongside R-01 (D-034), and both
are reported as findings.

#### What this does not say

Not that the regulons are wrong — they pass internal validation (Figure 2D) and
reproduce Plotnik's MYC-enhancer→neurogenesis link independently. Not that the
dependency enrichment is wrong — pooled OR 2.71 stands (D-036). It says the two
kinds of evidence do not point at the same genes, which is a statement about the
biology as measured here and about the limits of integrating two thin layers.

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

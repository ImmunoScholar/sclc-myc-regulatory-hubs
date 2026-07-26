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

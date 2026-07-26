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

### D-012 · ENV · 2026-07-26 — GSVA outstanding; needs one system library

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

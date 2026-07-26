# Research Journal

A dated, running record of what was actually done, what was found, what broke,
and what it cost. Written as work happens, not reconstructed afterwards.
Negative results and dead ends are recorded here deliberately — they are part of
the scientific record and they stop the same wall being walked into twice.

Entries are newest-last.

---

## 2026-07-26 · M0–M1 · Scoping and specification

**Done.** Literature map, dataset inventory, gap identification, novelty
validation, analysis architecture, dependency inventory, risk log, roadmap.

**The finding that reshaped the project.** Full-text review of Plotnik et al.
2024 (Mol Cancer Res, PMID 38747975) — the paper behind GSE230649 — showed the
authors had already published paralog-specific MYC/MYCN/MYCL1 enhancer
characterisation: E-box central-dinucleotide preferences (MYC `CAGATG`, MYCN
`CACATG`, MYCL `CACCTG`), 18,823 / 4,017 / 5,688 active regions, ~84%
MYCN-in-MYC overlap. The first analytical step of the original plan was already
done, by the people who generated the data.

The gap was **re-scoped rather than defended** (D-001). Competing with a
published foundation is weaker than building on it.

**Second scope hit.** Zhang et al. 2025 (Clin Cancer Res) already links ecDNA
MYC-paralog amplification to spatial immune exclusion in SCLC. The spatial aim
was demoted from discovery to a restricted coherence check.

**Net effect on novelty grade:** HIGH → MODERATE-HIGH. Recorded honestly rather
than defended.

---

## 2026-07-26 · M2 · Blocking decisions resolved

**Single-cell layer dropped** (D-005). Synapse/HTAN cannot be authorised in this
environment, so the Chan 2021 atlas is unreachable. The pre-agreed decision rule
— drop rather than substitute CDX data — was invoked as written. Cost: patient
translation now rests on two bulk cohorts alone. Benefit: the project is now pure
R with no Python dependency.

**Hardware measured, not assumed.** WSL2 Ubuntu 24.04 (noble), R 4.6.1
"Happy Hop", Bioconductor 3.23, 6 cores, 948 GB free on the native Linux
filesystem — and **10 GiB RAM, not the 16 GiB the plan had assumed.**

That last number matters. Chromosome-by-chromosome bedGraph processing moves from
"fallback if needed" to a **hard design requirement** (new risk R-13). A
whole-genome `rtracklayer::import()` over 8.6 GB would exhaust memory. Swap
exists but is treated as a failure signal, not a safety net.

**GitHub** SSH already working (`ImmunoScholar`); no key generation needed.

---

## 2026-07-26 · M3 · Repository and environment

**Two environment traps hit and worked around — worth recording because both
produce silent, confusing failures:**

1. **MSYS2 path mangling.** Invoking WSL through a Git-Bash-backed shell rewrote
   `/home/priya/...` into `C:/Program Files/Git/home/priya/...`, creating an
   empty directory tree in the wrong filesystem. Detected immediately (the
   scaffold listing showed Windows paths), inspected to confirm it contained no
   files, and removed. All WSL invocation now goes through PowerShell, and
   multi-step operations are written to script files rather than passed as
   inline strings — PowerShell also swallows `|`, `$()` and `!` inside arguments
   passed to native executables.

2. **`bioconductor.org` is unroutable from this machine.** DNS resolves fine
   (18.161.246.x, AWS CloudFront) but every TCP connection times out at 25 s,
   IPv4-forced or not. `BiocManager` reported 0 visible packages. Posit P3M does
   not mirror Bioconductor (404 at every candidate path). Probed four official
   mirrors: GWDG, TU Dortmund and TUNA all returned HTTP 200; AARNet 404'd.
   Adopted GWDG (D-006). This is a network workaround, not a scientific choice,
   and is flagged as such in the decision log and in `.Rprofile`.

**`.gitignore` verified adversarially before any download existed.** Eleven decoy
files mimicking real downloads (`GSE230649_H1048_MYC.bedGraph`, `GSE230649_RAW.tar`,
`IMfirst_DSP_rawcounts.xlsx`, `.rds` intermediates, `.Renviron`, a stray bedGraph
in the repo root) were planted and all eleven confirmed ignored. Nine deliverable
paths (README, `renv.lock`, figures, results tables, metadata manifest) confirmed
still trackable — an over-broad ignore file is a failure too. Backed by a
`pre-commit` hook rejecting any staged file over 5 MB (D-010).

**Apt audit before asking for sudo.** Of four apparently-missing packages, two
were phantom renames on noble (`libfreetype6-dev`, `libtiff5-dev` — the real
`libfreetype-dev` and `libtiff-dev` were already installed), one is unnecessary
(`libgit2-dev`), and one (`pandoc`) is genuinely needed but not until M10.
**Net sudo required at M3: none** (D-009).

---

**Environment built and locked.** R 4.6.1, Bioconductor 3.23, **34 of 35
packages installed**, `renv.lock` written with `type = "all"` so every installed
package is pinned rather than only those renv can infer from code.

Much of the install linked from an existing renv cache, so it took minutes rather
than the hours a cold source build would have needed.

**One genuine failure: GSVA.** It hard-depends on `SpatialExperiment` →
`magick` → the system library `libmagick++-dev`. Verified by walking the
dependency tree rather than guessing — `magick` is in GSVA's recursive *hard*
dependency set, not Suggests. GSVA is the secondary scoring method only, so
nothing on the critical path is blocked (D-012).

Two package names in the original inventory were also wrong for noble
(`libfreetype6-dev`, `libtiff5-dev` do not exist under those names) and the real
packages were already installed — worth checking before ever reporting a
dependency as missing.

**Environment tests written and passing** (`tests/test_environment.R`): 61
assertions covering directory structure, package loading, and — more usefully —
config invariants. The tests assert that `genome_build` is hg19, that
`chromosome_wise` processing stays enabled, and that **MOES has no weights
field**. That last one turns decision D-003 from a note in a document into
something that fails a test if someone later "improves" the framework by adding
weights.

**First commit:** 39 files, 662 KB, no data. Verified before committing that no
data payload, no data-format extension, no file over 5 MB, and no credential
file was staged.

**Clean-clone test run early, not saved for M11.** Cloned the pushed remote into
a throwaway directory, ran `renv::restore()` (196 packages, all linked from
cache), and ran the environment suite *inside the clone* — all passing. This is a
fair test of the lockfile rather than of a cold machine, since the renv cache is
shared, but it catches the failures that actually happen: an incomplete lockfile,
a file never committed, a hard-coded absolute path.

**Two reconciliation false alarms, both chased down and both benign.** Recording
them so neither is investigated twice:

1. `renv.lock` holds 213 packages but the project library has 198. The 15
   "missing" are all R **recommended** packages (MASS, Matrix, lattice, survival,
   mgcv…). renv deliberately does not duplicate base/recommended packages into
   the project library — they come from the R installation. All 15 verified
   loadable.
2. My own diagnostic then reported 13 of those 15 as VERSION DIFFERS. That was a
   bug in the diagnostic, not the environment: the lockfile stores `7.3-65` while
   `packageVersion()` returns `7.3.65`, and I was comparing strings. Fixed to
   compare `package_version()` objects. Worth flagging as a general trap — a
   naive string comparison of R version numbers will manufacture mismatches.

`renv::status()` also reports packages as `used = n`. Expected: the snapshot used
`type = "all"` and no analysis scripts exist yet to `library()` them. Resolves as
scripts are written.

**GSVA resolved the same day.** `libmagick++-dev` and `pandoc` installed; `magick`
2.9.1 then came down as a prebuilt P3M binary rather than compiling, and GSVA
2.6.3 built cleanly. Environment is **35/35, 231 packages locked**. `pandoc`
3.1.3 also means the M10 reporting dependency is satisfied early.

Rather than just note that GSVA is back, the test suite now asserts
`tumour_scoring.secondary_method == "gsva"`. The reasoning for keeping GSVA is
that `singscore` and `AUCell` are both rank-based — pairing them would make the
sensitivity analysis a restatement rather than an independent check. That
argument now fails a test if someone later swaps in AUCell for convenience.

**M3 closed:** 35/35 packages, 64 test assertions passing, clean-clone restore
verified, repository pushed.

---

## 2026-07-26 · M4 · Data acquisition — verification before download

**Goal.** Manifest first, then resumable checksum-aware downloads, then QC.

**The decision to verify everything against GEO's own records before fetching a
byte turned out to be the most valuable hour of the project so far.** Five things
were not what the frozen inventory implied. None of them would have thrown an
error; all of them would have quietly produced wrong results.

**1. Genome builds are split, badly.** Builds were read from each series'
`!Sample_data_processing` "Assembly" declaration rather than from file names. The
keystone, NEUROD1 and POU2F3 datasets are hg19. **All three supporting ATAC
datasets and the ASCL1 ChIP are hg38.** 23 of 66 files need lifting.

That matters more than a routine liftOver, because the consensus region universe
requires ≥2 independent ATAC datasets and only one of them is in the project
build. Risk R-05 escalated MED → HIGH. Policy set in D-014: **never lift signal;
call intervals in the native build and lift only intervals; report loss rate.**
Lifting a bigWig across builds silently misattributes coverage, which is exactly
the kind of error that produces a plausible figure.

**2. GSE269424 is not native ATAC.** Its titles are `H524_ASCL1`, `H524_EGFP`,
`Lu139_NEUROD1`, `Lu139_EGFP` — the TF names are **overexpressed transgenes**, not
antibodies. It is a chromatin-remodelling experiment. Had the "ATAC-seq of SCLC
cell lines" label been taken at face value, engineered accessibility would have
contaminated the region universe underpinning Aim 1. Only the four EGFP control
arms are used now.

**3. The lineage-TF controls barely share cell lines with the keystone.** This is
the finding with real scientific consequences (new risk R-14, HIGH). POU2F3 covers
NCIH1048, NCIH211 and NCIH526 — three keystone lines across two paralogs, in hg19,
with peak BEDs. ASCL1 covers SHP-77 and nothing else usable. **NEUROD1 covers only
H446, which has zero keystone overlap.**

So the R-01 confounding analysis cannot be uniformly within-line. Its resolution is
now stated per TF: within-line for POU2F3, n=1 descriptive for ASCL1, subtype-level
only for NEUROD1. The frozen design is unchanged; what changes is that a blanket
"we controlled for lineage TFs" claim would be false, and the report must say so.

**4. GSE249362 is both the best control and the cheapest.** Triaging its 185
supplementary files found POU2F3 peak BEDs at ~1.5 MB total. Taking the bigWigs
instead would have been ~20 GB for strictly less usable information.

**5. GSE60052 has no raw counts.** Only `normalized.log2` is deposited, confirmed
by reading the header. DESeq2 is therefore unusable on this cohort; limma on the
log2 matrix it is. Both were already in the dependency inventory, so this selects
between existing options rather than changing anything. Header whitespace and the
`.normal` suffix are both asserted in QC.

**My own errors, recorded because the pattern matters more than the instances:**

- Three representative sample IDs in an early probe were **guessed** rather than
  read from the series records, and returned unrelated data — mouse liver, a kidney
  biopsy, MCF7. A wrong accession still returns a valid-looking record, so the
  failure is silent. Re-derived from `!Series_sample_id`.
- A `SHP77` file pattern matched **nothing**, because the file names spell it
  `SHP-77` while the sample titles say `SHP77`. That would have silently dropped
  the only keystone-overlapping ASCL1 samples in the project. `expected_files` is
  now asserted per dataset so a zero match is a hard failure, never an empty
  success.
- The manifest builder wrote a **16-row manifest instead of 65** when NCBI
  rate-limited enumeration, and still printed "manifest built". Same silent-partial
  -success class as the WSL mangling at M3. It now caches filelists and HEAD sizes
  on disk, and refuses to write at all if any dataset fails to enumerate.
- `!isTRUE(vector)` reported all 70 endpoints as non-resumable. Second
  vectorisation bug of this project after the version-string comparison; worth
  remembering that `isTRUE`/`identical` silently collapse vectors.

**Also settled:** DepMap serves a Cloudflare bot check and explicitly asks users
not to scrape it, so it is a documented manual download (D-016) rather than
something to work around. cBioPortal's `sclc_ucologne_2015` re-confirmed to carry
**no copy-number profile** — only mutations, expression and structural variants.

**Integrity, honestly.** GEO publishes **no MD5 or SHA256** for any of these
supplementary files. So the manifest records `checksum_algorithm =
"none_published"` rather than an empty field, and integrity rests on GEO's
published byte size (the only independent check) plus SHA256 trust-on-first-use
plus `gzip -t`. TOFU proves nothing has changed since we fetched it; it does not
prove we fetched what was deposited. Stated in D-017 rather than glossed.

**Built:** `config/datasets.yml` (curated registry), `01_build_manifest.R`
(generates the manifest — sizes and URLs never hand-typed), `02_download.sh`
(resumable, size-verified, gzip-tested, deliberately serial), `03_verify.R` (the
gate: size + SHA256 + gzip + HTML-error-page detection), `04_qc_report.R`
(structural QC, build assertion, cell-line alias map).

**Verified working:** manifest builds to 70 rows / 66 automated files / 12.08 GB,
all endpoints confirmed resumable (HTTP 206). Verify gate correctly reports 66
missing and exits non-zero while separating manual-required items from failures.

**Next:** download completes (~12 GB, serial by choice), then verify, then QC.

---

<!-- Template for future entries:

## YYYY-MM-DD · Mn · Short title

**Goal.**
**Done.**
**Result / numbers.**
**Sanity checks run.**
**What broke.**
**Decisions logged.** D-0xx
**Next.**

-->

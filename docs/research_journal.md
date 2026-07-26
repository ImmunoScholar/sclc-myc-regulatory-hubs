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

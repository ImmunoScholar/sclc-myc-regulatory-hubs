# Dependency Inventory

Target platform: **Ubuntu (local)**, R-first with minimal Python. Nothing paid, nothing cloud-required.
Versions are pinned at Phase 3 via `renv.lock`; the list below is the intended set, not yet installed.

---

## System dependencies (apt)

Install before any R package compilation — several Bioconductor packages fail to build without these.

```bash
sudo apt update && sudo apt install -y build-essential gfortran libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev zlib1g-dev libbz2-dev liblzma-dev libhdf5-dev pandoc git curl wget
```

| Package | Needed by |
|---|---|
| `build-essential`, `gfortran` | compiling R packages from source |
| `libcurl4-openssl-dev`, `libssl-dev` | `curl`, `httr`, downloads |
| `libxml2-dev` | `XML`, `xml2` |
| `libfontconfig1-dev`, `libharfbuzz-dev`, `libfribidi-dev`, `libfreetype6-dev`, `libpng-dev`, `libtiff5-dev`, `libjpeg-dev` | `ragg`, `systemfonts`, `textshaping` — required for publication-quality figure rendering |
| `zlib1g-dev`, `libbz2-dev`, `liblzma-dev` | `Rsamtools`, `Rhtslib`, `rtracklayer` |
| `libhdf5-dev` | HDF5-backed single-cell objects (conditional layer) |
| `pandoc` | Quarto/rmarkdown rendering |

> Note: `bedtools` and `MACS2` are **not** required under the chosen architecture — bedGraph signal is quantified over ATAC-defined regions inside R via `rtracklayer`/`GenomicRanges`. This is a deliberate simplification (see architecture Part A / dataset trap on GSE230649).

---

## R (target: R ≥ 4.4)

### Infrastructure
`renv` · `here` · `sessioninfo` · `targets` (optional pipeline orchestration) · `quarto`

### Genomic ranges & signal
`GenomicRanges` · `rtracklayer` · `GenomeInfoDb` · `IRanges` · `S4Vectors` · `BSgenome.Hsapiens.UCSC.hg19` · `TxDb.Hsapiens.UCSC.hg19.knownGene` · `org.Hs.eg.db` · `annotatr`

> hg19 is mandated by GSE230649's source alignment. `liftOver` (via `rtracklayer`) only if a documented, justified lift to hg38 is required.

### Expression & differential analysis
`DESeq2` · `limma` · `edgeR` · `SummarizedExperiment` · `sva` (batch diagnostics)

### Signature scoring
`GSVA` · `singscore` · `AUCell` · `msigdbr` (Hallmark comparators)

### Network inference
`GENIE3` · `igraph`

### Evidence integration
`RobustRankAggreg` — the MOES engine

### Data access
`GEOquery` · `cBioPortalData` (or `cgdsr` fallback) · `readxl` (GeoMx .xlsx) · `data.table` · `arrow` (optional, large tables)

### Visualisation
`ggplot2` · `ComplexHeatmap` · `circlize` · `patchwork` · `ggrepel` · `scales` · `ragg` · `viridisLite` · `colorspace`

### ~~Conditional (single-cell layer only)~~ — **REMOVED 2026-07-26**
~~`Seurat` · `SingleCellExperiment` · `scater` · `zellkonverter`~~
The single-cell layer was formally dropped (risk log R-03). These packages are **not** installed. Do not add them back without reopening R-03.

---

## Python (minimal — only if specific optional paths are taken)

Target: Python ≥ 3.10 in a dedicated venv, **not** conda.

| Package | Needed for | Status |
|---|---|---|
| `arboreto` (GRNBoost2) | Speed fallback for GENIE3 | Optional — only if GENIE3 runtime is prohibitive |
| ~~`pyscenic`~~ | ~~MYC-program cell-state activity~~ | **Dropped** — single-cell layer excluded (R-03) |
| ~~`synapseclient`~~ | ~~HTAN/Synapse download~~ | **Dropped** — Synapse route abandoned (R-03) |

```bash
python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
```

**Status 2026-07-26: the single-cell layer is dropped, so the project has no required Python dependency.** `arboreto` remains the only conditional one, and only if GENIE3 proves too slow. Preferred outcome: a pure-R project.

---

## Git / GitHub

| Requirement | Detail |
|---|---|
| Repository | Create **manually on github.com** at the start of Phase 3 — Claude does not create it for you |
| Name | `sclc-myc-regulatory-hubs` |
| Visibility | Public (portfolio) — but create it **empty**, no auto-README, to avoid a merge conflict on first push |
| Protocol | **SSH recommended.** HTTPS requires a Personal Access Token on every push since password auth was removed |
| SSH setup | `ssh-keygen -t ed25519 -C "vpriyaangel@gmail.com"` → add `~/.ssh/id_ed25519.pub` to GitHub → Settings → SSH and GPG keys → verify with `ssh -T git@github.com` |
| Identity | `git config --global user.name` / `user.email` must be set before the first commit |
| Large files | **No data files in git.** `data/raw/` and `data/processed/` go in `.gitignore`. GSE230649 alone is 8.6 GB — far beyond GitHub limits. Git LFS is *not* recommended (free-tier quota is small); use a download script + checksum manifest instead |

**Credential warning:** you will be blocked at first push if SSH keys are not configured. Do this before Phase 3, not during it.

---

## Containers (optional)

**Recommendation: defer.** Docker/Apptainer adds real value only once the pipeline is stable. Premature containerisation slows iteration. `renv.lock` + a pinned R version + the apt list above is sufficient reproducibility for a portfolio project. Revisit at Phase 7 — if added, use `rocker/r-ver:<pinned>` as the base.

## GitHub Actions (optional)

**Recommendation: yes, but lightweight.** A single workflow that (a) checks `renv.lock` restores, (b) runs `R CMD check`-style lint on scripts, (c) confirms no data files were committed. **Do not attempt to run the analysis in CI** — the data are far too large for free runners.

---

## Hardware — **MEASURED 2026-07-26** (no longer an assumption)

Target machine: Windows 11 host, **WSL2 Ubuntu-24.04**, user `priya`. Verified directly:

| Resource | Needed | Actual | Verdict |
|---|---|---|---|
| Disk (Linux FS, `/dev/sdd`) | ~40 GB | **948 GB free** of 1007 GB | Ample |
| RAM | 16 GB comfortable | **10 GiB** + 8 GiB swap | **Binding constraint** → see R-13 |
| Cores | ≥2 | **6** | Fine; parallelise GENIE3 |
| R | 4.5+ | **4.6.1 "Happy Hop"** (2026-06-24) | Matches the MI project toolchain |
| git | any | **2.43.0**, identity configured | Ready |
| GitHub SSH | required | `id_ed25519` present, authenticates as `ImmunoScholar` | Ready |
| `gh` CLI | optional | not installed | Not required — SSH is sufficient |

**Consequence of the 10 GB measurement:** chromosome-by-chromosome bedGraph processing is mandatory, not the fallback plan. This is now a design requirement of every coordinate-handling script (R-13).

| Runtime | GENIE3 on ~50 lines × ~15k genes: minutes to low hours single-threaded; parallelise across the 6 available cores |
|---|---|

### Build location

Project root: **`/home/priya/projects/sclc-myc-regulatory-hubs`** on the native Linux filesystem.
Do **not** build under `/mnt/c/` — WSL2 I/O across the Windows filesystem boundary is slow enough to matter for an 8.6 GB bedGraph workload.

**Invocation gotcha:** calling WSL from PowerShell leaks the Windows `HOME` (`C:UsersPriya`) into the login shell, which will break `renv` and `~`-relative paths. Invoke WSL through the Bash tool (`wsl.exe -d Ubuntu-24.04 -- bash -lc '...'`), where `HOME=/home/priya` resolves correctly, or set `HOME` explicitly.

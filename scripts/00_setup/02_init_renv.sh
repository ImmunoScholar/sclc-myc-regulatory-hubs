#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 02_init_renv.sh — initialise renv and write the project .Rprofile.
#
# Run from the project root:  bash scripts/00_setup/02_init_renv.sh
# Idempotent.
#
# Repo choices (see docs/decision_log.md, D-006):
#   CRAN  -> Posit P3M noble binaries (fast; avoids compiling ~200 CRAN deps)
#   BioC  -> GWDG mirror, because bioconductor.org itself is unroutable from
#            this network (DNS resolves, TCP to 18.161.246.x times out).
# -----------------------------------------------------------------------------
set -euo pipefail

[ -f .gitignore ] || { echo "FATAL: run from the project root" >&2; exit 1; }

Rscript -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org")'
Rscript -e 'renv::init(bare = TRUE, restart = FALSE)'

# renv::init writes a minimal .Rprofile; replace it with ours (activate FIRST).
cat > .Rprofile <<'RPROFILE'
# --- renv must activate before anything else --------------------------------
source("renv/activate.R")

# --- Repositories ------------------------------------------------------------
# CRAN via Posit P3M: serves precompiled noble binaries when the HTTP user agent
# identifies the platform, which turns an hours-long source build into minutes.
# Bioconductor via the GWDG mirror: bioconductor.org is unroutable from this
# network (see docs/decision_log.md D-006). GWDG is an official mirror serving
# identical content.
local({
  osr <- tryCatch(readLines("/etc/os-release", warn = FALSE), error = function(e) character(0))
  is_noble_linux <- identical(Sys.info()[["sysname"]], "Linux") &&
    any(grepl("noble", osr, fixed = TRUE))

  cran <- if (is_noble_linux) {
    "https://packagemanager.posit.co/cran/__linux__/noble/latest"
  } else {
    "https://cloud.r-project.org"
  }

  bioc_mirror <- "https://ftp.gwdg.de/pub/misc/bioconductor"
  bioc_ver    <- "3.23"

  options(
    repos = c(
      CRAN     = cran,
      BioCsoft = file.path(bioc_mirror, "packages", bioc_ver, "bioc"),
      BioCann  = file.path(bioc_mirror, "packages", bioc_ver, "data/annotation"),
      BioCexp  = file.path(bioc_mirror, "packages", bioc_ver, "data/experiment")
    ),
    BioC_mirror = bioc_mirror,
    timeout = 600
  )

  # Required for P3M to return binaries rather than source tarballs.
  if (is_noble_linux) {
    options(HTTPUserAgent = sprintf(
      "R/%s R (%s)", getRversion(),
      paste(getRversion(), R.version$platform, R.version$arch, R.version$os)
    ))
  }
})

# --- Reproducibility hygiene --------------------------------------------------
options(
  warnPartialMatchArgs   = TRUE,
  warnPartialMatchDollar = TRUE,
  scipen                 = 10,
  timeout                = 600
)
RPROFILE

echo "=== .Rprofile written; verifying repos resolve ==="
Rscript -e '
cat("R", as.character(getRversion()), "\n")
print(getOption("repos"))
ap_cran <- nrow(available.packages(repos = getOption("repos")[["CRAN"]]))
ap_bioc <- nrow(available.packages(repos = getOption("repos")[["BioCsoft"]]))
cat("CRAN packages visible:    ", ap_cran, "\n")
cat("BioCsoft packages visible:", ap_bioc, "\n")
if (ap_cran < 1000 || ap_bioc < 1000) stop("repository resolution FAILED")
cat("repository check: PASS\n")
'

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

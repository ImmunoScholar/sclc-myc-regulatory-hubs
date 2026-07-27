# -----------------------------------------------------------------------------
# 06_install_m5_deps.R — packages needed for the M5 upgrades.
#
#   Rscript scripts/00_setup/06_install_m5_deps.R
#
# BSgenome.Hsapiens.UCSC.hg19 (~700 MB) REVERSES decision D-008, which deferred
# it because no planned step required genomic sequence. Motif validation is now
# a planned step (D-021), so the reversal is deliberate and logged as D-022.
#
# Why it matters: it supplies the only check of our derived regions that we
# cannot tune. Region COUNTS are set by our threshold and universe size; sequence
# content is not. If our paralog regions genuinely capture paralog-specific
# binding, they should recover Plotnik's E-box central dinucleotides
# (MYC CAGATG / MYCN CACATG / MYCL CACCTG) without us having asked them to.
# -----------------------------------------------------------------------------

pkgs <- c(
  "BSgenome.Hsapiens.UCSC.hg19",  # sequence, for motif validation (D-022)
  "Biostrings",                   # motif counting
  "colorspace",                   # colourblind simulation (already present)
  "boot",                         # bootstrap for MOES rank stability
  "gt",                           # publication tables
  "systemfonts"                   # font handling for embedded-font PDF output
)

cat("R           :", as.character(getRversion()), "\n")
cat("library     :", renv::paths$library(), "\n\n")

failures <- character(0)
for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [have]    %-32s %s\n", p, as.character(packageVersion(p))))
    next
  }
  cat(sprintf("  [install] %s ...\n", p))
  ok <- tryCatch({ renv::install(p, prompt = FALSE); TRUE },
                 error = function(e) { cat("  [FAIL]   ", p, ":", conditionMessage(e), "\n"); FALSE })
  if (!ok || !requireNamespace(p, quietly = TRUE)) {
    failures <- c(failures, p)
  } else {
    cat(sprintf("  [ok]      %-32s %s\n", p, as.character(packageVersion(p))))
  }
}

cat("\n")
if (length(failures)) {
  cat("FAILED:", paste(failures, collapse = ", "), "\n")
  cat("RESULT: FAIL — do not snapshot.\n")
  quit(status = 1)
}

# Confirm the genome actually loads and is the build we expect.
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
g <- BSgenome.Hsapiens.UCSC.hg19
cat("genome      :", BSgenome::providerVersion(g), "\n")
cat("chr1 length :", format(GenomeInfoDb::seqlengths(g)[["chr1"]], big.mark = ","), "\n")
stopifnot(GenomeInfoDb::seqlengths(g)[["chr1"]] == 249250621)  # hg19 chr1
cat("build check : PASS (chr1 matches hg19)\n")

cat("\nRESULT: PASS. Next: Rscript scripts/00_setup/04_snapshot.R\n")

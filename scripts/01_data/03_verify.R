# -----------------------------------------------------------------------------
# 03_verify.R — verify every downloaded file against the manifest, and maintain
# the project checksum ledger.
#
#   Rscript scripts/01_data/03_verify.R
#
# THIS IS THE GATE. No downstream analysis script may run against data that has
# not passed this check. It is deliberately pedantic.
#
# What it checks, per file:
#   1. presence
#   2. exact byte size against the manifest (GEO's published size)
#   3. SHA256 against data/metadata/checksums.sha256
#        - first successful run RECORDS the hash (trust-on-first-use)
#        - every later run ENFORCES it; a change is an error
#   4. gzip integrity for .gz files
#   5. non-trivial size (a 0-byte or HTML-error-page-sized file is not data)
#   6. that no file is secretly an NCBI/Cloudflare HTML error page
#
# Then it writes back into the manifest the two fields that can only be known
# after acquisition: download_date and sha256_observed.
#
# On the honesty of TOFU: recording our own hash proves a file has not changed
# or corrupted since we fetched it. It does NOT prove we fetched what the authors
# deposited, because GEO publishes no checksums (verified 2026-07-26). The
# published byte size is the only independent check available, and it is used.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(digest)
})

MAN_PATH   <- "data/metadata/dataset_manifest.csv"
CHECKS     <- "data/metadata/checksums.sha256"
stopifnot(file.exists(MAN_PATH))

man <- read.csv(MAN_PATH, stringsAsFactors = FALSE)
man$dest <- file.path("data/raw", man$dataset_id, man$file_name)

# Rows that are not fetchable files
skip_status <- c("api", "manual_required", "pending_triage")

ledger <- if (file.exists(CHECKS)) {
  read.delim(CHECKS, stringsAsFactors = FALSE)
} else {
  data.frame(file = character(0), sha256 = character(0),
             size = numeric(0), first_seen = character(0),
             stringsAsFactors = FALSE)
}

# A downloaded "file" that is really an error page. NCBI returns HTML on 403 and
# DepMap returns a Cloudflare challenge; both are small and start with markup.
looks_like_html <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  head_bytes <- readBin(con, "raw", n = 512)
  if (!length(head_bytes)) return(FALSE)
  txt <- rawToChar(head_bytes[head_bytes != as.raw(0)])
  grepl("^\\s*(<!DOCTYPE|<html|<\\?xml)", txt, ignore.case = TRUE)
}

gzip_ok <- function(path) {
  # gzip -t is the authoritative test; catches truncation a size check misses.
  status <- suppressWarnings(
    system2("gzip", c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE))
  identical(as.integer(status), 0L)
}

results <- data.frame()
n_pass <- 0; n_fail <- 0; n_new <- 0; n_missing <- 0; n_skipped <- 0
problems <- character(0)

cat("=================== verifying ===================\n")

for (i in seq_len(nrow(man))) {
  r  <- man[i, ]
  id <- paste0(r$dataset_id, "/", r$file_name)

  if (r$status %in% skip_status) {
    n_skipped <- n_skipped + 1
    next
  }

  if (!file.exists(r$dest)) {
    n_missing <- n_missing + 1
    problems <- c(problems, paste0("MISSING       ", id))
    next
  }

  sz  <- file.info(r$dest)$size
  bad <- character(0)

  # --- 2. size against manifest ----------------------------------------------
  if (!is.na(r$size_bytes_expected) && sz != r$size_bytes_expected) {
    bad <- c(bad, sprintf("size %d != expected %.0f", sz, r$size_bytes_expected))
  }

  # --- 5/6. is this actually data? -------------------------------------------
  if (sz == 0) {
    bad <- c(bad, "zero bytes")
  } else if (sz < 2048 && looks_like_html(r$dest)) {
    bad <- c(bad, "content is an HTML error page, not data")
  }

  # --- 4. gzip integrity ------------------------------------------------------
  if (grepl("\\.gz$", r$file_name) && sz > 0) {
    if (!gzip_ok(r$dest)) bad <- c(bad, "gzip -t failed (corrupt)")
  }

  # --- 3. SHA256, TOFU then enforced -----------------------------------------
  sha <- digest::digest(r$dest, algo = "sha256", file = TRUE)
  hit <- which(ledger$file == id)
  if (length(hit)) {
    if (!identical(ledger$sha256[hit[1]], sha)) {
      bad <- c(bad, sprintf("SHA256 CHANGED (ledger %s..., now %s...)",
                            substr(ledger$sha256[hit[1]], 1, 12), substr(sha, 1, 12)))
    }
  } else if (!length(bad)) {
    # Only record a first-use hash for a file that otherwise passed. Never
    # enshrine the hash of a file we already believe is broken.
    ledger <- rbind(ledger, data.frame(
      file = id, sha256 = sha, size = sz,
      first_seen = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      stringsAsFactors = FALSE))
    n_new <- n_new + 1
  }

  man$sha256_observed[i] <- sha
  man$download_date[i]   <- format(file.info(r$dest)$mtime, "%Y-%m-%d")

  if (length(bad)) {
    n_fail <- n_fail + 1
    problems <- c(problems, paste0("FAIL          ", id, " :: ",
                                   paste(bad, collapse = "; ")))
  } else {
    n_pass <- n_pass + 1
  }
}

# --- persist ------------------------------------------------------------------
ledger <- ledger[order(ledger$file), ]
write.table(ledger, CHECKS, sep = "\t", row.names = FALSE, quote = FALSE)
man$dest <- NULL
write.csv(man, MAN_PATH, row.names = FALSE)

# --- report -------------------------------------------------------------------
cat("passed              : ", n_pass, "\n", sep = "")
cat("failed              : ", n_fail, "\n", sep = "")
cat("missing             : ", n_missing, "\n", sep = "")
cat("new hashes recorded : ", n_new, "\n", sep = "")
cat("skipped (api/manual): ", n_skipped, "\n", sep = "")
cat("ledger entries      : ", nrow(ledger), "\n", sep = "")

if (length(problems)) {
  cat("\n=================== problems ===================\n")
  cat(paste0("  ", problems, collapse = "\n"), "\n")
}

# --- manual-acquisition reminder ---------------------------------------------
manual <- man[man$status == "manual_required", ]
if (nrow(manual)) {
  cat("\n============ MANUAL ACQUISITION REQUIRED ============\n")
  cat("These cannot be automated and are NOT counted as failures above.\n")
  for (i in seq_len(nrow(manual))) {
    p <- file.path("data/raw", manual$dataset_id[i], manual$file_name[i])
    cat(sprintf("  [%s] %s\n", if (file.exists(p)) "present" else "ABSENT ",
                manual$file_name[i]))
  }
  cat("\nDepMap serves a Cloudflare bot check and asks users not to scrape the\n")
  cat("portal, so it is downloaded by hand from https://depmap.org/portal/data_page/\n")
  cat("into data/raw/depmap/ . The release must then be pinned in\n")
  cat("config/params.yml (depmap.release), which is null as a deliberate tripwire.\n")
}

cat("\n")
if (n_fail > 0 || n_missing > 0) {
  cat("RESULT: FAIL — ", n_fail, " corrupt/mismatched, ", n_missing,
      " missing. Downstream analysis MUST NOT run.\n", sep = "")
  cat("Re-run scripts/01_data/02_download.sh to resume incomplete files.\n")
  quit(status = 1)
}
cat("RESULT: PASS — every automated file verified against the manifest.\n")
cat("Next: Rscript scripts/01_data/04_qc_report.R\n")

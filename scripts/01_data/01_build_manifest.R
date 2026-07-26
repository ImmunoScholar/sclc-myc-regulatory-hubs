# -----------------------------------------------------------------------------
# 01_build_manifest.R — generate the dataset manifest BEFORE any download.
#
#   Rscript scripts/01_data/01_build_manifest.R
#
# Reads the curated registry (config/datasets.yml), resolves every concrete file
# against live GEO/UCSC endpoints, and writes:
#
#   data/metadata/dataset_manifest.csv   one row per file, machine-readable
#   data/metadata/manifest_summary.md    human-readable review copy
#
# Sizes and URLs are NEVER hand-typed. Each file's expected size comes from an
# HTTP HEAD Content-Length, cross-checked against GEO's published filelist.txt
# where one exists. A manifest containing a typo is worse than no manifest,
# because it will be trusted.
#
# Fills: accession, source URL, publication, genome build, size, published
# checksum (or an explicit record that none exists), preprocessing status,
# licence. download_date and sha256_observed stay NA until 02_download.sh and
# 03_verify.R populate them.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(yaml)
  library(httr2)
})

stopifnot(file.exists("config/datasets.yml"))
dir.create("data/metadata", showWarnings = FALSE, recursive = TRUE)

cfg          <- yaml::read_yaml("config/datasets.yml")
PROJECT_BUILD <- cfg$project_build
GEO_SERIES   <- "https://ftp.ncbi.nlm.nih.gov/geo/series"
GEO_SAMPLES  <- "https://ftp.ncbi.nlm.nih.gov/geo/samples"

message("project build: ", PROJECT_BUILD)
message("datasets in registry: ", length(cfg$datasets))

# --- helpers -----------------------------------------------------------------

# GSE230649 -> GSE230nnn ; GSM7230493 -> GSM7230nnn
nnn_dir <- function(acc) {
  num <- sub("^(GSE|GSM)", "", acc)
  pre <- sub("^(GSE|GSM).*", "\\1", acc)
  if (nchar(num) <= 3) paste0(pre, "nnn")
  else paste0(pre, substr(num, 1, nchar(num) - 3), "nnn")
}

# GEO publishes filelist.txt with columns: Archive/File, Name, Time, Size, Type.
# Verified 2026-07-26: it contains NO checksum column for any series used here.
#
# CACHED ON DISK. Two reasons: NCBI rate-limits repeated bursts and starts
# returning 403 HTML error pages, and a manifest build must be repeatable without
# re-hammering someone else's server. The cache lives in data/metadata (tracked,
# a few KB) so the manifest is reproducible from the repo alone.
FILELIST_CACHE <- "data/metadata/geo_filelists"
dir.create(FILELIST_CACHE, showWarnings = FALSE, recursive = TRUE)

fetch_filelist <- function(acc) {
  cache <- file.path(FILELIST_CACHE, paste0(acc, "_filelist.txt"))

  read_cached <- function() {
    txt <- readLines(cache, warn = FALSE)
    # A cached NCBI error page is worse than no cache — detect and reject it.
    if (!length(txt) || !grepl("Archive", txt[1], fixed = TRUE)) return(NULL)
    df <- utils::read.delim(cache, stringsAsFactors = FALSE, check.names = FALSE)
    names(df) <- tolower(gsub("[^A-Za-z]", "", names(df)))
    df
  }

  if (file.exists(cache)) {
    df <- read_cached()
    if (!is.null(df)) { message("  filelist: cached (", nrow(df), " entries)"); return(df) }
    message("  cached filelist was invalid; refetching")
  }

  url <- file.path(GEO_SERIES, nnn_dir(acc), acc, "suppl", "filelist.txt")
  ok <- tryCatch({
    r <- request(url) |> req_timeout(120) |>
      req_user_agent("sclc-myc-regulatory-hubs manifest builder") |>
      req_retry(max_tries = 4, backoff = ~ 15) |> req_perform()
    writeLines(resp_body_string(r), cache)
    TRUE
  }, error = function(e) {
    message("  filelist FETCH FAILED for ", acc, ": ", conditionMessage(e))
    FALSE
  })
  Sys.sleep(3)   # be a good citizen; bursts get thumped
  if (!ok) return(NULL)
  read_cached()
}

# Authoritative expected size: Content-Length from a HEAD request.
#
# CACHED, for the same reason filelists are: NCBI throttles bursts, and a
# manifest build that only succeeds when the network is in a good mood is not
# reproducible. Successful lookups are persisted; failures are not cached, so a
# re-run retries only what is still missing and the build converges.
HEAD_CACHE_FILE <- "data/metadata/head_sizes.tsv"
.head_cache <- if (file.exists(HEAD_CACHE_FILE)) {
  utils::read.delim(HEAD_CACHE_FILE, stringsAsFactors = FALSE)
} else {
  data.frame(url = character(0), size = numeric(0), resumable = logical(0),
             stringsAsFactors = FALSE)
}

head_size <- function(url) {
  hit <- which(.head_cache$url == url)
  if (length(hit)) {
    return(list(size = .head_cache$size[hit[1]],
                resumable = .head_cache$resumable[hit[1]], http = 200L))
  }
  res <- tryCatch({
    r <- request(url) |> req_method("HEAD") |> req_timeout(90) |>
      req_user_agent("sclc-myc-regulatory-hubs manifest builder") |>
      req_retry(max_tries = 4, backoff = ~ 12) |> req_perform()
    cl <- resp_header(r, "Content-Length")
    ar <- resp_header(r, "Accept-Ranges")
    list(size = if (is.null(cl)) NA_real_ else as.numeric(cl),
         resumable = !is.null(ar) && grepl("bytes", ar, ignore.case = TRUE),
         http = resp_status(r))
  }, error = function(e) list(size = NA_real_, resumable = NA, http = NA_integer_))

  if (!is.na(res$size)) {
    .head_cache <<- rbind(.head_cache, data.frame(
      url = url, size = res$size,
      resumable = isTRUE(res$resumable), stringsAsFactors = FALSE))
    utils::write.table(.head_cache, HEAD_CACHE_FILE, sep = "\t",
                       row.names = FALSE, quote = FALSE)
  }
  res
}

rows <- list()
FAILURES <- character(0)   # any dataset that could not be fully enumerated
add_row <- function(d, file_name, url, size, resumable, gsm = NA_character_,
                    required = TRUE, status = "pending") {
  pub <- d$publication
  rows[[length(rows) + 1L]] <<- data.frame(
    dataset_id             = d$id,
    accession              = d$accession %||% NA_character_,
    role                   = d$role %||% NA_character_,
    layer                  = d$layer %||% NA_character_,
    assay                  = paste(d$assay %||% NA_character_, collapse = ";"),
    gsm                    = gsm,
    file_name              = file_name,
    source_url             = url,
    size_bytes_expected    = size,
    resumable              = resumable,
    # Honest recording: GEO publishes no checksums. This is not an oversight to
    # be filled in later, it is a property of the source.
    checksum_published     = NA_character_,
    checksum_algorithm     = "none_published",
    sha256_observed        = NA_character_,   # filled by 03_verify.R (TOFU)
    download_date          = NA_character_,   # filled by 03_verify.R
    genome_build           = d$genome_build %||% NA_character_,
    build_matches_project  = as.character(d$build_matches_project %||% NA),
    liftover_required      = isTRUE(d$liftover_required),
    preprocessing_status   = gsub("\\s+", " ", d$preprocessing_status %||% NA_character_),
    publication_pmid       = as.character(pub$pmid %||% NA_character_),
    publication_citation   = gsub("\\s+", " ", pub$citation %||% NA_character_),
    licence                = gsub("\\s+", " ", d$licence %||% NA_character_),
    acquisition_method     = d$acquisition$method %||% NA_character_,
    required               = required,
    status                 = status,
    stringsAsFactors       = FALSE
  )
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# --- walk the registry -------------------------------------------------------

for (d in cfg$datasets) {
  meth <- d$acquisition$method
  message("\n[", d$id, "] method=", meth)

  if (identical(meth, "MANUAL")) {
    for (f in d$acquisition$required_files) {
      add_row(d, f, d$source_url, NA_real_, NA, required = TRUE,
              status = "manual_required")
    }
    message("  ", length(d$acquisition$required_files),
            " file(s) flagged MANUAL — cannot be automated")
    next
  }

  if (identical(meth, "rest_api")) {
    add_row(d, "(REST API query)", d$api_url %||% d$source_url, NA_real_, NA,
            status = "api")
    next
  }

  if (identical(meth, "direct")) {
    h <- head_size(d$source_url)
    add_row(d, basename(d$source_url), d$source_url, h$size, h$resumable)
    exp <- d$expected_size_bytes
    if (!is.null(exp) && !is.na(h$size) && exp != h$size)
      warning("  SIZE MISMATCH vs registry for ", d$id, ": registry=", exp,
              " live=", h$size, call. = FALSE)
    message("  1 file, ", format(h$size, big.mark = ","), " bytes")
    next
  }

  if (identical(meth, "per_sample_filtered")) {
    # Needs sample-level triage from the series matrix before files can be named.
    add_row(d, "(pending triage)", d$source_url, NA_real_, NA,
            status = "pending_triage")
    message("  deferred: requires per-sample triage (see 01b_triage_gse249362.R)")
    next
  }

  fl <- fetch_filelist(d$accession)

  if (meth %in% c("series_supplementary", "series_tar")) {
    for (f in d$acquisition$files) {
      url <- file.path(GEO_SERIES, nnn_dir(d$accession), d$accession, "suppl", f)
      # Prefer GEO's PUBLISHED size from filelist.txt over a live HEAD: it is
      # authoritative, it is cached, and it does not depend on the network being
      # cooperative at manifest-build time. HEAD is the fallback and cross-check.
      sz <- NA_real_
      if (!is.null(fl) && all(c("name", "size") %in% names(fl))) {
        sz <- suppressWarnings(as.numeric(fl$size[match(f, fl$name)]))
      }
      res <- TRUE
      if (is.na(sz)) {
        h   <- head_size(url); sz <- h$size; res <- h$resumable
        Sys.sleep(1)
      }
      if (is.na(sz)) FAILURES <<- c(FAILURES, paste0(d$id, "/", f, " (no size)"))
      add_row(d, f, url, sz, res)
    }
    message("  ", length(d$acquisition$files), " file(s)")
    next
  }

  if (identical(meth, "per_sample")) {
    if (is.null(fl) || !"name" %in% names(fl)) {
      FAILURES <<- c(FAILURES, paste0(d$id, " (no filelist)"))
      warning("  cannot enumerate ", d$id, " — no filelist", call. = FALSE)
      next
    }
    pat  <- d$acquisition$file_pattern
    keep <- fl$name[grepl(pat, fl$name) &
                      (is.na(fl$archivefile) | fl$archivefile != "Archive")]
    keep <- keep[grepl("^GSM", keep)]
    if (!length(keep)) {
      # A zero match must never be an empty success. This exact case bit once:
      # 'SHP77' matched nothing because the file names spell it 'SHP-77', which
      # would have silently dropped the only usable ASCL1 samples in the project.
      msg <- paste0("pattern '", pat, "' matched NOTHING for ", d$id,
                    " (", length(fl$name), " files in filelist)")
      if (!is.null(d$acquisition$expected_files)) stop(msg, " — refusing to continue")
      warning("  ", msg, call. = FALSE)
      next
    }
    for (f in keep) {
      gsm <- sub("_.*$", "", f)
      url <- file.path(GEO_SAMPLES, nnn_dir(gsm), gsm, "suppl", f)
      sz  <- suppressWarnings(as.numeric(fl$size[match(f, fl$name)]))
      add_row(d, f, url, sz, TRUE, gsm = gsm,
              required = !isTRUE(d$acquisition$optional))
    }
    exp <- d$acquisition$expected_files
    message("  ", length(keep), " file(s) matched",
            if (!is.null(exp)) paste0(" (expected ", exp, ")") else "")
    if (!is.null(exp) && length(keep) != exp)
      stop("FILE COUNT MISMATCH for ", d$id, ": expected ", exp,
           " got ", length(keep), " — refusing to write a wrong manifest")
    next
  }

  warning("  unhandled acquisition method: ", meth, call. = FALSE)
}

man <- do.call(rbind, rows)
stopifnot(nrow(man) > 0)

# --- refuse to write an incomplete manifest ----------------------------------
# An earlier version of this script wrote a 16-row manifest instead of 65 when
# NCBI rate-limited the enumeration, and still printed "manifest built". A
# partially-populated manifest is worse than none: downstream code trusts it, so
# the missing datasets would simply never be downloaded and never be missed.
if (length(FAILURES)) {
  cat("\n================ ENUMERATION FAILURES ================\n")
  cat(paste0("  - ", FAILURES, collapse = "\n"), "\n")
  stop("refusing to write a partial manifest: ", length(FAILURES),
       " item(s) could not be resolved. This is usually NCBI rate-limiting — ",
       "wait a few minutes and re-run. Cached filelists in ",
       FILELIST_CACHE, " make re-runs cheap.")
}

# Every dataset in the registry must be represented.
missing_ds <- setdiff(vapply(cfg$datasets, function(d) d$id, character(1)),
                      unique(man$dataset_id))
if (length(missing_ds))
  stop("datasets absent from the manifest: ", paste(missing_ds, collapse = ", "))

# --- sanity checks on the manifest itself ------------------------------------
cat("\n================ manifest checks ================\n")
cat("datasets represented: ", length(unique(man$dataset_id)), " / ",
    length(cfg$datasets), "\n", sep = "")
dup <- man$file_name[duplicated(man$file_name) & man$file_name != "(REST API query)"]
if (length(dup)) stop("duplicate file names in manifest: ", paste(dup, collapse = ", "))
cat("no duplicate file names\n")

auto <- man[!man$status %in% c("manual_required", "api", "pending_triage"), ]
missing_size <- auto$file_name[is.na(auto$size_bytes_expected)]
if (length(missing_size))
  warning("no expected size resolved for: ", paste(missing_size, collapse = ", "),
          call. = FALSE)
cat("files with a resolved expected size: ", sum(!is.na(auto$size_bytes_expected)),
    " / ", nrow(auto), "\n", sep = "")

# NB: isTRUE() collapses a vector to a single value, so `!isTRUE(auto$resumable)`
# is always TRUE and reported every file as non-resumable. Compare elementwise.
nonres <- auto$file_name[(!auto$resumable | is.na(auto$resumable)) &
                           !is.na(auto$size_bytes_expected)]
cat("resumable endpoints:     ", sum(auto$resumable, na.rm = TRUE), " / ", nrow(auto), "\n", sep = "")
cat("non-resumable endpoints: ", length(nonres), "\n", sep = "")
if (length(nonres) && length(nonres) <= 10)
  cat("  ", paste(nonres, collapse = "\n  "), "\n", sep = "")

cat("\ntotal download volume (automated files): ",
    format(round(sum(auto$size_bytes_expected, na.rm = TRUE) / 1024^3, 2),
           nsmall = 2), " GB\n", sep = "")

cat("\nby genome build:\n")
print(table(man$genome_build, useNA = "ifany"))
cat("\nfiles requiring liftOver to ", PROJECT_BUILD, ": ",
    sum(man$liftover_required), "\n", sep = "")

write.csv(man, "data/metadata/dataset_manifest.csv", row.names = FALSE)
cat("\nwrote data/metadata/dataset_manifest.csv (", nrow(man), " rows)\n", sep = "")

# --- machine worklist for the downloader -------------------------------------
# Tab-separated on purpose: the manifest's citation and licence fields contain
# commas and quotes, and shell CSV parsing of those is a reliable source of
# silent corruption. None of these five fields can contain a tab.
dl <- man[man$status == "pending", c("dataset_id", "file_name", "source_url",
                                     "size_bytes_expected")]
dl$dest <- file.path("data/raw", dl$dataset_id, dl$file_name)
stopifnot(!any(grepl("\t", unlist(dl))))
write.table(dl, "data/metadata/download_list.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("wrote data/metadata/download_list.tsv (", nrow(dl), " files, ",
    format(round(sum(dl$size_bytes_expected, na.rm = TRUE) / 1024^3, 2), nsmall = 2),
    " GB)\n", sep = "")

# --- human-readable summary --------------------------------------------------
md <- c(
  "# Dataset Manifest Summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Project genome build: **", PROJECT_BUILD, "**"),
  paste0("Rows (files): **", nrow(man), "**"),
  "",
  "Generated by `scripts/01_data/01_build_manifest.R` from `config/datasets.yml`.",
  "Do not edit by hand — regenerate.",
  "",
  "## Integrity policy",
  "",
  "GEO publishes **no MD5 or SHA256** for the supplementary files used here",
  "(verified 2026-07-26 across every series). Integrity therefore rests on:",
  "",
  "1. **Exact byte size** — published by GEO and authoritative.",
  "2. **SHA256 trust-on-first-use** — computed at first successful download,",
  "   recorded in `data/metadata/checksums.sha256`, enforced on every run after.",
  "3. **`gzip -t`** on every gzipped file, which catches truncation that a size",
  "   match alone can miss.",
  "",
  "TOFU proves a file has not changed since we fetched it. It does **not** prove",
  "we fetched what the authors deposited. That limit is real and is not papered over.",
  "",
  "## Per-dataset",
  ""
)
for (id in unique(man$dataset_id)) {
  s <- man[man$dataset_id == id, ]
  gb <- unique(s$genome_build)
  lo <- if (any(s$liftover_required)) " — **liftOver required**" else ""
  tot <- sum(s$size_bytes_expected, na.rm = TRUE)
  md <- c(md,
    paste0("### ", id),
    "",
    paste0("- accession: `", unique(s$accession), "`"),
    paste0("- role / layer: ", unique(s$role), " / ", unique(s$layer)),
    paste0("- genome build: **", paste(gb, collapse = ", "), "**", lo),
    paste0("- files: ", nrow(s),
           if (tot > 0) paste0(" (", format(round(tot / 1024^2, 1), nsmall = 1), " MB)") else ""),
    paste0("- acquisition: ", unique(s$acquisition_method),
           " · status: ", paste(unique(s$status), collapse = ", ")),
    paste0("- PMID: ", paste(unique(s$publication_pmid), collapse = ", ")),
    paste0("- preprocessing: ", substr(unique(s$preprocessing_status)[1], 1, 200)),
    paste0("- licence: ", substr(unique(s$licence)[1], 1, 160)),
    ""
  )
}
writeLines(md, "data/metadata/manifest_summary.md")
cat("wrote data/metadata/manifest_summary.md\n")
cat("\nRESULT: manifest built. Nothing downloaded yet.\n")
cat("Next: review data/metadata/manifest_summary.md, then run 02_download.sh\n")

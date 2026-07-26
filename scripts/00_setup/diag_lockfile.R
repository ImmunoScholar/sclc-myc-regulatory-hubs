# Diagnostic: reconcile renv.lock against what is actually loadable.
# A lockfile entry that never restores would be a silent reproducibility hole —
# but R's base/recommended packages are supplied by the R installation itself and
# deliberately not duplicated into the project library. Check LOADABILITY, not
# just presence in the project library, or you will chase a phantom.
# Read-only. Run from the project root.

lock   <- jsonlite::fromJSON("renv.lock")
locked <- names(lock$Packages)
inst   <- rownames(installed.packages(lib.loc = renv::paths$library()))

cat("locked in renv.lock:            ", length(locked), "\n")
cat("present in project library:     ", length(inst), "\n")

not_in_proj <- setdiff(locked, inst)
cat("locked but not in project lib:  ", length(not_in_proj), "\n\n")

if (length(not_in_proj)) {
  base_rec <- rownames(installed.packages(priority = c("base", "recommended")))
  cat("Reconciling those", length(not_in_proj), "packages:\n")
  unresolved <- character(0)
  for (p in not_in_proj) {
    is_br     <- p %in% base_rec
    loadable  <- requireNamespace(p, quietly = TRUE)
    where     <- if (loadable) dirname(find.package(p, quiet = TRUE)[1]) else NA_character_
    lock_ver  <- lock$Packages[[p]]$Version
    have_ver  <- if (loadable) as.character(packageVersion(p)) else NA_character_
    flag <- if (!loadable) "UNRESOLVED" else if (is_br) "base/recommended" else "other lib"
    # Compare as version objects, NOT as strings: R normalises the separator, so
    # the lockfile's "7.3-65" and packageVersion()'s "7.3.65" are the same
    # version. A string comparison here reports 13 false mismatches.
    match_ver <- if (loadable) {
      isTRUE(package_version(lock_ver) == package_version(have_ver))
    } else FALSE
    cat(sprintf("  %-14s lock=%-10s have=%-10s %-17s %s\n",
                p, lock_ver, ifelse(is.na(have_ver), "-", have_ver), flag,
                if (!loadable) "" else if (match_ver) "version OK" else "VERSION DIFFERS"))
    if (!loadable) unresolved <- c(unresolved, p)
  }
  cat("\n")
  if (length(unresolved)) {
    cat("RESULT: FAIL —", length(unresolved), "locked package(s) are not loadable:",
        paste(unresolved, collapse = ", "), "\n")
    quit(status = 1)
  }
  cat("RESULT: PASS — every locked package is loadable.\n")
  cat("Base/recommended packages come from the R installation and are intentionally\n")
  cat("not copied into the project library. Their versions therefore track the\n")
  cat("installed R (4.6.1), which renv.lock also pins.\n")
} else {
  cat("RESULT: PASS — project library covers the entire lockfile.\n")
}

cat("\nNote on renv::status() reporting packages as 'used = n': expected at this\n")
cat("stage. The snapshot was taken with type = \"all\", and no analysis scripts\n")
cat("that library() these packages exist yet. This resolves as scripts are added.\n")

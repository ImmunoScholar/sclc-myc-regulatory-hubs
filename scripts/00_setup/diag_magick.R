# Diagnostic: why does GSVA pull in 'magick'? Determines whether a system
# library install (sudo) is genuinely required or avoidable.
# Read-only. Run from the project root.

repos <- getOption("repos")
ap <- available.packages(repos = c(repos[["CRAN"]], repos[["BioCsoft"]], repos[["BioCann"]]))
cat("packages visible:", nrow(ap), "\n\n")

# Hard dependencies only (what must be installed for GSVA to work).
hard <- tools::package_dependencies("GSVA", db = ap, recursive = TRUE,
                                    which = c("Depends", "Imports", "LinkingTo"))[[1]]
cat("GSVA recursive HARD dependency count:", length(hard), "\n")
cat("magick in HARD deps? ", "magick" %in% hard, "\n\n")

# Including Suggests (what renv may pull if configured to).
soft <- tools::package_dependencies("GSVA", db = ap, recursive = TRUE,
                                    which = "all")[[1]]
cat("magick in deps incl. Suggests? ", "magick" %in% soft, "\n\n")

# Which direct dependants of magick sit inside GSVA's tree?
if ("magick" %in% soft) {
  candidates <- union(hard, soft)
  who <- character(0)
  for (p in candidates) {
    d <- tryCatch(tools::package_dependencies(p, db = ap, recursive = FALSE,
                                              which = c("Depends","Imports","LinkingTo"))[[1]],
                  error = function(e) character(0))
    if ("magick" %in% d) who <- c(who, paste0(p, " [hard]"))
    d2 <- tryCatch(tools::package_dependencies(p, db = ap, recursive = FALSE,
                                               which = "Suggests")[[1]],
                   error = function(e) character(0))
    if ("magick" %in% d2) who <- c(who, paste0(p, " [suggests]"))
  }
  cat("packages in GSVA's tree that require magick:\n")
  if (length(who)) cat(paste0("  ", who, collapse = "\n"), "\n") else cat("  none found\n")
}

cat("\nrenv dependency setting:\n")
print(renv::settings$package.dependency.fields())

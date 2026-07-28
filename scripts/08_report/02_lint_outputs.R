# -----------------------------------------------------------------------------
# 02_lint_outputs.R — enforce the contract's interpretive limits on written output.
#
#   Rscript scripts/08_report/02_lint_outputs.R
#
# The project contract (section 7) prohibits specific vocabulary in all reports,
# figures and code comments, and requires the MOES interpretive-limit statement to
# appear in the README and in the legend of every figure that presents MOES.
#
# Those are testable claims, and this project has been caught three times by
# claims recorded in config that no code ever verified (D-038, D-039, D-040). So
# they are checked here rather than asserted.
#
# Occurrences inside an explicit DISCLAIMER — a sentence that denies the claim,
# e.g. "should not be read as identifying a causal mechanism" — are legitimate and
# are listed as reviewed exceptions below. Anything else fails.
# -----------------------------------------------------------------------------

BANNED <- c("predicts", "driver of", "therapeutic target", "causes", "validated")

# Reviewed exceptions: file + term + why the use is legitimate. Any occurrence
# NOT matching one of these fails the lint. Kept small and specific on purpose —
# a broad allowlist would defeat the check.
EXCEPT <- list(
  list(file = "report.qmd", term = "validated",
       why = "the contract's own disclaimer sentence, quoted verbatim"),
  list(file = "README.md", term = "validated",
       why = "the contract's own disclaimer sentence, quoted verbatim"),
  list(file = "fig04_moes_null_caption.md", term = "validated",
       why = "the contract's own disclaimer sentence, quoted verbatim")
)

targets <- c("report.qmd", "README.md",
             list.files("figures", pattern = "_caption[.]md$", full.names = TRUE))
targets <- targets[file.exists(targets)]

cat("=========== prohibited vocabulary ===========\n")
cat("scanning ", length(targets), " files for: ", paste(BANNED, collapse = ", "), "\n\n", sep = "")

hits <- data.frame()
for (f in targets) {
  ln <- readLines(f, warn = FALSE)
  for (b in BANNED) {
    ix <- grep(paste0("\\b", b, "\\b"), ln, ignore.case = TRUE)
    for (i in ix) hits <- rbind(hits, data.frame(
      file = f, line = i, term = b, text = trimws(ln[i]), stringsAsFactors = FALSE))
  }
}

allowed <- function(f, t) any(vapply(EXCEPT, function(e)
  basename(e$file) == basename(f) && e$term == t, logical(1)))

fails <- 0
if (!nrow(hits)) {
  cat("  no occurrences\n")
} else {
  for (i in seq_len(nrow(hits))) {
    ok <- allowed(hits$file[i], hits$term[i])
    if (!ok) fails <- fails + 1
    cat(sprintf("  [%s] %-28s :%-4d %-18s %s\n",
                if (ok) "reviewed" else "FAIL", basename(hits$file[i]), hits$line[i],
                hits$term[i], substr(hits$text[i], 1, 68)))
  }
}

# --- the MOES interpretive limit must appear where the contract says ----------
cat("\n=========== MOES interpretive limit ===========\n")
# Checks the statement's SUBSTANCE, not one exact string: both "heuristic" and
# "not a predictive" must appear. Matching a single fixed phrase would fail on
# harmless rewording, and this check exists to catch the claim going missing, not
# to police punctuation.
must_carry <- c("report.qmd", "README.md", "figures/fig04_moes_null_caption.md")
missing <- 0
for (f in must_carry) {
  if (!file.exists(f)) { cat(sprintf("  MISSING FILE %s\n", f)); missing <- missing + 1; next }
  txt <- paste(readLines(f, warn = FALSE), collapse = " ")
  present <- grepl("heuristic", txt, ignore.case = TRUE) &&
             grepl("not a predictive", txt, ignore.case = TRUE)
  cat(sprintf("  %-40s %s\n", f, if (present) "present" else "ABSENT"))
  if (!present) missing <- missing + 1
}

cat("\n==============================================================\n")
cat("prohibited-vocabulary failures : ", fails, "\n", sep = "")
cat("missing interpretive limits    : ", missing, "\n", sep = "")
if (fails == 0 && missing == 0) {
  cat("RESULT: outputs comply with the contract's interpretive limits\n")
} else {
  cat("RESULT: FAILED\n")
}
cat("==============================================================\n")
if (fails > 0 || missing > 0) quit(status = 1)

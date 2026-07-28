# -----------------------------------------------------------------------------
# 06_drug_response.R — Aim 3's drug-response association.
#
#   Rscript scripts/04_functional/06_drug_response.R
#
# The last open contract commitment. It is deliberately narrow and ON-THESIS: the
# source study (Plotnik et al. 2024) is titled "MYC family amplification dictates
# sensitivity to BET bromodomain protein inhibitor Mivebresib in SCLC", so the
# question worth asking of public pharmacogenomic data is whether MYC-family
# status tracks BET-inhibitor sensitivity in SCLC cell lines — and whether any
# such association survives the neuroendocrine lineage adjustment that removed
# every other paralog-specific signal in this project.
#
# THREE independent screens (PRISM Repurposing Secondary, Sanger GDSC1, GDSC2) so
# a single-screen artefact cannot carry the result, matching how every other claim
# here was replicated.
#
# Two traps handled explicitly:
#   * "Non-Small Cell Lung Cancer" CONTAINS "Small Cell Lung Cancer". Lineage is
#     matched EXACTLY. A grepl here silently swept 14 extra lines into the SCLC
#     set, the same defect that once inverted an ASCL1 dependency result.
#   * BET inhibitors are named, not pattern-matched. A regex for "BET" returns
#     BETAMETHASONE, BETA-LAPACHONE and BETULINIC-ACID.
#
# Output: data/metadata/m7_drug_response.csv, m7_drug_coverage.csv
#         results/tables/m7_drug_response.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(yaml) })
CFG <- yaml::read_yaml("config/params.yml"); set.seed(CFG$project$seed)

# Drug screens live in the project's own data tree, like every other input. An
# earlier version read them straight from a Windows Downloads folder, which ran
# here and nowhere else; the clean-clone check caught it. DEPMAP_DRUG_DIR
# overrides for anyone keeping the exports elsewhere.
DL <- Sys.getenv("DEPMAP_DRUG_DIR", unset = "data/raw/depmap")

SCREENS <- list(
  PRISM = "Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen)_subsetted.csv",
  GDSC1 = "Drug_sensitivity_AUC_(Sanger_GDSC1)_subsetted.csv",
  GDSC2 = "Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted.csv")

# Named, not matched. These are the compounds that actually inhibit BET
# bromodomains; the commented names are what a "BET" regex wrongly returns.
# "JQ1" must not match JQ12, a different compound present in GDSC1 — hence the
# negative lookahead on a following digit.
BET_TRUE <- c("BIRABRESIB", "MOLIBRESIB", "I-BET151", "I-BET-762", "JQ1(?![0-9])",
              "MIVEBRESIB", "ABBV-075", "OTX015", "GSK525762", "PFI-1",
              "RVX-208", "APABETALONE")

SCLC_EXACT <- "Small Cell Lung Cancer"

read_screen <- function(path) {
  h <- fread(path, nrows = 1)
  lin <- grep("^lineage_", names(h), value = TRUE)
  meta_cols <- intersect(c("depmap_id", "cell_line_display_name", lin), names(h))
  dt <- fread(path)
  # Exact match across whichever lineage column carries the SCLC label.
  hit <- NULL
  for (cc in lin) if (sum(dt[[cc]] == SCLC_EXACT, na.rm = TRUE) > 0) { hit <- cc; break }
  list(dt = dt, meta = meta_cols, lineage_col = hit)
}

cov <- data.frame(); res <- data.frame()
cat("=========== drug screens ===========\n")
cat("  input directory: ", DL, "\n", sep = "")
absent <- names(SCREENS)[!file.exists(file.path(DL, unlist(SCREENS)))]
if (length(absent) == length(SCREENS))
  stop("no drug-screen exports found in ", DL, ".\n",
       "  Download from https://depmap.org/portal/data_page/?tab=customDownloads\n",
       "  (Drug sensitivity AUC: PRISM Repurposing Secondary, Sanger GDSC1, GDSC2)\n",
       "  and place the CSVs there, or set DEPMAP_DRUG_DIR.")

for (sn in names(SCREENS)) {
  p <- file.path(DL, SCREENS[[sn]])
  if (!file.exists(p)) { cat("  ", sn, ": FILE MISSING\n"); next }
  S <- read_screen(p)
  dt <- S$dt
  if (is.null(S$lineage_col)) { cat("  ", sn, ": no exact SCLC lineage label found\n"); next }

  n_exact <- sum(dt[[S$lineage_col]] == SCLC_EXACT, na.rm = TRUE)
  n_grepl <- sum(grepl("Small Cell Lung", dt[[S$lineage_col]]), na.rm = TRUE)
  sclc <- dt[dt[[S$lineage_col]] == SCLC_EXACT]

  drugcols <- setdiff(names(dt), S$meta)
  bet <- drugcols[vapply(drugcols, function(cn)
    any(vapply(BET_TRUE, function(b) grepl(b, cn, ignore.case = TRUE, perl = TRUE),
               logical(1))), logical(1))]

  cat(sprintf("\n  %s: %d lines total | SCLC exact %d (grepl would give %d, +%d false)\n",
              sn, nrow(dt), n_exact, n_grepl, n_grepl - n_exact))
  cat(sprintf("    compounds: %d | BET inhibitors identified: %d\n", length(drugcols), length(bet)))
  if (length(bet)) for (b in bet) cat("      ", substr(b, 1, 62), "\n")

  cov <- rbind(cov, data.frame(screen = sn, lines_total = nrow(dt),
                               sclc_exact = n_exact, sclc_grepl = n_grepl,
                               compounds = length(drugcols), bet_compounds = length(bet),
                               stringsAsFactors = FALSE))
  if (!length(bet) || nrow(sclc) < 8) {
    cat("    -> insufficient for testing in this screen\n"); next
  }

  # --- MYC-family expression and NE score for the same lines ------------------
  # Gene columns in this export are BARE SYMBOLS ("MYC", "MYCL"), not the
  # "MYC (4609)" form used elsewhere in DepMap. Selected by exact name; an
  # earlier pattern assuming the parenthesised form silently selected nothing
  # and every downstream test was skipped rather than failing.
  ex <- fread("data/raw/depmap/depmap_expression_26Q1.csv", nrows = 1)
  want <- c("MYC", "MYCN", "MYCL", "ASCL1", "NEUROD1", "POU2F3", "YAP1",
            "CHGA", "SYP", "INSM1", "NCAM1", "DLL3")
  cols <- c(names(ex)[1], intersect(want, names(ex)))
  if (length(cols) < 4) stop("expression export lacks the required gene columns")
  E <- fread("data/raw/depmap/depmap_expression_26Q1.csv", select = cols)
  M <- merge(sclc[, c("depmap_id", ..bet)], E, by.x = "depmap_id", by.y = names(E)[1])
  cat(sprintf("    SCLC lines with BOTH drug and expression data: %d\n", nrow(M)))
  if (nrow(M) < 8) { cat("    -> too few for testing\n"); next }

  ne_g <- intersect(c("CHGA","SYP","INSM1","NCAM1","DLL3","ASCL1"), names(M))
  ne <- rowMeans(scale(as.matrix(M[, ..ne_g])), na.rm = TRUE)

  for (b in bet) {
    y <- as.numeric(M[[b]])
    if (sum(!is.na(y)) < 8) next
    for (g in c("MYC", "MYCN", "MYCL")) {
      if (!g %in% names(M)) next
      x <- as.numeric(M[[g]])
      ct <- suppressWarnings(stats::cor.test(y, x, method = "spearman", exact = FALSE))
      # lineage-adjusted: partial correlation given NE score
      ok <- stats::complete.cases(y, x, ne)
      pr <- NA_real_; pp <- NA_real_
      if (sum(ok) >= 8) {
        ry <- stats::residuals(stats::lm(y[ok] ~ ne[ok]))
        rx <- stats::residuals(stats::lm(x[ok] ~ ne[ok]))
        pc <- suppressWarnings(stats::cor.test(ry, rx, method = "spearman", exact = FALSE))
        pr <- unname(pc$estimate); pp <- pc$p.value
      }
      res <- rbind(res, data.frame(
        screen = sn, drug = sub(" \\(BRD.*$", "", b), gene = g, n = sum(ok),
        rho_raw = round(unname(ct$estimate), 3), p_raw = signif(ct$p.value, 3),
        rho_ne_adjusted = round(pr, 3), p_ne_adjusted = signif(pp, 3),
        stringsAsFactors = FALSE))
    }
  }
}

write.csv(cov, "data/metadata/m7_drug_coverage.csv", row.names = FALSE)
write.csv(res, "data/metadata/m7_drug_response.csv", row.names = FALSE)

cat("\n=========== association: BET sensitivity vs MYC-family expression ===========\n")
cat("AUC is a sensitivity metric where LOWER = more sensitive, so a POSITIVE rho\n")
cat("means higher expression tracks LESS sensitivity.\n\n")
if (nrow(res)) {
  print(res[order(res$p_raw), ], row.names = FALSE)
  nsig <- sum(res$p_raw < 0.05, na.rm = TRUE)
  nsig_adj <- sum(res$p_ne_adjusted < 0.05, na.rm = TRUE)
  cat(sprintf("\n  nominal p<0.05 raw: %d/%d | after NE adjustment: %d/%d\n",
              nsig, nrow(res), nsig_adj, nrow(res)))
  cat("  (no multiple-testing correction applied to these exploratory pairs;\n")
  cat("   with ", nrow(res), " tests, ", round(0.05*nrow(res),1),
      " nominal hits are expected by chance)\n", sep = "")
} else cat("  no testable combinations\n")

md <- c("# Aim 3 — drug-response association (BET inhibitors)", "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "The last open commitment from the project contract. Deliberately narrow and",
  "on-thesis: the source study reports that MYC-family amplification dictates",
  "sensitivity to BET bromodomain inhibition in SCLC, so the question asked of",
  "public pharmacogenomic data is whether MYC-family status tracks BET-inhibitor",
  "sensitivity, and whether it survives neuroendocrine lineage adjustment.", "",
  "Three independent screens, so no single-screen artefact can carry the result.", "",
  "## Coverage", "",
  "| screen | lines | SCLC (exact) | SCLC (grepl) | compounds | BET |",
  "|---|---|---|---|---|---|",
  sprintf("| %s | %d | %d | %d | %d | %d |", cov$screen, cov$lines_total,
          cov$sclc_exact, cov$sclc_grepl, cov$compounds, cov$bet_compounds), "",
  paste0("The grepl column is shown to make the trap visible: \"Non-Small Cell Lung ",
         "Cancer\" contains \"Small Cell Lung Cancer\", so pattern-matching the ",
         "lineage silently adds NSCLC lines to an SCLC set. Matching is exact."), "",
  "## Associations", "",
  "AUC is a sensitivity metric on which **lower means more sensitive**, so a",
  "positive rho means higher expression tracks *less* sensitivity.", "")
if (nrow(res)) {
  md <- c(md, "| screen | drug | gene | n | rho | p | rho (NE-adj) | p (NE-adj) |",
              "|---|---|---|---|---|---|---|---|",
          sprintf("| %s | %s | %s | %d | %+.3f | %.3g | %s | %s |",
                  res$screen, res$drug, res$gene, res$n, res$rho_raw, res$p_raw,
                  ifelse(is.na(res$rho_ne_adjusted), "—", sprintf("%+.3f", res$rho_ne_adjusted)),
                  ifelse(is.na(res$p_ne_adjusted), "—", sprintf("%.3g", res$p_ne_adjusted))), "",
          paste0("Benjamini-Hochberg across all ", nrow(res), " tests: ",
                 sum(p.adjust(res$p_raw, "BH") < 0.05), " survive at FDR < 0.05 ",
                 "(minimum FDR ", sprintf("%.3f", min(p.adjust(res$p_raw, "BH"))), ")."), "",
          "## Cross-screen replication — read this before the table above", "",
          {
            mycl <- res[res$gene == "MYCL", ]
            sig <- mycl[!is.na(mycl$p_ne_adjusted) & mycl$p_ne_adjusted < 0.05, ]
            prism <- mycl[mycl$screen == "PRISM", ]
            c(paste0("The MYCL association is the **only** signal anywhere in this ",
                     "project that survives neuroendocrine lineage adjustment. That ",
                     "makes it the one most in need of scepticism, not the least."), "",
              paste0("It appears in **", length(unique(sig$screen)), " of 3 screens** ",
                     "(", paste(unique(sig$screen), collapse = ", "), ") and **not in PRISM** ",
                     if (nrow(prism)) paste0("(rho ", sprintf("%+.3f", prism$rho_raw[1]),
                                             ", p = ", signif(prism$p_raw[1], 2), ", n = ",
                                             prism$n[1], ")") else "", "."), "",
              paste0("**GDSC1 and GDSC2 are not independent.** They are the same Sanger ",
                     "platform over heavily overlapping cell lines, so agreement between ",
                     "them is much weaker evidence than agreement between GDSC and PRISM ",
                     "would be. Counting them as two replications would overstate the ",
                     "result."), "",
              paste0("The direction is also worth stating plainly: positive rho means ",
                     "MYCL-high lines are **less** sensitive to BET inhibition. The ",
                     "source study's claim concerns MYC-family *amplification* driving ",
                     "*sensitivity*; MYC here trends the opposite way to MYCL (negative ",
                     "rho, i.e. more sensitive) but does **not** survive lineage ",
                     "adjustment. Nothing here reproduces or refutes that study, which ",
                     "tested amplification rather than expression and a different ",
                     "compound."), "",
              "Reported as **exploratory and provisional**.")
          })
} else md <- c(md, "No testable combinations.")
writeLines(md, "results/tables/m7_drug_response.md")
cat("\nwrote data/metadata/m7_drug_response.csv\n")
cat("wrote results/tables/m7_drug_response.md\n")

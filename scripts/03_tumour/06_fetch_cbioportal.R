# -----------------------------------------------------------------------------
# 06_fetch_cbioportal.R — fetch George et al. 2015 expression for replication.
#
#   Rscript scripts/03_tumour/06_fetch_cbioportal.R
#
# WHY. The gap statement promised replication in a second independent cohort, and
# a NEGATIVE result needs replication more than a positive one would: a single-
# cohort null is much weaker than a null that reproduces in independently
# generated data. M6's finding (lineage dominates, paralog signal absent) rests
# on GSE60052 alone until this runs.
#
# George et al. 2015 (Nature, PMID 26168399) via cBioPortal `sclc_ucologne_2015`.
# Independent of GSE60052: different patients, different sequencing, different
# processing pipeline. If lineage dominance reproduces there, the M6 conclusion is
# solid; if it does not, that must be known BEFORE it is written up as the
# headline.
#
# RE-VERIFIED CONSTRAINT (M4, config/datasets.yml): this study has NO copy-number
# profile — only MUTATION_EXTENDED, MRNA_EXPRESSION (continuous and z-score) and
# STRUCTURAL_VARIANT. Paralog amplification status cannot come from here.
#
# Uses httr2 directly rather than cBioPortalData, which pulls a heavy dependency
# tree (AnVIL) for one REST call.
#
# Output: data/processed/tumour/george2015.rds
#         data/metadata/george2015_fetch_report.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(httr2); library(jsonlite); library(org.Hs.eg.db); library(yaml)
})

CFG   <- yaml::read_yaml("config/params.yml")
STUDY <- "sclc_ucologne_2015"
BASE  <- "https://www.cbioportal.org/api"
OUT   <- "data/processed/tumour"
CACHE <- file.path(OUT, "_cbio_cache"); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

api_get <- function(path) {
  request(paste0(BASE, path)) |> req_timeout(120) |>
    req_retry(max_tries = 3) |> req_perform() |> resp_body_json(simplifyVector = TRUE)
}

# ---- confirm the study still looks as verified at M4 ------------------------
cat("=========== study profiles ===========\n")
prof <- api_get(paste0("/studies/", STUDY, "/molecular-profiles"))
print(prof[, c("molecularProfileId", "molecularAlterationType", "datatype")])
if (any(grepl("COPY_NUMBER", prof$molecularAlterationType)))
  cat("\nNOTE: a copy-number profile now exists — this differs from the M4 check.\n")
expr_prof <- prof$molecularProfileId[
  prof$molecularAlterationType == "MRNA_EXPRESSION" & prof$datatype == "CONTINUOUS"][1]
stopifnot(!is.na(expr_prof))
cat("\nusing expression profile: ", expr_prof, "\n", sep = "")

# ---- samples -----------------------------------------------------------------
sl <- api_get(paste0("/studies/", STUDY, "/sample-lists"))
pick <- sl$sampleListId[grepl("rna_seq", sl$sampleListId)][1]
if (is.na(pick)) pick <- sl$sampleListId[grepl("_all$", sl$sampleListId)][1]
samples <- api_get(paste0("/sample-lists/", pick, "/sample-ids"))
cat("sample list: ", pick, "  (", length(samples), " samples)\n", sep = "")

# ---- genes to fetch ----------------------------------------------------------
reg <- readRDS("data/processed/regions/regulons.rds")$regulons
NE_MARKERS <- c("ASCL1","INSM1","CHGA","SYP","NCAM1","DLL3","CALCA","GRP","UCHL1","SYT11")
want <- unique(c(unlist(reg), "MYC","MYCN","MYCL1","MYCL",
                 "ASCL1","NEUROD1","POU2F3","YAP1", NE_MARKERS))
cat("genes requested: ", length(want), "\n", sep = "")

sym2eg <- suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = want, keytype = "SYMBOL", column = "ENTREZID",
  multiVals = "first"))
sym2eg <- sym2eg[!is.na(sym2eg)]
cat("mapped to Entrez IDs: ", length(sym2eg), "\n\n", sep = "")

# ---- fetch in batches, cached ------------------------------------------------
eg <- unique(as.integer(sym2eg))
BATCH <- 400
chunks <- split(eg, ceiling(seq_along(eg) / BATCH))
cat("=========== fetching (", length(chunks), " batches) ===========\n", sep = "")

all_rows <- list()
for (i in seq_along(chunks)) {
  cf <- file.path(CACHE, sprintf("batch_%02d.rds", i))
  if (file.exists(cf)) {
    all_rows[[i]] <- readRDS(cf); cat(sprintf("  batch %2d [cached]\n", i)); next
  }
  body <- list(entrezGeneIds = chunks[[i]], sampleIds = samples)
  r <- request(paste0(BASE, "/molecular-profiles/", expr_prof, "/molecular-data/fetch")) |>
    req_url_query(projection = "SUMMARY") |>
    req_body_json(body) |> req_timeout(300) |> req_retry(max_tries = 3) |>
    req_perform() |> resp_body_json(simplifyVector = TRUE)
  all_rows[[i]] <- r
  saveRDS(r, cf)
  cat(sprintf("  batch %2d -> %s rows\n", i, format(nrow(r), big.mark = ",")))
  Sys.sleep(1)
}
d <- do.call(rbind, all_rows)
cat("\ntotal rows: ", format(nrow(d), big.mark = ","), "\n", sep = "")
stopifnot(nrow(d) > 0)

# ---- assemble the matrix -----------------------------------------------------
eg2sym <- stats::setNames(names(sym2eg), as.integer(sym2eg))
d$symbol <- eg2sym[as.character(d$entrezGeneId)]
d <- d[!is.na(d$symbol) & !is.na(d$value), ]
syms <- sort(unique(d$symbol)); samp <- sort(unique(d$sampleId))
M <- matrix(NA_real_, length(syms), length(samp), dimnames = list(syms, samp))
M[cbind(match(d$symbol, syms), match(d$sampleId, samp))] <- d$value

cat("matrix: ", nrow(M), " genes x ", ncol(M), " samples\n", sep = "")
cat("missing entries: ", sprintf("%.1f%%", 100 * mean(is.na(M))), "\n", sep = "")

# drop genes measured in too few samples to score
keep <- rowMeans(!is.na(M)) >= 0.8
M <- M[keep, , drop = FALSE]
cat("genes retained (measured in >=80% of samples): ", nrow(M), "\n\n", sep = "")

# ---- scale check: is this log-transformed like GSE60052? --------------------
cat("=========== scale check ===========\n")
cat("range: ", round(min(M, na.rm = TRUE), 2), " to ",
    round(max(M, na.rm = TRUE), 2), "\n", sep = "")
looks_log <- max(M, na.rm = TRUE) < 100
cat(if (looks_log) "consistent with log scale\n"
    else "NOT log scale — values exceed 100. singscore is rank-based so this is\ntolerable, but note the difference from GSE60052 in the write-up.\n")

# ---- coverage ----------------------------------------------------------------
cat("\n=========== regulon coverage in George 2015 ===========\n")
cov <- data.frame()
for (p in names(reg)) {
  n_present <- sum(reg[[p]] %in% rownames(M))
  cov <- rbind(cov, data.frame(paralog = p, regulon_size = length(reg[[p]]),
                               present = n_present,
                               pct = round(100 * n_present / length(reg[[p]]), 1)))
  cat(sprintf("  %-6s %3d of %3d (%.1f%%)\n", p, n_present, length(reg[[p]]),
              100 * n_present / length(reg[[p]])))
}

cat("\n=========== key genes present? ===========\n")
for (g in c("MYC","MYCN","MYCL1","MYCL","ASCL1","NEUROD1","POU2F3","YAP1"))
  cat(sprintf("  %-8s %s\n", g, if (g %in% rownames(M)) "yes" else "NO"))

saveRDS(list(expr = M, coverage = cov, profile = expr_prof,
             sample_list = pick, n_samples = ncol(M),
             note = "George et al. 2015 via cBioPortal; NO copy-number profile in this study"),
        file.path(OUT, "george2015.rds"))
write.csv(cov, "data/metadata/george2015_fetch_report.csv", row.names = FALSE)

cat("\nwrote data/processed/tumour/george2015.rds\n")
cat("\nRESULT: fetched. Next: Rscript scripts/03_tumour/07_replicate_m6.R\n")

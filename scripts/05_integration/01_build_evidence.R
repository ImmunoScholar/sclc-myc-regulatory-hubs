# -----------------------------------------------------------------------------
# 01_build_evidence.R — assemble the per-gene evidence MOES aggregates.
#
#   Rscript scripts/05_integration/01_build_evidence.R
#
# Two admitted domains (D-037), so two evidence vectors per gene:
#
#   cis_regulatory  aggregate peak-to-gene link score within a paralog's active
#                   regions. PARALOG-RESOLVED — a different vector per paralog.
#   functional      evidence of SCLC-selective CRISPR dependency. GENE-LEVEL ONLY
#                   — one vector, identical for all three paralogs (D-036).
#
# The cis score is the SAME quantity that defines the regulons in 21_regulons.R.
# Rather than trusting that this script reimplements it correctly, the
# reconstruction is CHECKED: the top-N genes it produces must equal the committed
# regulons exactly, or the script stops. A silently-drifted reimplementation would
# make every downstream MOES rank wrong in a way nothing else would catch.
#
# Output: data/processed/integration/moes_evidence.rds
#         data/metadata/moes_evidence_summary.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(yaml)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db)
})

CFG <- yaml::read_yaml("config/params.yml")
RG  <- CFG$regulons
set.seed(CFG$project$seed)
dir.create("data/processed/integration", showWarnings = FALSE, recursive = TRUE)

nrm  <- readRDS("data/processed/signal/region_signal_normalised.rds")
p2g  <- readRDS("data/processed/regions/peak_to_gene.rds")
REG  <- readRDS("data/processed/regions/regulons.rds")
DEP  <- readRDS("data/processed/functional/selective_dependency.rds")
u    <- nrm$regions
links <- p2g$retained

cat("universe regions   : ", format(length(u), big.mark = ","), "\n", sep = "")
cat("retained links     : ", format(nrow(links), big.mark = ","), "\n", sep = "")

# --- promoter-proximal assignment, identical to 21_regulons.R -----------------
tx   <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tssg <- promoters(tx, upstream = 0, downstream = 1)
tssg <- tssg[as.character(seqnames(tssg)) %in% paste0("chr", c(1:22, "X"))]
promg <- GenomicRanges::resize(tssg, 2000, fix = "center")
ovp <- GenomicRanges::findOverlaps(u, promg)
prom_sym <- suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = names(tssg)[subjectHits(ovp)], keytype = "ENTREZID",
  column = "SYMBOL", multiVals = "first"))
prom_link <- data.frame(region = queryHits(ovp), gene = unname(prom_sym),
                        score = 1, source = "promoter", stringsAsFactors = FALSE)
prom_link <- prom_link[!is.na(prom_link$gene), ]
dist_link <- data.frame(region = links$region, gene = links$gene,
                        score = links$score, source = "distal",
                        stringsAsFactors = FALSE)
all_link <- rbind(prom_link, dist_link)
cat("promoter links     : ", format(nrow(prom_link), big.mark = ","), "\n", sep = "")
cat("total assignments  : ", format(nrow(all_link), big.mark = ","), "\n\n", sep = "")

# --- cis-regulatory evidence, per paralog -------------------------------------
agg_for <- function(mask) {
  L <- all_link[all_link$region %in% which(mask), ]
  sort(tapply(L$score, L$gene, sum), decreasing = TRUE)
}

PARALOGS <- names(REG$regulons)
cis <- list()
cat("=========== cis-regulatory evidence ===========\n")
for (p in PARALOGS) {
  a <- agg_for(REG$regions[[p]])
  cis[[p]] <- a
  # THE CHECK. Not "did it run" but "is it the same object the project already
  # committed". Any divergence here means MOES would rank against a different
  # cis layer than the one the regulons, the M5 gate and Figure 2 describe.
  rebuilt   <- names(a)[seq_len(min(length(a), RG$max_size))]
  committed <- REG$regulons[[p]]
  if (!identical(rebuilt, committed)) {
    n_shared <- length(intersect(rebuilt, committed))
    stop("cis reconstruction does NOT reproduce the committed ", p, " regulon: ",
         n_shared, "/", length(committed), " genes shared, order ",
         if (identical(sort(rebuilt), sort(committed))) "differs" else "and membership differ",
         ". Refusing to run MOES on a drifted cis layer.")
  }
  cat(sprintf("  %-6s %s genes scored | top-%d reproduces committed regulon EXACTLY\n",
              p, format(length(a), big.mark = ","), RG$max_size))
}

# --- functional evidence, gene level ------------------------------------------
# Signed evidence of SCLC-selective dependency: -log10(p) carrying the sign of
# the effect, so a gene that is MORE essential in SCLC (delta < 0) with a small
# p-value scores highest, and a gene that is LESS essential in SCLC scores
# lowest. No hard essentiality gate is applied: gating would create thousands of
# ties at the boundary and RRA would then rank tied genes arbitrarily. Genes that
# are essential nowhere carry small |delta| and large p, so they fall low on this
# scale without needing to be excluded by hand.
res <- DEP$per_gene
res <- res[!is.na(res$p) & !is.na(res$delta), ]
res$func_score <- -log10(pmax(res$p, .Machine$double.xmin)) * ifelse(res$delta < 0, 1, -1)
func <- setNames(res$func_score, res$gene)
func <- sort(func, decreasing = TRUE)

cat("\n=========== functional evidence ===========\n")
cat("  genes with a dependency test: ", format(length(func), big.mark = ","), "\n", sep = "")
cat("  top 5 by SCLC-selective evidence: ",
    paste(names(func)[1:5], collapse = ", "), "\n", sep = "")
cat("  NOTE: one vector, identical across paralogs. This domain cannot\n")
cat("  contribute anything paralog-specific (D-036).\n")

# --- shared gene universe -----------------------------------------------------
# RRA compares ranks over a common universe. A gene scored in one domain but
# absent from the other has no comparable rank, and imputing it a middling rank
# would invent evidence. Restricted to the intersection, with the loss reported
# rather than passed over.
cis_genes <- Reduce(union, lapply(cis, names))
universe_genes <- sort(intersect(cis_genes, names(func)))
cat("\n=========== shared universe ===========\n")
cat("  genes with cis evidence      : ", format(length(cis_genes), big.mark = ","), "\n", sep = "")
cat("  genes with functional evidence: ", format(length(func), big.mark = ","), "\n", sep = "")
cat("  genes in BOTH (MOES universe) : ", format(length(universe_genes), big.mark = ","), "\n", sep = "")
cat("  dropped, cis-only            : ",
    format(length(setdiff(cis_genes, names(func))), big.mark = ","), "\n", sep = "")
cat("  dropped, functional-only     : ",
    format(length(setdiff(names(func), cis_genes)), big.mark = ","), "\n", sep = "")

ev <- list(cis = cis, func = func, universe = universe_genes,
           all_link = all_link, paralogs = PARALOGS,
           n_cis_only = length(setdiff(cis_genes, names(func))),
           n_func_only = length(setdiff(names(func), cis_genes)))
saveRDS(ev, "data/processed/integration/moes_evidence.rds")

summ <- data.frame(
  paralog = PARALOGS,
  n_genes_cis = vapply(cis, length, integer(1)),
  n_in_universe = vapply(cis, function(a) length(intersect(names(a), universe_genes)), integer(1)),
  regulon_reproduced = TRUE,
  stringsAsFactors = FALSE)
summ <- rbind(summ, data.frame(paralog = "functional (all paralogs)",
                               n_genes_cis = length(func),
                               n_in_universe = length(universe_genes),
                               regulon_reproduced = NA))
write.csv(summ, "data/metadata/moes_evidence_summary.csv", row.names = FALSE)
cat("\nwrote data/processed/integration/moes_evidence.rds\n")
cat("wrote data/metadata/moes_evidence_summary.csv\n")

# -----------------------------------------------------------------------------
# 20_peak_to_gene.R — link distal regions to genes, ABC-style.
#
#   Rscript scripts/02_regulatory/20_peak_to_gene.R
#
# Replaces "nearest gene", which is the assignment the prior work used and a known
# weak one (D-001). Score combines distance decay with activity correlation:
#
#     score = exp(-d / decay) * max(0, rho)
#
# WHAT THE CORRELATION ACTUALLY IS. The config originally said "expression
# correlated", but GSE230649 has no RNA-seq and our expression data are tumour
# cohorts, not these cell lines. Promoter H3K27ac is used as the activity proxy —
# an ABC-family model, NOT expression correlation. Every output says so. Real CCLE
# expression arrives with DepMap at M7 (D-027).
#
# WHY TIERS RATHER THAN ONE CUTOFF. n = 10 lines. Spearman |rho| = 0.3 gives
# p ~ 0.4; significance needs |rho| >= 0.64. A single 0.3 threshold would look like
# a filter and behave like noise, so links carry a confidence tier and the high
# tier is the only one treated as established.
#
# NULL VALIDATION. Distal and promoter H3K27ac both come from the same
# experiments, so lines with globally high signal could correlate everything.
# Fold-over-background removes per-track scaling, but a residual global component
# may survive. Candidate pairs are therefore compared against a null of random
# region-gene pairs >1 Mb apart. If real pairs do not out-correlate the null, the
# correlation term is measuring line-level covariance and is dropped.
#
# Output: data/processed/regions/peak_to_gene.rds
#         data/metadata/peak_to_gene_summary.csv
#         results/tables/m5_peak_to_gene.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(yaml)
})
source("R/genome_utils.R")

CFG <- yaml::read_yaml("config/params.yml")
P2G <- CFG$peak_to_gene
TIER <- P2G$correlation_tiers
DECAY <- P2G$distance_decay_kb * 1000
set.seed(CFG$project$seed)

nrm  <- readRDS("data/processed/signal/region_signal_normalised.rds")
M    <- nrm$mean_fob; meta <- nrm$meta; u <- nrm$regions
k27  <- which(meta$assay == "H3K27ac")
K    <- M[, k27, drop = FALSE]
lines <- meta$line[k27]
cat("H3K27ac matrix: ", nrow(K), " regions x ", ncol(K), " lines\n", sep = "")
cat("activity proxy : promoter H3K27ac (NOT expression) — D-027\n\n")

# --- TSS and promoter activity per gene --------------------------------------
tx  <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
tss <- promoters(tx, upstream = 0, downstream = 1)
tss <- tss[as.character(seqnames(tss)) %in% ANALYSIS_CHROMS_UCSC]
prom <- GenomicRanges::resize(tss, 2000, fix = "center")

# a gene's promoter activity = mean H3K27ac over universe regions in its promoter
ovp <- GenomicRanges::findOverlaps(prom, u)
if (!length(ovp)) stop("no promoter overlaps the universe")
gene_act <- matrix(NA_real_, nrow = length(prom), ncol = ncol(K),
                   dimnames = list(names(prom), lines))
sp <- split(subjectHits(ovp), queryHits(ovp))
for (gi in names(sp)) gene_act[as.integer(gi), ] <- colMeans(K[sp[[gi]], , drop = FALSE])
has_act <- rowSums(is.na(gene_act)) == 0
cat("genes with promoter activity measurable: ", sum(has_act), " of ", length(prom), "\n\n", sep = "")

# --- candidate links: distal regions to genes within max_distance ------------
is_prom_region <- IRanges::overlapsAny(u, prom)
distal <- which(!is_prom_region)
cat("distal regions: ", format(length(distal), big.mark = ","), "\n", sep = "")

# resize_trim(), not resize(): this is the widest window in the pipeline and the
# only one that produced out-of-bound ranges (54). Trimming is a proven no-op
# here — see scripts/00_setup/test_resize_trim.R — and the window is in any case
# only a prefilter, since candidates are cut below by an exact distance test. It
# is used so the warning stops masking a future out-of-bound range that WOULD
# matter.
cand <- GenomicRanges::findOverlaps(
  u[distal], resize_trim(tss, 2 * P2G$max_distance))
cand <- data.frame(region = distal[queryHits(cand)], gene_i = subjectHits(cand))
cand <- cand[has_act[cand$gene_i], ]
cand$distance <- abs(start(u)[cand$region] + width(u)[cand$region] / 2 -
                     start(tss)[cand$gene_i])
cand <- cand[cand$distance <= P2G$max_distance, ]
cat("candidate region-gene pairs: ", format(nrow(cand), big.mark = ","), "\n\n", sep = "")

rowcor <- function(A, B) {
  ra <- t(apply(A, 1, rank)); rb <- t(apply(B, 1, rank))
  ra <- ra - rowMeans(ra); rb <- rb - rowMeans(rb)
  num <- rowSums(ra * rb)
  den <- sqrt(rowSums(ra^2) * rowSums(rb^2))
  ifelse(den > 0, num / den, NA_real_)
}

cat("computing correlations...\n")
cand$rho <- rowcor(K[cand$region, , drop = FALSE], gene_act[cand$gene_i, , drop = FALSE])

# --- null: random pairs > 1 Mb apart -----------------------------------------
n_null <- min(200000, nrow(cand))
nl <- data.frame(region = sample(distal, n_null, replace = TRUE),
                 gene_i = sample(which(has_act), n_null, replace = TRUE))
far <- as.character(seqnames(u))[nl$region] != as.character(seqnames(tss))[nl$gene_i] |
       abs(start(u)[nl$region] - start(tss)[nl$gene_i]) > 1e6
nl <- nl[far, ]
nl$rho <- rowcor(K[nl$region, , drop = FALSE], gene_act[nl$gene_i, , drop = FALSE])

cat("\n=========== null validation ===========\n")
cat(sprintf("candidate pairs  n=%s  median rho %.3f  mean |rho| %.3f  frac >=0.64  %.3f\n",
            format(nrow(cand), big.mark = ","), median(cand$rho, na.rm = TRUE),
            mean(abs(cand$rho), na.rm = TRUE), mean(cand$rho >= TIER$high, na.rm = TRUE)))
cat(sprintf("null pairs       n=%s  median rho %.3f  mean |rho| %.3f  frac >=0.64  %.3f\n",
            format(nrow(nl), big.mark = ","), median(nl$rho, na.rm = TRUE),
            mean(abs(nl$rho), na.rm = TRUE), mean(nl$rho >= TIER$high, na.rm = TRUE)))
enr <- mean(cand$rho >= TIER$high, na.rm = TRUE) / max(mean(nl$rho >= TIER$high, na.rm = TRUE), 1e-9)
cat(sprintf("enrichment of high-tier links over null: %.2fx\n", enr))
null_ok <- enr > 1.5
cat("null validation: ", if (null_ok) "PASS — correlation carries real information"
    else "FAIL — correlation indistinguishable from line-level covariance", "\n", sep = "")

# --- tier and score -----------------------------------------------------------
cand$tier <- ifelse(cand$rho >= TIER$high, "high",
             ifelse(cand$rho >= TIER$moderate, "moderate",
             ifelse(cand$rho >= TIER$weak, "weak", "none")))
cand$dist_weight <- exp(-cand$distance / DECAY)
cand$score <- cand$dist_weight * pmax(0, cand$rho)

sym <- suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = names(tss)[cand$gene_i], keytype = "ENTREZID",
  column = "SYMBOL", multiVals = "first"))
cand$gene <- unname(sym)
cand <- cand[!is.na(cand$gene), ]

cat("\n=========== links by tier ===========\n")
tb <- table(factor(cand$tier, levels = c("high","moderate","weak","none")))
for (k in names(tb))
  cat(sprintf("  %-9s %8s pairs (%5.1f%%)\n", k, format(as.integer(tb[[k]]), big.mark = ","),
              100 * tb[[k]] / nrow(cand)))

kept <- cand[cand$tier %in% c("high","moderate"), ]
cat("\nretained (high + moderate): ", format(nrow(kept), big.mark = ","), " links, ",
    format(length(unique(kept$region)), big.mark = ","), " regions, ",
    format(length(unique(kept$gene)), big.mark = ","), " genes\n", sep = "")
cat("median links per region: ", stats::median(table(kept$region)), "\n", sep = "")
cat("median links per gene  : ", stats::median(table(kept$gene)), "\n", sep = "")

cat("\ncomparison with nearest-gene assignment:\n")
nn <- GenomicRanges::distanceToNearest(u[distal], tss)
nearest_gene <- rep(NA_character_, length(u))
nearest_gene[distal[queryHits(nn)]] <- unname(suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = names(tss)[subjectHits(nn)], keytype = "ENTREZID",
  column = "SYMBOL", multiVals = "first")))
agree <- mean(kept$gene == nearest_gene[kept$region], na.rm = TRUE)
cat(sprintf("  %.1f%% of retained links point to the nearest gene\n", 100 * agree))
cat("  (100%% would mean the model adds nothing over nearest-gene)\n")

saveRDS(list(links = cand, retained = kept, null = nl,
             enrichment_over_null = enr, null_pass = null_ok,
             activity_proxy = "promoter_H3K27ac"),
        "data/processed/regions/peak_to_gene.rds")

summ <- data.frame(
  metric = c("candidate_pairs","high_tier","moderate_tier","weak_tier",
             "retained_links","unique_regions","unique_genes",
             "null_enrichment","pct_agree_nearest_gene"),
  value = c(nrow(cand), sum(cand$tier=="high"), sum(cand$tier=="moderate"),
            sum(cand$tier=="weak"), nrow(kept), length(unique(kept$region)),
            length(unique(kept$gene)), round(enr,2), round(100*agree,1)),
  stringsAsFactors = FALSE)
write.csv(summ, "data/metadata/peak_to_gene_summary.csv", row.names = FALSE)

md <- c("# M5 — peak-to-gene linking", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        "**Activity proxy is promoter H3K27ac, not expression.** GSE230649 contains",
        "no RNA-seq; CCLE expression arrives with DepMap at M7 (D-027). This is an",
        "ABC-family activity model and should be described as such.", "",
        paste0("**n = 10 lines.** Spearman |rho| = 0.3 gives p ~ 0.4; significance needs ",
               ">= 0.64. Links therefore carry tiers and only the high tier is treated ",
               "as established."), "",
        "| metric | value |", "|---|---|")
for (i in seq_len(nrow(summ)))
  md <- c(md, paste0("| ", summ$metric[i], " | ", format(summ$value[i], big.mark = ","), " |"))
writeLines(md, "results/tables/m5_peak_to_gene.md")

cat("\nwrote data/processed/regions/peak_to_gene.rds\n")
if (!null_ok) {
  cat("\nRESULT: FAIL — null validation did not pass. Do not use the correlation term.\n")
  quit(status = 1)
}
cat("\nRESULT: PASS\n")

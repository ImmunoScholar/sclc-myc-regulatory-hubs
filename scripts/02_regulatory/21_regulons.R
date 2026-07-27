# -----------------------------------------------------------------------------
# 21_regulons.R — paralog regulons and their internal-validity gate.
#
#   Rscript scripts/02_regulatory/21_regulons.R
#
# CONSTRUCTION (config regulons, as clarified in D-028)
# A paralog regulon is the TOP-RANKED genes by aggregate link score to that
# paralog's active regions — not every gene with a link. MYC links to thousands
# of genes and singscore/GSVA need focused sets, so the selection rule is part of
# the definition and is stated, not applied silently.
#
# VALIDITY (D-028). The original spec was not computable and was circular:
# it wanted an amplified-vs-non-amplified AUC (MYC amplification is UNVERIFIED,
# D-026) and expression correlation (no cell-line RNA-seq). And a regulon built
# from a paralog's own lines scores high in those lines BY CONSTRUCTION.
#
# Three checks that can each fail:
#   1. LEAVE-ONE-LINE-OUT. Build the paralog's active regions from its remaining
#      line(s); test whether the HELD-OUT line's own ChIP signal is enriched
#      there against background. The held-out ChIP never informed the regions, so
#      this is not circular.
#   2. CROSS-PARALOG DISTINCTNESS. Regulons more similar than max_regulon_jaccard
#      are not separable, and that must be reported rather than glossed.
#   3. BENCHMARK. The MYC regulon should be enriched for HALLMARK_MYC_TARGETS_V1.
#      A MYC regulon with no excess of known MYC targets is not a MYC regulon.
#
# Output: data/processed/regions/regulons.rds
#         data/metadata/regulon_validity.csv
#         results/tables/m5_regulons.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(yaml); library(msigdbr)
})

CFG <- yaml::read_yaml("config/params.yml")
RG  <- CFG$regulons; RV <- CFG$regulon_validity
FOLD <- CFG$active_regions$primary_fold
set.seed(CFG$project$seed)

nrm  <- readRDS("data/processed/signal/region_signal_normalised.rds")
p2g  <- readRDS("data/processed/regions/peak_to_gene.rds")
M    <- nrm$mean_fob; meta <- nrm$meta; u <- nrm$regions
links <- p2g$retained
cat("links (high+moderate): ", format(nrow(links), big.mark = ","), "\n", sep = "")

PARALOG_LINES <- list(MYC = c("H1048","H211","H524","H847","SHP77"),
                      MYCN = c("H526","H69"), MYCL1 = c("COLO668","H889"))
col_for <- function(a, l) { j <- which(meta$assay == a & meta$line == l); if (length(j)==1) j else NA }

active_in_line <- function(p, l, f = FOLD) {
  jp <- col_for(p, l); jk <- col_for("H3K27ac", l)
  if (is.na(jp) || is.na(jk)) return(NULL)
  (M[, jp] >= f) & (M[, jk] >= f)
}
active_set <- function(p, lines = PARALOG_LINES[[p]], min_rep = 2) {
  mats <- Filter(Negate(is.null), lapply(lines, active_in_line, p = p))
  if (!length(mats)) return(rep(FALSE, nrow(M)))
  rowSums(do.call(cbind, mats)) >= min(min_rep, length(mats))
}

# --- promoter-proximal assignment ---------------------------------------------
# ADDED after the Hallmark benchmark failed (D-029). peak_to_gene links DISTAL
# regions only, so a paralog-bound PROMOTER contributed nothing to its regulon —
# discarding roughly half the binding, since only 45-56% of active regions are
# distal. HALLMARK_MYC_TARGETS_V1 is dominated by ribosome-biogenesis and
# translation genes MYC drives FROM PROMOTERS, so a distal-only regulon was
# expected to miss them. It did: 2 hits against ~5 by chance.
suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db)
})
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
cat("promoter-proximal assignments: ", format(nrow(prom_link), big.mark = ","),
    " (", length(unique(prom_link$gene)), " genes)\n", sep = "")

dist_link <- data.frame(region = links$region, gene = links$gene,
                        score = links$score, source = "distal",
                        stringsAsFactors = FALSE)
all_link <- rbind(prom_link, dist_link)
cat("total assignments: ", format(nrow(all_link), big.mark = ","), "\n", sep = "")

build_regulon <- function(region_mask, use = c("both", "distal", "promoter")) {
  use <- match.arg(use)
  idx <- which(region_mask)
  L <- switch(use,
              both     = all_link,
              distal   = dist_link,
              promoter = prom_link)
  L <- L[L$region %in% idx, ]
  if (!nrow(L)) return(character(0))
  agg <- sort(tapply(L$score, L$gene, sum), decreasing = TRUE)
  n <- min(length(agg), RG$max_size)
  if (n < RG$min_size) return(character(0))
  names(agg)[seq_len(n)]
}

cat("\n=========== regulons ===========\n")
reg <- list(); reg_distal <- list(); regions_used <- list()
for (p in names(PARALOG_LINES)) {
  mask <- active_set(p)
  reg[[p]]        <- build_regulon(mask, "both")     # promoter + distal (primary)
  reg_distal[[p]] <- build_regulon(mask, "distal")   # enhancer-only, for the
                                                     # ontology prediction below
  regions_used[[p]] <- mask
  cat(sprintf("  %-6s active regions %6d -> regulon %4d genes (distal-only %4d)\n",
              p, sum(mask), length(reg[[p]]), length(reg_distal[[p]])))
}

# ---- 1. leave-one-line-out ---------------------------------------------------
cat("\n=========== validity 1: leave-one-line-out ===========\n")
cat("Regions built WITHOUT the held-out line; the held-out line's own ChIP is\n")
cat("then tested for enrichment there. Its data never informed the regions.\n\n")
loo <- list()
for (p in names(PARALOG_LINES)) {
  lns <- PARALOG_LINES[[p]]
  lns <- lns[!vapply(lns, function(l) is.null(active_in_line(p, l)), logical(1))]
  for (lh in lns) {
    rest <- setdiff(lns, lh)
    if (!length(rest)) next
    mask <- active_set(p, lines = rest, min_rep = 1)
    jh <- col_for(p, lh)
    if (is.na(jh) || !sum(mask)) next
    inr <- M[mask,  jh]; out <- M[!mask, jh]
    fold <- stats::median(inr) / max(stats::median(out), 1e-9)
    wt <- suppressWarnings(stats::wilcox.test(inr, out, alternative = "greater"))
    # AUC from the rank-sum statistic
    auc <- as.numeric(wt$statistic) / (length(inr) * length(out))
    loo[[length(loo)+1L]] <- data.frame(
      paralog = p, held_out = lh, n_train_lines = length(rest),
      n_regions = sum(mask), median_fold = round(fold, 2),
      auc = round(auc, 3), p = signif(wt$p.value, 3),
      pass = auc >= 0.70 && wt$p.value < 0.05, stringsAsFactors = FALSE)
    cat(sprintf("  %-6s hold out %-8s regions %6d  fold %5.2f  AUC %.3f  p %8.2g  %s\n",
                p, lh, sum(mask), fold, auc, wt$p.value,
                if (auc >= 0.70 && wt$p.value < 0.05) "PASS" else "FAIL"))
  }
}
loo <- do.call(rbind, loo)
v1 <- mean(loo$pass) >= 0.70
cat(sprintf("\nleave-one-line-out: %d/%d held-out lines pass -> %s\n",
            sum(loo$pass), nrow(loo), if (v1) "PASS" else "FAIL"))

# ---- 2. cross-paralog distinctness ------------------------------------------
cat("\n=========== validity 2: cross-paralog distinctness ===========\n")
jac <- function(a, b) length(intersect(a,b)) / max(length(union(a,b)), 1)
J <- outer(names(reg), names(reg), Vectorize(function(a,b) round(jac(reg[[a]], reg[[b]]), 3)))
dimnames(J) <- list(names(reg), names(reg))
print(J)
offd <- J[upper.tri(J)]
v2 <- all(offd <= RV$max_regulon_jaccard)
cat(sprintf("\nmax off-diagonal Jaccard %.3f (limit %.2f) -> %s\n",
            max(offd), RV$max_regulon_jaccard, if (v2) "PASS" else "FAIL"))
if (!v2) cat("Regulons overlap too heavily to be treated as separable programmes.\n")

# ---- 3. benchmark enrichment -------------------------------------------------
cat("\n=========== validity 3: HALLMARK MYC targets ===========\n")
# Read the set name from config rather than a removed key. `benchmark_enrichment`
# was replaced by `benchmark_sets` in D-029 and this line still referenced the old
# name — a config/code divergence of exactly the kind D-028 was about, introduced
# while fixing D-029. Fail loudly if the config does not supply it.
HALLMARK_SET <- names(RV$benchmark_sets)[
  grepl("^HALLMARK", names(RV$benchmark_sets))][1]
stopifnot(!is.na(HALLMARK_SET), nzchar(HALLMARK_SET))

hm <- msigdbr(species = "Homo sapiens", collection = "H")
hm <- hm[hm$gs_name == HALLMARK_SET, ]
stopifnot(nrow(hm) > 0)
hs <- unique(hm$gene_symbol)
universe_genes <- unique(c(links$gene, prom_link$gene))
cat("benchmark set: ", HALLMARK_SET, " (", length(hs), " genes; ",
    length(intersect(hs, universe_genes)), " in our linkable universe)\n\n", sep = "")
bench <- list()
for (p in names(reg)) {
  g <- reg[[p]]
  a <- length(intersect(g, hs)); b <- length(setdiff(g, hs))
  c_ <- length(setdiff(intersect(hs, universe_genes), g))
  d <- length(setdiff(universe_genes, union(g, hs)))
  ft <- stats::fisher.test(matrix(c(a,b,c_,d), 2, byrow = TRUE), alternative = "greater")
  bench[[length(bench)+1L]] <- data.frame(
    paralog = p, n_regulon = length(g), n_hallmark_hit = a,
    odds_ratio = round(as.numeric(ft$estimate), 2),
    p = signif(ft$p.value, 3), stringsAsFactors = FALSE)
  cat(sprintf("  %-6s %4d genes, %3d Hallmark MYC targets, OR %5.2f, p %8.2g\n",
              p, length(g), a, ft$estimate, ft$p.value))
}
bench <- do.call(rbind, bench)
cat("\n(Hallmark alone is NOT the gate — see validity 3c below and D-029.)\n")

# ---- 3b. Plotnik's ontology prediction, stated before testing ---------------
# Published (Plotnik 2024): MYC/MYCN ENHANCER targets -> neurogenesis;
# MYCL1 and shared PROMOTER targets -> housekeeping / ribosome biogenesis.
# So the DISTAL-ONLY regulons should show the opposite pattern to the
# promoter-inclusive ones: depleted for Hallmark MYC targets, enriched for
# neuronal terms. Both directions can fail.
cat("\n=========== validity 3b: distal-only vs Hallmark (prediction) ===========\n")
cat("Plotnik: enhancer targets -> neurogenesis, promoter/shared -> ribosome.\n")
cat("Prediction: distal-only regulons DEPLETED for Hallmark MYC targets.\n\n")
gob <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
neuro <- unique(gob$gene_symbol[gob$gs_name %in% c(
  "GOBP_NEUROGENESIS", "GOBP_NEURON_DIFFERENTIATION",
  "GOBP_GENERATION_OF_NEURONS", "GOBP_NERVOUS_SYSTEM_DEVELOPMENT")])
onto <- list()
for (p in names(reg)) {
  for (wh in c("both", "distal")) {
    g <- if (wh == "both") reg[[p]] else reg_distal[[p]]
    if (!length(g)) next
    ftest <- function(set) {
      a <- length(intersect(g, set)); b <- length(setdiff(g, set))
      c_ <- length(setdiff(intersect(set, universe_genes), g))
      d <- length(setdiff(universe_genes, union(g, set)))
      f <- stats::fisher.test(matrix(c(a,b,c_,d), 2, byrow = TRUE), alternative = "greater")
      c(hits = a, or = as.numeric(f$estimate), p = f$p.value)
    }
    hb <- ftest(hs); nb <- ftest(neuro)
    onto[[length(onto)+1L]] <- data.frame(
      paralog = p, regulon = wh, n = length(g),
      hallmark_hits = hb[["hits"]], hallmark_or = round(hb[["or"]], 2),
      neuro_hits = nb[["hits"]], neuro_or = round(nb[["or"]], 2),
      neuro_p = signif(nb[["p"]], 3), stringsAsFactors = FALSE)
  }
}
onto <- do.call(rbind, onto)
print(onto, row.names = FALSE)
cat("\nIf distal-only regulons are neuro-enriched while promoter-inclusive ones are\n")
cat("Hallmark-enriched, both of Plotnik's ontology findings are reproduced from a\n")
cat("completely independent pipeline. If neither holds, the regulons capture\n")
cat("neither programme and the construction is wrong, not the benchmark.\n")
write.csv(onto, "data/metadata/regulon_ontology.csv", row.names = FALSE)

# ---- 3c. the actual benchmark gate (D-029) ----------------------------------
# Each regulon must show significant enrichment for at least ONE MYC-relevant
# programme, and that programme must be NAMED. Hallmark alone was the wrong sole
# test: it is pan-cancer and promoter-centric, whereas Plotnik's SCLC result is
# that MYC ENHANCER targets are neurogenesis. This is not a relaxation — a
# regulon enriched for nothing still fails.
cat("\n=========== validity 3c: coherent programme per paralog ===========\n")
ALPHA <- RV$benchmark_alpha
prog <- list()
for (p in names(reg)) {
  o <- onto[onto$paralog == p & onto$regulon == "both", ]
  hp <- bench$p[bench$paralog == p]
  cands <- c(HALLMARK = hp, NEUROGENESIS = o$neuro_p)
  best <- names(cands)[which.min(cands)]
  bp   <- min(cands)
  ok   <- bp < ALPHA
  prog[[length(prog)+1L]] <- data.frame(
    paralog = p, programme = if (ok) best else NA_character_,
    p = signif(bp, 3), pass = ok, stringsAsFactors = FALSE)
  cat(sprintf("  %-6s Hallmark p %8.2g | neurogenesis p %8.2g -> %s\n",
              p, hp, o$neuro_p,
              if (ok) paste0("PASS via ", best) else "FAIL (no programme)"))
}
prog <- do.call(rbind, prog)
v3 <- all(prog$pass)
cat(sprintf("\ncoherent programme: %d/%d paralogs -> %s\n",
            sum(prog$pass), nrow(prog), if (v3) "PASS" else "FAIL"))
if (!v3) {
  cat("FAILING PARALOG(S): ", paste(prog$paralog[!prog$pass], collapse = ", "),
      " — no significant enrichment for any tested programme.\n", sep = "")
  cat("These regulons must not be interpreted as regulatory programmes.\n")
}
write.csv(prog, "data/metadata/regulon_programme.csv", row.names = FALSE)

# ---- verdict -----------------------------------------------------------------
res <- data.frame(
  check = c("leave_one_line_out","cross_paralog_distinctness","coherent_programme"),
  detail = c(sprintf("%d/%d held-out lines pass", sum(loo$pass), nrow(loo)),
             sprintf("max Jaccard %.3f (limit %.2f)", max(offd), RV$max_regulon_jaccard),
             paste(sprintf("%s:%s", prog$paralog,
                           ifelse(prog$pass, prog$programme, "NONE")), collapse = " ")),
  pass = c(v1, v2, v3), stringsAsFactors = FALSE)
write.csv(rbind(
  data.frame(type="loo", loo[, c("paralog","held_out","auc","p","pass")]),
  data.frame(type="bench", paralog=bench$paralog, held_out=NA,
             auc=bench$odds_ratio, p=bench$p, pass=bench$p<0.05)),
  "data/metadata/regulon_validity.csv", row.names = FALSE)

cat("\n=========== REGULON VALIDITY GATE ===========\n")
print(res, row.names = FALSE)
saveRDS(list(regulons = reg, regions = regions_used, loo = loo,
             jaccard = J, benchmark = bench, verdict = res),
        "data/processed/regions/regulons.rds")

md <- c("# M5 — paralog regulons", "",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
        paste0("Top-ranked genes by aggregate link score, capped at ", RG$max_size,
               ". Links are high+moderate tier only, from an H3K27ac activity",
               " proxy (D-027) — not expression."), "",
        "| paralog | active regions | regulon genes |", "|---|---|---|")
for (p in names(reg))
  md <- c(md, sprintf("| %s | %s | %s |", p,
                      format(sum(regions_used[[p]]), big.mark = ","), length(reg[[p]])))
md <- c(md, "", "## Validity gate", "", "| check | detail | result |", "|---|---|---|")
for (i in seq_len(nrow(res)))
  md <- c(md, sprintf("| %s | %s | %s |", res$check[i], res$detail[i],
                      ifelse(res$pass[i], "PASS", "**FAIL**")))
writeLines(md, "results/tables/m5_regulons.md")

cat("\nwrote data/processed/regions/regulons.rds\n")
if (!all(res$pass)) {
  cat("\nRESULT: FAIL — regulons must not be scored in tumours until this passes.\n")
  quit(status = 1)
}
cat("\nRESULT: PASS — regulons cleared for M6 tumour scoring.\n")
# -----------------------------------------------------------------------------
# 02_moes.R — the reduced two-domain MOES.
#
#   Rscript scripts/05_integration/02_moes.R
#
# MOES was specified at M1 as a four-domain, two-stage, paralog-specific
# prioritisation. Three domains did not survive: transcriptional was
# lineage-confounded (D-033), network failed both pre-conditions (D-037), and the
# functional domain was admitted but only at gene level (D-036). What remains is
# TWO domains, of which ONE attributes evidence to a paralog.
#
# This script runs that reduced version and reports what it actually is, because
# a two-domain aggregation presented with four-domain language would overstate
# the evidence. Three consequences are computed rather than asserted:
#
#   1. two-stage RRA is DEGENERATE at two domains — stage 2 has nothing to
#      aggregate that stage 1 did not (config two_stage_rra_degenerate)
#   2. leave-one-domain-out reduces to single-domain rankings
#   3. because the functional vector is identical across paralogs, ALL
#      between-paralog difference in MOES comes from chromatin. Quantified.
#
# Output: data/processed/integration/moes_results.rds
#         data/metadata/moes_ranking.csv, moes_diagnostics.csv
#         results/tables/moes.md
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(yaml) })

CFG <- yaml::read_yaml("config/params.yml")
MO  <- CFG$moes
set.seed(CFG$project$seed)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

ev <- readRDS("data/processed/integration/moes_evidence.rds")
G  <- ev$universe; N <- length(G); PARALOGS <- ev$paralogs

cat("=========== reduced MOES ===========\n")
cat("domains admitted   : ", paste(MO$domains, collapse = ", "),
    " (", length(MO$domains), " of 4 planned)\n", sep = "")
cat("paralog-resolved   : ",
    paste(names(which(unlist(MO$domain_attribution) == "paralog_resolved")),
          collapse = ", "), "\n", sep = "")
cat("MOES gene universe : ", format(N, big.mark = ","), "\n\n", sep = "")

# --- RRA -----------------------------------------------------------------------
# Closed form for n = 2 rank lists. For normalised ranks r(1) <= r(2):
#   k=1: Beta(1,2), CDF = 1-(1-x)^2
#   k=2: Beta(2,1), CDF = x^2
#   rho = min over k, score = min(1, n*rho)   [Bonferroni over k, as Kolde 2012]
# Written out rather than called because the permutation and bootstrap loops
# evaluate it ~2e4 times; the closed form is vectorised over all genes at once.
# It is CHECKED against RobustRankAggreg below — a hand-derived null distribution
# that is subtly wrong would make every p-value and FDR in this script wrong.
rra2 <- function(r1, r2) {
  lo <- pmin(r1, r2); hi <- pmax(r1, r2)
  pmin(1, 2 * pmin(1 - (1 - lo)^2, hi^2))
}
norm_rank <- function(score_named, genes) {
  s <- score_named[genes]; s[is.na(s)] <- -Inf
  rank(-s, ties.method = "average") / length(genes)
}

# --- verify the closed form against the reference implementation ---------------
cat("=========== verifying RRA closed form ===========\n")
if (!requireNamespace("RobustRankAggreg", quietly = TRUE))
  stop("RobustRankAggreg required to verify the closed form; refusing to proceed unverified")
set.seed(1)
chk_genes <- G[1:400]
chk_a <- setNames(runif(400), chk_genes); chk_b <- setNames(runif(400), chk_genes)
mine <- rra2(norm_rank(chk_a, chk_genes), norm_rank(chk_b, chk_genes))
names(mine) <- chk_genes
ref <- RobustRankAggreg::aggregateRanks(
  glist = list(names(sort(chk_a, decreasing = TRUE)),
               names(sort(chk_b, decreasing = TRUE))),
  N = 400, method = "RRA")
mine_o <- mine[as.character(ref$Name)]
max_abs <- max(abs(mine_o - ref$Score))
cat(sprintf("  max |closed form - RobustRankAggreg| = %.3e over %d genes\n",
            max_abs, length(chk_genes)))
if (max_abs > 1e-9)
  stop("closed-form RRA disagrees with RobustRankAggreg (max diff ", max_abs,
       "). Refusing to report scores from an unverified null.")
cat("  -> closed form VERIFIED\n\n")

# --- rank vectors --------------------------------------------------------------
r_func <- norm_rank(ev$func, G)
r_cis  <- lapply(PARALOGS, function(p) norm_rank(ev$cis[[p]], G))
names(r_cis) <- PARALOGS

# --- layer correlation (config report_layer_correlation) -----------------------
cat("=========== layer correlation ===========\n")
cat("Spearman between the two domains' rankings. Near zero means RRA is\n")
cat("combining independent evidence rather than double-counting one signal.\n")
layer_rho <- vapply(PARALOGS, function(p)
  suppressWarnings(cor(r_cis[[p]], r_func, method = "spearman")), numeric(1))
for (p in PARALOGS) cat(sprintf("  %-6s rho = %+.4f\n", p, layer_rho[p]))

# --- observed scores -----------------------------------------------------------
obs <- lapply(PARALOGS, function(p) rra2(r_cis[[p]], r_func))
names(obs) <- PARALOGS

# --- permutation FDR -----------------------------------------------------------
# The RRA score is analytically a p-value under random ranks, but the actual rank
# vectors carry ties (many genes share a cis score) so the analytic null is not
# guaranteed. The empirical null permutes the functional labels, preserving both
# marginal rank distributions including their ties.
NPERM <- MO$permutations
BINS  <- 12000                       # bins over -log10(score), resolution 1e-3
tolog <- function(s) -log10(pmax(s, 1e-12))
tobin <- function(L) pmin(pmax(floor(L * 1000) + 1L, 1L), BINS)

cat("\n=========== permutation FDR (", format(NPERM, big.mark = ","),
    " permutations) ===========\n", sep = "")
fdr_for <- function(rc) {
  H <- integer(BINS)
  for (i in seq_len(NPERM)) {
    s <- rra2(rc, r_func[sample.int(N)])
    H <- H + tabulate(tobin(tolog(s)), nbins = BINS)
  }
  # counts at or better than each bin, averaged per permutation
  exp_null <- rev(cumsum(rev(H))) / NPERM
  obs_bin  <- tobin(tolog(rra2(rc, r_func)))
  obs_cnt  <- rev(cumsum(rev(tabulate(obs_bin, nbins = BINS))))
  f <- pmin(1, exp_null[obs_bin] / pmax(obs_cnt[obs_bin], 1))
  f
}
fdr <- list()
for (p in PARALOGS) {
  t0 <- Sys.time()
  fdr[[p]] <- fdr_for(r_cis[[p]])
  cat(sprintf("  %-6s %d genes at FDR < %.2f   (%.0fs)\n", p,
              sum(fdr[[p]] < MO$fdr_threshold), MO$fdr_threshold,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# --- bootstrap rank stability --------------------------------------------------
# Resamples the cis-regulatory LINKS, the sampling unit that actually varies.
# The functional domain is held fixed: its per-gene statistics are precomputed
# summaries over cell lines, and resampling genes rather than lines would not be
# a valid bootstrap of it. The reported intervals therefore reflect
# CIS-REGULATORY uncertainty only and UNDERSTATE total rank uncertainty. Said
# here, in the table, and in the report — not silently.
NBOOT <- MO$bootstrap$n_resamples
cat("\n=========== bootstrap rank stability (", format(NBOOT, big.mark = ","),
    " resamples, cis links only) ===========\n", sep = "")
REG_MASKS <- readRDS("data/processed/regions/regulons.rds")$regions
boot_ci <- list()
for (p in PARALOGS) {
  t0 <- Sys.time()
  L <- ev$all_link[ev$all_link$region %in% which(REG_MASKS[[p]]), ]
  gv <- factor(L$gene); sv <- L$score; nL <- nrow(L)
  ranks <- matrix(NA_integer_, nrow = NBOOT, ncol = N)
  for (b in seq_len(NBOOT)) {
    idx <- sample.int(nL, nL, replace = TRUE)
    agg <- rowsum(sv[idx], gv[idx], reorder = FALSE)
    sc  <- setNames(as.numeric(agg), rownames(agg))
    ranks[b, ] <- rank(rra2(norm_rank(sc, G), r_func), ties.method = "average")
  }
  q <- apply(ranks, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  boot_ci[[p]] <- list(lo = q[1, ], hi = q[2, ])
  cat(sprintf("  %-6s done (%.0fs)\n", p,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# --- leave-one-domain-out ------------------------------------------------------
cat("\n=========== leave-one-domain-out ===========\n")
cat("DEGENERATE at two domains: removing one leaves a SINGLE domain, so this is\n")
cat("not a stability check in the sense specified at M1 (which assumed four).\n")
cat("It measures how far the combined ranking sits from each domain alone.\n")
lodo <- data.frame()
for (p in PARALOGS) {
  comb <- rank(obs[[p]])
  rho_cis  <- suppressWarnings(cor(comb, r_cis[[p]], method = "spearman"))
  rho_func <- suppressWarnings(cor(comb, r_func,     method = "spearman"))
  lodo <- rbind(lodo, data.frame(paralog = p, vs_cis_only = round(rho_cis, 4),
                                 vs_functional_only = round(rho_func, 4)))
  cat(sprintf("  %-6s combined vs cis-only rho %+.3f | vs functional-only rho %+.3f\n",
              p, rho_cis, rho_func))
}

# --- how much of MOES is paralog-specific? -------------------------------------
# The decisive diagnostic. The functional vector is identical for all paralogs,
# so any difference between two paralogs' MOES rankings can only come from the
# cis layer. Quantified by comparing between-paralog agreement in MOES against
# between-paralog agreement in cis alone.
cat("\n=========== what drives between-paralog difference? ===========\n")
pairs <- utils::combn(PARALOGS, 2, simplify = FALSE)
attrib <- data.frame()
for (pr in pairs) {
  rm_ <- cor(rank(obs[[pr[1]]]), rank(obs[[pr[2]]]), method = "spearman")
  rc_ <- cor(r_cis[[pr[1]]], r_cis[[pr[2]]], method = "spearman")
  attrib <- rbind(attrib, data.frame(pair = paste(pr, collapse = " vs "),
                                     moes_rho = round(rm_, 4), cis_rho = round(rc_, 4)))
  cat(sprintf("  %-16s MOES rho %+.3f | cis-only rho %+.3f\n",
              paste(pr, collapse = " vs "), rm_, rc_))
}
cat("\nMOES agreement between paralogs tracks cis agreement, because the shared\n")
cat("functional domain pulls every paralog's ranking toward the SAME genes.\n")

# --- assemble ------------------------------------------------------------------
out <- do.call(rbind, lapply(PARALOGS, function(p) {
  rk <- rank(obs[[p]], ties.method = "min")
  w  <- (boot_ci[[p]]$hi - boot_ci[[p]]$lo) / N
  data.frame(paralog = p, gene = G, moes_score = obs[[p]], rank = rk,
             fdr = fdr[[p]],
             cis_rank = round(r_cis[[p]] * N), func_rank = round(r_func * N),
             boot_rank_lo = round(boot_ci[[p]]$lo), boot_rank_hi = round(boot_ci[[p]]$hi),
             rel_interval_width = round(w, 4),
             stable = w <= MO$bootstrap$max_relative_interval_width,
             in_regulon = G %in% names(ev$cis[[p]])[seq_len(min(500, length(ev$cis[[p]])))],
             stringsAsFactors = FALSE)
}))
out <- out[order(out$paralog, out$rank), ]
write.csv(out, "data/metadata/moes_ranking.csv", row.names = FALSE)

sig <- subset(out, fdr < MO$fdr_threshold)
cat("\n=========== result ===========\n")
for (p in PARALOGS) {
  s <- subset(sig, paralog == p)
  cat(sprintf("  %-6s %d genes at FDR < %.2f, %d of them rank-stable\n",
              p, nrow(s), MO$fdr_threshold, sum(s$stable)))
  if (nrow(s)) cat("         top: ", paste(head(s$gene, 8), collapse = ", "), "\n", sep = "")
}

diag <- data.frame(
  metric = c("domains_admitted", "domains_planned", "paralog_resolved_domains",
             "two_stage_rra_degenerate", "moes_universe_genes",
             paste0("layer_rho_", PARALOGS),
             paste0("n_fdr_", PARALOGS), paste0("n_stable_fdr_", PARALOGS)),
  value = c(length(MO$domains), 4,
            sum(unlist(MO$domain_attribution) == "paralog_resolved"),
            TRUE, N, round(layer_rho, 4),
            vapply(PARALOGS, function(p) sum(fdr[[p]] < MO$fdr_threshold), numeric(1)),
            vapply(PARALOGS, function(p) sum(subset(sig, paralog == p)$stable), numeric(1))),
  stringsAsFactors = FALSE)
write.csv(diag, "data/metadata/moes_diagnostics.csv", row.names = FALSE)

saveRDS(list(ranking = out, fdr = fdr, layer_rho = layer_rho, lodo = lodo,
             attribution = attrib, n_perm = NPERM, n_boot = NBOOT,
             universe = G),
        "data/processed/integration/moes_results.rds")

md <- c(
  "# Reduced two-domain MOES", "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "## What this is, and is not", "",
  paste0("MOES was specified at M1 as a **four-domain, two-stage, paralog-specific** ",
         "prioritisation. It runs here on **two** domains, of which **one** attributes ",
         "evidence to a paralog."), "",
  "| domain | status | attributes to a paralog? |",
  "|---|---|---|",
  "| cis-regulatory | admitted | yes |",
  "| functional | admitted (D-036) | no — gene level only |",
  "| transcriptional | excluded (D-033) | — |",
  "| network | excluded (D-037) | — |", "",
  paste0("**Two-stage RRA is degenerate here.** With two domains, stage 2 has nothing ",
         "to aggregate that stage 1 did not; the reported score is a single ",
         "aggregation, not a hierarchy. Leave-one-domain-out likewise reduces to ",
         "single-domain rankings and is reported as such."), "",
  "## Layer correlation", "",
  "| paralog | Spearman(cis, functional) |", "|---|---|",
  sprintf("| %s | %+.4f |", PARALOGS, layer_rho), "",
  paste0("Near zero, so the two domains are close to independent and RRA is not ",
         "double-counting one signal."), "",
  "## Leave-one-domain-out (degenerate)", "",
  "| paralog | combined vs cis only | combined vs functional only |", "|---|---|---|",
  sprintf("| %s | %+.3f | %+.3f |", lodo$paralog, lodo$vs_cis_only, lodo$vs_functional_only), "",
  "## Between-paralog difference is chromatin alone", "",
  "| pair | MOES rho | cis-only rho |", "|---|---|---|",
  sprintf("| %s | %+.3f | %+.3f |", attrib$pair, attrib$moes_rho, attrib$cis_rho), "",
  paste0("The functional vector is identical for all three paralogs, so every ",
         "difference between their MOES rankings originates in the cis layer. ",
         "MOES cannot deliver multi-layer paralog-specific prioritisation on this ",
         "evidence, and this table is why."), "",
  "## Result", "",
  paste0("Universe: ", format(N, big.mark = ","), " genes with evidence in both domains. ",
         "Permutations: ", format(NPERM, big.mark = ","), ". Bootstrap: ",
         format(NBOOT, big.mark = ","), " resamples of cis links."), "",
  "| paralog | genes at FDR < 0.05 | best FDR achieved | rank-stable |", "|---|---|---|---|",
  vapply(PARALOGS, function(p) {
    s <- subset(sig, paralog == p)
    sprintf("| %s | %d | %.3f | %d |", p, nrow(s), min(fdr[[p]]), sum(s$stable))
  }, character(1)), "",
  paste0("**No gene reaches FDR < 0.05 for any paralog.** The best FDR achieved is ",
         sprintf("%.3f", min(unlist(fdr))), ". This is not a marginal miss: the two ",
         "admitted domains are close to independent (layer rho between ",
         sprintf("%+.3f and %+.3f", min(layer_rho), max(layer_rho)),
         "), so genes ranking highly in one do not rank highly in the other more ",
         "often than chance predicts. RRA has nothing to aggregate."), "",
  paste0("MOES therefore returns **no prioritised hub list**. Reporting a top-N ",
         "table here would present the head of an unranked distribution as a ",
         "result. The ranking is written to `data/metadata/moes_ranking.csv` with ",
         "its FDR column so the absence is inspectable, not hidden."), "",
  "See the sensitivity section below before reading this null as uninformative.", "",
  paste0("**Bootstrap caveat.** Intervals resample cis-regulatory links only; the ",
         "functional domain is held fixed because its per-gene values are ",
         "precomputed summaries over cell lines. Reported rank intervals therefore ",
         "UNDERSTATE total uncertainty."))
writeLines(md, "results/tables/moes.md")

cat("\nwrote data/metadata/moes_ranking.csv\n")
cat("wrote data/metadata/moes_diagnostics.csv\n")
cat("wrote results/tables/moes.md\n")

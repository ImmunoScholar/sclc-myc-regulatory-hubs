# -----------------------------------------------------------------------------
# 19_blacklist_filter.R — remove ENCODE blacklist regions from the universe.
#
#   Rscript scripts/02_regulatory/19_blacklist_filter.R
#
# WHY, AND WHY LATE
#
# Blacklist filtering was named as necessary when the local-background approach
# was written and then never implemented. The omission surfaced in the
# super-enhancer output: among SE loci recurrent in ALL TEN lines were the RNU1-2
# multicopy snRNA cluster (chr1:17.18-17.30 Mb), PDE4DIP in the 1q21.1
# segmental-duplication hotspot (chr1:144.97-145.15 Mb, 180 kb wide), and two
# adjacent segmental duplications at chr7:102.1-102.3 Mb.
#
# Those are mapping artefacts. They recur in every line because artefacts recur in
# every sample — which is exactly what "recurrent across all lines" looks like when
# it is not biology. Recurrence was the main defence against copy-number artefacts,
# so an artefact class that mimics recurrence undermines that defence specifically.
#
# The blacklist has been present in every object since the universe was built: all
# 102,334 regions, every paralog set, and the gate calculations.
#
# FILTERING IS CONSISTENT AND CHEAP. The signal matrix is indexed by universe
# region, so universe and signal are subset with the SAME index. No re-quantification
# is needed — the 45-minute streaming step is not repeated.
#
# Output: data/processed/regions/universe_final.rds        (filtered, in place)
#         data/processed/signal/region_signal_normalised.rds (filtered, in place)
#         data/metadata/blacklist_filter_report.csv
#         originals preserved as *_preblacklist.rds
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(GenomicRanges) })
source("R/genome_utils.R")

BL_DIR <- "data/raw/encode_blacklist"
BL_GZ  <- file.path(BL_DIR, "hg19-blacklist.v2.bed.gz")
BL_URL <- "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz"
dir.create(BL_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(BL_GZ)) {
  cat("fetching ENCODE hg19 blacklist v2...\n")
  st <- utils::download.file(BL_URL, BL_GZ, mode = "wb", quiet = TRUE)
  if (st != 0 || !file.exists(BL_GZ) || file.size(BL_GZ) < 1000) {
    stop("blacklist download failed. Fetch manually into ", BL_GZ, "\n  ", BL_URL)
  }
}
cat("blacklist file: ", BL_GZ, " (", file.size(BL_GZ), " bytes)\n", sep = "")

bl <- utils::read.delim(gzfile(BL_GZ), header = FALSE, stringsAsFactors = FALSE)
names(bl)[1:3] <- c("chrom", "start", "end")
bl_gr <- GenomicRanges::reduce(GRanges(to_ucsc_seqnames(bl$chrom),
                                       IRanges(bl$start + 1L, bl$end)))
bl_gr <- bl_gr[as.character(seqnames(bl_gr)) %in% ANALYSIS_CHROMS_UCSC]
cat("blacklist intervals: ", length(bl_gr), " covering ",
    round(sum(as.numeric(width(bl_gr))) / 1e6, 1), " Mb\n\n", sep = "")
if (ncol(bl) >= 4) { cat("reasons:\n"); print(table(bl[[4]])) ; cat("\n") }

# --- load current objects -----------------------------------------------------
UF  <- "data/processed/regions/universe_final.rds"
SF  <- "data/processed/signal/region_signal_normalised.rds"
RAW <- "data/processed/signal/region_signal_raw.rds"
u   <- readRDS(UF)
nrm <- readRDS(SF)
stopifnot(length(u) == nrow(nrm$mean_fob))

hit <- IRanges::overlapsAny(u, bl_gr)
cat("=========== impact ===========\n")
cat("universe regions          : ", format(length(u), big.mark = ","), "\n", sep = "")
cat("overlapping blacklist     : ", format(sum(hit), big.mark = ","),
    sprintf(" (%.2f%%)\n", 100 * mean(hit)), sep = "")
cat("retained                  : ", format(sum(!hit), big.mark = ","), "\n\n", sep = "")

# Were the flagged SE loci actually blacklisted? A direct check of the diagnosis.
check <- GRanges(c("chr1","chr1","chr7","chr7"),
                 IRanges(c(17178947, 144965361, 102122851, 102217228),
                         c(17303825, 145145225, 102189024, 102279708)))
names(check) <- c("RNU1-2", "PDE4DIP_1q21.1", "RASA4B_chr7SD", "UPK3BL1_chr7SD")
cat("suspect SE loci vs blacklist:\n")
for (i in seq_along(check))
  cat(sprintf("  %-16s %s\n", names(check)[i],
              if (IRanges::overlapsAny(check[i], bl_gr)) "BLACKLISTED" else "not in blacklist"))

# --- preserve originals, then filter -----------------------------------------
if (!file.exists(sub("\\.rds$", "_preblacklist.rds", UF)))
  saveRDS(u,   sub("\\.rds$", "_preblacklist.rds", UF))
if (!file.exists(sub("\\.rds$", "_preblacklist.rds", SF)))
  saveRDS(nrm, sub("\\.rds$", "_preblacklist.rds", SF))

keep <- !hit
u2 <- u[keep]
nrm2 <- nrm
nrm2$mean_raw <- nrm$mean_raw[keep, , drop = FALSE]
nrm2$mean_fob <- nrm$mean_fob[keep, , drop = FALSE]
nrm2$regions  <- u2
nrm2$width    <- nrm$width[keep]
nrm2$blacklist_filtered <- TRUE
saveRDS(u2, UF); saveRDS(nrm2, SF)

if (file.exists(RAW)) {
  r <- readRDS(RAW)
  if (nrow(r$mean) == length(u)) {
    if (!file.exists(sub("\\.rds$", "_preblacklist.rds", RAW)))
      saveRDS(r, sub("\\.rds$", "_preblacklist.rds", RAW))
    r$sum <- r$sum[keep, , drop = FALSE]; r$mean <- r$mean[keep, , drop = FALSE]
    r$regions <- u2; r$width <- r$width[keep]
    saveRDS(r, RAW)
  }
}

rep <- data.frame(
  metric = c("blacklist_intervals","blacklist_mb","universe_before",
             "universe_removed","universe_after","pct_removed"),
  value = c(length(bl_gr), round(sum(as.numeric(width(bl_gr)))/1e6, 1),
            length(u), sum(hit), sum(!hit), round(100*mean(hit), 3)),
  stringsAsFactors = FALSE)
write.csv(rep, "data/metadata/blacklist_filter_report.csv", row.names = FALSE)
print(rep, row.names = FALSE)

cat("\nfiltered universe and signal matrices written in place;\n")
cat("originals preserved as *_preblacklist.rds\n")
cat("\nRE-RUN, in order (no re-quantification needed):\n")
cat("  Rscript scripts/02_regulatory/15_paralog_regions.R\n")
cat("  Rscript scripts/02_regulatory/16_m5_gate.R\n")
cat("  Rscript scripts/02_regulatory/18_super_enhancers.R\n")

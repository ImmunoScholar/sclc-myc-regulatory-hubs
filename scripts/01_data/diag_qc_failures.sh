#!/usr/bin/env bash
# Diagnose the two QC failure clusters. Read-only.
cd "$(dirname "$0")/../.." || exit 1

F=data/raw/GSE230649/GSM7230503_MYC_H1048_R1_treat_pileup.bedgraph.gz

echo "=========== 1. bedGraph: first 5 lines (raw) ==========="
zcat "$F" | head -5 | cat -A | cut -c1-120

echo
echo "=========== 2. bedGraph: unique chrom names in first 500k lines ==========="
zcat "$F" | head -500000 | cut -f1 | sort -u | tr '\n' ' '
echo

echo
echo "=========== 3. bedGraph: ALL unique chrom names (full stream, one file) ==========="
zcat "$F" | cut -f1 | sort -u | tr '\n' ' '
echo

echo
echo "=========== 4. bedGraph: line count and chrom order as encountered ==========="
zcat "$F" | cut -f1 | uniq | head -30 | tr '\n' ' '
echo
echo "total lines: $(zcat "$F" | wc -l)"

echo
echo "=========== 5. GeoMx workbooks: sheet inventory ==========="
Rscript -e '
suppressPackageStartupMessages(library(readxl))
for (f in c("data/raw/GSE261348/GSE261348_IMfirst_DSP_rawcounts.xlsx",
            "data/raw/GSE261345/GSE261345_CANTABRICO_DSP_rawcounts.xlsx")) {
  cat("\n---", basename(f), "---\n")
  sh <- excel_sheets(f)
  cat("sheets:", length(sh), "\n")
  for (s in sh) {
    d <- tryCatch(read_excel(f, sheet = s, n_max = 3), error = function(e) NULL)
    if (is.null(d)) { cat(sprintf("  %-28s <unreadable>\n", s)); next }
    full <- tryCatch(read_excel(f, sheet = s), error = function(e) NULL)
    nr <- if (is.null(full)) NA else nrow(full)
    cat(sprintf("  %-28s cols=%-5d rows=%-6s first cols: %s\n",
                s, ncol(d), nr,
                paste(head(names(d), 4), collapse=" | ")))
  }
}
' 2>&1 | grep -v "out-of-sync"

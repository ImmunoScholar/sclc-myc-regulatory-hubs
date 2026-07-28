#!/usr/bin/env bash
cd /home/priya/projects/sclc-myc-regulatory-hubs || exit 1
rm -f .fix_renv.R

echo "=========== renv gap closed? ==========="
Rscript -e '
  lock <- names(renv::lockfile_read()$Packages)
  dep  <- unique(renv::dependencies(progress=FALSE, dev=FALSE)$Package)
  base <- rownames(installed.packages(priority=c("base","recommended")))
  gap  <- setdiff(setdiff(dep, lock), base)
  gap  <- gap[vapply(gap, function(p) requireNamespace(p, quietly=TRUE), logical(1))]
  cat("used-but-unrecorded packages:", if(length(gap)) paste(gap, collapse=", ") else "NONE", "\n")
' 2>&1 | grep -viE "out-of-sync|inconsistent|^ *package|^ *[A-Za-z0-9._]+ +y +|See .renv::status|^$"

echo
echo "=========== repos configured (restore depends on these) ==========="
Rscript -e 'print(getOption("repos"))' 2>&1 | grep -v out-of-sync

echo
echo "=========== the TSS resize in 20_peak_to_gene.R ==========="
grep -n "resize\|promoters(\|trim\|seqlengths\|out-of-bound" scripts/02_regulatory/20_peak_to_gene.R

echo
echo "=========== same pattern elsewhere? ==========="
grep -rn "resize(" scripts/ --include=*.R | grep -v "^scripts/02_regulatory/20_peak_to_gene.R"

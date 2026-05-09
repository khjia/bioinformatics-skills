#!/usr/bin/env Rscript
# 01_prepare_data.R — Convert PLINK BED to rMVP big.matrix (run once).
# Usage: Rscript 01_prepare_data.R
#   Assumes the working dir is <proj>/05.gwas/rMVP/.
#   Edit `bed_prefix` if your BED lives elsewhere.

suppressPackageStartupMessages(library(rMVP))

gwas_dir   <- ".."
bed_prefix <- file.path(gwas_dir, "genotype_full_bed", "full_snps_140")  # EDIT
out_dir    <- "."

MVP.Data(fileBed = bed_prefix,
         fileOut = file.path(out_dir, "mvp"),
         verbose = TRUE)

cat("Data conversion complete.\n")
print(list.files(out_dir, pattern = "^mvp\\."))

#!/usr/bin/env Rscript
# 02_run_rMVP.R — run GLM + MLM + FarmCPU for a slice of traits.
# Usage: Rscript 02_run_rMVP.R <trait_start_idx> <trait_end_idx>
#   Indices are 1-based into the phenotype columns (excluding FID/IID).
#
# Inputs (relative to this script's dir = <proj>/05.gwas/rMVP):
#   mvp.geno.{bin,desc,ind,map}        # from 01_prepare_data.R
#   ../pheno_filtered.tsv              # FID IID trait_01 ...
#   ../regenie_covariates/quant_pc3.tsv# FID IID PC1 PC2 PC3
# Outputs:
#   results/trait_XX.{GLM,MLM,FarmCPU}.csv
#   kinship.rds (cached)

suppressPackageStartupMessages({
    library(rMVP); library(bigmemory)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 02_run_rMVP.R <start> <end>")
trait_start <- as.integer(args[1]); trait_end <- as.integer(args[2])

# ---- genotype -------------------------------------------------------------
geno <- attach.big.matrix("mvp.geno.desc")
map  <- read.table("mvp.geno.map", header = TRUE, stringsAsFactors = FALSE)
ind_raw <- readLines("mvp.geno.ind")

# rMVP strips leading zeros; restore 3-digit zero-padding to match fam/pheno.
# Adjust width if your project uses a different format.
ind <- sprintf("%03d", as.integer(ind_raw))

# ---- phenotype ------------------------------------------------------------
pheno_all <- read.delim("../pheno_filtered.tsv", stringsAsFactors = FALSE,
                        colClasses = c("FID" = "character", "IID" = "character"))
pheno_all <- pheno_all[match(ind, pheno_all$IID), ]
if (any(is.na(pheno_all$IID))) stop("IID mismatch between geno and pheno")
trait_names <- colnames(pheno_all)[-(1:2)]

# ---- PC covariates --------------------------------------------------------
pc_data <- read.delim("../regenie_covariates/quant_pc3.tsv", stringsAsFactors = FALSE,
                      colClasses = c("FID" = "character", "IID" = "character"))
pc_data <- pc_data[match(ind, pc_data$IID), ]
if (any(is.na(pc_data$IID))) stop("IID mismatch between geno and PC")
CV <- as.matrix(pc_data[, c("PC1", "PC2", "PC3")])

# ---- kinship --------------------------------------------------------------
K_path <- "kinship.rds"
if (file.exists(K_path)) {
    cat("Loading cached kinship ...\n"); K <- readRDS(K_path)
} else {
    cat("Computing VanRaden kinship ...\n")
    K <- MVP.K.VanRaden(geno, verbose = TRUE); saveRDS(K, K_path)
}

# ---- run GWAS -------------------------------------------------------------
dir.create("results", showWarnings = FALSE)
setwd("results")   # rMVP writes to cwd
bonf <- 0.05 / nrow(map)
ncpu <- as.integer(Sys.getenv("OMP_NUM_THREADS", "16"))

for (i in trait_start:min(trait_end, length(trait_names))) {
    tname <- trait_names[i]
    y <- as.numeric(pheno_all[[tname]])
    n_valid <- sum(!is.na(y))
    cat(sprintf("\n[%d/%d] %s  n_valid=%d\n", i, trait_end, tname, n_valid))
    if (n_valid < 50) { cat("  skip (n_valid < 50)\n"); next }

    pheno_df <- data.frame(IID = ind, y = y, stringsAsFactors = FALSE)
    colnames(pheno_df)[2] <- tname

    tryCatch({
        MVP(
            phe        = pheno_df,
            geno       = geno,
            map        = map,
            K          = K,
            CV.MLM     = CV,
            CV.FarmCPU = CV,
            nPC.GLM    = 3,
            maxLoop    = 10,
            method.bin = "static",
            threshold  = bonf,
            method     = c("GLM", "MLM", "FarmCPU"),
            ncpus      = ncpu,
            file.output = TRUE,
            file.type  = "csv",
            outpath    = ".",
            plot.type  = NULL,          # Python post-processor makes plots
            verbose    = FALSE
        )
    }, error = function(e) cat("  ERR:", conditionMessage(e), "\n"))
}
cat("\nDone range", trait_start, "-", trait_end, "\n")

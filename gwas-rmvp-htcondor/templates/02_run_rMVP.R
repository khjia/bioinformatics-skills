#!/usr/bin/env Rscript
# 02_run_rMVP.R — main GWAS driver
# Usage: Rscript 02_run_rMVP.R <trait_start> <trait_end>
#
# Features:
#   - Kinship cache with genotype fingerprint validation (auto-rebuild on BED change)
#   - Per-trait retry (3 attempts) — survives transient IO / memory hiccups
#   - 4 models in one pass: GLM, MLM, FarmCPU, BLINK
#   - PC count via env var RMVP_NPC (default 3) for adaptive tuning sweeps
#   - English-only trait labels in output filenames (CJK-safe)

suppressPackageStartupMessages({
  library(rMVP)
  library(bigmemory)
  library(tools)
})

args        <- commandArgs(trailingOnly = TRUE)
trait_start <- as.integer(args[1])
trait_end   <- as.integer(args[2])
npc         <- as.integer(Sys.getenv("RMVP_NPC", "3"))
ncpus       <- as.integer(Sys.getenv("OMP_NUM_THREADS", "16"))
max_retry   <- as.integer(Sys.getenv("RMVP_RETRY", "3"))
outdir      <- Sys.getenv("RMVP_OUTDIR", "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[CFG] trait %d-%d | nPC=%d | ncpus=%d | retry=%d | out=%s\n",
            trait_start, trait_end, npc, ncpus, max_retry, outdir))

# ---- genotype --------------------------------------------------------------
geno <- attach.big.matrix("mvp.geno.desc")
map  <- read.table("mvp.geno.map", header = TRUE, stringsAsFactors = FALSE)

# rMVP strips leading zeros from fam IIDs — restore zero-padding to 3 digits
ind_raw <- readLines("mvp.geno.ind")
ind     <- sprintf("%03d", as.integer(ind_raw))

# ---- phenotype + PC covariates --------------------------------------------
pheno_all <- read.delim("../pheno_filtered.tsv", stringsAsFactors = FALSE,
                        colClasses = c("FID" = "character", "IID" = "character"))
pheno_all <- pheno_all[match(ind, pheno_all$IID), ]
trait_names <- colnames(pheno_all)[-(1:2)]

pc_data <- read.delim("../regenie_covariates/quant_pc3.tsv", stringsAsFactors = FALSE,
                      colClasses = c("FID" = "character", "IID" = "character"))
pc_data <- pc_data[match(ind, pc_data$IID), ]
CV <- if (npc == 0) NULL else as.matrix(pc_data[, paste0("PC", 1:npc)])

# ---- kinship cache with fingerprint ----------------------------------------
# Hash the genotype backing file so a BED change auto-invalidates K.
geno_bin <- "mvp.geno.bin"
fp_now <- if (file.exists(geno_bin)) {
  paste(file.info(geno_bin)$size,
        format(file.info(geno_bin)$mtime, "%s"),
        md5sum(geno_bin), sep = "_")
} else {
  stop("mvp.geno.bin missing — run 01_prepare_data.R first")
}

K_path  <- "kinship.rds"
fp_path <- "kinship.fp"
fp_old  <- if (file.exists(fp_path)) readLines(fp_path, n = 1) else ""

if (file.exists(K_path) && identical(fp_old, fp_now)) {
  cat("[K] cache hit, loading kinship.rds\n")
  K <- readRDS(K_path)
} else {
  if (file.exists(K_path)) {
    cat("[K] genotype changed (fingerprint mismatch) — rebuilding\n")
  } else {
    cat("[K] no cache — building VanRaden kinship\n")
  }
  K <- MVP.K.VanRaden(geno, verbose = TRUE)
  saveRDS(K, K_path)
  writeLines(fp_now, fp_path)
}

# ---- run loop with per-trait retry -----------------------------------------
bonf   <- 0.05 / nrow(map)
models <- c("GLM", "MLM", "FarmCPU", "BLINK")

run_one <- function(tname, y) {
  pheno_df <- data.frame(IID = ind, y, stringsAsFactors = FALSE)
  colnames(pheno_df)[2] <- tname
  MVP(
    phe        = pheno_df,
    geno       = geno,
    map        = map,
    K          = K,
    CV.GLM     = CV,
    CV.MLM     = CV,
    CV.FarmCPU = CV,
    nPC.GLM    = npc,
    nPC.MLM    = npc,
    nPC.FarmCPU = npc,
    maxLoop    = 10,
    method.bin = "static",
    threshold  = bonf,
    method     = models,
    ncpus      = ncpus,
    file.output = TRUE,
    file.type  = "csv",
    outpath    = outdir,
    verbose    = FALSE
  )
}

for (i in trait_start:min(trait_end, length(trait_names))) {
  tname <- trait_names[i]
  y     <- as.numeric(pheno_all[[tname]])
  n_ok  <- sum(!is.na(y))

  if (n_ok < 50) {
    cat(sprintf("[SKIP] %s n_valid=%d (<50)\n", tname, n_ok))
    next
  }

  done <- FALSE
  for (attempt in 1:max_retry) {
    t0 <- Sys.time()
    res <- tryCatch(
      { run_one(tname, y); TRUE },
      error = function(e) {
        cat(sprintf("[ERR][try %d/%d] %s: %s\n",
                    attempt, max_retry, tname, conditionMessage(e)))
        FALSE
      }
    )
    if (isTRUE(res)) {
      cat(sprintf("[OK] %s n=%d in %.1fs\n",
                  tname, n_ok, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      done <- TRUE
      break
    }
    Sys.sleep(30)
    gc(verbose = FALSE)
  }
  if (!done) {
    cat(sprintf("[FAIL] %s after %d attempts — moving on\n", tname, max_retry))
  }
}
cat("[DONE] batch finished\n")

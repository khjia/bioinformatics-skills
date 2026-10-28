## 02_make_inputs.R ---- build Fast3VmrMLM input files (trait.csv, PCA.csv)
## Env: PROJECT, PHENO_TSV, FAM, EIGENVEC, TRAIT_CSV, PCA_CSV, NPC
suppressMessages(library(data.table))

PROJECT   <- Sys.getenv("PROJECT")
stopifnot(nzchar(PROJECT))
PHENO_TSV <- Sys.getenv("PHENO_TSV", file.path(PROJECT, "01.phenotype/pheno_all.tsv"))
FAM       <- Sys.getenv("FAM",       file.path(PROJECT, "02.ld_filter/beds/full_snps_renamed.fam"))
EIGENVEC  <- Sys.getenv("EIGENVEC",  file.path(PROJECT, "03.pca/pca.eigenvec"))
TRAIT_CSV <- Sys.getenv("TRAIT_CSV", file.path(PROJECT, "06.Fast3VmrMLM/inputs/trait.csv"))
PCA_CSV   <- Sys.getenv("PCA_CSV",   file.path(PROJECT, "06.Fast3VmrMLM/inputs/PCA.csv"))
NPC       <- as.integer(Sys.getenv("NPC", "3"))

dir.create(dirname(TRAIT_CSV), showWarnings = FALSE, recursive = TRUE)

## ---- phenotype aligned to genotype order ----
phe <- fread(PHENO_TSV)
setnames(phe, 1:2, c("FID", "IID"))
phe[, IID := as.character(IID)]
fam <- fread(FAM, header = FALSE)
ids <- as.character(fam$V2)

miss <- setdiff(ids, phe$IID)
if (length(miss)) stop("samples missing from phenotype: ", paste(head(miss, 5), collapse = ", "))
phe <- phe[match(ids, phe$IID)]
trait_cols <- setdiff(names(phe), c("FID", "IID"))

out <- data.table(`<Phenotype>` = ids)
out <- cbind(out, phe[, ..trait_cols])
fwrite(out, TRAIT_CSV, sep = ",", na = "NA")
cat(sprintf("trait.csv: %d samples x %d traits\n", nrow(out), length(trait_cols)))

## ---- PCA (first NPC PCs) ----
ev <- fread(EIGENVEC, skip = 1, header = FALSE)
pc <- ev[, 1:(2 + NPC)]
setnames(pc, c("FID", "IID", paste0("pc_", seq_len(NPC))))
pc[, IID := as.character(IID)]
if (!identical(pc$IID, ids)) {
  pc <- pc[match(ids, pc$IID)]
  if (anyNA(pc$IID)) stop("PCA samples do not match genotype samples")
}
pc_out <- data.table(`<pca>` = pc$IID)
pc_out <- cbind(pc_out, pc[, paste0("pc_", seq_len(NPC)), with = FALSE])
fwrite(pc_out, PCA_CSV, sep = ",", na = "NA")
cat(sprintf("PCA.csv: %d samples x %d PCs\n", nrow(pc_out), NPC))
cat("DONE 02\n")

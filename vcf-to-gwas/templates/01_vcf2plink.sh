#!/usr/bin/env bash
## 01_vcf2plink.sh ---- VCF -> PLINK bed + LD-pruned set + PCA
## Env: PROJECT (required), VCF, THREADS, GENO, MAF, PRUNE_WIN, PRUNE_STEP, PRUNE_R2, NPCS
set -euo pipefail

: "${PROJECT:?set PROJECT to the analysis root}"
VCF=${VCF:-$PROJECT/dataset3.recode.vcf}
THREADS=${THREADS:-16}
GENO=${GENO:-0.2}
MAF=${MAF:-0.05}
PRUNE_WIN=${PRUNE_WIN:-50}
PRUNE_STEP=${PRUNE_STEP:-10}
PRUNE_R2=${PRUNE_R2:-0.2}
NPCS=${NPCS:-10}

BEDS=$PROJECT/02.ld_filter/beds
PCAD=$PROJECT/03.pca
mkdir -p "$BEDS" "$PCAD"

echo "[1/5] VCF -> bed (geno=$GENO, maf=$MAF)"
plink2 --vcf "$VCF" --allow-extra-chr --geno "$GENO" --maf "$MAF" \
       --make-bed --threads "$THREADS" --out "$BEDS/full_snps"

echo "[2/5] standardise variant IDs (chr:pos:ref:alt)"
plink2 --bfile "$BEDS/full_snps" --allow-extra-chr \
       --set-all-var-ids '@:#:$r:$a' --make-bed --threads "$THREADS" \
       --out "$BEDS/full_snps_renamed"

echo "[3/5] LD pruning (indep-pairwise $PRUNE_WIN $PRUNE_STEP $PRUNE_R2)"
plink2 --bfile "$BEDS/full_snps" --allow-extra-chr \
       --set-all-var-ids '@:#:$r:$a' \
       --indep-pairwise "$PRUNE_WIN" "$PRUNE_STEP" "$PRUNE_R2" \
       --threads "$THREADS" --out "$BEDS/ld_pruned"

echo "[4/5] build pruned bed for PCA"
plink2 --bfile "$BEDS/full_snps_renamed" --allow-extra-chr \
       --extract "$BEDS/ld_pruned.prune.in" --make-bed --threads "$THREADS" \
       --out "$PCAD/ld_pruned_snps"

echo "[5/5] PCA ($NPCS PCs)"
plink2 --bfile "$PCAD/ld_pruned_snps" --allow-extra-chr \
       --pca "$NPCS" --threads "$THREADS" --out "$PCAD/pca"

echo
echo "SNPs: $(wc -l < "$BEDS/full_snps_renamed.bim")  |  pruned: $(wc -l < "$PCAD/ld_pruned_snps.bim")  |  samples: $(wc -l < "$BEDS/full_snps_renamed.fam")"
echo "Chromosomes: $(cut -f1 "$BEDS/full_snps_renamed.bim" | sort -u | tr '\n' ' ')"
echo "DONE 01"

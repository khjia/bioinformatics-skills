#!/usr/bin/env bash
## 34_ldblock_trait.sh ---- LDBlockShow locus panel for one trait
## Usage: bash 34_ldblock_trait.sh fa|vee <trait>
set -euo pipefail

PREFIX=$1; TR=$2
WDIR=${GWAS_WDIR:-/media/nfs1/hermes/project/yuanbaofeng/09.postgwas2}
PROJ=/media/nfs1/hermes/project/yuanbaofeng
PLINK2=/media/nfs1/hermes/miniforge3/bin/plink2
LDB=${LDBSHOW:-$WDIR/bin/LDBlockShow}
[ -x "$LDB" ] || LDB=/tmp/LDBlockShow-main/bin/LDBlockShow
BFILE=$PROJ/02.ld_filter/beds/full_snps_renamed
SUMMARY=$WDIR/results/${PREFIX}_${TR}_report/tab_summary.tsv
if [ "$PREFIX" = "fa" ]; then MID=$PROJ/06.Fast3VmrMLM/results/${TR}_midresult.csv; else MID=$PROJ/08.vee/results/${TR}_midresult.csv; fi
GFF=$PROJ/final.gene.mRNA.longest.gff3
FIG=$WDIR/figures/${PREFIX}_${TR}_report
mkdir -p "$FIG"

[ -s "$SUMMARY" ] || { echo "no summary: $TR"; exit 0; }
read -r CHR APOS < <(awk -F'\t' 'NR==2{print $3, int($4)}' "$SUMMARY")
[ -n "${CHR:-}" ] || { echo "no anchor: $TR"; exit 0; }

FROM=$(( APOS > 60000 ? APOS - 60000 : 1 ))
TO=$(( APOS + 60000 ))
WORK=$(mktemp -d /tmp/ldb.XXXXXX)
trap 'rm -rf "$WORK"' EXIT

"$PLINK2" --bfile "$BFILE" --allow-extra-chr --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
  --maf 0.1 --max-alleles 2 --snps-only --geno 0.25 --export vcf --out "$WORK/region" >/dev/null 2>&1 || true
NSNP=$(grep -vc '^#' "$WORK/region.vcf" 2>/dev/null || echo 0)
if [ "$NSNP" -lt 2 ]; then
  "$PLINK2" --bfile "$BFILE" --allow-extra-chr --chr "$CHR" --from-bp "$FROM" --to-bp "$TO" \
    --maf 0.05 --max-alleles 2 --snps-only --geno 0.25 --export vcf --out "$WORK/region" >/dev/null 2>&1 || true
  NSNP=$(grep -vc '^#' "$WORK/region.vcf" 2>/dev/null || echo 0)
fi
if [ "$NSNP" -lt 2 ]; then echo "too few SNP in $CHR:$FROM-$TO ($NSNP)"; exit 0; fi

awk -F, -v c="$CHR" -v f="$FROM" -v t="$TO" 'NR>1 && $2==c && $3>=f && $3<=t {print $2"\t"$3"\t"$5}' "$MID" > "$WORK/gwas.txt"
awk -F'\t' -v c="$CHR" -v f="$FROM" -v t="$TO" '$1==c && $4>=f && $5<=t && ($3=="gene" || $3=="exon" || $3=="mRNA")' "$GFF" > "$WORK/genes.gff3"

"$LDB" -InVCF "$WORK/region.vcf" -OutPut "$WORK/ld" -Region "${CHR}:${FROM}-${TO}" \
  -InGWAS "$WORK/gwas.txt" -InGFF "$WORK/genes.gff3" -SeleVar 2 -BlockType 5 > "$WORK/ld.log" 2>&1 || true
if [ -s "$WORK/ld.svg" ]; then
  cp "$WORK/ld.svg" "$FIG/ldblock.svg"
  echo "OK $TR $CHR:$FROM-$TO SNP=$NSNP"
else
  echo "LDB failed $TR"; tail -3 "$WORK/ld.log"
fi

#!/usr/bin/env bash
## 40_compute_all.sh ---- run every per-trait step for all 22 traits
## Usage: bash 40_compute_all.sh
set -u
WDIR=${GWAS_WDIR:-/media/nfs1/hermes/project/yuanbaofeng/09.postgwas2}
export GWAS_WDIR=$WDIR
export PATH=$WDIR/bin:/media/nfs1/hermes/miniforge3/bin:$PATH
R=/media/nfs1/hermes/miniforge3/bin/Rscript
LOG=$WDIR/logs/traits_compute_$(date +%Y%m%d_%H%M).log
mkdir -p "$WDIR/logs"

FA="oil_g100g protein_g100g seed_weight_g fa_C16_0 fa_C16_1 fa_C18_0 fa_C18_1 fa_C18_2 fa_C18_3n6 fa_C18_3n3 fa_C20_0 fa_C20_1 fa_C20_2 fa_C22_0 fa_C22_1 fa_C24_0 fa_C24_1"
VEE="vee_alpha vee_beta vee_gamma vee_delta vee_total"

run_one() {
  local p=$1 tr=$2
  echo "===== [$p $tr] $(date +%H:%M:%S)" | tee -a "$LOG"
  $R "$WDIR/scripts/20_trait_report_figs.R" "$p" "$tr" >>"$LOG" 2>&1 || echo "FAIL figs $tr" | tee -a "$LOG"
  $R "$WDIR/scripts/31_susie_trait.R"      "$p" "$tr" >>"$LOG" 2>&1 || echo "FAIL susie $tr" | tee -a "$LOG"
  $R "$WDIR/scripts/32_variant_annotation_trait.R" "$p" "$tr" >>"$LOG" 2>&1 || echo "FAIL annot $tr" | tee -a "$LOG"
  $R "$WDIR/scripts/33_known_loci.R"       "$p" "$tr" >>"$LOG" 2>&1 || echo "FAIL known $tr" | tee -a "$LOG"
  bash "$WDIR/scripts/34_ldblock_trait.sh" "$p" "$tr" >>"$LOG" 2>&1 || echo "FAIL ld $tr" | tee -a "$LOG"
}

for tr in $FA;  do run_one fa  "$tr"; done
for tr in $VEE; do run_one vee "$tr"; done
echo "ALL COMPUTE DONE $(date)" | tee -a "$LOG"

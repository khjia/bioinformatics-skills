#!/usr/bin/env bash
## run_all.sh ---- VCF -> PLINK/QC -> PCA -> Fast3VmrMLM -> QTN annotation
## Usage:
##   export PROJECT=/path/to/project
##   bash run_all.sh
## Optional env: VCF, THREADS, GENO, MAF, PRUNE_WIN, PRUNE_STEP, PRUNE_R2,
##               NPC, TRAITS, SVRAD, SVPAL, SVMLOD, RSCRIPT, PY,
##               RUN_REPORTS=1 (render per-trait HTML; needs 09.postgwas2/scripts/00_config.R)
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
: "${PROJECT:?set PROJECT to the analysis root}"
RSCRIPT=${RSCRIPT:-Rscript}
PY=${PY:-python3}
THREADS=${THREADS:-16}

echo "=== 01 VCF -> PLINK + PCA ==="
bash "$HERE/01_vcf2plink.sh"

echo "=== 02 Fast3VmrMLM inputs ==="
"$RSCRIPT" "$HERE/02_make_inputs.R"

echo "=== 02b link genotype into 06.Fast3VmrMLM/inputs ==="
IN=$PROJECT/06.Fast3VmrMLM/inputs
mkdir -p "$IN"
for ext in bed bim fam; do
  ln -sf "../../02.ld_filter/beds/full_snps_renamed.$ext" "$IN/genotype.$ext"
done

echo "=== 03 Fast3VmrMLM ==="
"$RSCRIPT" "$HERE/03_run_fast3vmrmlm.R"

echo "=== 04 annotate QTN ==="
"$PY" "$HERE/04_annotate_qtl.py"

echo "=== 05 gene coordinates ==="
"$RSCRIPT" "$HERE/05_gene_coords.R"

## ---- optional: per-trait HTML reports ----
if [ "${RUN_REPORTS:-0}" = "1" ]; then
  echo "=== 06 install report scripts + render per-trait HTML ==="
  W=$PROJECT/09.postgwas2
  mkdir -p "$W/scripts" "$W/report"
  cp "$HERE/20_trait_report_figs.R" "$W/scripts/"
  cp "$HERE/22_render_all.R"        "$W/scripts/"
  cp "$HERE/21_trait_report.Rmd"    "$W/report/"
  if [ ! -f "$W/scripts/00_config.R" ]; then
    echo "ERROR: create $W/scripts/00_config.R from templates/00_config.example.R first" >&2
    exit 1
  fi
  GWAS_WDIR=$W RSCRIPT=$RSCRIPT "$RSCRIPT" "$W/scripts/22_render_all.R"
fi

echo
echo "pipeline done."
if [ "${RUN_REPORTS:-0}" != "1" ]; then
  echo "reports: copy templates/20…22 into <PROJECT>/09.postgwas2 and run 22_render_all.R"
fi

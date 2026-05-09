#!/bin/bash
# run_batch.sh — HTCondor wrapper
# Usage: run_batch.sh <trait_start> <trait_end>
set -eo pipefail

export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16

# EDIT absolute path
cd /abs/path/to/project/05.gwas/rMVP
/media/nfs1/hermes/miniforge3/bin/Rscript 02_run_rMVP.R "$1" "$2"

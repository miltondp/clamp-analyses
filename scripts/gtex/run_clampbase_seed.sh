#!/usr/bin/env bash
# Run an independent GTEx CLAMPbase replicate.
set -euo pipefail
seed=${1:?rng seed required}
index=${2:?replicate index required}
out=${3:?output directory required}
inputs=${GTEX_INPUTS:-output/01_model_building/01_gtex}
threads=${CLAMP_THREADS:-8}
export OPENBLAS_NUM_THREADS=$threads MKL_NUM_THREADS=$threads OMP_NUM_THREADS=$threads
Rscript scripts/gtex/fit_clampbase_seed.R \
  --fbm-rds "$inputs/gtex_fbm_filt.rds" --backing "$inputs/FBMgtex_preproc_filtered.bk" \
  --genes-rds "$inputs/gtex_genes.rds" --samples-rds "$inputs/gtex_samples.rds" \
  --out-dir "$out" --seed "$seed" --seed-index "$index" \
  --expected-k "$(Rscript -e 'cat(as.integer(readRDS(commandArgs(TRUE)[1])))' "$inputs/CLAMP_K_gtex.rds")"

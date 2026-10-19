#!/bin/bash
# Run the CLAMPbase coverage trait workflow.
#SBATCH -J cb_cov
#SBATCH -c 12
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo"
exec conda run --no-capture-output -n snakemake snakemake \
  --snakefile workflow/Snakefile --use-conda --cores "${SLURM_CPUS_PER_TASK:-12}" phenoplier_coverage "$@"

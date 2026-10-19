#!/bin/bash
# Run the fixed-rank CLAMPbase trait workflow.
#SBATCH -J cb_sat
#SBATCH -c 12
#SBATCH --mem=64G
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo"
exec conda run --no-capture-output -n snakemake snakemake \
  --snakefile workflow/Snakefile --use-conda --cores "${SLURM_CPUS_PER_TASK:-12}" phenoplier_saturation "$@"

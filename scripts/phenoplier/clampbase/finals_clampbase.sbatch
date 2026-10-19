#!/bin/bash
# Run the full-compendium CLAMPbase trait workflow.
#SBATCH -J cb_final
#SBATCH -c 12
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo"
exec conda run --no-capture-output -n snakemake snakemake \
  --snakefile workflow/Snakefile --use-conda --cores "${SLURM_CPUS_PER_TASK:-12}" phenoplier_finals "$@"

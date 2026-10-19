#!/bin/bash
# Run PhenoPLIER trait associations for metaLVs.
#SBATCH -J metalv_gls
#SBATCH -c 32
#SBATCH --mem=150G
#SBATCH -o %x_%j.out
set -euo pipefail

INPUT="${1:?usage: run_gls_metalv.sbatch <metalv_z.pkl> [name]}"
NAME="${2:-metalv}"
CONDA_SH="${PHENOPLIER_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
ENV_NAME="${PHENOPLIER_ENV:-phenoplier-cli-neo}"
WS="${PHENOPLIER_HOME:?set PHENOPLIER_HOME to a writable PhenoPLIER workspace}"

source "$CONDA_SH"
conda activate "$ENV_NAME"
export PHENOPLIER_HOME="$WS" PHENOPLIER_ROOT_DIR="$WS"
export XDG_CACHE_HOME="$WS/../cache" TMPDIR="$WS/../tmp" MPLCONFIGDIR="$WS/../cache/mpl" NUMBA_CACHE_DIR="$WS/../cache/numba"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$TMPDIR"
cd "$WS/.."

echo "[start] $(date -Is) host=$(hostname) version=$(phenoplier -v 2>&1 | tail -1)"
phenoplier shortcut gls \
  --input "$INPUT" \
  --name "$NAME" \
  --model-namespace clamp \
  --lv-percentile 0.01 \
  --trait-filter biomedical \
  --executor local \
  --n-jobs "${SLURM_CPUS_PER_TASK:-32}"
echo "[done] $(date -Is)"

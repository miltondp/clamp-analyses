#!/usr/bin/env bash
# Snakemake activates the declared environment. --conda-env is a manual override.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
args=("$@")
conda_env=""
while [[ $# -gt 0 ]]; do
  if [[ "$1" == --conda-env ]]; then conda_env="$2"; shift; fi
  shift
done
if [[ -n "$conda_env" ]]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$conda_env" || exit 1
fi
if [[ -n "${CONDA_PREFIX:-}" ]]; then
  export PATH="$CONDA_PREFIX/bin:$PATH"
fi
exec python "$here/workflow.py" store "${args[@]}"

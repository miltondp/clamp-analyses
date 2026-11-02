#!/usr/bin/env bash
# Activate the isolated CLI environment before invoking the workflow adapter.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
args=("$@")
conda_env=phenoplier-cli-neo
while [[ $# -gt 0 ]]; do
  if [[ "$1" == --conda-env ]]; then conda_env="$2"; shift; fi
  shift
done
if [[ -d "$conda_env" ]]; then
  export PATH="$(cd "$conda_env" && pwd)/bin:$PATH"
else
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$conda_env"
  export PATH="$CONDA_PREFIX/bin:$PATH"
fi
exec python "$here/workflow.py" gls "${args[@]}"

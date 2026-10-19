#!/usr/bin/env bash
# Run trait aggregation in the PhenoPLIER environment.
set -euo pipefail

conda_env="$1"
shift

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$script_dir/aggregate_traits.py" "$@"

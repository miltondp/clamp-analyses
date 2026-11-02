#!/usr/bin/env bash
# Optional named environment for manual CLI use; Snakemake manages its own prefix
# from the same YAML. Existing environments are never silently updated.
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
target_env="${PHENOPLIER_TARGET_ENV:-phenoplier-cli-neo}"
conda env create --name "$target_env" --file "$repo_root/envs/phenoplier-cli-neo.yaml"
conda run --no-capture-output -n "$target_env" \
  python "$repo_root/scripts/phenoplier/check_env.py" --require-version 0.7.3
printf '\nEnvironment ready. See scripts/phenoplier/README.md for workspace/reference setup.\n'

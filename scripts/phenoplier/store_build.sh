#!/usr/bin/env bash
# Combine a fitted model and its traits into a study store.
set -euo pipefail

usage() {
  cat >&2 <<USAGE
usage: $0 --rds <model.rds> --name <project> --store-out <out.h5> [--conda-env NAME] [--require-version X.Y.Z] [--cohort NAME] [--workspace PATH]
Workspace: --workspace, else \$PHENOPLIER_HOME, else \$HOME/phenoplier.
USAGE
  exit 2
}

conda_env=phenoplier-cli-neo
require_version=""
cohort=phenomexcan_rapid_gwas
rds=""
name=""
store_out=""
workspace_arg=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rds) rds="$2"; shift 2 ;;
    --name) name="$2"; shift 2 ;;
    --store-out) store_out="$2"; shift 2 ;;
    --workspace) workspace_arg="$2"; shift 2 ;;
    --conda-env) conda_env="$2"; shift 2 ;;
    --require-version) require_version="$2"; shift 2 ;;
    --cohort) cohort="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done
[[ -n "$rds" && -n "$name" && -n "$store_out" ]] || usage
[[ -f "$rds" ]] || { echo "model not found: $rds" >&2; exit 1; }

workspace="${workspace_arg:-${PHENOPLIER_HOME:-$HOME/phenoplier}}"
export PHENOPLIER_HOME="$workspace" PHENOPLIER_ROOT_DIR="$workspace"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"
export PATH="${CONDA_PREFIX:-}/bin:${PATH}"

installed="$(phenoplier --version 2>&1 | tail -1 | sed -E 's/.*v([0-9][0-9A-Za-z.+-]*).*/\1/')"
if [[ -n "$require_version" && "$installed" != "$require_version" ]]; then
  echo "[ERROR] phenoplier-cli $installed in env '$conda_env' but $require_version is required" \
    "(workflow/config/phenoplier.yaml: version)" >&2
  exit 1
fi

project="$workspace/projects/$name"
gls_dir="$project/results/gls/phenoplier/${cohort}"
if [[ ! -d "$gls_dir" ]]; then
  echo "[ERROR] $name: no GLS results dir at $gls_dir" >&2
  exit 1
fi

mkdir -p "$(dirname "$store_out")"
tmp="${store_out}.tmp.$$"
rm -f "$tmp"
echo "[store] $(date -Is)  $name -> $store_out"
phenoplier store build \
  --output "$tmp" \
  --cohort "$cohort" \
  --clamp-rds "$rds" \
  --gls-dir "$gls_dir" \
  --method-name CLAMP
python -c "import h5py, sys; h5py.File(sys.argv[1], 'r').close()" "$tmp"
mv -f "$tmp" "$store_out"
echo "[done] $(date -Is)  $name"

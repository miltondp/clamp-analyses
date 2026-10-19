#!/usr/bin/env bash
# Create analysis environments and install the pinned CLAMP package.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

MINICONDA_PREFIX="${MINICONDA_PREFIX:-$HOME/miniconda3}"

log() { printf '\n[setup.sh] %s\n' "$1"; }

if command -v conda >/dev/null 2>&1; then
    log "conda already installed: $(command -v conda)"
    CONDA_BASE="$(conda info --base)"
else
    log "conda not found, installing Miniconda to $MINICONDA_PREFIX"
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64) installer="Miniconda3-latest-Linux-x86_64.sh" ;;
        *)
            echo "[setup.sh] unsupported platform $(uname -s)-$(uname -m)." >&2
            echo "[setup.sh] this repo's env lock files are pinned to linux-64;" >&2
            echo "[setup.sh] install Conda/Mambaforge yourself and re-run this script." >&2
            exit 1
            ;;
    esac
    tmp_installer="$(mktemp -d)/miniconda.sh"
    curl -fsSL "https://repo.anaconda.com/miniconda/${installer}" -o "$tmp_installer"
    bash "$tmp_installer" -b -p "$MINICONDA_PREFIX"
    rm -f "$tmp_installer"
    CONDA_BASE="$MINICONDA_PREFIX"
fi

source "$CONDA_BASE/etc/profile.d/conda.sh"

env_exists() {
    conda env list | awk '{print $1}' | grep -qx "$1"
}

create_env_from_lock() {
    local name="$1" lock="$2"
    if env_exists "$name"; then
        log "env '$name' already exists, skipping (remove it manually with 'conda env remove -n $name' to recreate)"
    else
        log "creating env '$name' from $lock"
        conda create --name "$name" --file "$lock" -y
    fi
}

create_env_from_yaml() {
    local name="$1" yaml="$2"
    if env_exists "$name"; then
        log "env '$name' already exists, skipping (remove it manually with 'conda env remove -n $name' to recreate)"
    else
        log "creating env '$name' from $yaml"
        conda env create -n "$name" -f "$yaml"
    fi
}

create_env_from_lock clamp-analyses envs/clamp-analyses.lock
create_env_from_lock gpu-kmeans envs/gpu-kmeans.lock
create_env_from_yaml snakemake envs/snakemake.yaml
conda run -n clamp-analyses python -c '
import subprocess
import sys
import yaml
with open(sys.argv[1]) as handle:
    dependencies = yaml.safe_load(handle)["dependencies"]
requirements = [requirement for dependency in dependencies if isinstance(dependency, dict)
                for requirement in dependency.get("pip", [])]
if not requirements:
    raise SystemExit("No pip requirements found in the environment YAML")
subprocess.check_call([sys.executable, "-m", "pip", "install", *requirements])
' envs/clamp-analyses.yaml

log "WARNING: using the repository-pinned CLAMP revision; do not update CLAMP independently"
conda run -n clamp-analyses Rscript scripts/install_clamp.R

log "setup complete. Stage the external inputs listed in README.md before running the workflow. Next steps:"
cat <<'EOF'

  conda activate snakemake
  snakemake --cores 4 --use-conda --snakefile workflow/Snakefile all

To verify everything (envs, Snakemake, data files), run:
  nbs/00_setup/00_check_setup.ipynb
EOF

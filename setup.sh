#!/usr/bin/env bash
# Create analysis environments and install the pinned CLAMP package.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

if [[ "$(uname -s)-$(uname -m)" != "Linux-x86_64" ]]; then
    echo "[setup.sh] explicit environment locks require Linux x86-64." >&2
    exit 1
fi

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
        conda create --name "$name" --file "$lock" --override-channels -c conda-forge -y
    fi
}

create_env_from_lock clamp-analyses envs/clamp-analyses.lock
create_env_from_lock gpu-kmeans envs/gpu-kmeans.lock
create_env_from_lock snakemake envs/snakemake.lock
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

# Register kernels inside their environments, without changing user-level kernels.
conda run -n clamp-analyses python -m ipykernel install --sys-prefix --name python3
conda run -n clamp-analyses Rscript -e 'IRkernel::installspec(user = FALSE, prefix = Sys.getenv("CONDA_PREFIX"))'
conda run -n gpu-kmeans python -m ipykernel install --sys-prefix --name python3

conda run -n clamp-analyses Rscript scripts/install_analysis_packages.R

log "WARNING: using the repository-pinned CLAMP revision; do not update CLAMP independently"
conda run -n clamp-analyses Rscript scripts/install_clamp.R

log "setup complete. Inputs are prepared by Snakemake; see README.md for data access and execution commands. Next steps:"
cat <<'EOF'

  conda activate snakemake
  snakemake --snakefile workflow/Snakefile --use-conda --cores 4 --dry-run all

To verify installed tools without starting analyses:
  python scripts/check_setup.py

Then follow README.md for inputs, hardware budgets and execution profiles.
EOF

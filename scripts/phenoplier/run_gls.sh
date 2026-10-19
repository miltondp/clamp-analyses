#!/usr/bin/env bash
# Run and validate trait associations for one model.
set -euo pipefail

usage() {
  cat >&2 <<USAGE
usage: $0 --rds <model.rds> --name <project> --summary-out <out.tsv.gz> [options]
  --conda-env NAME            (phenoplier-cli-neo)
  --require-version X.Y.Z     fail unless the installed phenoplier-cli is this release
  --namespace NS              workspace model namespace (clamp)
  --cohort NAME               (phenomexcan_rapid_gwas)
  --lv-percentile FLOAT       (0.01)
  --trait-filter all|biomedical  (biomedical)
  --expected-phenotypes N     completeness gate (2366 = biomedical)
  --executor local|slurm      phenoplier-cli's own executor (local)
  --cluster NAME              phenoplier-cli cluster profile, slurm only ('')
  --n-jobs N                  cores for the local executor (4)
  --workers-step6 N --blas-step6 N --workers-step7 N --blas-step7 N  (1 8 6 2)
  --workspace PATH            phenoplier workspace (default \$PHENOPLIER_HOME, then \$HOME/phenoplier)
USAGE
  exit 2
}

conda_env=phenoplier-cli-neo
require_version=""
namespace=clamp
cohort=phenomexcan_rapid_gwas
lv_percentile=0.01
trait_filter=biomedical
expected_phenotypes=2366
executor=local
cluster=""
n_jobs=4
workers_step6=1
blas_step6=8
workers_step7=6
blas_step7=2
rds=""
name=""
summary_out=""
workspace_arg=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rds) rds="$2"; shift 2 ;;
    --name) name="$2"; shift 2 ;;
    --summary-out) summary_out="$2"; shift 2 ;;
    --workspace) workspace_arg="$2"; shift 2 ;;
    --conda-env) conda_env="$2"; shift 2 ;;
    --require-version) require_version="$2"; shift 2 ;;
    --namespace) namespace="$2"; shift 2 ;;
    --cohort) cohort="$2"; shift 2 ;;
    --lv-percentile) lv_percentile="$2"; shift 2 ;;
    --trait-filter) trait_filter="$2"; shift 2 ;;
    --expected-phenotypes) expected_phenotypes="$2"; shift 2 ;;
    --executor) executor="$2"; shift 2 ;;
    --cluster) cluster="$2"; shift 2 ;;
    --n-jobs) n_jobs="$2"; shift 2 ;;
    --workers-step6) workers_step6="$2"; shift 2 ;;
    --blas-step6) blas_step6="$2"; shift 2 ;;
    --workers-step7) workers_step7="$2"; shift 2 ;;
    --blas-step7) blas_step7="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done
[[ -n "$rds" && -n "$name" && -n "$summary_out" ]] || usage
[[ -f "$rds" ]] || { echo "model not found: $rds" >&2; exit 1; }

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace="${workspace_arg:-${PHENOPLIER_HOME:-$HOME/phenoplier}}"
export PHENOPLIER_HOME="$workspace" PHENOPLIER_ROOT_DIR="$workspace"
if [[ ! -f "$workspace/config.toml" ]]; then
  echo "[ERROR] no phenoplier workspace at $workspace (missing config.toml):" \
    "run 'phenoplier workspace init' + 'workspace link' there, or set phenoplier.yaml: workspace" >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"
export PATH="${CONDA_PREFIX:-}/bin:${PATH}"

project="$workspace/projects/$name"
cfg="$project/pipeline_config.yaml"
results="$project/results/gls/phenoplier"
summary="$results/gls-summary-${cohort}.tsv.gz"
fingerprint_file="$project/clamp_source.sha256"

installed="$(phenoplier --version 2>&1 | tail -1 | sed -E 's/.*v([0-9][0-9A-Za-z.+-]*).*/\1/')"
echo "[start] $(date -Is)  name=$name  rds=$rds  executor=$executor  n_jobs=$n_jobs"
echo "[version] phenoplier-cli $installed  workspace=$workspace"
if [[ -n "$require_version" && "$installed" != "$require_version" ]]; then
  echo "[ERROR] phenoplier-cli $installed in env '$conda_env' but $require_version is required" \
    "(workflow/config/phenoplier.yaml: version); rebuild the env with scripts/phenoplier/setup_env.sh" >&2
  exit 1
fi

fingerprint="$(sha256sum "$rds" | cut -d' ' -f1)"
if [[ -d "$project" ]]; then
  if [[ -s "$fingerprint_file" ]]; then
    recorded="$(cut -d' ' -f1 < "$fingerprint_file")"
    if [[ "$recorded" != "$fingerprint" ]]; then
      stale="${project}.stale.$(date +%Y%m%dT%H%M%S)"
      echo "[stale] $(date -Is)  model changed since this project was initialised:"
      echo "[stale]   recorded $recorded"
      echo "[stale]   current  $fingerprint"
      echo "[stale]   moving $project -> $stale and re-registering $namespace/$name"
      exec 9>"$workspace/.clamp_registry.lock"; flock 9
      phenoplier model remove --force "$namespace/$name" || true
      flock -u 9; exec 9>&-
      mv "$project" "$stale"
    fi
  elif [[ -s "$cfg" ]]; then
    echo "[warn] $name: project predates fingerprinting; reusing it as-is." \
      "Delete $project to recompute from the current model." >&2
  fi
fi

if [[ ! -s "$summary" ]]; then
  if [[ ! -s "$cfg" ]]; then
    cluster_args=()
    if [[ -n "$cluster" ]]; then
      cluster_args=(--cluster "$cluster")
    fi
    lock="$workspace/.clamp_registry.lock"
    exec 9>"$lock"
    echo "[init] $(date -Is)  waiting for registry lock $lock"
    flock 9
    phenoplier model remove --force "$namespace/$name" 2>/dev/null || true
    echo "[init] $(date -Is)  registering model + writing $cfg"
    phenoplier shortcut gls \
      --input "$rds" \
      --name "$name" \
      --model-namespace "$namespace" \
      --cohort "$cohort" \
      --lv-percentile "$lv_percentile" \
      --trait-filter "$trait_filter" \
      --executor "$executor" \
      --n-jobs "$n_jobs" \
      ${cluster_args[@]+"${cluster_args[@]}"} \
      --unlock
    printf '%s  %s\n' "$fingerprint" "$rds" > "$fingerprint_file"
    flock -u 9
    exec 9>&-
  else
    echo "[init] $(date -Is)  resuming: $cfg exists"
  fi

  python - "$cfg" "$workers_step6" "$blas_step6" "$workers_step7" "$blas_step7" <<'PY'
import sys, yaml
path, nw6, bl6, nw7, bl7 = sys.argv[1], *map(int, sys.argv[2:6])
with open(path) as fh:
    cfg = yaml.safe_load(fh)
cfg["n_workers_step6"] = nw6
cfg["blas_threads_step6"] = bl6
cfg["n_workers_step7"] = nw7
cfg["blas_threads_step7"] = bl7
with open(path, "w") as fh:
    yaml.safe_dump(cfg, fh, sort_keys=False)
print(f"[pin] step6={nw6}x{bl6} step7={nw7}x{bl7} -> {path}")
PY

  rm -f "$project/.snakemake/locks/"*.lock 2>/dev/null || true
  echo "[run] $(date -Is)"
  phenoplier workflow gls run --config-file "$cfg"
else
  echo "[run] skip: $summary exists"
fi

if [[ ! -s "$summary" ]]; then
  summary="$results/gls-summary-phenomexcan.tsv.gz"
fi
if [[ ! -s "$summary" ]]; then
  summary="$(ls -1 "$results"/gls-summary-*.tsv.gz 2>/dev/null | head -n1 || true)"
fi
if [[ -z "$summary" || ! -s "$summary" ]]; then
  echo "[ERROR] $name: no gls-summary-*.tsv.gz under $results (run failed?)" >&2
  exit 1
fi

if ! python "$here/verify_summary.py" "$summary" "$expected_phenotypes"; then
  echo "[ERROR] $name: summary failed the completeness gate; removing it so a requeue resumes step 7" >&2
  rm -f "$summary" "${summary%.tsv.gz}.pkl.gz"
  exit 1
fi

mkdir -p "$(dirname "$summary_out")"
tmp="${summary_out}.tmp.$$"
cp -f "$summary" "$tmp"
mv -f "$tmp" "$summary_out"
stem="${summary_out%.tsv.gz}"
if [[ -f "${summary%.tsv.gz}.pkl.gz" ]]; then
  cp -f "${summary%.tsv.gz}.pkl.gz" "${stem}.pkl.gz"
fi
if [[ -f "$project/trait_filter_excluded.tsv" ]]; then
  cp -f "$project/trait_filter_excluded.tsv" "${stem}_trait_filter_excluded.tsv"
fi
echo "[done] $(date -Is)  $name -> $summary_out"

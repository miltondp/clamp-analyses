#!/usr/bin/env bash
# Run CLAMPbase trait associations and publish the study store.
set -euo pipefail
name="$1"; rds="$2"

WS="${PHENOPLIER_HOME:-$HOME/phenoplier}"
COHORT=phenomexcan_rapid_gwas
CENTRAL="${CLAMPBASE_CENTRAL:?set CLAMPBASE_CENTRAL to the study collection dir}"
NJOBS="${CLAMPBASE_NJOBS:-${SLURM_CPUS_PER_TASK:-12}}"
NW_S6="${CLAMPBASE_NW_STEP6:-1}"; BL_S6="${CLAMPBASE_BLAS_STEP6:-8}"
NW_S7="${CLAMPBASE_NW_STEP7:-6}"; BL_S7="${CLAMPBASE_BLAS_STEP7:-2}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PHENOPLIER_HOME="$WS" PHENOPLIER_ROOT_DIR="$WS"

proj="$WS/projects/$name"
cfg="$proj/pipeline_config.yaml"
summary="$proj/results/gls/phenoplier/gls-summary-${COHORT}.tsv.gz"
store="$CENTRAL/stores/${name}.h5"
mkdir -p "$CENTRAL/stores" "$CENTRAL/summaries"

echo "[start] $(date -Is)  $name  njobs=$NJOBS  mem=${SLURM_MEM_PER_NODE:-?}MB  ver=$(phenoplier -v 2>&1 | tail -1)"

if [ ! -s "$summary" ]; then
  if [ ! -s "$cfg" ]; then
    phenoplier shortcut gls --input "$rds" --name "$name" --model-namespace "${CLAMPBASE_NS:-clampbase}" \
      --trait-filter biomedical --lv-percentile 0.01 --executor local --n-jobs "$NJOBS" --unlock
  fi

  python - "$cfg" "$NW_S6" "$BL_S6" "$NW_S7" "$BL_S7" <<'PY'
import sys, yaml
p, nw6, bl6, nw7, bl7 = sys.argv[1], *map(int, sys.argv[2:6])
c = yaml.safe_load(open(p))
c["n_workers_step6"] = nw6; c["blas_threads_step6"] = bl6
c["n_workers_step7"] = nw7; c["blas_threads_step7"] = bl7
yaml.safe_dump(c, open(p, "w"), sort_keys=False)
print(f"patched step6={nw6}x{bl6} step7={nw7}x{bl7} ->", p)
PY

  rm -f "$proj/.snakemake/locks/"*.lock 2>/dev/null || true
  phenoplier workflow gls run --config-file "$cfg"
else
  echo "[gls] skip (summary exists)"
fi

if [ ! -s "$summary" ]; then
  echo "[ERROR] $name: GLS produced no combined summary (run failed?). Refusing to build a store." >&2
  exit 1
fi
if ! python "$HERE/../verify_summary.py" "$summary"; then
  echo "[ERROR] $name: summary incomplete/invalid. Refusing to build a store; requeue to finish step7." >&2
  exit 1
fi

if [ ! -s "$store" ]; then
  echo "[store] $(date -Is)"
  tmp="${store}.tmp.$$"
  rm -f "$tmp"
  if ! phenoplier store build --output "$tmp" --cohort "$COHORT" --clamp-rds "$rds" \
        --gls-dir "$proj/results/gls/phenoplier/${COHORT}" --method-name CLAMP; then
    echo "[ERROR] $name: store build failed." >&2; rm -f "$tmp"; exit 1
  fi
  if ! python -c "import h5py,sys; h5py.File(sys.argv[1],'r').close()" "$tmp"; then
    echo "[ERROR] $name: built store failed to open." >&2; rm -f "$tmp"; exit 1
  fi
  mv -f "$tmp" "$store"
else
  echo "[store] skip (exists)"
fi

if ! cp -f "$summary" "$CENTRAL/summaries/${name}.tsv.gz"; then
  echo "[ERROR] $name: could not publish summary to $CENTRAL/summaries" >&2; exit 1
fi
echo "[done] $(date -Is)  $name"

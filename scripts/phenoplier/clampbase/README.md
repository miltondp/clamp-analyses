# CLAMPbase — GLS run (pico)

Runs `phenoplier shortcut gls` for the **CLAMPbase** model family — the same
finals / saturation / coverage design wired in `workflow/rules/phenoplier.smk`.
**41 models:** finals 3 (archs4/gtex/recount2) + saturation-k1728 17 + coverage 21.
Every run uses `--trait-filter biomedical` (2,366 of 4,049 traits).

## The hardened per-model runner

All three sbatch scripts call `_run_clampbase_model.sh <name> <rds>`, which runs
one model end-to-end (register → GLS → store) with the fixes this pipeline needs:

- **Serial step 6** (`n_workers_step6=1`, in-process) — the forked
  `ProcessPoolExecutor` used to deadlock (workers parked in `futex`, 0 % CPU).
- **Pinned step 7** (`n_workers_step7 × blas_threads_step7`; default 6×2 for a
  SLURM cgroup) — an un-pinned run over-detects the node's cores and spawns ~29
  threads/job, thrashing the box ~10×. No-cgroup hosts export
  `CLAMPBASE_NW_STEP7=12 CLAMPBASE_BLAS_STEP7=1 OMP_NUM_THREADS=1`.
- **Register once, then resume.** `shortcut gls --unlock` registers the model +
  writes the config *without* computing; on requeue it does **not** re-init, so
  the expensive step1–6 sentinels survive and the run resumes.
- **Two gates before a store is trusted:** (1) the combined summary must exist and
  be *complete* — `verify_summary.py` checks 2,366 unique phenotypes, no
  zero-p-values, NaN only where `lv_degenerate=True`; `summarize` will happily
  concatenate a *partial* step-7, so this catches an interrupted run. (2) the
  store is built to a `.tmp`, verified to open, then atomically renamed — a
  failed build never leaves a partial `.h5` that a later `-s` check would reuse.

## Reproduce

Use the Snakemake targets from the repository root:

```bash
snakemake --snakefile workflow/Snakefile --use-conda --cores 12 phenoplier_coverage phenoplier_saturation phenoplier_finals
```

The existing `.sbatch` files delegate to these targets and require the named
`snakemake` environment. They use the same configured paths and producers as
the main workflow. Configure accounts and partitions through your site profile.

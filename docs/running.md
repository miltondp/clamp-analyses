# Hardware, environments and execution profiles

## Hardware

Memory values below are configured job requests, **not measured peak RSS**.
The full pipeline uses a single node per job; it does not distribute one matrix
fit across nodes. Disk-backed matrices still need substantial RAM and fast disk.

| Stage | Request (decimal MB) |
| --- | ---: |
| ARCHS4 preprocessing | 200000 |
| Coverage/full-data SVD and CLAMPbase | 320000 in config; SLURM coverage fit override 500000 |
| Canonical ARCHS4 CLAMPfull | 320000 |
| Stability CLAMPfull (also used by CRISPR) | 400000 |
| Saturation | Up to 264000 initially, 528000 on the third attempt |
| PhenoPLIER GLS | 100000 |
| GTEx timing (Flashier) | 110000 |
| GTEx GPU clustering (host RAM) | 128000 |

Coverage declares a retry policy capped at 800000 MB; saturation is capped at
675000 MB, with its current largest configured job reaching 528000 MB after two
retries. Not every fitting rule uses the retry multiplier. The supplied SLURM
profiles expose 800000 MB as the **aggregate scheduling pool**, not as each
job's RAM request. Recommend a roughly 1 TB node for the complete workflow with
headroom and retry capacity; the existing settings do not prove a lower physical
minimum. Increasing `--cores` does not reduce per-job memory.

The default local profile has only 96000 MB total and 64000 MB fit/ORA pools;
it is for smaller selected jobs. Snakemake can cap a request at a global resource
limit, so lowering those pools is not a safe way to make large fits fit locally.
Use budgets that actually satisfy the jobs you select.

A compatible NVIDIA GPU is required for GTEx k-means and projection clustering.
For SLURM, configure the GPU rules for your site's GPU partition and GRES syntax:
`kmeans_clustering_gtex`, `projection_benchmark_cluster`, and
`projection_benchmark_null_cluster` (confirm names with `--list-rules`).
Use your site's plugin-supported `slurm_extra`/GPU resource settings; the pinned plugin maps `gpu=1` to `--gpus=1`. Route to a GPU partition
on clusters that separate CPU and GPU nodes.

## Environments

`./setup.sh` creates three named environments. Explicit `.lock` files pin Linux
Conda package builds; YAML files are the editable specifications. Pip requirements
come from `envs/clamp-analyses.yaml`; these are version pins, not a transitive
hash lock. Non-Conda comparison methods have Git commit pins in
`scripts/install_analysis_packages.R`. CLAMP has its separate commit pin in
`scripts/install_clamp.R`. Do not update a single method independently.

For manual installation, with Conda available:

```bash
conda create -n clamp-analyses --file envs/clamp-analyses.lock
conda create -n gpu-kmeans --file envs/gpu-kmeans.lock
conda create -n snakemake --file envs/snakemake.lock
./setup.sh  # reuses these environments; installs pip, comparison methods and CLAMP
conda activate snakemake
python scripts/check_setup.py
```

Existing environments are reused, not reconciled to their locks. For an exact
rebuild, use a separate Conda installation/prefix or explicitly recreate your own
environments; setup never deletes them. Install Inkscape separately for Figure 4.
PhenoPLIER uses another environment and its reference bundle; see its linked guide.
Anonymous source access returned HTTP 404 during this audit. Obtain repository
access from the authors; for SSH authentication, set
`PHENOPLIER_REPO_URL=git@github.com:pivlab/phenoplier-cli.git` when running its setup.
To update the runner lock, solve `envs/snakemake.yaml` in a disposable environment,
verify the runner/plugin, then export with `conda list --explicit --md5`.

## Pico / server_cu

1. Log in to Pico using your own CU access details. Keep the clone, data, outputs,
   and Conda environments on storage visible to compute nodes, not login-node `/tmp`.
2. Run setup on a node where software installation is permitted. Activate
   `snakemake`; confirm `sbatch`, `squeue`, `sinfo`, and Conda are available.
3. `workflow/profiles/server_cu` uses `defq`, as recorded in the existing Pico
   site configuration. Check `sinfo -o '%P %m %c %l %G'` for current memory,
   wall-time and GPU limits and your eligibility. The pinned SLURM plugin maps `gpu=1` to `--gpus=1`; all GPU rules
   now declare this resource.
   Verified read-only with `sinfo` on 2026-10-19: `defq` reports 900000 MB,
   168 CPUs, infinite time and `gpu:1` (NVIDIA L40S, approximately 46 GB VRAM).
   No jobs were submitted.
4. Stage the [inputs](inputs.md). For PhenoPLIER, the recorded bundle is
   `/pividori_lab/data/phenoplier_full_data.tar.gz` and workspace
   `/pividori_lab/phenoplier_workspace`. Confirm availability/permissions. Use
   your own writable workspace/copy if sharing would modify someone else's data.
   Set `PHENOPLIER_HOME` or a private config override with `phenoplier.workspace`.
5. Start `tmux`, then run from the repository root:

```bash
conda activate snakemake
snakemake --profile workflow/profiles/server_cu --dry-run all
snakemake --profile workflow/profiles/server_cu all
# Substitute benchmark or archs4_biology to run one phase.
squeue -u "$USER"
```

For personal account/partition/GPU overrides, copy this profile to a user-owned
location and edit it. Avoid committing personal accounts or absolute paths.
The top-level Pico profile handles SLURM submission. Keep
`phenoplier.target: local` so PhenoPLIER runs inside each allocated job;
`phenoplier_target=server_cu` instead enables its separate nested scheduler.

## Another SLURM server

```bash
mkdir -p ~/.config/snakemake/clamp-slurm
cp workflow/profiles/slurm/config.yaml ~/.config/snakemake/clamp-slurm/config.yaml
```

Edit that copy: add `slurm_partition` and, where required, `slurm_account` and
`slurm_qos` to `default-resources`; set GPU partition/GRES in `set-resources` for
GPU rules. Check each partition's per-node RAM, CPU and time limits. Some fits
request weeks; route them to a suitable partition rather than shortening requests
and silently accepting incomplete fits. Set `jobs`, `cores`, `mem_mb`,
`fit_mem_mb`, and slot pools to the capacity available to you. Keep `tasks=1`,
`use-conda: true`, and rule-specific high-memory overrides. Environment names
must resolve on compute nodes. The profile selects the Snakefile relative to
the repository root, so launch there.

```bash
conda activate snakemake
snakemake --profile ~/.config/snakemake/clamp-slurm --dry-run all
snakemake --profile ~/.config/snakemake/clamp-slurm all
```

## Local machine

On a 192 GB machine with a compatible GPU, a benchmark example is:

```bash
snakemake --profile workflow/profiles/local --cores 16 \
  --resources mem_mb=160000 fit_mem_mb=160000 ora_mem_mb=64000 timing_slot=1 \
  --dry-run benchmark
```

On a roughly 1 TB machine with a compatible GPU, a full-workflow example is:

```bash
snakemake --profile workflow/profiles/local --cores 32 \
  --resources mem_mb=800000 fit_mem_mb=800000 ora_mem_mb=64000 large_fit_slot=1 timing_slot=1 \
  --dry-run all
```

After inspecting the DAG, remove `--dry-run` to execute. Adapt budgets to actual
available hardware and avoid concurrent timing fits when comparing runtimes.

# CLAMP analyses

## 1. Overview

This repository contains the code supporting **“A genetically informed atlas of
transcriptomic programs from massive human RNA-seq compendia.”**
[Manuscript](https://pivlab.github.io/clamp-manuscript/).

Pseudobulk and GTEx benchmarks evaluate recovery of cell types, tissues and
pathways. ARCHS4 analyses examine biological programs, drug responses and CRISPR
screens. Consensus models combine independent fits across ARCHS4, GTEx and recount2.
The workflow also produces manuscript figures and supplementary tables.

The analyses are organized as a **Snakemake workflow**. Each rule defines a step's
inputs, command and outputs. Selecting a target tells Snakemake to prepare its
missing prerequisites in order and reuse completed outputs. The main entrypoint
is [`workflow/Snakefile`](workflow/Snakefile); profiles configure local or server
execution. A `--dry-run` previews the required steps without executing them.

## 2. Repository structure

| Location | Contents |
| --- | --- |
| `setup.sh`, `envs/` | Environment installation, specifications and locks |
| `workflow/` | Main `Snakefile`, analysis rules, configuration and execution profiles |
| `scripts/` | Data preparation, model fitting and analysis code |
| `nbs/` | Analysis and figure notebooks |
| `references/` | Dataset manifests, fixed cohorts and recorded timing results |
| `data/` | Downloaded inputs and bundled reference datasets |
| `output/` | Generated models, analyses and logs |
| `figures/` | Publication figures and tables |

Commands below assume the repository root.

## 3. Installation and requirements

The locks support **Linux x86-64**. Setup uses Conda to install analysis Python
3.11/R 4.4, Snakemake 8.30.0, the SLURM executor and a separate GPU environment
(Python 3.13, CUDA 12.9). Git, curl and internet access are required.

```bash
git clone https://github.com/pivlab/clamp-analyses.git
cd clamp-analyses
./setup.sh
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
python scripts/check_setup.py
```

With existing Conda, use its installation path in the `source` command. Setup
installs the pinned CLAMP revision and notebook kernels.

## 4. Data availability

Inputs download and process automatically when their analysis is requested,
except for the HeartMap download below.

| Dataset | Source and preparation |
| --- | --- |
| Brain: Mathys, Xiong | [UCSC Aging Brain](https://cells.ucsc.edu/ad-aging-brain/); pseudobulk matrices from the [versioned archive](https://doi.org/10.5281/zenodo.16649011) |
| Blood: OneK1K, Perez | CELLxGENE versions `1e44db10-b572-46cc-adae-dcc7acd44ca6`, `c55dc602-d168-4d15-acc1-5de4f2f5d551`; archived pseudobulk matrices |
| Lung: HLCA | CELLxGENE version `688185ad-11c2-4172-a53a-f4f1f4076860`; donor aggregation |
| Heart: HeartMap | [SCP3689](https://singlecell.broadinstitute.org/single_cell/study/SCP3689), [study DOI](https://doi.org/10.1038/s44161-026-00831-5); donor aggregation |
| GTEx | [GTEx v8](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression); expression and annotations |
| ARCHS4 | [Human gene expression v2.5](https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.5.h5); filtering and normalization |
| recount2 | [Archived expression and metadata](https://ndownloader.figshare.com/files/10881866); extraction and preprocessing |
| Projection datasets | GSE133218, GSE193336 and E-MTAB-6701; acquisition and processing in `workflow/config/projections.yaml` |

For HeartMap, open [study SCP3689 in the Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP3689).
Create an account if needed and sign in. In the study's **Download** tab, download
`HeartMap_V1.0_raw_counts.h5ad` and place it at the following path relative to the
repository root (create the directories if needed):

```text
data/pseudobulk/Heart_Datar2026/HeartMap_V1.0_raw_counts.h5ad
```

Pathway gene sets, cell-type markers and CRISPR screen data are included.
Source URLs, versions and checksums are recorded in `workflow/config/*inputs.json`.

## 5. Quick start

This synthetic example needs no external datasets or GPU. It processes four
cells using the production pseudobulk code and executes a notebook that checks
expression values and cell-type fractions.

```bash
snakemake --profile workflow/profiles/local --cores 1 quickstart
```

Expected: `output/00_quickstart/report.json` contains `"status": "PASS"`, alongside
expression/fraction CSVs and `quickstart.ipynb`. The validation run took about
17 seconds including scheduler/environment startup; processing took 1.1 seconds.
It requests one CPU and 2 GB RAM.

## 6. Reproducing manuscript results

| Target | What it requests |
| --- | --- |
| `benchmark` | Pseudobulk and GTEx model benchmarks and biological evaluations |
| `archs4_biology` | ARCHS4 fits, biological analyses, cross-compendium comparisons and consensus models |
| `publish_figures` | Manuscript figures and supplementary tables, exported to `figures/` |
| `all` | Both analysis phases and publication outputs |

Snakemake reuses existing outputs and builds missing prerequisites. In particular,
`publish_figures` renders publication outputs from analysis results, but missing
results can trigger upstream analyses. Its dry-run shows whether any fits would run.

### Phase 1: benchmarks on the host

The `local` profile runs on the current machine. Its defaults are defined in
[`workflow/profiles/local/config.yaml`](workflow/profiles/local/config.yaml).
The CPU and RAM budgets can be set directly in the command:

```bash
conda activate snakemake
CPUS=8
RAM_MB=128000  # example scheduling budget in MB; adjust to the host

snakemake --profile workflow/profiles/local --cores "$CPUS" \
  --resources mem_mb="$RAM_MB" gpu=1 --dry-run benchmark
snakemake --profile workflow/profiles/local --cores "$CPUS" \
  --resources mem_mb="$RAM_MB" gpu=1 benchmark
```

`--cores` limits CPU scheduling and `--resources mem_mb` sets the RAM budget
Snakemake schedules against across concurrent local jobs. Choose a budget below
the machine's physical RAM, leaving headroom for the operating system. These are
scheduling limits, not operating-system memory caps.

Each rule also declares its own RAM request. Those defaults are in
`workflow/rules/` and `workflow/config/`. To change one request, add, for example,
`--set-resources kmeans_clustering_gtex:mem_mb=64000` to both commands. Increasing
or decreasing this value changes the declared allocation, not the computation's
actual memory needs; an allocation that is too small can cause failure. The GTEx
clustering rule defaults to 128,000 MB and requires a GPU; `gpu=1` allows one GPU
job at a time. CPU-only options are in Section 7.

The targets below select smaller parts of the workflow. They are aggregate rules
in [`workflow/Snakefile`](workflow/Snakefile); their steps are defined in
[`pseudobulk.smk`](workflow/rules/pseudobulk.smk) and
[`gtex.smk`](workflow/rules/gtex.smk). The same CPU/RAM budgets apply:

```bash
# Pseudobulk analyses, donor-bulk evaluations and associated panels
snakemake --profile workflow/profiles/local --cores "$CPUS" \
  --resources mem_mb="$RAM_MB" gpu=1 pseudobulk_complete

# GTEx model fitting and biological evaluations
snakemake --profile workflow/profiles/local --cores "$CPUS" \
  --resources mem_mb="$RAM_MB" gpu=1 gtex_complete

# Pseudobulk input preparation only
snakemake --profile workflow/profiles/local --cores "$CPUS" \
  --resources mem_mb="$RAM_MB" gpu=1 pseudobulk_inputs
```

Shared dependencies are still included: for example, some pseudobulk publication
panels also read GTEx results. Adding `--dry-run` previews the jobs for any target.

If the local machine has enough memory for the full workflow (roughly 1 TB RAM
for the largest configured requests), a compatible NVIDIA GPU and sufficient
storage, both phases and publication outputs can run locally:

```bash
snakemake --profile workflow/profiles/local --dry-run all
snakemake --profile workflow/profiles/local all
```

### Phase 2: ARCHS4 biology on a server

Snakemake runs on the server, where it schedules jobs and reuses completed
results. The server session can be opened from a workstation through SSH;
profiles configure execution and resources, not the SSH connection.

**Server preparation.** The server needs its own checkout of this repository,
the environments from Section 3 and access to the inputs in Section 4. On SLURM,
the project and environments must be accessible to the compute nodes. Available
inputs download as needed; existing model/results files can also be reused.

**SLURM server.** Connect using your server account and hostname, enter the
checkout and activate its environment. Replace the placeholders and adjust the
Conda path if needed:

```bash
ssh USER@HOST
cd /SHARED/PROJECT/clamp-analyses
source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
```

The SLURM profile is **`workflow/profiles/slurm/config.yaml`**. Adapt its settings
to the server; the administrator supplies partition/account names and node limits.

```bash
nano workflow/profiles/slurm/config.yaml
```

| Setting in the profile | What to configure |
| --- | --- |
| `jobs`, `cores` | Submission concurrency and CPU scheduling budget |
| `resources` | RAM/resource budgets in MB and fit/GPU concurrency limits |
| `default-resources` | Add `slurm_partition=YOUR_PARTITION`; add `slurm_account=YOUR_ACCOUNT` only when required |
| `default-resources` entry `runtime` | Default wall time in minutes |
| `set-resources`, `set-threads` | Per-rule RAM, wall time, partition and CPU requests; these override defaults |
| Top-level `slurm-qos` | Add `slurm-qos: YOUR_QOS` only if required |

Keep per-job RAM requests within one node's capacity, and allow enough wall time
for each fit. For separate GPU partitions, add this entry under the **existing**
`set-resources` section:

```yaml
  kmeans_clustering_gtex:
    slurm_partition: YOUR_GPU_PARTITION
    gpu: 1
```

### Publication outputs or complete reproduction

Publication outputs can use the local profile:

```bash
# Publication figures and supplementary tables
snakemake --profile workflow/profiles/local --dry-run publish_figures
snakemake --profile workflow/profiles/local publish_figures

# Both analysis phases and publication outputs
snakemake --profile workflow/profiles/local --dry-run all
snakemake --profile workflow/profiles/local all

# One figure or the supplementary tables
snakemake --profile workflow/profiles/local output/99_panels/fig2/fig2.pdf
snakemake --profile workflow/profiles/local supplementary_tables
```

Figure rules and their input dependencies are in
[`workflow/rules/panels.smk`](workflow/rules/panels.smk).
`consensus_models` requests consensus models, ORA and LINCS only.

### Rerunning one notebook

Notebook rules reuse existing upstream results. For example, this reruns only the
GTEx clustering report notebook; missing prerequisites are prepared if needed:

```bash
snakemake --profile workflow/profiles/local --dry-run kmeans_clustering_report_gtex --forcerun kmeans_clustering_report_gtex
snakemake --profile workflow/profiles/local kmeans_clustering_report_gtex --forcerun kmeans_clustering_report_gtex
```

Source: `nbs/03_model_biology/01_gtex/00_kmeans_clustering.ipynb`.
The executed copy is `00_kmeans_clustering.executed.ipynb` in the same directory;
report tables are under `output/03_model_biology/01_gtex/00_kmeans_clustering/`.
Other notebooks have corresponding rules in `workflow/rules/`.

## 7. Computational resources

RAM requests are predefined in [`workflow/config/`](workflow/config/), the
`resources` sections of [`workflow/rules/`](workflow/rules/) and profile overrides
in [`workflow/profiles/`](workflow/profiles/). These settings are adjustable.

| Workflow part | RAM per job (configured) |
| --- | --- |
| Pseudobulk | Up to 64 GB |
| GTEx | Up to 96 GB |
| `archs4_biology` | Up to 528 GB |
| Figures | Default 8 GB |

### Without a GPU

The current GTEx and projection clustering code uses CuPy/cuML and has **no
automatic CPU fallback**. `quickstart` and CPU-only targets such as
`full_models_pseudobulk` and `full_models_gtex` can run without a GPU. Complete
`benchmark`/`all` runs require GPU access for those clustering steps, or their
already-generated outputs. A GPU server can perform those steps and its outputs
can be transferred back to the host. Setting `gpu=0` does not convert GPU code
to CPU execution.

## 8. Outputs

`output/01_model_building/` contains models; `02_model_performance/` benchmarks;
`03_model_biology/` biological analyses; `98_final_models/` final and consensus
bundles; and `99_panels/` figures/tables. Publication copies go to `figures/`.

Logs accompany outputs; scheduler logs are under `.snakemake/`. Repeating the
same command resumes work and reuses completed outputs.

## 9. Citation

Cite the [manuscript](https://pivlab.github.io/clamp-manuscript/) and the
[software repository](https://github.com/pivlab/clamp-analyses), including the
commit used. No archived software DOI is currently recorded here.

## 10. License and contact

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Bundled datasets retain
their source-specific license notices. Support: [GitHub issues](https://github.com/pivlab/clamp-analyses/issues).

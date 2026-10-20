# CLAMP analyses

Reproduce the CLAMP benchmarks, ARCHS4 biological analyses, and publication figures.
Run commands from the repository root. Large scientific inputs are downloaded by
Snakemake or supplied separately; cloning alone does not include them.

## 1. Check your machine

Locked environments support **Linux x86-64**. You need Git, curl, internet access
for setup/downloads, and a shared filesystem when using multiple compute nodes.
GTEx clustering and projection benchmarks require a compatible NVIDIA GPU/driver
(the `gpu-kmeans` environment uses CUDA). Figure 4 requires `inkscape` on PATH.

| Work | Configured RAM per job, not measured minimum |
| --- | --- |
| Benchmark (pseudobulk + GTEx) | Some jobs request 96–128 GB; plan for at least 192 GB RAM |
| ARCHS4 preprocessing | 200 GB |
| Full ARCHS4 models | 320–400 GB; supplied SLURM profiles request 500 GB for coverage fits |
| Retry allocations | Up to 528 GB in saturation; coverage retry policy is capped at 800 GB |

For the complete ARCHS4 workflow, use a **node/machine with approximately 1 TB
RAM** to accommodate the profiles' 800 GB scheduling budget and system overhead.
This is a conservative allocation recommendation, not a measured minimum.
A sufficiently large workstation works; otherwise use a high-memory server or
cluster. RAM must fit on **one node**, not be summed across nodes. Smaller
subsamples have smaller fits but still require preprocessing and shared models.
See [hardware and profiles](docs/running.md) for resource details.

Reserve multi-terabyte storage for a full rebuild: raw inputs, file-backed
matrices, replicate models, and scratch files coexist. Exact usage depends on
which analyses you run; check available space before starting large downloads.

## 2. Install

```bash
git clone https://github.com/pivlab/clamp-analyses.git
cd clamp-analyses
./setup.sh
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snakemake
python scripts/check_setup.py
```

If setup installed Conda and your shell cannot find it, first run
`source ~/miniconda3/etc/profile.d/conda.sh` (or your `MINICONDA_PREFIX`).
Setup creates `clamp-analyses`, `gpu-kmeans`, and `snakemake` from explicit Conda
locks, installs pinned pip/comparison packages, and installs the required CLAMP
revision. Existing environments are reused; verification reports mismatches.
Keep the pinned CLAMP revision. [Manual setup and locks](docs/running.md#environments).

## 3. Prepare inputs

Small reference inputs (MSigDB 2026.1, CellMarker, cell/tissue references and
ARCHS4 gene lengths) are included in Git. Download GTEx annotations and verify
the pinned references before starting:

```bash
snakemake --profile workflow/profiles/local reference_inputs
python scripts/inputs/reference_inputs.py --check
```

[Sources, DOIs and checksums](docs/reference_inputs.md) document each reference.
The [input checklist](docs/inputs.md) lists remaining cohort/reference inputs.
ARCHS4 biology additionally needs PhenoPLIER/PhenomeXcan and drug–disease/CRISPR
inputs. Some projection analyses build shared benchmark prerequisites.

For ARCHS4 trait analyses, also install and initialize **phenoplier-cli v0.5.2**
using the [environment/workspace instructions](scripts/phenoplier/README.md#3-phenoplier-cli-environment-and-workspace).
This separate environment is not installed by `setup.sh`. The source repository
currently requires access from the authors; arrange access before running phase 2.

## 4. Run

Choose a profile below. Always dry-run first; full fits and repeated timing runs
can take days. Use the same profile for the dry-run and execution.

```bash
# Complete workflow, including figures and supplementary tables
snakemake --profile workflow/profiles/server_cu --dry-run all
snakemake --profile workflow/profiles/server_cu all

# Phase 1: pseudobulk + GTEx benchmarks, donor-bulk, timing and figures
snakemake --profile workflow/profiles/server_cu benchmark

# Phase 2: ARCHS4 models and biology, comparisons, figures and tables
snakemake --profile workflow/profiles/server_cu archs4_biology
```

| Where | Profile/instructions |
| --- | --- |
| **Pico / server_cu (CU)** | `workflow/profiles/server_cu`; uses SLURM partition `defq`. Follow [Pico setup](docs/running.md#pico--server_cu) before submission. |
| Another SLURM server | Copy `workflow/profiles/slurm` to a user-owned profile; set partition/account, node memory and time limits. [Instructions](docs/running.md#another-slurm-server). |
| Local machine | `workflow/profiles/local` defaults to 24 cores / 96 GB for individual smaller targets. For a complete phase, set budgets matching your hardware; [examples](docs/running.md#local-machine). |

The profiles select `workflow/Snakefile` and activate each rule's named Conda
environment. `all` is the default target. Existing individual targets remain
available with `snakemake --profile <profile> --list-target-rules`.

## 5. Results and resuming

Models: `output/01_model_building/`; benchmarks: `output/02_model_performance/`;
biology: `output/03_model_biology/`; final models: `output/98_final_models/`;
figures/tables: `output/99_panels/`. `all` exports publication copies to `figures/`.
Executed notebook copies live alongside source notebooks; job logs live beside
analysis outputs and scheduler logs under `.snakemake/`.

Resume by rerunning the same command: completed outputs are reused and profiles
retry incomplete jobs. Keep cluster schedulers running in `tmux`/`screen`.
Inspect errors before forcing jobs; never force an existing final-model bundle.
For one notebook rule, use `snakemake --profile <profile> <rule>`; force only that
rule with `--forcerun <rule>` when intentionally regenerating its outputs.

[Analysis/provenance reference](docs/analysis-reference.md) ·
[PhenoPLIER details](scripts/phenoplier/README.md) ·
[Validation and limitations](docs/validation.md)

## License and acknowledgments

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Supported by NHGRI, NICHD, NSF, and NEI.

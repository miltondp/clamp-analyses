# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Analysis code (not a package) that consumes the external `CLAMP` R package (https://github.com/chikinalab/CLAMP). Everything here is notebooks that build/evaluate CLAMP latent-variable models on bulk RNA-seq datasets (ARCHS4, GTEx, recount2) and compare to baselines (PLIER, NMF/PCA, MOFA, CoGAPS, flashier, GenomicSuperSignature).

The `CLAMP` package itself lives in a separate repo and must be cloned + `devtools::install_local()`'d into the active conda env before notebooks will run (see README for the exact incantation).

## Environments

Two conda envs, intentionally separated to avoid R/Bioconductor vs RAPIDS conflicts:

- `envs/clamp-analyses.yaml` — default. R 4.4 + Bioconductor + Python data stack. Use for everything unless GPU is needed.
- `envs/gpu-kmeans.yaml` — RAPIDS/cuML/cupy. Only for GPU-accelerated clustering/benchmark notebooks.

Each notebook's first markdown cell states which env it requires. Honor that — switching envs mid-pipeline is expected.

Create with `conda env create -f envs/<file>.yaml`, then install CLAMP into the activated env via `Rscript -e "devtools::install_local('<path-to-CLAMP-repo>', force=TRUE, dependencies=FALSE)"`.

## Running notebooks

Interactive: `jupyter lab` from the project root (the `nb_conda_kernels` package exposes both envs as kernels).

Headless / batch (SLURM): scripts in `jobs/*.batch` execute notebooks via `jupyter nbconvert --to notebook --execute <path> --output <name>_executed.ipynb`. The batch files are written for the lab's SLURM cluster (`/pividori_lab/...` paths, plierv2 conda env) and will need path edits to run elsewhere. They set `NUMBA/MKL/OPENBLAS/NUMEXPR/OMP_NUM_THREADS=30` to cap thread fan-out — preserve this pattern when adding new batch jobs.

## Config and paths

`config.R` is the single source of truth for dataset URLs, file paths, and per-dataset CLAMP hyperparameters (`MAX_ITER`, `MAX_U_UPDATES`, `RANDOM_SVD_SEED`, etc.). Every R notebook starts with `source(here("config.R"))` and then references `config$ARCHS4$...`, `config$GTEx$...`, `config$recount2$...`.

Paths use `here::here()` rooted at the project (anchored by `clamp.Rproj`). `data/` (inputs) and `output/` (model artifacts) are both gitignored — they're large and reproducible from the notebooks. When adding a new dataset, add a `config$<NAME>` block following the existing shape rather than hardcoding paths in the notebook.

## Notebook layout

The `nbs/` tree is a pipeline staged by number — directories run in order, files within run in order:

- `00_setup/` — env bootstrapping
- `01_model_building/{archs4,gtex,recount2}/` — download → preprocess → SVD → CLAMPbase → CLAMPfull → multi-model variants (`_M` suffix = with M matrix; `multimodel` = ensemble of MULTICLAMP runs)
- `02_model_comparisons/` — Rmd comparisons across hyperparameters/methods
- `03_model_performance/` — held-out evaluation
- `03_model_biology/` and `04_model_biology/` — biological interpretation. Both exist side-by-side: `04_model_biology/{archs4,gtex}/` holds the older projection/clustering/CRISPR work; `03_model_biology/00_archs4/` holds the newer drug-disease association pipeline (which has its own scoped `CLAUDE.md`). The numbering collision is historical — neither subtree supersedes the other; check both when looking for biology-stage analyses.

`old/` subdirs contain superseded `.Rmd`/`.html` versions — read for context, don't extend.

## Conventions

- R: 2-space indent (per `clamp.Rproj`).
- Set `set.seed(123)` after library loads; CLAMP's internal seeds come from `config$<DATASET>$...$RANDOM_SVD_SEED`.
- For new analyses on an existing dataset, add a numbered notebook in the appropriate `nbs/<stage>/<dataset>/` dir; don't reorganize the numbering of existing files.

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

Analysis code (not a package) for the CLAMP manuscript. It fits CLAMP (external R package `chikinalab/CLAMP`) and baseline latent-variable models on pseudobulk scRNA-seq, GTEx, ARCHS4 and recount2, evaluates them, and renders the manuscript figures and tables. README.md has the full data-access, resource and server/SLURM guide.

## Orchestration: Snakemake

Don't run notebooks or scripts by hand when a rule exists.

- `workflow/Snakefile` loads `workflow/config/*.yaml` and includes one `workflow/rules/*.smk` per analysis area. Each area has its own aggregate rule. The public targets are the rules at the bottom of the Snakefile (e.g. `quickstart`, `benchmark`, `archs4_biology`, `publish_figures`; `all` is the default; see README §6).
- Profiles live in `workflow/profiles/`; site-specific ones go in the gitignored `workflow/profiles/site-*/`. Profiles use `rerun-triggers: mtime` on purpose, so editing a rule does **not** invalidate multi-day ARCHS4 results. Use `--forcerun <rule>` to rerun deliberately.
- Rules declare `conda: "<env-name>"` and custom resources. The valid resource names are the `resources:` keys in `workflow/profiles/*/config.yaml`. Reuse them, and add any new one to every profile.

Common commands (run from the repo root with the `snakemake` env active):

```bash
snakemake --profile workflow/profiles/local --cores 1 quickstart            # smoke test, no data/GPU
snakemake --profile workflow/profiles/local --dry-run <target>              # preview jobs
snakemake --profile workflow/profiles/local <rule> --forcerun <rule>        # rerun one notebook/step
snakemake --profile workflow/profiles/local output/99_panels/fig2/fig2.pdf  # one figure
python scripts/check_setup.py                                               # all install/pin checks; no data or analyses
```

Dry-run before any real run: missing upstream outputs can trigger days-long fits and ~1 TB-RAM jobs.

## Environments

`./setup.sh` creates these from `envs/*.lock` (versions are in README §3):
- `clamp-analyses`: R + Python. The default for nearly all rules.
- `gpu-kmeans`: RAPIDS/cuML for GPU clustering. It has no CPU fallback.
- `snakemake`: the orchestrator.
- PhenoPLIER rules use `envs/phenoplier-cli-neo.yaml`, which Snakemake builds itself (see `scripts/phenoplier/README.md`).

CLAMP is pinned in `scripts/install_clamp.R`, and the comparison-method R packages are pinned in `scripts/install_analysis_packages.R`. Both accept `--check` (run in `clamp-analyses`). Don't install or upgrade them any other way; change the pin in the script.

## Code layout and data flow

- `scripts/<area>/`: CLI R/Python scripts (data prep, fits, ORA) called from rule `shell:` blocks with explicit `--flag` arguments. The areas roughly mirror `workflow/rules/*.smk`.
- `nbs/`: report and figure notebooks in numbered stage directories (`00_setup` … `99_panels`). They are executed by Snakemake's `notebook:` directive, so inputs come from the `snakemake` object, not hardcoded paths. Each run writes a sibling `<name>.executed.ipynb`, which is committed. The first cells name the required env. To find the rule that runs a notebook, grep `workflow/rules/` for its path. Some panels (e.g. fig1) are built by scripts in `scripts/panels/`, not notebooks.
- `workflow/config/`: per-area parameters (paths, hyperparameters, seeds, resource sizes) plus `*inputs.json` manifests with the source URLs, versions and sha256 of downloaded inputs. New parameters and seeds go here (seeds reach scripts as `--seed`), not in notebooks or scripts. Don't add to the legacy root `config.R`.
- Root `config.yaml`: shared figure color palettes, read by `scripts/panels/` and the panel notebooks.
- `output/`: all generated artifacts (gitignored; layout in README §8). `publish_figures` copies the final figures and tables into the tracked `figures/`.
- `data/`: downloaded inputs (gitignored), except a small allowlist of bundled references (see `.gitignore`). A few inputs need a manual download (README §4).
- `references/`: fixed cohorts, dataset manifests and recorded timing results (the `timing_mode: recorded|recompute` config switch).

## Tests

Offline `unittest` tests live in `tests/` and `scripts/<area>/tests/`. Run them in the `clamp-analyses` env:

```bash
python -m unittest discover -s tests
python -m unittest tests.test_reference_inputs.ReferenceInputTests.test_download_and_reuse  # single test
```

## Conventions

- R: paths via `here::here()`.
- New analysis: add a script in `scripts/<area>/` and/or a notebook in the matching `nbs/` stage, plus a rule in the matching `.smk`, and hook the rule into that area's aggregate rule. Don't renumber existing notebooks.

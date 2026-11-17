# PhenoPLIER GLS workflow

The workflow pins **phenoplier-cli v0.7.3**, the release used for the October
2026 campaign. Current comparisons use `phenomexcan-curated`: 2,824 available PhenomeXcan rapid-GWAS
traits from the CLI's 2,866-code allow-list. Unfiltered runs contain 4,049 traits.
Changing the trait family changes the cross-trait BH-FDR correction.

To use completed pico results, start with **Collect completed results** below.
The reference setup is needed only for new GLS computations.

## Dedicated software environment

All PhenoPLIER rules, including collection, declare
[`envs/phenoplier-cli-neo.yaml`](../../envs/phenoplier-cli-neo.yaml).
With `--use-conda`, Snakemake creates and activates its own environment under
`.snakemake/conda/`; an existing environment named `phenoplier-cli-neo` is not
required or silently reused. The wrappers keep that activated environment.
The CLI is a hard software dependency for these rules; collection checks it but
does not require workspace initialization or reference data.

The environment pins v0.7.3 by its full Git commit, installed non-editably, plus
NumPy, SciPy, pandas, statsmodels, HDF5 bindings and the inner workflow engine to
versions in that release's `poetry.lock`. PhenoPLIER's Python 3.12 / Snakemake 9
stack stays separate from this repository's Python 3.11 / Snakemake 8 orchestrator.
Git authentication for the private `pivlab/phenoplier-cli` repository is required.

Prepare the environment on a node with package-download access, without running
any analysis or requiring campaign files:

```bash
conda activate snakemake
snakemake -s workflow/Snakefile phenoplier_environment --cores 1 --use-conda \
  --conda-create-envs-only --conda-frontend conda
# Import/version checks only; writes environment.json, never launches GLS:
snakemake -s workflow/Snakefile phenoplier_environment --cores 1 --use-conda
```

The environment report is `output/03_model_biology/03_phenoplier/environment.json`.
These pins improve repeatability but are not a complete binary lock: some
transitive packages, Conda builds and BLAS libraries can still vary. Exact
cross-host numerical identity is not guaranteed (see
[upstream #242](https://github.com/pivlab/phenoplier-cli/issues/242)). A validated
platform lock or container digest would be the next stronger deployment option,
tracked in [upstream #252](https://github.com/pivlab/phenoplier-cli/issues/252).
Snakemake's [environment deployment guide](https://snakemake.readthedocs.io/en/v8.30.0/snakefiles/deployment.html)
describes the YAML and platform pin-file approaches.

For interactive CLI work, `bash scripts/phenoplier/setup_env.sh` creates a named
environment from the same YAML; set `PHENOPLIER_TARGET_ENV` to choose another name.
It refuses to overwrite an existing environment. Activate it before manual CLI
commands. `run_gls.sh` and `store_build.sh` also accept an explicit `--conda-env`
name or prefix for manual use; Snakemake never passes that override.

## Workspace and references (recomputation only)

Using the CLI environment above:

```bash
export PHENOPLIER_HOME=/absolute/path/to/new-phenoplier-workspace
export PHENOPLIER_ROOT_DIR="$PHENOPLIER_HOME"
phenoplier workspace init --home "$PHENOPLIER_HOME"
phenoplier workspace link /path/to/extracted/phenoplier_full_data
phenoplier doctor
```

Use a new workspace for recomputation; old v0.5.x model registries require the
CLI's migration procedure. Set `phenoplier.workspace` in
`workflow/config/phenoplier.yaml` for batch jobs. Set `conda_env_file` only when
intentionally changing the checked-in environment specification.

**A complete public reference download is still unavailable in v0.7.3.** The
`reference_subset` download mode (formerly `full_data`) does not assemble all the
inputs needed for GLS. This is tracked in
[phenoplier-cli #185](https://github.com/pivlab/phenoplier-cli/issues/185).
The lab archives are `/pividori_lab/data/phenoplier_full_data.tar.gz` on pico and
`/pl/active/pivlab/projects/hzhang/data/phenoplier_full_data.tar.gz` on Alpine.
Extract an archive before linking it. `workspace link` also links `models/`,
which must be writable; follow the pinned release's
[shared-installation instructions](https://github.com/pivlab/phenoplier-cli/blob/v0.7.3/docs/guides/shared-installations.md).
Do not claim a successful `reference_subset` download is a complete setup.

The required reference family is **GTEx v8 / MASHR**, not an interchangeable
1000G default. It includes:

- Template GWAS: `data/phenomexcan/gwas_parsing/full/22617_7112.txt.gz`.
- The 49 raw tissue CSVs under
  `data/phenomexcan/gene_assoc/spredixcan/rapid_gwas_project/22617_7112/`,
  named `22617_7112-gtex_v8-{tissue}-2018_10.csv`.
- Template `smultixcan_22617_7112_ccn30.tsv.gz` under
  `data/phenomexcan/gene_assoc/smultixcan/rapid_gwas_project/`.
- The 4,049 phenotype files under `data/smultixcan/rapid_gwas/`, MASHR
  prediction models, GTEx genotype/parquet references, gene metadata and trait
  descriptions. The processed tissue pickle matrices do not replace the raw CSVs.

Repository figure inputs are separate from the CLI workspace. Their acquisition
is declared in `workflow/config/phenoplier_inputs.json`; existing copies can be
staged with `PHENOPLIER_DATA_ROOT` or `phenoplier_input_root`. The exact historical
`smultixcan-mashr-pvalues.pkl` and `fastenloc-torus-rcp.pkl` still lack working
public acquisition routes in that manifest. Public TSVs differ in precision.
See [the input audit](../../references/phenoplier/input_equivalence_audit.json).
The 49 processed S-PrediXcan matrices already have repository download/conversion
rules. These remaining data-access constraints cannot be fixed by a CLI upgrade.

## The October campaign: 192 explicit units

`workflow/config/phenoplier_models.tsv` lists all 192 units delivered under:

```text
/pividori_lab/phenoplier_workspace/projects/marc/10-15-2026-trait-filtering-fix/
```

On November 2, 2026, a read-only pico inventory confirmed every listed unit has
its curated TSV, full-precision pickle, trait-exclusion sidecar and study store.
The 43 previously pending saturation units are now delivered. All 102 new units
have `PASS` entries in `INDEX-new.tsv`; the 90 historical units are described in
`INDEX-curated.tsv`. The preliminary ARCHS4 canonical summary contains 89 LVs;
the earlier handoff incorrectly listed 1,728 for this result.
Directory replicate indices and fitted RNG seeds are recorded separately.

| Batch | Units | Fit traits | Delivered comparison |
|---|---:|---|---|
| Historical reruns, including CLAMPbase and CLAMPfull models | 90 | all 4,049 | curated summary/store |
| Full consensus models | 3 | curated 2,824 | curated summary/store |
| Coverage fractions 1, 5, 10, 25, 50, 75, indices 4–10 | 42 | curated 2,824 | curated summary/store |
| ARCHS4 full-data stability RNG seeds 1001–1007 | 7 | curated 2,824 | curated summary/store |
| GTEx full-data comparator indices 4–10 | 7 | curated 2,824 | curated summary/store |
| Additional saturation K=1728, including 1%/index 1 | 43 | curated 2,824 | curated summary/store |

### Dry runs

The dedicated campaign target covers all 192 units. On pico, or with the campaign
mounted locally:

```bash
snakemake -s workflow/Snakefile phenoplier_campaign -n --cores 2 --use-conda
# Optional recomputation DAG for the 102 new units:
snakemake -s workflow/Snakefile phenoplier_campaign_new -n --cores 28 --use-conda \
  --config phenoplier_campaign_mode=recompute
# Environment check DAG, requiring no scientific inputs:
snakemake -s workflow/Snakefile phenoplier_environment -n --cores 1 --use-conda
```

A dry run resolves inputs and jobs; it does not install software or execute GLS.
Collection needs the delivered files visible at `phenoplier_results_root`.
Recomputing all 192 requires recovery of the historical inputs described below;
the full recompute DAG correctly stops if they are missing. The campaign target
does not include downstream plotting or change the existing figure selection.

### Collect completed results (default)

Activate the repository's `snakemake` environment (8.30.0). On pico:

```bash
conda activate snakemake
snakemake -s workflow/Snakefile phenoplier_campaign --cores 2 --use-conda
```

This target validates each full-precision summary's phenotype × LV grid, P values
and FDR, then creates absolute symlinks to the delivered summaries and stores.
It does not require model matrices, CLI reference data, registration or GLS jobs.
The source tree must remain available. Each unit exposes these files under
`output/03_model_biology/03_phenoplier/campaign/<unit>/`:

- `curated.tsv.gz` and `curated.pkl.gz`;
- `study-phenomexcan-curated.h5` and `curated_excluded.tsv`;
- `source.json` identifying the manifest entry and original files.

The collector checks that all source files exist and are nonempty; it revalidates
the pickle summary, but does not rerun the delivered store's numerical QA.
The original `VERIFY` reports and campaign indices remain in the source tree.
Historical all-trait summaries also remain there under `summary/`.

For a local copy or mounted campaign directory, override its root:

```bash
snakemake -s workflow/Snakefile phenoplier_campaign --cores 2 --use-conda \
  --config phenoplier_results_root=/absolute/path/to/copied-campaign
```

`phenoplier_campaign_new` selects just the 102 new units. Individual output paths
are also valid targets. Missing recorded files cause missing-input errors; there
is no fallback to fitting models. This explicit campaign target is separate from
`all`. Existing figure targets retain their own model selection and paths.

### Recompute explicitly

Only use this mode when new computations are intended, after completing the CLI
setup above. Recomputed outputs go to `campaign-recomputed/`, so delivered and
recomputed files cannot silently satisfy each other's targets.

```bash
snakemake -s workflow/Snakefile phenoplier_campaign_new --cores 28 \
  --resources mem_mb=200000 --use-conda --config phenoplier_campaign_mode=recompute
```

`campaign.source_root: .` connects the 102 new units to this repository's upstream
model rules. Set it to another checkout root to consume existing fitted inputs.
Consensus uses `traits/phenoplier_Z.pkl` and preserves the `lv_mapping.csv` needed
by downstream metaLV consumers.

Recomputing historical units additionally requires their exact matrices under
`<campaign.historical_root>/<unit>/model/<unit>_z.pkl`. The handoff located 70 of
90 under `/pividori_lab/haoyu_projects/trait-filtering-fix/run`; recovery of the
remaining inputs has not been verified. The `canon_archs4_prelim89` input must
contain the original 89 selected LVs: the 1,728-column matrix in the earlier
inventory is not interchangeable. Missing inputs or wrong dimensions fail
explicitly. These limitations do not affect collection of the delivered results.
Historical recomputation fits all 4,049 traits, then recomputes the curated
summary's BH-FDR without a second fit. It also publishes all-trait sidecars.

The existing `phenoplier_coverage`, `phenoplier_saturation`, `phenoplier_finals`,
`phenoplier_canonical`, `phenoplier_final_stores` and `phenoplier_long_table`
targets keep their downstream paths. The historical recount2 final model in the
campaign is not a substitute for the current model selected by the repository.
Do not mix the old biomedical outputs with curated results.

## Store contents and model locations

The delivered and newly built study stores contain association results, trait
metadata and LV identifiers; they do **not** embed the Z matrix, fitted RDS or
hyperparameters. This also applies to the existing final-store targets. Use the
separate upstream model files for loadings and model interpretation. The missing
`/clamp` group is intentional here; see
[phenoplier-cli #235](https://github.com/pivlab/phenoplier-cli/issues/235).

- **New campaign models:** the manifest's `source` column gives each upstream
  path relative to `campaign.source_root` (the checkout by default). Coverage,
  saturation, stability and comparator rules produce `CLAMPbase/Z.csv`; fitted
  model RDS files live beside those exports where the producer supplies them.
- **Consensus models:** `output/98_final_models/consensus/<dataset>/` (CLAMPfull) and `output/98_final_models/consensus_base/<dataset>/` (CLAMPbase)
  contains the model products. `traits/phenoplier_Z.pkl` is the workflow's
  LV-renamed GLS input, and `traits/lv_mapping.csv` maps back to metaLV names.
- **Delivered copies of the 102 new inputs:** `<results_root>/<unit>/model/`
  contains the registered Z pickle and model metadata. The upstream producer
  paths are retained in the manifest; this delivery is useful for checking the
  exact input used by the recorded run.
- **Historical models:** use `campaign.historical_root` plus the manifest's
  `source` path. These old fits are not all reproducible from the current model
  rules; the recovery caveat above still applies. Do not substitute a current
  refit when interpreting a historical store.
- **Existing final-model targets:** the model input is the `CLAMPbase.rds`
  selected by `a4_final_base_rds` in `workflow/rules/phenoplier.smk`, including
  the linked ARCHS4 coverage model.

## Matrix names and downstream mapping

v0.7.3 accepts genes-by-modules CSV/TSV, pandas pickle and fitted RDS inputs in
`phenoplier model register`. It normalizes arbitrary module identifiers to LV
names and saves their original names beside the registered matrix. No B matrix
or fitted RDS is required for matrix-only GLS.

The existing consensus/MGM conversion rules remain because their downstream
consumers use `lv_mapping.csv` (`phenoplier_lv`, `metaLV`) and
`metalv/collect_summary.py`. Consensus uses every module; the CRISPR selection
consumes the consensus ARCHS4 summary and recomputes BH-FDR over its
GSEA-selected modules, as on `main`. The historical
`--rds` wrapper argument accepts matrices too (`--model-input` is an alias).

## Execution and validation

For recomputation, the wrappers activate the CLI environment, then call
`workflow.py`. It registers an input and initializes `run gls` with all scientific and execution options
before launching the generated config. Source content, CLI version, trait family,
cohort and LV percentile determine a separate project/model key. Changed inputs
therefore cannot reuse an old project's results. It does not remove registrations,
projects, Snakemake locks or completion records. Registration and same-project
calls are serialized. The CLI's own run-state checks remain in effect.

For an interrupted recomputation, the CLI project is
`<workspace>/projects/<namespace>/<model>/<fingerprint>/` (campaign namespace:
`clamp/oct2026`). Its `clamp_source.json` records the source fingerprint; the
failed command in the GLS log names `pipeline_config.yaml` in that directory.
From the same CLI environment, inspect and preview recovery first:

```bash
phenoplier run status -p /absolute/path/to/project --workflow gls --verify
phenoplier run resume gls -p /absolute/path/to/project --dry-run
# When ready to continue computations, under the intended allocation/profile:
phenoplier run resume gls -p /absolute/path/to/project
```

Then rerun the original outer Snakemake target to validate and publish its outputs.
A settled previous attempt can be launched again normally. If workers still run
or their ownership cannot be established, resolve that first and follow the CLI's
recovery instructions. The wrapper never deletes locks or forces takeover. For a
SLURM profile, run recovery from the same scheduler context or select the intended
`--executor`/`--profile`; do not launch heavy recovery on a login node.

Defaults: GTEx v8/MASHR, condition number 30, distance 5 Mb, LV percentile 0.01,
`keep-first` duplicate genes, correlation mode `sub`, recomputed LV membership,
and `gene_size gene_size_log gene_density gene_density_log` covariates.
The local executor uses 14 cores per model, with BLAS/OpenMP pinned to one thread.
Prefer two or three concurrent models on an appropriately sized machine. Keep the
CLI's step-3 memory requests; do not lower them to increase chromosome concurrency.
The nested local scheduler is capped at the outer job's CPU and memory allocation;
SLURM fan-out is capped at half the CPU budget because chromosome jobs use two
CPUs each. Each outer GLS job requests 100 GB and 72 hours. Budget roughly 19 GB temporary
disk per model and 2.5 GB delivered output; a full campaign is substantial compute.

The CLI verifies its recorded steps before publication. The summary gate checks
phenotype count, a complete LV × trait grid, duplicates, P values and FDR, and
campaign LV counts. Stores explicitly select the same summary and trait family.
Historical all-trait fits are summarized again with the curated filter; this
recomputes BH-FDR without refitting step 7. Comparing hosts requires numerical
tolerances on the pickle summary (the TSV is rounded), not byte-identical output.

Steps 1–5 are intentionally computed through the supported per-model workflow.
Sharing them by gene set requires CLI cache/provenance support (upstream #148);
this integration does not copy records between models or claim shared caching.

The obsolete hand-run `sbatch` launchers, duplicate Phenoplier READMEs and unused
Phenoplier coverage/QA notebook were removed. Workflow-driven coverage/saturation
reports and the scripts used by consensus and CRISPR analyses remain.

# PhenoPLIER GLS workflow

The workflow pins **phenoplier-cli v0.7.3**, the release used for the October
2026 campaign. Current comparisons use `phenomexcan-curated`: 2,824 available PhenomeXcan rapid-GWAS
traits from the CLI's 2,866-code allow-list. Unfiltered runs contain 4,049 traits.
Changing the trait family changes the cross-trait BH-FDR correction.

To use completed pico results, start with **Collect completed results** below.
The reference setup is needed only for new GLS computations.

## Fresh environment and references

From the repository root:

```bash
bash scripts/phenoplier/setup_env.sh
conda activate phenoplier-cli-neo
export PHENOPLIER_HOME=/absolute/path/to/new-phenoplier-workspace
export PHENOPLIER_ROOT_DIR="$PHENOPLIER_HOME"
phenoplier workspace init --home "$PHENOPLIER_HOME"
phenoplier workspace link /path/to/extracted/phenoplier_full_data
phenoplier doctor
```

The installer needs authenticated access to `pivlab/phenoplier-cli`. It uses the
release's installer and environment, including R support for existing RDS inputs.
Use a new workspace for this campaign; old v0.5.x model registries require the
CLI's migration procedure. Set `phenoplier.workspace` and `conda_env` in
`workflow/config/phenoplier.yaml` for batch jobs. Both wrappers accept an environment
name or its absolute directory.

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

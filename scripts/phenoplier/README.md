# PhenoPLIER GLS workflow

The workflow pins **phenoplier-cli v0.7.3**, the release used for the October
2026 campaign. It replaces v0.5.2 and its retired `biomedical` filter. Current
comparisons use `phenomexcan-curated`: 2,824 available PhenomeXcan rapid-GWAS
traits from the CLI's 2,866-code allow-list. Unfiltered runs contain 4,049 traits.
Changing the trait family changes the cross-trait BH-FDR correction.

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

`workflow/config/phenoplier_models.tsv` records every unit from the October 27
handoff: 149 delivered plus 43 saturation units in flight at that date. Status is
an inventory snapshot, not a completion claim or a reason to exclude a target.
Directory replicate indices and fitted RNG seeds are recorded separately.

| Batch | Units | Fit traits | Published comparison |
|---|---:|---|---|
| Historical reruns, including CLAMPbase and CLAMPfull models | 90 | all 4,049 | all-trait sidecars and curated summary/store |
| Full consensus models | 3 | curated 2,824 | curated summary/store |
| Coverage fractions 1, 5, 10, 25, 50, 75, indices 4–10 | 42 | curated 2,824 | curated summary/store |
| ARCHS4 full-data stability RNG seeds 1001–1007 | 7 | curated 2,824 | curated summary/store |
| GTEx full-data comparator indices 4–10 | 7 | curated 2,824 | curated summary/store |
| Saturation K=1728, including the formerly excluded 1%/index 1 cell | 43 | curated 2,824 | curated summary/store |

Activate the repository's `snakemake` environment (pinned to 8.30.0) first.

```bash
conda activate snakemake
# All 102 new units, using their upstream model producers when sources are absent:
snakemake -s workflow/Snakefile phenoplier_campaign_new --cores 28 --resources mem_mb=200000 --use-conda

# All 192, once the historical matrices are recovered:
snakemake -s workflow/Snakefile phenoplier_campaign --cores 28 --resources mem_mb=200000 --use-conda
```

The campaign config has independent `source_root`, `historical_root` and
`output_root` paths. `source_root: .` connects new units to this repository's model
rules. Set it to another checkout root to consume already fitted inputs there.
Consensus inputs use the workflow's equivalent `traits/phenoplier_Z.pkl`, with its
`lv_mapping.csv`; the campaign's original pickle lived one directory higher.
The converted matrices retain the same values and order.

Historical fits must be staged as
`<historical_root>/<unit>/model/<unit>_z.pkl`. The handoff located 70 of 90 under
`/pividori_lab/haoyu_projects/trait-filtering-fix/run`; the other 20 still need
recovery. This folder can be used directly as `historical_root` on pico. Missing
matrices produce explicit missing-input errors. The campaign never substitutes
current refits for missing historical models. Pickle bytes can differ after
serialization; compare gene/LV labels and values, not only file hashes.

Each unit publishes to `<output_root>/<unit>/`:

- `curated.tsv.gz` and the full-precision `curated.pkl.gz`;
- `study-phenomexcan-curated.h5`, built with that exact summary and trait filter;
- historical units additionally publish `all-traits.{tsv,pkl}.gz`;
- available trait-exclusion sidecars and GLS/store logs.

Individual summary or store paths are also valid Snakemake targets. The campaign
is explicit and separate from `all` so the historical input recovery and the
additional expensive models do not silently expand the existing figure workflow.
Existing `phenoplier_coverage`, `phenoplier_saturation`, `phenoplier_finals`,
`phenoplier_canonical`, `phenoplier_final_stores` and `phenoplier_long_table`
targets retain their downstream paths, but now use v0.7.3 and curated traits.
The old long table retains its existing model selection; it is not the 192-unit
campaign index. Do not mix old biomedical outputs with current curated results;
force the affected GLS targets to rerun when migrating an existing output tree.

## Matrix names and downstream mapping

v0.7.3 accepts genes-by-modules CSV/TSV, pandas pickle and fitted RDS inputs in
`phenoplier model register`. It normalizes arbitrary module identifiers to LV
names and saves their original names beside the registered matrix. No B matrix
or fitted RDS is required for matrix-only GLS.

The existing consensus/MGM conversion rules remain because their downstream
consumers use `lv_mapping.csv` (`phenoplier_lv`, `metaLV`) and
`metalv/collect_summary.py`. Consensus uses every module; CRISPR MGM uses the
GSEA-selected subset. Their mappings and projects are distinct. The historical
`--rds` wrapper argument accepts matrices too (`--model-input` is an alias).

## Execution and validation

The wrappers activate the CLI environment, then call `workflow.py`. It registers
an input and initializes `run gls` with all scientific and execution options
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

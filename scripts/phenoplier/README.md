# PhenoPLIER setup TODO and integration handoff

This is the handoff for the phenoplier-cli maintainer. The CLAMP workflow already
uses **phenoplier-cli v0.7.3** with the expert-reviewed `phenomexcan-curated` trait
filter (2,824 of 4,049 phenotypes) to compute GLS. Existing results were generated
with that CLI and an initialized reference workspace (see "As-run results"). The unresolved task is reproducing
that environment and its inputs for a new user. This document describes the current
workflow after PR #53, not the older CLAMPfull_bp launchers.

## Setup problems to resolve

| Item | Current evidence | Needed for automatic setup |
| --- | --- | --- |
| CLI environment | `scripts/phenoplier/setup_env.sh` installs tag v0.7.3 into `phenoplier-cli-neo`. Repository setup does not call it. Any other release fails the explicit version check. | A reproducible v0.7.3 installation and confirmed reviewer repository access; do not silently upgrade scientific runs. |
| MGM names and matrix input | MGM/consensus columns can be named `metaLV...` or other module IDs, while the current integration requires `LV1`, `LV2`, ... . Their source is `metaZ.csv` / `Z.csv`, not a fitted-model RDS. | Support matrix-only inputs and arbitrary module names, or document a supported conversion with lossless identifier mapping. |
| Workspace | `run_gls.sh` requires `config.toml`. Existing author workspace links the reference bundle's `data/` and `software/`. | Initialization plus acquisition of the actual references below in a new writable workspace. |
| Downloader coverage | At v0.5.2, `full_data` includes 1000G covariance and MASHR models; its S-MultiXcan download entry is explicitly unimplemented. `t21` includes rapid-GWAS results and GTEx covariance. | A supported download route for this exact GTEx/rapid-GWAS GLS workflow. The 1000G panel is not an interchangeable default. |
| GTEx references | Prior audit: Zenodo 11851741 was restricted. Author covariance MD5 `0d805415b43049c0bcb2d59d3579db65` differs from CLI's `0d7895b07665d5d3afab1ba26d445901`. This alone does not prove numerical differences. Step 1 computes covariance and may require the original genotype/parquet references. | Identify the references actually consumed by v0.5.2, whether cached covariance can be reused, and the authorized acquisition route. |
| Template inputs | The original workspace contains the template GWAS, 49 raw tissue CSVs and template S-MultiXcan result. They are not the processed 49 tissue pickle matrices. | Download or reproduce these exact templates, including SNP-count fields consumed by GLS. |
| Downstream exact matrices | Two historical Box links failed during the previous audit; local pickles match upstream historical checksums. Public TSVs differ in precision. | Restore acquisition of the exact P-value and fastENLOC matrices, or validate a scientifically equivalent reconstruction. These are figure/CRISPR inputs, separate from CLI GLS. |

The prior input audit is in [input_equivalence_audit.json](../../references/phenoplier/input_equivalence_audit.json).
It records checks performed previously; access conditions have not been rechecked here.
There is no new preflight system. PhenoPLIER remains TODO; full `all` and
`archs4_biology` execution require this setup even when their DAGs resolve.

## Required inputs

Paths below are relative to the CLI workspace's data directory unless stated otherwise.

| Input | Location / acquisition status |
| --- | --- |
| CLAMP model | Workflow-generated `.rds`; consensus/MGM inputs are genes-by-LV pandas pickles. No external model download is needed. |
| Template GWAS | `phenomexcan/gwas_parsing/full/22617_7112.txt.gz`; exists in author bundle; public preparation remains TODO. |
| Template S-PrediXcan | `phenomexcan/gene_assoc/spredixcan/rapid_gwas_project/22617_7112/22617_7112-gtex_v8-{tissue}-2018_10.csv`, 49 tissues; author bundle. |
| Template S-MultiXcan | `phenomexcan/gene_assoc/smultixcan/rapid_gwas_project/smultixcan_22617_7112_ccn30.tsv.gz`; author bundle. |
| Phenotype results | `smultixcan/rapid_gwas/*.tsv.gz`, 4,049 files; public Zenodo 13274496. Repository helper `scripts/inputs/rapid_gwas_inputs.py` exists; CLI workspace wiring remains TODO. |
| MASHR prediction models | Zenodo 3518299, `mashr_eqtl.tar`; prior audit matched all 98 members to author copies. |
| GTEx panel and gene metadata | Exact workspace/config paths resolved by CLI; genotype/parquet dependencies and covariance reuse need maintainer confirmation. Keep GTEx v8/MASHR. |
| Trait descriptions, coding tables, gene mappings | CLI common download resources. Needed for biomedical filtering and gene identifiers. |
| External software | Bundle's `software/`; maintainer should specify the minimal packages/tools needed by GLS alone. Avoid requiring unused GWAS/TWAS regeneration tools. |

Separate repository inputs are declared in `workflow/config/phenoplier_inputs.json`:
LINCS, gold-standard drug-disease associations, ontology mappings, Z scores and gene
name mappings have pinned public sources. The CRISPR `lipid_DEG.csv` (50,836 bytes)
and phenotype metadata are bundled with provenance. All 49 processed S-PrediXcan
matrices can now be downloaded/converted automatically, with exact numerical
agreement to author copies. They are no longer an unresolved input.

The two unresolved repository matrices are:

- `data/phenomexcan/gene_assoc/smultixcan-mashr-pvalues.pkl`, SHA256
  `f2480ff11d33e815e330579393644a16c64ee7f3c9c8301677a814702d9e274a`.
- `data/phenomexcan/gene_assoc/fastenloc-torus-rcp.pkl`, SHA256
  `c9da922040769cb6f3169edeafe6871bf8ae3817ba466dc5a9e712b5d3a3b341`.

Existing-copy staging supports `PHENOPLIER_DATA_ROOT` or `phenoplier_input_root`.
That is a reuse option, not a complete public acquisition solution.

## Models and expected GLS outputs

Path abbreviations: `A=output/01_model_building/02_archs4`,
`F=output/98_final_models/clampfull/canonical`,
`C=output/98_final_models/consensus_models`,
`M=output/03_model_biology/02_archs4/04_crispercas/00_mgm`.
`ds` is `archs4`, `gtex` or `recount2`.

| Family | Models / project name | Published GLS output |
| --- | --- | --- |
| Coverage CLAMPbase | 70: fractions 1,5,10,25,50,75,100 x replicates 1-10 (`archs4.coverage.pathway` seeds); `cov_rs{f}_seed{s}_CLAMPbase`, except 100% replicates 2-10 named `archs4_rs100_seed{2,3,1001..1007}_CLAMPbase`. Replicates 1-3 start from `CLAMPbase.rds`, 4-10 from `CLAMPbase/Z.csv` | `A/clampbase_coverage/summaries/cov_rs{f}_seed{s}.tsv.gz` |
| GTEx coverage comparators | 9: replicates 2-10, RNG seeds 456,789,1001-1007; `gtex_rs100_seed{rep}_CLAMPbase` | `A/clampbase_coverage/comparators/gtex/summaries/gtex_rs100_seed{rep}.tsv.gz` |
| Saturation CLAMPbase GLS | 17: fractions 1,5,10,25,50,75 x replicates 1,2,3 at K=1728, excluding (1,1728,1); `sat_rs{f}_k1728_seed{s}_CLAMPbase`. **Pending: 43 more** (replicates 4-10 at every fraction plus (1,1728,1)); their models are `A/01_saturation/study_saturation_rs{f}_k1728_seed_{s}/CLAMPbase/Z.csv` | `A/clampbase_saturation/summaries/sat_rs{f}_k1728_seed{s}.tsv.gz` |
| Final CLAMPbase | `phenoplier.trait_finals`: GTEx `final_gtex_CLAMPbase`; ARCHS4 reuses coverage 100%, replicate 1. **Pending:** recount2 (the current Gavish-Donoho model, 364 LVs; the as-run summary is the earlier 544-LV model) | `A/clampbase_finals/summaries/final_{ds}.tsv.gz` |
| Canonical CLAMPfull | 3: `canon_{ds}_CLAMPfull_canonical`; input `F/{ds}/CLAMPfull_canonical.rds` | `F/{ds}/traits/canon_{ds}.tsv.gz` |
| Consensus (MGM) | 3: `consensus_{ds}_CLAMPfull`, all metaLVs; input `C/{ds}/traits/phenoplier_Z.pkl` | `C/{ds}/traits/gls_summary.tsv.gz`; mapped `consensus_{ds}.tsv.gz` |

`consensus_models` requests the consensus GLS. The CRISPR MGM selection
(`crispercas.gls_summary`) reads `C/archs4/traits/consensus_archs4.tsv.gz` and
recomputes Benjamini-Hochberg over the GSEA-significant metaLVs; there is no
separate CRISPR MGM GLS run.

## Snakemake integration: model outputs to downstream results

PhenoPLIER is a stage of the same Snakemake DAG. Each GLS job takes a model or
matrix produced by an upstream rule as a declared `input`, then publishes a
summary as a declared `output`. Downstream analysis rules consume those summaries
as inputs. No manual model copying, GLS invocation or result placement should be
required once the separate environment/workspace setup is resolved.

| Upstream model / matrix | GLS rule | Downstream consumer |
| --- | --- | --- |
| Coverage CLAMPbase model outputs | `gls_bp_coverage_base` | `aggregate_archs4_coverage_traits`, `aggregate_archs4_saturation_traits` (100% references), and `link_final_base` |
| Saturation CLAMPbase outputs at **K=1728 only** | `gls_bp_saturation_base` | `aggregate_archs4_saturation_traits`, then `saturation_traits_report_archs4` |
| GTEx coverage comparator models | `gls_gtex_coverage_base` | `aggregate_archs4_coverage_traits` |
| Final GTEx/recount2 CLAMPbase models | `gls_final_base` | Coverage/final trait aggregation and `store_final_base` |
| Canonical CLAMPfull outputs from `fit_canonical_clampfull` | `gls_canonical` | `aggregate_phenoplier_traits` |
| Coverage / GTEx replicate `CLAMPbase/Z.csv` (replicates without an `.rds`) | `phenoplier_z_input` -> `gls_bp_coverage_base` / `gls_gtex_coverage_base` | as above |
| Consensus `Z.csv` from `consensus_model` | `consensus_gls_input` -> `consensus_gls` | `consensus_traits` maps LV names back to module identifiers; `crispercas_metalv_select` reads the ARCHS4 result |

The current saturation selection is `saturation_k_values: [1728]` in
`workflow/config/phenoplier.yaml`. `phenoplier_saturation` requests exactly those
17 GLS outputs; the report uses those summaries plus full-data coverage references.
Other K values are fitted for the saturation model analysis, not for PhenoPLIER.

Snakemake resolves upstream producers from the declared model paths (or explicit
`rules.<name>.output` references). Missing models are generated first, followed by
GLS and its consumers. Requesting a downstream trait/figure target therefore
includes the required GLS jobs. Updated saturation model files invalidate their
GLS outputs; the wrapper fingerprints model content to prevent stale project reuse.
The environment and reference acquisition issues above remain
TODO even though these model-to-summary dependencies are already connected.

Additional model families generated before these analyses:

- Coverage: 70 CLAMPbase cells under `A/02_coverage_study/<configured-level>/study_coverage_rs{f}_seed_{s}/` (replicates 4-10 hold `CLAMPbase/Z.csv` only).
- Saturation model fitting covers 107 cells across K=86,173,432,864,1296,1728 under `A/01_saturation/`. **PhenoPLIER runs only on K=1728: 17 models**, excluding fraction 1%, replicate 1. No GLS runs are requested for K=86,173,432,864 or 1296.
- Canonical finals: `F/{ds}/` contains `CLAMPfull_canonical.rds`, `Z.csv`, `B.csv`, `summary.csv`, `manifest.json`.
- Stability seeds: 456,789,1001,1002,1003,1004,1005,1006,1007. ARCHS4 uses `A/05_stability`; GTEx/recount2 use `output/01_model_building/{01_gtex,03_recount2}/stability/seed{seed}/CLAMPfull_canonical/`.
- Each compendium contributes its reference plus nine seed CLAMPfull models to consensus. ARCHS4 reference is S123; GTEx/recount2 references are S124. GTEx reuses two coverage base fits, adding seven; recount2 adds nine base fits. GTEx/recount2 together add 18 CLAMPfull fits.
- Consensus: `C/{ds}/` contains `Z.csv`, `B.csv`, `manifest.json`, `metaZ.csv`, cluster membership/summary, `models/` links, `noise/`, `ora.csv`, `ora/`, and `lincs-projection.tsv`. ARCHS4 reuses the CRISPR MGM clustering.
- Optional final CLAMPbase stores: `A/clampbase_finals/stores/final_{ds}_CLAMPbase.h5`, target `phenoplier_final_stores`.
- Combined traits: `output/03_model_biology/02_archs4/05_traits/phenoplier_traits_long.csv.gz`, target `phenoplier_long_table`.

The authoritative lists are `workflow/config/{archs4,consensus,phenoplier}.yaml`
and the corresponding rule files. Pseudobulk benchmark models are not GLS inputs.

## MGM input compatibility: module names and Z matrices

MGM and consensus models have module identifiers such as `metaLV1`, `metaLV2`,
or other names. These names are not accepted by the current PhenoPLIER integration,
which expects `LV1`, `LV2`, ... . The original identifiers must survive in the
published results so traits can be joined back to clustering, ORA and CRISPR outputs.

These analyses start from **gene-loading matrices (`metaZ.csv` or `Z.csv`), not a
fitted CLAMP `.rds` object**. The current workaround is:

1. Read the CSV with gene identifiers as rows and module identifiers as columns.
2. Rename columns to consecutive `LV1`, `LV2`, ... and save the correspondence in
   `lv_mapping.csv` (`phenoplier_lv`, `metaLV`).
3. Serialize the matrix as a pandas pickle and pass that file to the CLI's `--input`.
4. Map GLS results back to the original module identifiers with `collect_summary.py`.

For consensus, `consensus_gls_input` converts the full `Z.csv` and keeps
`lv_mapping.csv`; the CRISPR MGM analysis uses the ARCHS4 consensus result.
The wrapper's argument is still called `--rds`, even though it passes these pickles
through to `phenoplier shortcut gls --input`.

**Maintainer TODO:** confirm and document matrix-only GLS support without requiring
an RDS, B matrix or fitted-model metadata; define the accepted orientation and gene
identifiers; support arbitrary module names or preserve a supplied name mapping
through registration, GLS summaries and exports. A direct `Z.csv` input should have
a documented conversion or native support. The repository currently supplies a
converted pickle, not the CSV itself. This compatibility path needs a fresh-workspace
test; existing RDS-based runs alone do not validate it.

## CLI contract and runtime behavior

Workspace selection: `phenoplier.workspace` config, then `PHENOPLIER_HOME`, then
`~/phenoplier`. Both `PHENOPLIER_HOME` and `PHENOPLIER_ROOT_DIR` are exported.
A writable `projects/` and model registry are required. The wrapper serializes
registration with `flock`; each model otherwise has an independent project.

`run_gls.sh` calls `shortcut gls --unlock` to initialize, then `workflow gls run
--config-file <workspace>/projects/<name>/pipeline_config.yaml`. It fingerprints
inputs, preserves sentinels for resumption, and moves stale projects aside when
an input fingerprint changes. Legacy projects without fingerprints are reused
with a warning; use a fresh workspace for reproduction.

Parameters: namespace `clamp`, cohort `phenomexcan_rapid_gwas`, GTEx v8/MASHR,
LV percentile 0.01, `phenomexcan-curated` trait filter, expected 2,824 phenotypes
(4,049 before filtering). Step 6 uses one worker x eight BLAS threads; step 7 uses six
workers x two BLAS threads. Outer GLS jobs request 12 CPUs, 100 GB, 4,320 minutes.
The nested executor defaults to `local`, including when the outer Snakemake job
uses Pico's `server_cu` submission profile.

CLI summary source:
`<workspace>/projects/<name>/results/gls/phenoplier/gls-summary-phenomexcan_rapid_gwas.tsv.gz`.
The wrapper publishes it to the paths above, together with available pickle and
trait-filter exclusion sidecars. Logs sit alongside each requested summary.
Summaries must contain `phenotype`, `lv`, `pvalue`, a complete phenotype x LV
grid and no duplicate pairs or zero/nonfinite P values; NaNs are allowed only for
flagged degenerate LVs. Downstream scripts also consume FDR and trait annotations.
MGM/consensus exports rename columns to LV1, LV2, ... and preserve a mapping CSV
back to metaLV IDs.

## Maintainer acceptance goal

From a new environment and empty workspace: install the pinned CLI, obtain its
references through documented commands, initialize a small synthetic project and
validate input resolution without reusing author outputs. Then validate one small
GLS case against the pinned release. Full model fits and a full GLS campaign are
outside the repository's lightweight validation. No files have been sent externally.

## As-run results

The summaries in `output/` come from the pico campaign
`server_cu:/pividori_lab/phenoplier_workspace/projects/marc/10-15-2026-trait-filtering-fix/`
(one folder per project name above; `INDEX-curated.tsv` and `INDEX-new.tsv` list
models, LVs and hits). Each published summary is that project's
`summary/gls-summary-phenomexcan_rapid_gwas-phenomexcan-curated.tsv.gz`, copied to
the output path its rule declares. The earlier biomedical (2,366-trait) summaries
are kept in `output/archive/traits_biomedical_20261026/`.

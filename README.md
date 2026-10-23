# CLAMP analyses

## Setup

CLAMP uses three Conda environments: `clamp-analyses` (core CPU modeling/analysis),
`gpu-kmeans` (GPU clustering), and `snakemake` (pipeline orchestrator).

```bash
./setup.sh
```

Installs Conda if needed, creates all three environments, and installs the pinned
[`chikinalab/CLAMP`](https://github.com/chikinalab/CLAMP) R package into
`clamp-analyses`. Safe to re-run. Verify with `nbs/00_setup/00_check_setup.ipynb`.

> [!WARNING]
> Requires CLAMP commit
> [`748559ae6d4a9e2982d8f4488bae616fe390195a`](https://github.com/chikinalab/CLAMP/commit/748559ae6d4a9e2982d8f4488bae616fe390195a) —
> do not update it independently.

<details>
<summary>Manual setup, step by step</summary>

1. Install Conda (Miniconda or Mambaforge).
2. (Optional, for `gpu-kmeans`) verify GPU drivers with `nvidia-smi`.
3. Create the environments:
   ```bash
   conda create --name clamp-analyses --file envs/clamp-analyses.lock
   conda run -n clamp-analyses python -c '
   import subprocess
   import sys
   import yaml
   with open(sys.argv[1]) as handle:
       dependencies = yaml.safe_load(handle)["dependencies"]
   requirements = [requirement for dependency in dependencies if isinstance(dependency, dict)
                   for requirement in dependency.get("pip", [])]
   if not requirements:
       raise SystemExit("No pip requirements found in the environment YAML")
   subprocess.check_call([sys.executable, "-m", "pip", "install", *requirements])
   ' envs/clamp-analyses.yaml
   conda create --name gpu-kmeans --file envs/gpu-kmeans.lock
   conda env create -n snakemake -f envs/snakemake.yaml
   ```
4. Install the CLAMP R package:
   ```bash
   conda activate clamp-analyses
   Rscript scripts/install_clamp.R
   Rscript scripts/install_clamp.R --check   # verify only
   ```
</details>

`envs/clamp-analyses.lock` pins the Conda packages and builds. Pip package
versions are read from the `pip:` section of `envs/clamp-analyses.yaml`; there
is no separate pip lock. `envs/snakemake.yaml` defines the separate workflow
runner environment.

Each notebook states its required environment in its first cell. Run all `snakemake`
commands from the repo root with `conda activate snakemake` and `--use-conda`:

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile <target>
```

## Fresh checkout

After setup and staging the external inputs below, run from the repository root:

```bash
conda activate snakemake
snakemake --snakefile workflow/Snakefile --use-conda --cores 16 all
```

`all` is also the default target. It includes publication figures and eight
supplementary tables, canonical model ORA, coverage/saturation reports, trait
recovery, drug–disease analysis, and the CRISPR metaLV analysis. Full ARCHS4 fits
require a high-memory cluster; use a site profile as described below. GPU
clustering requires the `gpu-kmeans` environment and a compatible NVIDIA GPU.
Figure 4 requires Inkscape on PATH.

The repository does not contain the large scientific inputs. Stage these before
starting; paths and column mappings are defined in `workflow/config/*.yaml`.
Snakemake downloads ARCHS4 v2.5, GTEx v8 expression, recount2, the configured GEO
and EBI projection datasets, DrugBank, and the configured Ensembl reference files.
It generates the two distinct GO:BP files through their existing rules.

| External input | Required local location and version |
| --- | --- |
| MSigDB pathways | `data/pathways/{c2.all,c2.cp,c2.cp.reactome,c5.go.cc,c8.all,h.all}.v2026.1.Hs.symbols.gmt`, obtained from MSigDB using its access requirements. Keep this release fixed. |
| Cell and tissue references | `data/pathways/Cell_marker_Human.xlsx` (sheet `human`, columns `cell_name`, `Symbol`), `Allen_Brain_Atlas_10x_scRNA_2021.gmt`, and `data/archs4/GTEx_Tissues_pathMat.rds`. Azimuth 2023 and ARCHS4 Cell-lines have download rules. |
| ARCHS4 gene lengths | `data/archs4/gene_lengths.rds`: the named gene-symbol length vector for the Ensembl 107 preprocessing input. This is an external reference input, not a fitted model. |
| Pseudobulk cohorts | Raw matrices and cell metadata at the exact paths in `workflow/config/pseudobulk.yaml`. For datasets with `build_pseudobulk: false`, also supply `bulk_expr.csv`, `patient_info.csv`, `truthFrac_v0.csv`, and, where used, `truthFrac_v1.csv` in their configured `pseudobulk_dir`. These fixed cohort inputs preserve the original sample selections. |
| Cell labels | `data/pseudobulk/cell_type_labels.csv`, used by Supplementary Figure 1. |
| GTEx annotations | `data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt` and `GTEx_Analysis_v8_xCell_scores_7_celltypes.txt`. |
| PhenoPLIER reference bundle | Set up phenoplier-cli v0.5.2 with `bash scripts/phenoplier/setup_env.sh`, initialize its workspace, and link the `phenoplier_full_data` bundle as described in `scripts/phenoplier/README.md`. Set `phenoplier.workspace` or `PHENOPLIER_HOME` for your machine. |
| PhenomeXcan | Stage its `phenomexcan` directory at `data/phenomexcan`, including phenotype metadata, S-MultiXcan P values/z-scores, fastENLOC RCPs, and gene mappings. Set `archs4.crispercas.phenomexcan_dir` to override it. S-PrediXcan needs the 49 tissue pickles under `data/phenomexcan/gene_assoc/spredixcan/pkl`. |
| Drug–disease inputs | `data/drug_disease_associations/{lincs-data.pkl,gold_standard.pkl,phenomexcan_traits_fullcode_to_efo.tsv,term_id_xrefs.tsv.gz,xrefs-prop-slim.tsv}`. |
| CRISPR screen | `data/archs4/crispr_cas9/lipid_DEG.csv` and `data/archs4/drug_diseases_associations/smultixcan-mashr-zscores.pkl`, as configured in `archs4.crispercas`. |
| Timing evidence | The recorded CoGAPS evidence files specified by `workflow/config/runtime_benchmark*.yaml`. |

Keep site-specific accounts, partitions, and storage paths in your own profile
or configuration override. Existing outputs and local symlinks are conveniences;
a fresh checkout must stage the listed external inputs itself.

## Pseudobulk

Benchmarks CLAMP against other latent-variable/matrix-decomposition methods on pseudobulk
profiles built from single-cell/single-nucleus RNA-seq cohorts (brain, heart, PBMC, lung),
then projects the fitted latent variables back onto individual cells to check biological
interpretability. A donor-bulk extension collapses the same cells into one real bulk-like
library per donor, as an independent check that results aren't an artifact of pseudobulk
construction. Configuration: `workflow/config/pseudobulk.yaml`.

```bash
# End to end
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile biology_pseudobulk panels_pseudobulk

# Donor-bulk extension
snakemake --cores 8 --resources donor_bulk_io=1 --use-conda \
  --snakefile workflow/Snakefile donor_bulk_report fig2_panel
```

## GTEx

Fits CLAMP and comparison methods on GTEx v8 bulk tissue expression, then checks whether
the latent variables recover known tissue/subtissue structure and specific biology (e.g.
liver cell-type composition via xCell). Configuration: `workflow/config/gtex.yaml`.

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile biology_gtex
```

## ARCHS4

CLAMP fit across the full ARCHS4 human RNA-seq compendium (~605k samples). Its largest
fits need ~500 GB RAM and run on Slurm, not a workstation.

`download_archs4` fetches the raw HDF5 and `preprocess_archs4` builds the filtered,
z-scored matrix under `output/01_model_building/02_archs4/00_preprocess/` that every
ARCHS4 fit reads. The full-compendium SVD and CLAMPbase are the 100% cells of the
coverage campaign.

The coverage and saturation campaigns fit CLAMPbase only.

The GO:BP ORA database is `data/pathways/go_bp.Hs.symbols.gmt` (SHA-256
`33e559df968d0d8c2028c0f735624555dc2edbebe955648c3351dfa31cfd2e2a`).
Snakemake generates this file when missing from the GOALL Biological Process
mapping in `org.Hs.eg.db` 3.20.0 and term names in `GO.db` 3.20.0, then verifies
its checksum. The `pathway_prior` rule separately downloads the Enrichr file
`GO_Biological_Process_2025.gmt`, which is not interchangeable with this input.

Canonical-prior CLAMPfull models are published under
`output/98_final_models/clampfull/canonical/`, one each for ARCHS4, GTEx, and
recount2. They are evaluated by ORA against GO:BP, canonical, Reactome, and
CellMarker. Drug--disease and projection analyses consume these canonical artifacts
from their dedicated workflows. Missing models are fitted by Snakemake. See [Canonical model + ORA](#canonical-model--ora).

Configuration: `workflow/config/archs4.yaml` and `workflow/config/recount2.yaml`.

### Trait recovery (PhenoPLIER GLS)

LV–trait association runs `pivlab/phenoplier-cli` GLS on the coverage, saturation
and final CLAMP models, and the trait-recovery reports aggregate the per-model
summaries. phenoplier-cli is pinned to **v0.5.2** and runs in its own conda env
(`scripts/phenoplier/setup_env.sh`, once per machine; every job checks the installed
version); `workflow/rules/phenoplier.smk` drives it one model at a time from this
DAG. Configuration: `workflow/config/phenoplier.yaml`; details in
`scripts/phenoplier/README.md`.

```bash
snakemake -n --snakefile workflow/Snakefile archs4_traits      # models -> GLS -> reports
snakemake -n --snakefile workflow/Snakefile phenoplier_finals  # GLS on the published finals
```

Dry-run first: a GLS job takes hours to a day per model, and any model that is not
already published under `output/98_final_models` is fitted first (the ARCHS4 fits
need a cluster, see above).

### Canonical model + ORA

`canonical_models` builds the shared C2-all + Reactome + CellMarker prior, then
fits one canonical CLAMPfull model for ARCHS4, GTEx, and recount2 using their
workflow-produced expression matrices, SVDs, ranks, and CLAMPbase fits. The
`archs4_canonical_ora` target evaluates those models. Coverage and saturation
remain CLAMPbase-only. Drug–disease analyses use only canonical CLAMPfull models,
plus the gene-space prediction baseline.

```bash
snakemake --cores 16 --use-conda --snakefile workflow/Snakefile canonical_models
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile archs4_canonical_ora
```

Complete existing canonical and stability CLAMPfull bundles are reused without
requiring their training inputs. A bundle contains the model RDS, `B.csv`,
`Z.csv`, `summary.csv`, `manifest.json`, and `complete`. Partial bundles cause an
error: restore their missing files or move the bundle aside before rebuilding.
Do not force these fitting rules on an existing bundle.

The canonical prior has 11,107 sets (7,670 C2-all, 1,839 Reactome, 1,598
CellMarker), MD5 `8ea2e09825753a90273eaa2dc03176f6`. Historical fit manifests
record `e668595a34ec9da67b69c87fc95a5c29`: that GMT has identical terms and gene
memberships, with `CellMarker` instead of `CellMarker_Human` in its description
field. Historical adoption records remain as records; downstream rules use the
fit manifest. ARCHS4 uses seed 123 and 5,000 iterations; GTEx/recount2 use seed
124 and 500 iterations, all with multiplier 100. Rebuilds record their actual
software versions and are not guaranteed to produce byte-identical models.

### Pathway replicates (Figure 3a, c, d)

Pathway recovery uses ten CLAMPbase replicates per group (`coverage.pathway` in
`workflow/config/archs4.yaml`). Study subsamples at 1-75% use sampling seeds 123, 124, 125,
33 and 1001-1006; each is fitted at its estimated K (SVD start) and at each fixed K (random
start). ARCHS4 100% (seeded SVD and random start for replicates 2-10) and GTEx use seeds
123, 456, 789 and 1001-1007. Replicates 1-3 are the coverage and
saturation cells; replicates 4-10 are fitted by `clampbase_pathway_cell`,
`clampbase_full_pathway_cell`, `fit_saturation_pathway_base` and `clampbase_gtex_pathway_replicate`, which keep only
`CLAMPbase/Z.csv` and `manifest.json`. Where a subsample's SVD has fewer components than the
requested K (1% at K = 1,728, `extended_svd`), `svd_saturation_extended` computes an
additional SVD with K components; the original SVD and estimated-K fit are unchanged.

```bash
snakemake --profile workflow/profiles/local coverage_pathway_ora saturation_pathway_ora
snakemake --profile workflow/profiles/local fig3_panel
```

## Running on a Slurm cluster

`workflow/profiles/slurm/config.yaml` holds the portable part of a cluster setup
(executor, job limits, memory pools, per-rule overrides for the ARCHS4 fits). Account,
partition and QOS are site-specific and stay out of the repo: derive a user-owned
profile from it, e.g. `~/.config/snakemake/clamp-slurm/config.yaml`, and add

```yaml
default-resources:
  - slurm_account=<account>
  - slurm_partition=<partition>
```

plus any `set-resources` a partition needs. The `snakemake` env already includes
`snakemake-executor-plugin-slurm`. Run the scheduler from a login node inside
`tmux`/`screen` (it lives for the whole DAG), always dry-running first:

```bash
conda activate snakemake
snakemake --profile ~/.config/snakemake/clamp-slurm -n <target>
snakemake --profile ~/.config/snakemake/clamp-slurm <target>
```

The conda envs must be visible from the compute nodes; `use-conda` looks them up by
name. `scripts/phenoplier/README.md` §6 has the complete site-profile example
(including the PhenoPLIER GLS jobs) and the monitoring/resume notes.



## Running a single notebook

Each notebook-backed rule runs one specific notebook and writes the executed copy
back to it. Target it by rule name like any other rule:

```bash
snakemake --cores 1 --use-conda --snakefile workflow/Snakefile liver_disentangle_xcell_rf_true_labels_gtex
```

Snakemake skips a rule whose outputs are already newer than its inputs, so editing
a notebook's cells alone won't trigger a re-run. Force one with `-f`/`--forcerun`
(add `-R`/`--forceall` to also force everything downstream of it):

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile -f <target>
```

## Rebuilding publication figures from notebooks

Run these from the repository root in the `clamp-analyses` environment after the
upstream analyses have produced their input files:

```bash
for name in fig2 fig3 supp1 supp2 supp3 supp4 supp5 supp6 supp7; do
  papermill "nbs/99_panels/${name}.ipynb" "nbs/99_panels/${name}.executed.ipynb" -k ir
done
```

Each notebook writes PDF, PNG, and SVG files to `output/99_panels/<name>/`.
The committed `figures/` folder is a copy of these outputs and of
`output/99_panels/tables/`; do not edit it by hand. Refresh it with the
`publish_figures` rule after regenerating any figure or table, then copy
`figures/<name>/<name>.*` into `clamp-manuscript/content/images/final_figures/`:

```bash
snakemake publish_figures --snakefile workflow/Snakefile --cores 1 --allowed-rules publish_figures
```

Figure 2 can also run through its Snakemake rule:

```bash
snakemake output/99_panels/fig2/fig2.pdf --snakefile workflow/Snakefile --cores 4 --use-conda -R fig2_panel
```

| Figure | Required upstream data |
| --- | --- |
| Figure 2 | Pseudobulk benchmark, grouped cross-validation, eight-by-eight hard-pair marker ORA grid, donor-bulk recovery, Perez single-cell projection, and GTEx xCell tables listed in its input-path cell and `workflow/rules/panels.smk`. |
| Figure 3 | ARCHS4 coverage and saturation tables under `output/03_model_biology/02_archs4/00_coverage` and `01_saturation`. |
| Figure 4 | `scripts/panels/assets/projection_scheme.svg` (rasterized with Inkscape), cytokine projection aggregates under `03_projections/aggregate`, and drug-disease tables under `02_drug_diseases_canonical/three_model_comparison/data`. |
| Figure 5 | CRISPR-Cas9 metaLV tables under `output/03_model_biology/02_archs4/04_crispercas/00_mgm/fig5/` (rule `crispercas_fig5_data`): metaLVs clustered across the ten CLAMP-full canonical models, PhenoPLIER GLS traits, ARCHS4 activity and PhenomeXcan S-MultiXcan P values and fastENLOC RCPs. |
| Supplements 1 and 3 | Pseudobulk benchmark and marker-recovery tables, and GTEx clustering, held-out subtissue confusion, and composition tables, respectively. |
| Supplement 2 | Donor-bulk recovery and UMAP tables under `output/03_model_biology/00_pseudobulk/06_donor_bulk_recovery/`. |
| Supplement 4 | Pseudobulk and GTEx computational timing outputs under `output/02_model_performance/`. |
| Supplement 5 | ARCHS4 coverage and saturation tables under `output/03_model_biology/02_archs4/`. |
| Supplement 6 | ARCHS4 projection aggregate results, model loadings under `output/98_final_models/`, projection benchmark results, and pathway GMT/marker files under `data/pathways/` (the required local input directory). The notebook regenerates its own `source_data` files. |
| Supplement 7 | Drug-disease prediction pickle and per-tissue comparison tables under `output/03_model_biology/02_archs4/02_drug_diseases_canonical/`; the notebook converts the pickle with the environment's Python. |
| Supplement 8 | metaLV clustering (`Z_all.npz`, `cluster_members.csv`) and CRISPR GSEA under `output/03_model_biology/02_archs4/04_crispercas/00_mgm/`. |
| Supplementary Tables 1-8 | `scripts/tables/build_supplementary_tables.py` (rule `supplementary_tables`), writing to `output/99_panels/tables/`. |

Figures 3-5 share the helpers in `scripts/panels/main_figure.R` (180 mm maximum width,
170 mm maximum height, text clamped to 5-7 pt).

## Citation

## License

This project is licensed under the [CC-BY 4.0 License](http://creativecommons.org/licenses/by/4.0/).

## Acknowledgments

Supported by the National Human Genome Research Institute,
The Eunice Kennedy Shriver National Institute of Child Health and Human Development,
the National Science Foundation, and the National Eye Institute.

## Revised Supplementary Figure 6 projection workflow

The Supplementary Figure 6 benchmarks (permutation-adjusted pseudobulk recovery,
stratified gene nulls and GTEx ARI) are built by `projection_benchmark_figure` in
`workflow/rules/projection_benchmarks.smk`; `supp6_panel` renders the figure.

The benchmark uses `benchmark_project.R` for model preparation and projection,
`benchmark.py` for scoring and clustering, and `benchmark_panels.R` for plotting.
CRISPR metaLV scripts live directly under `scripts/archs4/crispercas/`; clustering,
selection, and figure-data scripts expose subcommands for their smaller steps.

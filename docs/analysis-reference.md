# Analysis and provenance reference

For installation, inputs, profiles and execution, start with the [README](../README.md).

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

CLAMP fit across the full ARCHS4 human RNA-seq compendium (~605k samples). Use the hardware requirements and execution profiles in the main README.

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

## Running a single notebook

Each notebook-backed rule runs one specific notebook and writes the executed copy
back to it. Target it by rule name like any other rule:

```bash
snakemake --cores 1 --use-conda --snakefile workflow/Snakefile liver_disentangle_xcell_rf_true_labels_gtex
```

The profiles use mtime-based rerun triggers; inspect the dry-run after edits. Force one with `--forcerun <rule>`
(`--forceall` also forces dependencies and can launch expensive fits):

```bash
snakemake --cores 4 --use-conda --snakefile workflow/Snakefile --forcerun <target>
```

## Rebuilding publication figures from notebooks

Use the Snakemake panel rules from the repository root; this activates the
correct environment and builds declared upstream inputs:

```bash
snakemake --profile workflow/profiles/server_cu fig2_panel fig3_panel fig4_panel fig5_panel
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
| Figure 5 | CRISPR-Cas9 metaLV tables under `output/03_model_biology/02_archs4/04_crispercas/03_metalv/fig5/` (rule `crispercas_fig5_data`): metaLVs clustered across the ten CLAMP-full canonical models, PhenoPLIER GLS traits, ARCHS4 activity and PhenomeXcan S-MultiXcan P values and fastENLOC RCPs. |
| Supplements 1 and 3 | Pseudobulk benchmark and marker-recovery tables, and GTEx clustering, held-out subtissue confusion, and composition tables, respectively. |
| Supplement 2 | Donor-bulk recovery and UMAP tables under `output/03_model_biology/00_pseudobulk/06_donor_bulk_recovery/`. |
| Supplement 4 | Pseudobulk and GTEx computational timing outputs under `output/02_model_performance/`. |
| Supplement 5 | ARCHS4 coverage and saturation tables under `output/03_model_biology/02_archs4/`. |
| Supplement 6 | ARCHS4 projection aggregate results, model loadings under `output/98_final_models/`, projection benchmark results, and pathway GMT/marker files under `data/pathways/` (the required local input directory). The notebook regenerates its own `source_data` files. |
| Supplement 7 | Drug-disease prediction pickle and per-tissue comparison tables under `output/03_model_biology/02_archs4/02_drug_diseases_canonical/`; the notebook converts the pickle with the environment's Python. |
| Supplement 8 | metaLV clustering (`Z_all.npz`, `cluster_members.csv`) and CRISPR GSEA under `output/03_model_biology/02_archs4/04_crispercas/03_metalv/`. |
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

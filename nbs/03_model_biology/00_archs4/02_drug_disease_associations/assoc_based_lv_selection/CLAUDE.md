# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working in this directory.

Read the **project-root** `CLAUDE.md` and the **parent**
`../CLAUDE.md` (the drug-disease pipeline) first. This file only covers what is specific to
this isolated sub-analysis.

## Running commands

Use the `clamp-analyses` conda env (e.g. `conda run -n clamp-analyses <cmd>`), same as the
parent pipeline. All paths use `here()` rooted at the project, so notebooks run from anywhere.

## What this directory is

An **isolated variant** of the ARCHS4 module-based drug→disease method that changes how the
top-`ntc` latent variables (LVs) are selected per disease.

- Parent `07_prediction_module_based_archs4.ipynb` selects the top-`ntc` LVs by
  **`|projected value|`** (magnitude of the CLAMP-projected S-PrediXcan z-scores).
- This variant selects the top-`ntc` LVs by their **GLS LV→disease association p-value**
  (smallest p = strongest association), read from
  `data/archs4/traits/hall_coverage_rs100_seed_1_CLAMPfull_hall/gls-summary-phenomexcan.tsv.gz`.

Motivation: in the gene-based method (`06`) the disease vector *is* the gene→disease
association (S-PrediXcan z-scores), so selection is driven by real disease associations. The
parent LV method's projected value is **not** a real LV→disease association, just a linear
transform of the gene z-scores. This variant makes the LV method use the proper GLS
LV→disease associations for selection.

**Important:** the GLS file has only `pvalue`/`fdr` — **no signed effect statistic** — so the
p-value can only drive *selection*. The dot-product score still uses the projected (signed)
values: `score = -1 * drug^T @ disease` in LV space (the `predict_dotprod_neg` convention).

## Pipeline (run in order)

1. `01_prediction_module_based_archs4_assoc.ipynb` — produces the prediction HDF5s
   (method name `module_based_archs4_assoc`), 49 tissues × 5 thresholds `{all, 5, 10, 25, 50}`.
2. `02_prediction_performance.ipynb` — isolated copy of parent `10`; aggregates AUROC/AUPRC
   across the 4 parent methods (06–09) **plus** this variant.
3. `03_prediction_performance_plots.ipynb` — isolated copy of parent `11`; ROC/PR figures and
   Fig-4 CSVs for the 5 methods.

## Isolation guarantees

- The parent pipeline (notebooks `06`–`13`) and the shared `output/99_panels/fig4` are
  **never modified**. `02`/`03` are self-contained copies that read parent prediction outputs
  but write only under
  `output/03_model_biology/00_archs4/02_drug_disease_associations/assoc_based_lv_selection/`.
- The one shared edit is additive and backward-compatible: `libs/drug_disease_utils.py` gained
  `_zero_nontop_by_ranking` and an optional `selection_ranking=None` kwarg on
  `predict_dotprod_neg`. With the default `None`, behaviour for `06`–`09` and `12` is unchanged.

## Conventions specific to this directory

- **GLS trait-name reconciliation.** The GLS file keys traits by a *base* code (e.g.
  `100001_raw`), while the projection columns / `ukb_efo` index use *full* codes with a
  description suffix (e.g. `100001_raw_Food_weight`). `01` reconciles them with a
  longest-underscore-boundary prefix match (`fullcode_to_base`): ~4,048/4,091 columns map,
  with zero collisions.
- **Missing-GLS fallback.** ~43 projection columns are external consortium GWAS (Astle blood
  counts, BCAC breast cancer, etc.) with no GLS data. For those, selection falls back to the
  parent `|projected value|` ranking. This keeps the evaluated (drug, disease) universe
  identical to parent `07`, so the AUROC comparison in `02` stays apples-to-apples.
- Because `ntc=None` (the `all_genes` case) applies no thresholding, this variant's
  `all_genes` predictions are **identical** to parent `07`'s — a useful sanity anchor.

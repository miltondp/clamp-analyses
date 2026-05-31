# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Read the **parent** `CLAUDE.md` (`../CLAUDE.md`, the drug-disease pipeline) and the **project-root**
`CLAUDE.md` first. This file only covers what is specific to this `tissue_agg_test/` directory.

## What this directory is

The **max-over-49-tissues ablation** (`NEXT_STEPS.md` Active concern #1). The pipeline collapses each
(drug, disease) pair's 49 tissue scores with a **max** before pooling into one AUROC
(`../signif_test/00_aggregate_predictions.ipynb`, the `_reduce_max` step). The memo flagged this as a
possible winner's-curse confound and as an unjustified asymmetry (mean over the 5 `n_top_genes`
thresholds, but max over tissues).

Two findings reshaped the concern and set the approach (recorded so they are not relitigated):
- **The "order statistic grows with n" mechanism is absent.** `signif_test/00` hard-asserts every
  pair has **exactly 49 tissues** — coverage is constant, so the memo's "more tissues → higher max"
  sub-mechanism and its "does coverage correlate with `true_class`" check are moot.
- **Mean/median-of-scores is not a clean ablation.** `max` encodes "one relevant tissue carries the
  signal for this pair"; mean assumes broad sharing. Swapping max→mean tests a *different model*, so
  it cannot adjudicate whether `max` is fishing noise.

So the instrument here is the alternative the memo itself suggested: compute **per-tissue AUROC/AUPRC,
then aggregate across tissues**. This removes the per-pair tissue-selection step entirely and asks the
decisive question: **does the cross-method ordering survive removing the max?**
- Ordering survives per-tissue → it is present in the average single tissue, **not a max artifact**.
- Edge present only under max → it lives in the selection step (real tissue-specificity *or* winner's
  curse), which only the deferred **analysis B** (tissue-selection biological validation) can separate.

It is a **read-only consumer** of the NB06–09 prediction HDF5s and `../signif_test/predictions_paired.pkl`.
**Do not** modify NB 06–13, `../signif_test/`, `libs/drug_disease_utils.py`, `config.R`, or the
scoring-sign convention from here.

## Notebooks (run in order, `clamp-analyses` env)

```
00_per_tissue_metrics.ipynb  →  per_tissue_metrics.csv (+ per_tissue_scores.pkl,
                                  max_aggregate_reference.csv)  →  01_per_tissue_summary.ipynb
```

- **`00_per_tissue_metrics.ipynb`** — slow I/O step (~3 min, re-reads the 980 NB06–09 HDF5s). Copies
  the settings / loading loop / fail-loud completeness asserts **verbatim from `signif_test/00`**
  (rank `score` over the full DOID distribution → inner-merge gold standard), then runs **only step 1**
  of that aggregation — `groupby([trait, drug, method, tissue]).apply(_reduce_mean)` (mean of ranks
  across the 5 `n_top_genes` thresholds) — and **deliberately stops before the `_reduce_max` over
  tissues**, keeping the tissue axis. Computes AUROC / AUPRC / `auprc_log2_enrich` per `(method, tissue)`.
  Outputs `per_tissue_scores.pkl` (134,260 rows = 685×4×49, reusable for analysis B), `per_tissue_metrics.csv`
  (196 rows = 4×49), and `max_aggregate_reference.csv` (recomputed from `predictions_paired.pkl`).

- **`01_per_tissue_summary.ipynb`** — fast statistics. Distribution summary per method (macro-mean /
  median / IQR of per-tissue AUROC and log2-AUPRC-enrichment) beside the max-aggregate, with a
  `selection_gain` column; a descriptive per-tissue ordering (B-wins-in-k/49 + Wilcoxon signed-rank
  across tissues); a **tissue-cluster paired bootstrap** (resample the 49 tissues, recompute each
  method's mean-per-tissue metric on the same resample — paired, `seed=42`, `N_BOOT=10000`); 6
  pairwise differences with CI + two-sided bootstrap p + BH within each metric family. Plots +
  `FINDINGS.md`.

## Load-bearing details

- **Coverage is constant at 49 tissues** (asserted in NB00). Every tissue scores the *identical* 685
  pairs, so the per-tissue label split is constant (531 pos / 154 neg, base_rate ≈ 0.775) and **every**
  per-tissue AUROC/AUPRC is defined — no eligibility filtering (unlike `per_disease_test`, where the
  *disease* node varied).
- **The cluster unit is the tissue.** `signif_test` resamples *pairs*; `per_disease_test` resamples
  *diseases*; here the estimand is "per-tissue metric, then aggregate," so we resample the 49 **tissues**
  (same paired-across-methods property: one tissue resample per iteration, applied to every method,
  diffs by subtraction).
- **CIs may be OPTIMISTIC.** The 49 tissues are **not independent** (shared 685 pairs, shared genes/LVs,
  GTEx inter-tissue correlation). The tissue-cluster bootstrap treats them as exchangeable, so it is
  anti-conservative — state this and read significance in the repo's "directional / consistent" framing.
  A pair-resample sensitivity is the honest robustness follow-up (vectorizable for AUROC, costly for
  AUPRC) but is deferred.
- **Mean per-tissue AUROC and the max-aggregate AUROC are DIFFERENT ESTIMANDS** (mean single-tissue
  discrimination vs best-tissue-selected discrimination). The contrast sizes the value of selection; it
  is not a like-for-like recomputation. Keep this framing in any prose.
- Helpers (`boot_ci`, `boot_pvalue_two_sided`, `bh_correct`) and `_get_tissue` / the loading loop are
  **copied verbatim** from `../signif_test/`, not imported, to keep the directory isolated. If the
  upstream aggregation changes, mirror it here.

## Paths

- Inputs (resolved via `pyprojroot.here()`):
  - NB06–09 HDF5s: `output/.../0{6,7,8,9}_prediction_*/lincs/predictions/dotprod_neg/*.h5`
  - Gold standard: `data/drug_disease_associations/gold_standard.pkl`
  - Status-quo frame: `output/.../signif_test/predictions_paired.pkl`
- Outputs (gitignored):
  `output/03_model_biology/00_archs4/02_drug_disease_associations/tissue_agg_test/`
  — `per_tissue_scores.pkl`, `per_tissue_metrics.csv`, `max_aggregate_reference.csv`,
  `per_tissue_macro_summary.csv`, `per_tissue_ordering.csv`, `per_tissue_bootstrap_results.csv`/`.pkl`,
  `figures/*.png`, `FINDINGS.md`.

## Out of scope this round (flagged, not done)

- **Analysis B — tissue-selection biological validation** (which tissue wins per correct prediction;
  concentrated/sensible vs uniform-random). `per_tissue_scores.pkl` retains the tissue axis precisely
  so this can be done later without re-reading the HDF5s.
- **The second max** — the UKB-trait→DOID collapse in `map_traits_to_doid` (`libs/drug_disease_utils.py`
  line 53), which the memo notes "compounds." Would be re-aggregated from the HDF5 `full_prediction` key.
- The standalone variance-vs-label diagnostic was **dropped**: across-tissue variance is mechanistically
  entangled with the legitimate "one relevant tissue" signal, so it cannot separate confound from signal.

## Verification

```
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  00_per_tissue_metrics.ipynb --output 00_per_tissue_metrics_executed.ipynb
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  01_per_tissue_summary.ipynb --output 01_per_tissue_summary_executed.ipynb
```

Confirm: `per_tissue_scores.pkl` has 134,260 rows; `per_tissue_metrics.csv` has 196 rows (4×49), no
NaNs, constant 531/154 label split; `max_aggregate_reference.csv` AUROC reproduces 0.583 / 0.625 /
0.602 / 0.612 (the NB00 sanity assert enforces this); `per_tissue_bootstrap_results.csv` has 12 rows
(2 metrics × 6 comparisons) with `p_value` + per-metric `p_value_bh`; figures render; `FINDINGS.md`
written with the per-tissue-vs-max-aggregate ordering comparison as the headline.

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Read the **parent** `CLAUDE.md` (`../CLAUDE.md`, the drug-disease pipeline) and the **project-root**
`CLAUDE.md` first. This file only covers what is specific to this `per_disease_test/` directory.

## What this directory is

The **per-disease performance** analysis (`NEXT_STEPS.md`, "Suggested order of attack" #1). NB10/11
report **pooled** AUROC/AUPRC across all (drug, disease) pairs; the paired bootstrap in
`../signif_test/` found no pairwise AUROC difference significant. Pooled AUROC also cannot separate
real signal from **disease-popularity** confounds.

This directory computes AUROC/AUPRC **within each disease** (across that disease's candidate drugs),
then summarizes the distribution across diseases. Fixing the disease node removes disease-level
popularity **by construction** — the rigorous version of "restrict to popular diseases."

It is a **read-only consumer** of `../signif_test/predictions_paired.pkl` (the already-aligned
685-pair × 4-method frame). **Do not** modify NB 06–13, `../signif_test/`, `libs/drug_disease_utils.py`,
`config.R`, or the scoring-sign convention from here.

## Notebooks (run in order, `clamp-analyses` env)

```
00_per_disease_metrics.ipynb  →  per_disease_metrics.csv  →  01_per_disease_summary.ipynb
```

- **`00_per_disease_metrics.ipynb`** — load `predictions_paired.pkl`, pivot to wide (shared
  `(trait, drug)` index), and for each disease × method compute `auroc`, `auprc`, `n_pos`, `n_neg`,
  `n_total`, `base_rate = n_pos/n_total`, and `auprc_log2_enrich = log2(auprc / base_rate)` (AUPRC
  fold-enrichment over the disease prior; 0 = no better than prior). AUROC defined only for diseases
  with ≥1 pos & ≥1 neg. `eligible_strict` flags diseases with ≥3 pos & ≥3 neg. Saves the long frame
  `per_disease_metrics.csv` (one row per `(trait, method)`).

- **`01_per_disease_summary.ipynb`** — distribution summary per method (macro-mean AUROC, median,
  IQR, frac>0.5; macro-mean **and** median `auprc_log2_enrich`), a **disease-level cluster bootstrap**
  (resample *diseases* with replacement, recompute each method's macro-mean on the same resampled
  disease set — paired across methods, `seed=42`, `N_BOOT=10000`), and the 6 pairwise macro-mean
  differences (CI + two-sided bootstrap p + BH within each metric family). Plots + `FINDINGS.md`.

## Load-bearing details

- **Macro-mean weights diseases EQUALLY** (one vote per disease). This is the estimand that removes
  disease popularity; n-weighting would re-pool it and is only ever a sensitivity footnote.
- **The cluster unit is the disease, not the pair.** The existing `signif_test` bootstrap resamples
  *pairs*; here we resample *diseases*, the independent unit once the disease node is fixed. Same
  paired-across-methods property (one resampled disease set per iteration, applied to every method),
  so differences are by subtraction.
- **No single-class skip needed.** Every AUROC-eligible disease has its per-disease metric defined,
  so resampling diseases never produces an undefined value (unlike the pair-level bootstrap).
- **AUPRC headline is the log2 fold-enrichment over the per-disease prior**, not raw AUPRC: base
  rates span ~0.36–0.95, so raw AUPRC mostly reads prevalence. Report macro-mean **and median**
  (the log ratio can produce large negative outliers for poorly-ranked diseases).
- Helpers (`boot_ci`, `boot_pvalue_two_sided`, `bh_correct`) are **copied verbatim from
  `../signif_test/01_paired_bootstrap_test.ipynb`**, not imported, to keep the directory isolated.

## Paths

- Input: `../signif_test/predictions_paired.pkl` (cols `trait` DOID, `drug` DB#, `method`, `score`,
  `true_class` ∈ {0,1}; 2740 rows = 685 pairs × 4 methods, 57 diseases), resolved via
  `pyprojroot.here()`.
- Outputs (gitignored):
  `output/03_model_biology/00_archs4/02_drug_disease_associations/per_disease_test/`
  — `per_disease_metrics.csv`, `per_disease_macro_summary.csv`, `per_disease_bootstrap_results.csv`,
  `figures/*.png`, `FINDINGS.md`.

## Known structure (verify, don't assume)

Of the 57 diseases in the 685-pair universe: **33** are AUROC-eligible (≥1 pos & ≥1 neg), **17** are
strict (≥3 & ≥3); 24 are AUROC-undefined (19 all-positive, 5 all-negative). The undefined diseases
are the most label-imbalanced — their exclusion is non-random (ties to `NEXT_STEPS.md` Active #3) and
must be stated as a limitation. Per-disease n is small (median ~14), so individual AUROCs are
near-categorical — only distributions / macro-means are interpretable, never per-disease rankings.

## Verification

```
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=python3 \
  00_per_disease_metrics.ipynb --output 00_per_disease_metrics_executed.ipynb
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=python3 \
  01_per_disease_summary.ipynb --output 01_per_disease_summary_executed.ipynb
```

Confirm: `per_disease_metrics.csv` has 33 AUROC-eligible diseases (`auroc` non-null) and 17 with
`eligible_strict=True`; **strict-17** macro-mean per-disease AUROC ≈ archs4 0.649 / gene 0.598
(matches the `NEXT_STEPS.md` preliminary — note the ordering **flips** on the full 33-disease set,
where gene ≈ 0.658 ≥ archs4 ≈ 0.636); `per_disease_bootstrap_results.csv` has 6 comparisons × 2
metric families × 2 eligibility sets (24 rows) with `p_value_bh`; figures render; `FINDINGS.md`
written. The one BH-significant difference is ARCHS4 > GTEx on the strict subset (both metrics);
ARCHS4-vs-gene is directional but n.s.

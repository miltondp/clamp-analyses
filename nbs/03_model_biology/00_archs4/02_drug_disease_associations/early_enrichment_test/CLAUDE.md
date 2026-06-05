# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Read the **parent** `CLAUDE.md` (`../CLAUDE.md`, the drug-disease pipeline) and the **project-root**
`CLAUDE.md` first. This file only covers what is specific to this `early_enrichment_test/` directory.

## What this directory is

The **early-enrichment** analysis (`NEXT_STEPS.md`, "Switch headline metric to early-enrichment").
NB10/11 and `../signif_test` report **AUROC/AUPRC**, which integrate over the *whole* ranked list.
Drug repurposing is a **top-of-the-list** problem — you only act on the handful of candidates you
would take into validation. This directory scores the **pooled global** ranking (all 685 pairs
together) with metrics that weight the top:

- **BEDROC** (Truchon & Bajorath 2007), α=20 — the headline single number (top ~5%).
- **Enrichment factor** EF@{1, 5, 10%} and **precision@k** at k ∈ {10, 20, 50} — interpretability.

It is a **read-only consumer** of `../signif_test/predictions_paired.pkl` (the aligned 685-pair ×
4-method frame). **Do not** modify NB 06–13, `../signif_test/`, `config.R`, or the scoring-sign
convention from here.

## Notebooks (run in order, `clamp-analyses` env)

```
00_early_enrichment_metrics.ipynb  →  early_enrichment_observed.csv  →  01_early_enrichment_summary.ipynb
```

- **`00_early_enrichment_metrics.ipynb`** — load `predictions_paired.pkl`, pivot to a shared
  `(trait, drug)` index, and compute pooled BEDROC / EF@{1,5,10%} / P@{10,20,50} per method. Saves
  the long frame `early_enrichment_observed.csv` (one row per `(metric, method)`, with the integer
  cut and the null/ceiling recorded).
- **`01_early_enrichment_summary.ipynb`** — **paired bootstrap** (resample the 685 *pairs*, recompute
  every method's metric on the *same* resample, diffs by subtraction; `seed=42`, `N_BOOT=10000`), the
  6 pairwise differences per metric (CI + two-sided bootstrap p + **BH within each metric family**),
  four figures, and `FINDINGS.md`.

## Load-bearing details

- **HIGH PREVALENCE → LOW HEADROOM is the central caveat.** Pooled base rate = 531/685 = **77.5%
  positive** (the *inverse* of the rare-positive regime these metrics were built for). So
  **BEDROC's null = base rate ≈ 0.775, NOT 0.5**; **EF is capped at 1/0.775 ≈ 1.29**; precision@k's
  null = 0.775. Read **differences** between methods, not absolute enrichment. State this first in any
  writeup.
- **`rdkit` is not in the `clamp-analyses` env** — BEDROC/EF/precision@k are hand-rolled in numpy
  (verified: random ≈ 0.775, perfect = 1.0, worst = 0.0). Re-pasted into both notebooks (not
  imported) to keep the directory isolated.
- **Ties** are broken deterministically by **input (pivot) order** via `np.argsort(-scores,
  kind='stable')` — fully reproducible, no random tie-breaking.
- **EF χ→K via `ceil`**: K = 7 / 35 / 69 at N=685 (recorded in the observed CSV).
- **BH within each of the 7 metric families** (each metric = one family of 6 comparisons → 42 rows),
  mirroring how `../signif_test` keeps AUROC and AUPRC separate.
- **Single-class skip guard kept** (`len(np.unique(yb)) < 2`) for parity with `../signif_test`;
  at 77.5% prevalence it essentially never triggers.
- **The cluster unit is the pair** (same as `../signif_test`), not the disease — this is the
  pooled-global view. `../per_disease_test` is the complementary per-disease (use-case) framing.
- Helpers `boot_ci`/`boot_pvalue_two_sided`/`bh_correct` are **copied verbatim from
  `../signif_test/01_paired_bootstrap_test.ipynb`**.

## Paths

- Input: `../signif_test/predictions_paired.pkl` (cols `trait` DOID, `drug` DB#, `method`, `score`,
  `true_class` ∈ {0,1}; 2740 rows = 685 pairs × 4 methods), resolved via `pyprojroot.here()`.
- Outputs (gitignored):
  `output/03_model_biology/00_archs4/02_drug_disease_associations/early_enrichment_test/`
  — `early_enrichment_observed.csv`, `early_enrichment_bootstrap_results.csv`/`.pkl`,
  `figures/early_enrichment_{bedroc,ef,pk}_dist.png`, `figures/early_enrichment_diff_forest.png`,
  `FINDINGS.md`.

## Known result (verify, don't assume)

No pairwise difference is BH-significant for any of the 7 metrics. Observed BEDROC ordering
**gene 0.944 > recount2 0.910 > ARCHS4 0.902 > GTEx 0.890** — the gene baseline concentrates true
treatments at the very top at least as densely as the LV models, i.e. switching the headline to the
early region does **not** hand the LV models an advantage. EF@1% = the 1.29 ceiling for *every*
method (the top-7 is all-positive everywhere) — a vivid illustration of the low-headroom caveat.
This reinforces the `../signif_test` H4 null. All four methods sit well above the 0.775 BEDROC null,
so the top of the list is genuinely enriched; the high prevalence caps cross-method daylight.

## Verification

```
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=python3 \
  00_early_enrichment_metrics.ipynb --output 00_early_enrichment_metrics.ipynb
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=python3 \
  01_early_enrichment_summary.ipynb --output 01_early_enrichment_summary.ipynb
```

Confirm: `early_enrichment_observed.csv` has 28 rows (4 methods × 7 metrics) with `value ≤ ceiling`;
`early_enrichment_bootstrap_results.csv` has **42 rows** (7 metrics × 6 comparisons) with
`p_value_bh`; the four figures render; `FINDINGS.md` leads with the high-prevalence/low-headroom
caveat. Sanity: a random-score control gives BEDROC ≈ 0.775 and EF ≤ 1.29.

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Read the **parent** `CLAUDE.md` (`../CLAUDE.md`, the drug-disease pipeline) and the **project-root**
`CLAUDE.md` first. This file only covers what is specific to this `signif_test/` directory.

## What this directory is

The **REVIEW H4 fix** — notebooks 10/11/13 report cross-method AUROC (gene-based ~0.583,
ARCHS4 ~0.625, GTEx ~0.602, recount2 ~0.612) with **no paired significance test**, so the ~4-pp gap
the analysis hinges on was unsupported.

This directory implements a **paired bootstrap** across the **full 4-method grid** — single-gene
baseline (NB06, `gene_based`), ARCHS4 (NB07, `module_based_archs4`), GTEx (NB08, `module_based_gtex`),
recount2 (NB09, `module_based_recount2`) — for **all 6 pairwise comparisons**, in both AUROC and AUPRC.

> History: it began as an isolated 2-method prototype (the "smallest meaningful comparison",
> gene vs ARCHS4) to validate the methodology before generalizing. It now covers the full grid and is
> the intended H4 deliverable. It remains a **read-only consumer** of the 06–09 prediction outputs.

**Do not** modify NB 06–13, `libs/drug_disease_utils.py`, `config.R`, or the scoring-sign convention
from here.

## Notebooks (run in order, `clamp-analyses` env)

```
00_aggregate_predictions.ipynb  →  predictions_paired.pkl  →  01_paired_bootstrap_test.ipynb
```

- **`00_aggregate_predictions.ipynb`** — slow I/O step. Reads the raw 06/07/08/09 prediction HDF5
  files and reproduces **NB10's exact aggregation** so all four methods are scored on an *identical*
  set of (drug, disease) pairs:
  1. Per file: `score = score.rank()` over the full DOID distribution, **then** inner-merge with the
     gold standard (this order matters — copied from NB10).
  2. `groupby([trait,drug,method,tissue]).apply(_reduce_mean)` — mean of ranks across the 5
     `n_top_genes` thresholds.
  3. `groupby([trait,drug,method]).apply(_reduce_max)` — max across the 49 tissues.

  Saves the long frame `predictions_paired.pkl` (`[trait, drug, method, score, true_class]`,
  shape `4 × N_PREDICTIONS` = `2740` rows). The set of methods is driven by `METHOD_THRESHOLDS` /
  `PREDICTIONS_DIRS`; `EXPECTED_METHODS = tuple(METHOD_THRESHOLDS)` and every assert loops over it, so
  adding/removing a method is a config-only change.

- **`01_paired_bootstrap_test.ipynb`** — fast statistics; the H4 fix. Pivots to a wide frame whose
  shared `(trait, drug)` index **is** the pairing, then runs **one** bootstrap (`seed=42`,
  `N_BOOT=10000`): each iteration draws a single resample and recomputes *every* method's metric on
  it (per-method bootstrap arrays). The 6 pairwise differences are derived by **subtracting** those
  shared-resample arrays (`metric_b − metric_a`), which keeps each comparison paired and makes all
  comparisons mutually consistent. Reports, per (comparison, metric): per-method 95% CIs, the
  difference's 95% CI, an `ci_excludes_0` flag, a two-sided bootstrap p-value, and a BH-corrected
  column. Saves `paired_bootstrap_results.csv`/`.pkl`, a **forest plot** `paired_bootstrap_diff.png`,
  and an appendix histogram grid `paired_bootstrap_hist.png`.

  - **Sign convention:** for a comparison `(A, B)` taken in `METHOD_ORDER =
    [gene_based, module_based_archs4, module_based_gtex, module_based_recount2]`, `diff = B − A`
    (later minus earlier). `gene_based` is first, so gene-vs-module rows read "module − gene"
    (> 0 ⇒ the module model beats the baseline).
  - **Results schema** (one row per `(metric, comparison)`, 12 rows): `metric, comparison, method_a,
    method_b, a_obs, a_lo, a_hi, b_obs, b_lo, b_hi, diff_obs, diff_boot_mean, diff_ci_lo, diff_ci_hi,
    ci_excludes_0, p_value, p_value_bh`. Column names `a_*`/`b_*` are metric-agnostic (the metric is
    named in the `metric` column).
  - **BH is applied within each metric family** (the 6 AUROC comparisons together, the 6 AUPRC
    separately) via `results.groupby('metric')['p_value'].transform(bh_correct)` — AUROC and AUPRC are
    distinct hypothesis families, so they are not pooled into one correction set.

## Paths

- Helpers (`_get_tissue`, `_reduce_mean`, `_reduce_max`) are **copied verbatim from NB10**, not
  imported, to keep the directory isolated. If you change the aggregation in NB10, mirror it here
  (or vice versa). `_get_tissue` already strips the `-projection-{archs4,gtex,recount2}` suffixes.
- Inputs resolved via `pyprojroot.here()`:
  - Gold standard: `data/drug_disease_associations/gold_standard.pkl`
    (cols `trait` DOID, `drug` DB#, `true_class` ∈ {0,1}).
  - NB06: `output/.../06_prediction_single_gene_based/lincs/predictions/dotprod_neg/*.h5`
  - NB07: `output/.../07_prediction_module_based_archs4/lincs/predictions/dotprod_neg/*.h5`
  - NB08: `output/.../08_prediction_module_based_gtex/lincs/predictions/dotprod_neg/*.h5`
  - NB09: `output/.../09_prediction_module_based_recount2/lincs/predictions/dotprod_neg/*.h5`
- Outputs (gitignored): `output/03_model_biology/00_archs4/02_drug_disease_associations/signif_test/`.

## Load-bearing details

- **The pairing is the whole point.** All four methods must be defined on the *identical* `(trait,
  drug)` set (NB01 asserts this), and each bootstrap iteration applies the *same* resample indices to
  every method. Resampling methods independently would break the paired test; the shared-resample-once
  loop guarantees it across all comparisons.
- **Completeness is asserted (fail-loud).** NB00 requires 49 tissues × 5 thresholds per method
  (`actual == N_TISSUES * len(thresholds) * N_PREDICTIONS`) for **all of NB06–NB09**, and asserts each
  `PREDICTIONS_DIRS` entry exists with a message naming the source notebook to run. Do not weaken these
  to work around partial data — the AUROC's "max across tissues" step needs all 49 tissues to reproduce
  the published numbers (a 3-tissue subset gives ~0.555 / ~0.598 for gene / ARCHS4).
- **NB08/NB09 must be run separately first.** At the time of writing they had not been executed (no
  output on disk), so NB00 **fails loud** at the `assert d.exists()` for the GTEx/recount2 dirs — by
  design, until those prediction jobs finish. NB10 was also unrun.
- `N_PREDICTIONS` is expected to be **685** unique (drug, disease) pairs (matches NB10).

## Verification

**Tier 1 — regression check, runnable without NB08/NB09.** Confirms the generalized bootstrap is
behavior-preserving against the original 2-method result. Restrict `METHOD_ORDER` to the two methods
present in an existing 2-method `predictions_paired.pkl` and run NB01's logic — it must reproduce:
AUROC `diff_obs = +0.042037`, 95% CI `[−0.000199, 0.084255]`, `p = 0.052195`; AUPRC `p = 0.650335`.
(The shared-resample-once loop is identical to the original per-pair loop under the same seed.)

**Tier 2 — full grid, once NB08 + NB09 are on disk.** Run headless in order:

```
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  00_aggregate_predictions.ipynb --output 00_aggregate_predictions_executed.ipynb
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  01_paired_bootstrap_test.ipynb --output 01_paired_bootstrap_test_executed.ipynb
```

Confirm: `predictions_paired.pkl` has `2740` rows; per-method sanity AUROC ≈ 0.583 / 0.625 / 0.602 /
0.612 (gene / archs4 / gtex / recount2); the results CSV has 12 rows with `p_value` + per-metric
`p_value_bh`; the forest plot renders. **The headline H4 answer is the
`ci_excludes_0` / `p_value` / `p_value_bh` on each AUROC comparison row** (gene-vs-module rows in
particular).

## p-value convention

`boot_pvalue_two_sided` uses the `(count+1)/(N+1)` convention (Davison & Hinkley 1997), so the
two-sided p-value is floored at `2/(N+1)` and never reads exactly `0`. The 95% CI on the difference
remains the primary inference.

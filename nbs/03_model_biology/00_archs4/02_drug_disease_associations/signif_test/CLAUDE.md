# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Read the **parent** `CLAUDE.md` (`../CLAUDE.md`, the drug-disease pipeline) and the **project-root**
`CLAUDE.md` first. This file only covers what is specific to this `signif_test/` sandbox.

## What this directory is

An **isolated prototype** for **REVIEW H4** — the top-priority open issue from the parent directory:
notebooks 10/11/13 report cross-method AUROC (gene-based ~0.583, ARCHS4 ~0.625, GTEx ~0.602,
recount2 ~0.612) with **no paired significance test**, so the ~4-pp gap the analysis hinges on is
currently unsupported.

This sandbox implements a **paired bootstrap** on the smallest meaningful comparison — single-gene
baseline (NB06, `gene_based`) vs ARCHS4 module model (NB07, `module_based_archs4`) — for both AUROC
and AUPRC. It is kept separate so the methodology can be validated before being generalized to the
full 4-method grid in NB 10/12/13.

**Do not** modify NB 06–13, `libs/drug_disease_utils.py`, `config.R`, or the scoring-sign convention
from here — this is a read-only consumer of their outputs.

## Notebooks (run in order, `clamp-analyses` env)

```
00_aggregate_predictions.ipynb  →  predictions_paired.pkl  →  01_paired_bootstrap_test.ipynb
```

- **`00_aggregate_predictions.ipynb`** — slow I/O step. Reads the raw 06/07 prediction HDF5 files
  and reproduces **NB10's exact aggregation** so the two methods are scored on an *identical* set of
  (drug, disease) pairs:
  1. Per file: `score = score.rank()` over the full DOID distribution, **then** inner-merge with the
     gold standard (this order matters — copied from NB10).
  2. `groupby([trait,drug,method,tissue]).apply(_reduce_mean)` — mean of ranks across the 5
     `n_top_genes` thresholds.
  3. `groupby([trait,drug,method]).apply(_reduce_max)` — max across the 49 tissues.

  Saves the long frame `predictions_paired.pkl` (`[trait, drug, method, score, true_class]`,
  shape `2 × N_PREDICTIONS`).

- **`01_paired_bootstrap_test.ipynb`** — fast statistics; the actual H4 fix. Pivots to a wide frame
  whose shared `(trait, drug)` index **is** the pairing, then resamples the *same* pairs for both
  methods (`seed=42`, `N_BOOT=10000`) and records the per-resample difference (ARCHS4 − gene). Reports
  per-method 95% CIs, the difference's 95% CI, an `ci_excludes_0` flag, a two-sided bootstrap p-value,
  and a BH-corrected column. Saves `paired_bootstrap_results.csv`/`.pkl` and `paired_bootstrap_diff.png`.

## Paths

- Helpers (`_get_tissue`, `_reduce_mean`, `_reduce_max`) are **copied verbatim from NB10**, not imported,
  to keep the sandbox isolated. If you change the aggregation in NB10, mirror it here (or vice versa).
- Inputs resolved via `pyprojroot.here()`:
  - Gold standard: `data/drug_disease_associations/gold_standard.pkl`
    (cols `trait` DOID, `drug` DB#, `true_class` ∈ {0,1}).
  - NB06: `output/.../06_prediction_single_gene_based/lincs/predictions/dotprod_neg/*.h5`
  - NB07: `output/.../07_prediction_module_based_archs4/lincs/predictions/dotprod_neg/*.h5`
- Outputs (gitignored): `output/03_model_biology/00_archs4/02_drug_disease_associations/signif_test/`.

## Load-bearing details

- **The pairing is the whole point.** Both methods must be defined on the *identical* `(trait, drug)`
  set (NB01 asserts this), and each bootstrap iteration applies the *same* resample indices to both
  methods. Resampling the methods independently would break the paired test.
- **Completeness is asserted.** NB00 requires 49 tissues × 5 thresholds per method
  (`actual == N_TISSUES * len(thresholds) * N_PREDICTIONS`). NB07 was incomplete (42/49 tissues) at the
  time this was written, so NB00 fails at that assert **by design** until NB07 finishes all 49 tissues.
  Do not weaken the assert to work around partial data — the AUROC's "max across tissues" step needs all
  49 tissues to reproduce the published 0.583 / 0.625 (a 3-tissue subset gives ~0.555 / ~0.598).
- `N_PREDICTIONS` is expected to be **685** unique (drug, disease) pairs (matches NB10).

## Verification

Once NB07 has all 49 tissues on disk, run headless in order:

```
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  00_aggregate_predictions.ipynb --output 00_aggregate_predictions_executed.ipynb
conda run -n clamp-analyses jupyter nbconvert --to notebook --execute \
  01_paired_bootstrap_test.ipynb --output 01_paired_bootstrap_test_executed.ipynb
```

Confirm: `predictions_paired.pkl` has `2 × 685` rows; observed AUROC ≈ 0.583 (gene) / ≈ 0.625 (archs4);
the results CSV reports a 95% CI and p-value on the AUROC difference; the histogram renders.
**The headline H4 answer is the `ci_excludes_0` / `p_value` on the AUROC-difference row.**

## p-value convention

`boot_pvalue_two_sided` uses the `(count+1)/(N+1)` convention (Davison & Hinkley 1997), so the
two-sided p-value is floored at `2/(N+1)` and never reads exactly `0`. The 95% CI on the difference
remains the primary inference.

## Generalizing to the full grid

The results table is structured as one row per method-pair comparison so it extends to the 4×3 grid
(gene/ARCHS4/GTEx/recount2). When you generalize: feed all four methods through the same aggregation,
loop the bootstrap over each method pair using the shared resample indices, and let the existing
`bh_correct` apply Benjamini-Hochberg across the full set of comparisons (it is a no-op for the single
comparison here).

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Read the **parent** `CLAUDE.md` (`../CLAUDE.md`, the drug-disease pipeline) and the **project-root**
`CLAUDE.md` first. This file only covers what is specific to this `null_adjust_test/` sandbox.

## What this directory is

An **isolated prototype** (a sibling of `signif_test/`) for a *different* question: when drug-disease
scores are built as `score = -1 * drugᵀ · disease`, can we **adjust each score against a null
distribution** so that scores become comparable across drugs/traits/methods and rankings become more
sensible? It is exploratory / feasibility — **not** integrated into NB06–13.

It is a **read-only consumer** of the inputs of NB00/01 (projections) and NB06/07 (gene-level data).
**Do not** modify NB06–13, `libs/drug_disease_utils.py`, `config.R`, or the scoring-sign convention
from here.

## Scope of this prototype (deliberately minimal)

- **Methods:** single-gene (NB06 logic) and module-based ARCHS4 (NB07 logic).
- **Tissue / threshold:** **Liver only**, **all-genes/all-LVs** (no top-N masking). The most
  lipid-relevant tissue; the simplest scoring case.
- **Target diseases:** the only two lipid DOIDs present in the gold standard:
  `DOID:1936` (atherosclerosis; 1 contributing UKB trait) and
  `DOID:3393` (coronary artery disease; 3 contributing UKB traits — the DOID score is a `max`
  over them, so the null aggregates the same way).

## The three nulls (what each tests)

For a disease vector `x_t` (genes or LVs) and the drug matrix `L`, the observed score column is
`s = -Lᵀ x_t`.

1. **Permute-disease (primary).** Permute the entries of `x_t` across genes/LVs (`B=1000`, seeded),
   recompute `s_b = -Lᵀ x_b`. Per (drug): `NES = (s-μ_b)/σ_b`, one-sided `p = (1+#{s_b≥s})/(B+1)`.
   Preserves the trait's magnitude distribution and the full drug matrix.
2. **Analytic Pearson.** Closed-form equivalent: `adj_pearson(d,t) = -corr(L[:,d], x_t)`. The
   permutation NES is *proportional* to this (≈ `sqrt(n-1)·adj_pearson`), so it is a no-simulation
   cross-check — NB01 asserts `corr(NES_permute, adj_pearson) > 0.95` at the trait level.
3. **Background-trait empirical null.** For drug d and lipid DOID, null = that drug's DOID-level
   scores across all *other* mapped DOIDs. Preserves gene-gene correlation (uses real trait
   vectors), unlike label permutation which breaks it.

4. **Z-shuffle null (module-only; NB03).** Shuffle the gene→LV assignment of CLAMP's loadings `Z`
   and re-project both LINCS and the disease data, then re-score — tests whether the *learned*
   modules carry signal vs. random gene groupings. See the efficiency note below.

Trait→DOID aggregation: adjusted scores are computed per contributing UKB trait, then aggregated to
the DOID by `max` (mirroring `map_traits_to_doid`).

## Notebooks (run in order, `clamp-analyses` env)

```
# Phase 1 — single tissue (Liver), all-genes, 2 lipid DOIDs
00_observed_scores.ipynb      →  observed_lipid_scores.pkl (+ background DOID matrix)
01_null_adjustment.ipynb      →  adjusted_lipid_scores.pkl (+ showcase null arrays)
02_compare_and_plots.ipynb    →  figures + FINDINGS.md
03_zshuffle_null_module.ipynb →  module_zshuffle_scores.pkl + FINDINGS_zshuffle.md
# Phase 2 — full scope: 49 tissues × 5 top-N thresholds, 685-pair universe (gene + ARCHS4)
04_multitissue_aggregate.ipynb → predictions_multitissue_aggregated.pkl
05_multitissue_eval.ipynb      → multitissue_bootstrap_results.csv + FINDINGS_multitissue.md
```

### NB04/05 — multi-tissue scope (Phase 2)

NB04 recomputes **raw**, **cosine**, **pearson**, and **background** scores per tissue × threshold
(reusing the NB10 rank→merge→mean→max aggregation), for gene-based + module ARCHS4 only. They form a
normalization ladder: **raw** keeps both vector norms and the mean; **cosine** `−(L_d·x)/(‖L_d‖‖x‖)`
removes only the norms (the magnitude); **pearson** `−corr(L_d, x_masked)` removes norms *and* the
mean (= cosine of centered vectors). The per-cell **permutation** null is infeasible at 49-tissue
scale, so the scalable adjustment is **analytic Pearson** (centered over all genes/LVs, zeros
included) — which is *exactly* the masked permute-disease NES; NB04 asserts this with a B=200
permutation spot-check, and asserts recomputed **raw** reproduces the published AUROCs (gene 0.583 /
ARCHS4 0.625) on the same 685-pair universe. NB05 runs the `signif_test` paired bootstrap on the
AUROC/AUPRC differences.

**Phase-2 result (load-bearing): removing the dot product's magnitude HURTS.** cosine, pearson, and
background all significantly *lower* AUROC (gene 0.583→0.527/0.526/0.485; ARCHS4 0.625→0.525/0.525/
0.518; all BH-sig) and AUPRC, and erase the ARCHS4-over-gene gap (raw +0.042 borderline →
cosine/pearson ≈ 0, n.s.). The raw dot product's magnitude carries real signal here; normalizing it
away removes it. Don't "fix" this by re-tuning — it is the finding. **cosine ≈ pearson** (diff
≈0.001/−0.0001, n.s.): the masked gene/LV vectors are effectively mean-zero, so the two
normalizations are interchangeable and cosine adds no new behavior beyond pearson — it is the clean,
literal test of NEXT_STEPS concern #2 ("is the signal just magnitude?"), answered: no, magnitude *is*
the signal. Z-shuffle (NB03) was **not** scaled to 49 tissues (needs B re-projections/tissue).

### NB03 efficiency / correctness note (load-bearing)

`projectCLAMP` is ridge regression `B = (ZᵀZ + L2·I)⁻¹ Zᵀ Y`, and **`ZᵀZ` is invariant under
permuting Z's gene rows**. So shuffling Z's rows by `τ` is *algebraically identical* to permuting
the gene rows of the data `Y` by `τ⁻¹` and reusing the real Z. NB03 precomputes the projector
`W = (ZᵀZ + L2·I)⁻¹ Zᵀ` once (via rpy2) and runs the whole null in numpy (`W · Y[σ]`, ~0.3 s/shuffle).
The same `σ` is applied to **both** LINCS and disease each iteration. NB03 asserts this fast path is
exact: numpy reprojection == `projectCLAMP` (~1e-14), one shuffle == R `projectCLAMP` on a
row-permuted Z (~1e-14), and identity scores == NB00. Do not "optimize" by dropping these checks.

- Trait→DOID mapping reuses `map_traits_to_doid` from `libs/drug_disease_utils.py` with the same
  three mapping files NB06/07 load. The scoring sign convention is copied, not re-derived.
- Inputs (read-only): LINCS gene-level `data/drug_disease_associations/lincs-data.pkl`; LINCS
  projection `output/.../01_lincs_projection_archs4/lincs/lincs-projection.pkl`; S-PrediXcan Liver
  raw `output/.../00_spredixcan_projection_archs4/spredixcan/raw/spredixcan-mashr-zscores-Liver-data.pkl`
  and projection `.../proj/spredixcan-mashr-zscores-Liver-projection-archs4.pkl`.
- Outputs (gitignored): `output/03_model_biology/00_archs4/02_drug_disease_associations/null_adjust_test/`.

## Deferred (TODO — not yet done)

- GTEx / recount2 module models (Phase 2 covered gene-based + ARCHS4 only).
- Z-shuffle null at 49-tissue scale, and beyond the 2 lipid DOIDs (NB03 is a Liver/lipid prototype).
- Full per-cell *permutation* null at scale (Pearson is its validated analytic stand-in).
- Integrating adjusted scores back into NB06–13 — currently unmotivated, since Phase 2 shows
  adjustment lowers AUROC.

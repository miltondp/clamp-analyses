# NEXT_STEPS.md

Methodological assessment of the drug→disease association pipeline and a
prioritized list of fixes / explorations. Written as a working memo so we don't
lose the reasoning. See `REVIEW.md` for the earlier write-up; this extends it.

## What's already sound (keep)

- **No label leakage in the representation.** LVs are learned unsupervised on
  RNA-seq (ARCHS4/GTEx/recount2), independent of PharmacotherapyDB. The latent
  space is not fit to the evaluation.
- **Scoring principle is well-grounded.** Negative dot product = "drug reverses
  disease signature" (Connectivity Map / LINCS reversal hypothesis);
  PharmacotherapyDB is the standard gold standard (Himmelstein/Hetionet lineage).
- **Single-gene baseline (NB 06) is the right control** — isolates whether the
  latent representation adds value beyond gene-level signal. Keep it central.
- **Three models in parallel give internal replication** (ARCHS4/GTEx/recount2).
- **Significance testing (REVIEW H4) is done and correct.** `signif_test/`
  implements a properly paired bootstrap (see "Previous concerns" below).

## What we're doing now: bootstrap ordering-stability

The paired bootstrap (Previous #1) finds **no** pairwise AUROC difference
significant at α=0.05. But a strict p<0.05 gate discards two real patterns in the
effect sizes, and we want to characterize them *without* dichotomizing:
1. **All three LV methods point positive vs the gene baseline** (Δ = +0.019 to
   +0.042) — consistent direction.
2. **AUROC is ordered by model size**: ARCHS4 (0.625) > recount2 (0.612) >
   GTEx (0.602) > gene (0.583), matching training-corpus size exactly.

**What we compute** (appended to `signif_test/01_paired_bootstrap_test.ipynb`,
reusing the already-computed, cross-method-aligned per-resample metric arrays —
no re-bootstrapping):
- **Pairwise dominance probabilities** `P(AUROC_A > AUROC_B)` across resamples
  (a threshold-free complement to the existing two-sided p-values).
- **Ordering-preservation fraction**: fraction of resamples where the full
  size-ordered chain holds (4-method `archs4>recount2>gtex>gene` and the LV-only
  `archs4>recount2>gtex`), plus the adjacent-pair decomposition and the
  "all three LV > gene simultaneously" fraction.
- **Spearman size-rank vs AUROC** per resample → observed value, bootstrap mean,
  95% CI, fraction == 1.0 (perfect concordance), fraction > 0.

Output: `output/.../signif_test/ordering_stability.csv` + a dominance-probability
heatmap `ordering_stability.png`.

**Interpretation guardrails** (stated in the notebook too, so the result isn't
over-read):
- The three models share the same gold-standard pairs and the same S-PrediXcan /
  LINCS input signatures — this is **correlated** corroboration, *not* three
  independent replications.
- "Size" is confounded with (a) dataset diversity/heterogeneity and (b) latent
  dimensionality (LV count is **not** fixed in `config.R`, so larger datasets
  likely yield more LVs / more capacity).
- n=3 models is weak for a trend (random ordering matches by chance with
  p≈1/6). This is **descriptive / hypothesis-generating**, not confirmatory; the
  definitive test is a within-dataset size gradient (subsample ARCHS4) — see
  "Worth exploring".

## Active concerns, ranked by impact on the conclusions

### 1. "max across 49 tissues" aggregation is a winner's-curse confound
Max of 49 values is an order statistic whose expectation grows with the number
of samples. Any pair with more available tissues / larger-variance scores gets a
systematically higher score for non-biological reasons.
- Check: does tissue coverage (or score variance) correlate with `true_class`?
- Asymmetry is unjustified: **mean** over `n_top_genes` but **max** over tissues.
  Ablate. Alternatives: per-tissue AUROC then aggregate; tissue-relevance-weighted
  combination.
- Same max-aggregation reappears in the **UKB→DOID collapse**
  (`map_traits_to_doid` keeps max when many UKB traits hit one DOID) — compounds.

### 2. Degree bias in the gold standard
Popular drugs / well-studied diseases have more curated indications; a predictor
tracking popularity scores well with zero biology. Add a **degree-only baseline**
(score = f(drug degree, disease degree) in the gold standard) and/or a
**degree-preserving permutation null** (Himmelstein 2017). If degree-only
approaches the model's AUROC, the latent model isn't adding biological signal.
Most likely "result is real but not for the claimed reason" failure mode.
Elevated in importance: now that the methods are statistically indistinguishable
on AUROC (see Previous #1), "does *any* method beat a degree null" is the live
question.

### 3. Raw dot product confounds direction with magnitude
Large-norm drug/disease vectors yield extreme dot products regardless of biology.
With `use_abs=True` on top, partly scoring "vector magnitude." Test **cosine
similarity** (normalize both vectors); if the gene→latent advantage disappears,
the signal was magnitude, not learned structure.

### 4. `use_abs=True` is in tension with the reversal narrative
Reversal is a *negative, sign-preserved* dot product. Taking absolute values at
top-LV selection discards the direction being claimed — it tests "shared loading
magnitude in the same LVs," not "opposite direction." Ablate `use_abs` True vs
False (one-line change, conceptual stakes).

### 5. Report the evaluation set after the inner join; confirm identical across methods
The inner join drops every pair lacking a LINCS drug signature *or* a UKB disease
trait — non-random dropout that tracks how well-studied a drug/disease is, which
tracks label. Report final N (pos/neg) out of 998, the DOID coverage of the
UKB→DOID mapping, and **assert** that 06/07/08/09 evaluate the exact same pair
universe (currently only a convention). NOTE: `signif_test/00_aggregate_predictions`
already enforces a shared `(trait, drug)` index across methods — partially
covers this; the reporting of N and DOID coverage is the remaining gap.

### 6. Soften "model X beats Y" prose in NB 10/11/13
Given the null significance result (Previous #1), audit NB 10/11/13 for any
sentence/figure that states or implies one method beats another. Reword to
"comparable, differences within bootstrap CI." (Not yet checked — flagged.)

## Worth exploring (higher-value science, not just fixes)

- **Compare latent spaces, not just latent-vs-gene.** The repo already builds
  PLIER, NMF/PCA, MOFA, etc. The compelling test is **CLAMP > PCA/NMF/PLIER latent
  space of equal dimension**, not just CLAMP > single gene. Shows CLAMP's
  *structure* matters, not merely dimensionality reduction. Low-hanging given
  existing models.
- **Switch headline metric to early-enrichment.** AUROC weights all thresholds
  equally; repurposing cares about the top of the list. With class imbalance
  (755/243 globally, smaller after join), report **AUPRC** and an early-enrichment
  metric (precision@k, EF, or BEDROC). (AUPRC is already in `signif_test`; extend
  to early-enrichment.)
- **Per-disease (or per-drug) AUROC, not just pooled.** Pooled AUROC can be
  dominated by a few easy diseases. Compute AUROC within each disease across its
  candidate drugs, then look at the distribution (+ a mixed-effects summary).
  Now more important: a per-disease view may reveal real separation that pooled
  AUROC (statistically null between methods) hides.
- **Interpret tissue selection as validation.** Record which tissue wins for each
  correct prediction. Mechanistically sensible winners (e.g., cardiac drugs →
  heart tissue) corroborate; random winners argue the max is noise-fishing (ties
  to Active #1).
- **Degree-preserving permutation null** (Active #2) doubles as the principled
  significance framework for the whole pipeline, not just pairwise AUROC diffs.

## Suggested order of attack

The controls that would most change confidence, in order:
1. **Degree-only baseline / permutation null** (Active #2) — now the central
   "is there any signal at all" test.
2. **Cosine vs dot product** (Active #3).
3. **Ablate max-over-tissues aggregation** (Active #1).
4. **Per-disease AUROC + early-enrichment metrics** (Worth exploring).

If a method survives the degree null and beats a matched PCA/NMF latent space,
the claim is defensible. Absent that, the honest headline is "latent and
gene-based methods are comparable; differences are within CI."

## Previous concerns (resolved / superseded)

### [RESOLVED] No significance test on AUROC differences (REVIEW H4)
**This was originally listed as the top priority.** It is already addressed by
`signif_test/01_paired_bootstrap_test.ipynb`, which is correctly designed:
- A **properly paired** bootstrap — one resample per iteration, *reused across
  all four methods*; pairwise diffs taken by subtraction (`boot_b - boot_a`),
  which cancels cross-method correlation (the part that matters for a difference
  test).
- Pivots on a shared `(trait, drug)` index, asserts identical pair coverage,
  reports CI + bootstrap p-value on each difference, BH-corrects within each
  metric family (AUROC and AUPRC separately).

**Finding (from `output/.../signif_test/paired_bootstrap_results.csv`): none of
the 6 pairwise AUROC differences are significant.** The headline ARCHS4-vs-gene
gap is +0.042, 95% CI [−0.0002, 0.084], raw p = 0.052, BH p = 0.31. All
module-vs-module comparisons clearly cross 0. So the ~4pp gap the narrative
leaned on is **within bootstrap noise** — exactly the worry H4 raised, now
confirmed.

**Why a two-way clustered bootstrap is NOT needed here.** A clustered bootstrap
(resampling drugs/diseases to account for pair non-independence) is *strictly
more conservative* — it can only **widen** these CIs. The current i.i.d.-pair
bootstrap is the one most prone to *over-claiming* significance, and it already
finds nothing significant. Clustering can only push an already-null result
further into null; it cannot overturn or rescue the conclusion. It would be, at
most, a low-value robustness footnote on a null. Skip it. (Clustering would only
matter if we had *found* significance and wanted to defend it against
non-independence — not the situation here.)

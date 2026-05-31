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
- **Significance test (H4) and ordering-stability are done.** `signif_test/`
  implements a properly paired bootstrap and a threshold-free ordering-stability
  analysis (see "Previous concerns"). A degree/popularity null was also explored
  and **set aside as not relevant** (see "Previous concerns").

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

### 2. Raw dot product confounds direction with magnitude
Large-norm drug/disease vectors yield extreme dot products regardless of biology.
With `use_abs=True` on top, partly scoring "vector magnitude." Test **cosine
similarity** (normalize both vectors); if the gene→latent advantage disappears,
the signal was magnitude, not learned structure.

### 3. `use_abs=True` is in tension with the reversal narrative
Reversal is a *negative, sign-preserved* dot product. Taking absolute values at
top-LV selection discards the direction being claimed — it tests "shared loading
magnitude in the same LVs," not "opposite direction." Ablate `use_abs` True vs
False (one-line change, conceptual stakes).

### 4. Report the evaluation set after the inner join; confirm identical across methods
The inner join drops every pair lacking a LINCS drug signature *or* a UKB disease
trait — non-random dropout that tracks how well-studied a drug/disease is, which
tracks label. Report final N (pos/neg) out of 998, the DOID coverage of the
UKB→DOID mapping, and **assert** that 06/07/08/09 evaluate the exact same pair
universe (currently only a convention). NOTE: `signif_test/00_aggregate_predictions`
already enforces a shared `(trait, drug)` index across methods — partially
covers this; the reporting of N and DOID coverage is the remaining gap.

### 5. Soften method-comparison prose in NB 10/11/13
Given the null significance result (Previous: H4), audit NB 10/11/13 for any
sentence/figure that states or implies one method *significantly* beats another.
Reword to "comparable; differences within bootstrap CI" (the consistent-but-not-
significant ordering is the honest framing — see Previous: ordering-stability).
(Not yet edited — flagged.)

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
- **Per-disease AUROC *and AUPRC*, not just pooled** — *planned next experiment.*
  Computing AUROC/AUPRC within each disease (across its candidate drugs) removes
  disease-level popularity *by construction* (the disease node is fixed) — this is
  the rigorous version of "restrict to popular diseases." Report the distribution
  across diseases (+ a mixed-effects / macro-mean summary), per method and vs a
  drug-popularity baseline. **Preliminary (quick check):** macro-mean per-disease
  AUROC archs4 ≈ 0.649 (>0.5 in 71% of diseases), gene ≈ 0.598 — i.e. signal
  survives once disease-popularity is gone. **Limitation:** only **17 of 57**
  diseases have ≥3 positives *and* ≥3 negatives (PharmacotherapyDB is ~78%
  positive), so the per-disease view is thin; AUPRC will be near the high base
  rate and must be read as Δ vs that base rate.
- **Interpret tissue selection as validation.** Record which tissue wins for each
  correct prediction. Mechanistically sensible winners (e.g., cardiac drugs →
  heart tissue) corroborate; random winners argue the max is noise-fishing (ties
  to Active #1).

## Suggested order of attack

The controls that would most change confidence, in order:
1. **Cosine vs dot product** (Active #2) — does a magnitude-free score change the
   gene→latent picture?
2. **Per-disease AUROC + AUPRC** (Worth exploring) — pooled AUROC may hide where,
   if anywhere, biology helps; removes disease-popularity by construction.
3. **Ablate max-over-tissues aggregation** (Active #1).
4. **Compare latent spaces (CLAMP vs PCA/NMF/PLIER)** (Worth exploring) — the test
   that CLAMP's *structure* matters, not just having a latent space.

The honest current headline is "latent and gene-based methods are comparable;
the consistent, size-ordered LV>gene advantage is suggestive but not significant."

## Previous concerns (resolved / superseded)

### [EXPLORED — SET ASIDE, not relevant] Degree / popularity null
We built a degree/popularity null (gold-standard degree-only baselines + an
incremental-value test). **The notebook and its outputs were removed; this entry
is the sole record.** It was set aside because, on investigation, the degree
confound is **not a relevant concern** for this pipeline:

- **The initial headline was misleading.** A label-derived popularity baseline
  scored *higher* in absolute AUROC (count-degree 0.66, prior-rate 0.79 vs methods
  0.58–0.62) — but that baseline peeks at the gold standard's own label graph,
  which the methods never see, so it is a *null model*, not a competing predictor.
- **Methods do NOT track *disease* popularity.** ρ(method score, disease-degree)
  ≈ 0.01–0.07; a popular disease does not get systematically higher scores. The
  "popular-disease" confound is empirically absent.
- **The cross-method ordering is popularity-independent.** gene and archs4 have
  ~equal *drug*-degree correlation (0.27 vs 0.26) but archs4 wins on *label*
  correlation (0.18 vs 0.12) → **archs4 > gene is real signal, not popularity**
  (the most degree-correlated method, gene, is the *worst* performer).
- **Per-disease holds** (disease fixed ⇒ disease-popularity removed): macro-mean
  per-disease AUROC archs4 ≈ 0.649, >0.5 in 71% of the 17 eligible diseases.

A residual *drug*-degree correlation (ρ≈0.26) inflates the *absolute* pooled
AUROC, **but MP does not consider this drug-degree aspect a relevant concern.**
The fair, popularity-controlled question is taken up properly by the planned
per-disease AUROC + AUPRC analysis (see "Worth exploring").

### [DONE] Bootstrap ordering-stability (beyond p-values)
Follow-up to the H4 null: instead of dichotomizing at p<0.05, characterize the
two effect-size patterns the gate discards — (1) all three LV methods point
positive vs gene, (2) AUROC is ordered by model size (ARCHS4 0.625 > recount2
0.612 > GTEx 0.602 > gene 0.583). Implemented as cells appended to
`signif_test/01_paired_bootstrap_test.ipynb`, reusing the cross-method-aligned
per-resample `boot_metric` arrays (no re-bootstrapping). Outputs:
`output/.../signif_test/ordering_stability.csv` + dominance heatmap
`ordering_stability.png`.

**Findings (AUROC, 10k resamples):**
- **Directional consistency is strong.** P(each LV > gene): ARCHS4 0.97,
  recount2 0.88, GTEx 0.80; all three LV > gene *simultaneously* in **0.75** of
  resamples.
- **Size gradient holds link-by-link.** P(archs4>recount2)=0.74,
  P(recount2>gtex)=0.70 — both clearly >0.5.
- **The strict 4-way chain is fragile.** Full ordering preserved in only **0.33**
  of resamples (= frac_rho==1, by construction). Spearman(size, AUROC): observed
  ρ=1.0, bootstrap mean 0.74, **frac ρ>0 = 0.97**, but 95% CI lower bound = 0.0
  (not bounded away from zero — consistent with n=3 being weak for a trend).

**Verdict:** supports the defensible framing "LV methods directionally and
consistently beat the gene baseline, and performance is size-ordered" — *without*
overclaiming, since no single difference clears p<0.05 and the strict ordering is
not resample-stable. Guardrails (correlated not independent replication; size
confounded with diversity + LV dimensionality; n=3 weak) are stated in the
notebook. **The definitive confirmation remains the within-dataset size-gradient
experiment (subsample ARCHS4)** — still listed under "Worth exploring".

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

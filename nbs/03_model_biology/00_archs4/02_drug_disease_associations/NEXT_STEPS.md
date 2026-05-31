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

### 1. `use_abs=True` is in tension with the reversal narrative
Reversal is a *negative, sign-preserved* dot product. Taking absolute values at
top-LV selection discards the direction being claimed — it tests "shared loading
magnitude in the same LVs," not "opposite direction." Ablate `use_abs` True vs
False (one-line change, conceptual stakes).

### 2. Report the evaluation set after the inner join; confirm identical across methods
The inner join drops every pair lacking a LINCS drug signature *or* a UKB disease
trait — non-random dropout that tracks how well-studied a drug/disease is, which
tracks label. Report final N (pos/neg) out of 998, the DOID coverage of the
UKB→DOID mapping, and **assert** that 06/07/08/09 evaluate the exact same pair
universe (currently only a convention). NOTE: `signif_test/00_aggregate_predictions`
already enforces a shared `(trait, drug)` index across methods — partially
covers this; the reporting of N and DOID coverage is the remaining gap.

### 3. Soften method-comparison prose in NB 10/11/13
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
- **Interpret tissue selection as validation.** Record which tissue wins for each
  correct prediction. Mechanistically sensible winners (e.g., cardiac drugs →
  heart tissue) corroborate; random winners argue the max is noise-fishing (ties
  to Active #1).

## Suggested order of attack

The controls that would most change confidence, in order:
1. **Compare latent spaces (CLAMP vs PCA/NMF/PLIER)** (Worth exploring) — the test
   that CLAMP's *structure* matters, not just having a latent space.


## Previous concerns (resolved / superseded)

### [DONE] "max across 49 tissues" aggregation (winner's-curse / asymmetry)
Implemented in `tissue_agg_test/` (NB00 per-tissue metrics → NB01 summary + tissue-cluster
bootstrap). Two reframings set the approach (do not relitigate):
- **The "order statistic grows with n" mechanism is absent.** Coverage is a *constant* 49 tissues
  per pair (asserted in `signif_test/00`), so "more tissues → higher max" and the proposed
  "coverage vs `true_class`" check are moot.
- **Mean/median-of-scores is not a clean ablation** — `max` encodes "one relevant tissue per pair";
  mean assumes broad sharing, a *different* model. The standalone variance-vs-label diagnostic was
  also **dropped**: across-tissue variance is mechanistically entangled with the legitimate
  one-relevant-tissue signal, so it cannot separate confound from signal.

So we used the memo's own suggested alternative: **per-tissue AUROC/AUPRC, then aggregate** (removes
per-pair tissue selection entirely), with a tissue-cluster paired bootstrap (resample the 49
tissues; `seed=42`, N=10k; BH within metric). AUPRC headline = `log2(AUPRC/base_rate)`.

**Findings:**
- **The max extracts genuine per-pair signal — it is not noise-fishing.** Max-aggregate AUROC sits
  well above mean-per-tissue for *every* method (selection gain +0.042 gene, +0.071 ARCHS4, +0.076
  GTEx, +0.065 recount2). If the max were a pure order-statistic artifact on noise, mean-per-tissue
  would be ~0.5 with no method ordering; instead single-tissue AUROC is 0.527–0.555 and the ordering
  is preserved. **Verdict: max is a legitimate, signal-extracting aggregation.** Notably the *LV*
  methods gain more from selection than the gene baseline (consistent with LVs capturing
  tissue-specific structure the max can exploit).
- **ARCHS4 ≥ gene is NOT a max artifact.** Per-tissue ordering ARCHS4 > recount2 > gene > GTEx;
  max-aggregate ordering ARCHS4 > recount2 > GTEx > gene. ARCHS4 tops both. ARCHS4 beats gene
  per-tissue (+0.0132; wins in **30/49** tissues; Wilcoxon p=0.002), i.e. the LV-vs-gene direction
  is present in the *average single tissue*, before any selection.
- **The one selection-manufactured effect is GTEx-vs-gene.** GTEx is *below* gene per-tissue
  (−0.0148) but *above* gene under max-aggregate — its apparent edge over the baseline comes entirely
  from tissue selection.
- **Significance is bootstrap-optimistic.** The tissue-cluster bootstrap flags ARCHS4>gene (BH≈0),
  GTEx<gene, GTEx<ARCHS4, recount2>GTEx as "significant," but the 49 tissues are **correlated**
  (shared 685 pairs, shared genes/LVs, GTEx inter-tissue correlation), so the bootstrap is
  anti-conservative. The robust, distribution-free claim is the consistent **direction** (30/49,
  Wilcoxon), not the p-value.

**Verdict:** the max-over-tissues step is defensible (it recovers real per-pair tissue-specific
signal, not a winner's curse), and the cross-method ordering — including ARCHS4 ≥ gene — is **robust
to removing it**. This *nuances* the H4 null: the pooled max-aggregate paired test was underpowered
(ARCHS4-vs-gene +0.042, BH 0.31), whereas the per-tissue decomposition shows ARCHS4 consistently
beats gene across tissues — though we do not over-claim significance given the correlated-tissue
caveat. **Guardrails:** anti-conservative cluster bootstrap; n=49 correlated clusters; methods share
inputs (not independent replication); mean-per-tissue and max-aggregate are *different estimands*.
**Out of scope (flagged):** the tissue-selection biological validation (analysis B — which tissue
wins per correct prediction; `per_tissue_scores.pkl` retains the tissue axis for it) and the second
max (UKB-trait→DOID collapse in `map_traits_to_doid`, which "compounds"). Outputs:
`output/.../tissue_agg_test/` (`per_tissue_metrics.csv`, `per_tissue_macro_summary.csv`,
`per_tissue_ordering.csv`, `per_tissue_bootstrap_results.csv`, `figures/`, `FINDINGS.md`).

### [DONE] Per-disease AUROC + AUPRC
Implemented in `per_disease_test/` (NB00 metrics → NB01 summary + disease-level
cluster bootstrap). Computes AUROC/AUPRC *within* each disease (across its candidate
drugs), removing **disease**-level popularity by construction. Reuses the aligned
685-pair × 4-method frame `signif_test/predictions_paired.pkl`. AUPRC headline =
`log2(AUPRC / base_rate)` (fold-enrichment over the per-disease prior; raw AUPRC just
reads the ~0.36–0.95 base rate). Significance via a **disease-level cluster bootstrap**
(resample *diseases* — the independent unit once the disease node is fixed — recompute
each method's **equal-weight** macro-mean on the same resampled set; paired across
methods, N=10k, seed 42; BH within each metric × eligibility family). No drug-popularity
baseline (dropped — disease popularity is already removed by construction; residual
*drug*-degree is not considered a relevant concern). A mixed-effects model was rejected
as over-parameterized for 17–33 heteroscedastic AUROC clusters.

**Disease universe:** of 57 diseases, **33** are AUROC-eligible (≥1 pos & ≥1 neg), **17**
strict (≥3 & ≥3); 24 are AUROC-undefined (19 all-positive, 5 all-negative). The undefined
diseases are the most label-imbalanced, so the per-disease view is conditioned on a
non-random subset (ties to Active #3).

**Findings (macro-mean per-disease AUROC):**
- **Reproduces the preliminary on the strict-17 subset:** ARCHS4 ≈ 0.649, gene ≈ 0.598,
  ARCHS4 >0.5 in 71%. Pooled AUROC also reproduces NB10 (gene 0.583 / archs4 0.625) from
  the same frame — the input is correct.
- **The ordering is threshold-sensitive.** On the full **33**-disease set the ordering
  *flips*: gene ≈ 0.658 ≥ ARCHS4 ≈ 0.636 ≈ GTEx 0.640 > recount2 0.622. The 16
  thinly-sampled diseases (1 pos or 1 neg) add noise that erases ARCHS4's edge.
- **Significance:** the **only** BH-significant differences are **ARCHS4 > GTEx on the
  strict-17 subset** — AUROC diff +0.084 (BH p=0.022) and log2-AUPRC-enrichment +0.104
  (BH p=0.006). **ARCHS4-vs-gene is directional but n.s.** (strict AUROC +0.051, BH 0.274;
  null on the full set). So the headline "LV beats single-gene" claim stays within noise,
  consistent with the H4 null; the one robust per-disease signal is *intra-LV*
  (ARCHS4 > GTEx among well-powered diseases), and it too is absent on the full 33-disease
  set.

**Verdict:** per-disease analysis **does not rescue** a significant LV-vs-gene advantage
once disease popularity is removed — it confirms the H4 / ordering-stability framing
(directional, consistent, not significant). The contribution is the *construction*
(disease popularity removed by construction + AUPRC enrichment + a correctly-clustered
bootstrap), not a new p-value. **Guardrails:** small per-disease n (median ~14, min 4)
makes individual AUROCs near-categorical — only distributions/macro-means are
interpretable; n=17–33 clusters is weak power (wide CIs); methods are correlated (shared
inputs), not independent replications. Outputs: `output/.../per_disease_test/`
(`per_disease_metrics.csv`, `per_disease_macro_summary.csv`,
`per_disease_bootstrap_results.csv`, `figures/`, `FINDINGS.md`).

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
- **Per-disease** (disease fixed ⇒ disease-popularity removed): macro-mean
  per-disease AUROC archs4 ≈ 0.649 on the strict-17 subset, >0.5 in 71%. NOTE: the
  full analysis (see [DONE] Per-disease AUROC + AUPRC) shows this is threshold-
  sensitive — on all 33 AUROC-eligible diseases the ordering flips (gene ≈ archs4)
  and archs4-vs-gene is not significant; the per-disease view confirms rather than
  overturns the H4 null.

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

### [RESOLVED] Raw dot product confounds direction with magnitude (cosine test)
The worry: large-norm drug/disease vectors yield extreme dot products regardless
of biology, so the score partly measures "vector magnitude." Proposed test: cosine
similarity (normalize both vectors); if the gene→latent advantage disappears, the
signal was magnitude, not learned structure.

**Done in `null_adjust_test/` (NB04/05, 49 tissues × 5 top-N, 685-pair universe,
paired bootstrap).** Added **cosine** (`−(L_d·x)/‖L_d‖‖x‖`, removes the norms and
nothing else) alongside raw / pearson / background. Result: removing the magnitude
*lowers* AUROC for both methods and erases the gap — cosine gene 0.527, ARCHS4
0.525 (vs raw 0.583 / 0.625; cosine-vs-raw BH-sig), ARCHS4-vs-gene under cosine
−0.003 (n.s.). And **cosine ≈ pearson** (diff ≈+0.001 / −0.0001, n.s.), so the
masked vectors are effectively mean-zero and the earlier Pearson result already
answered this. **Conclusion (reframed):** magnitude is *signal*, not a confound —
the dot-product norm encodes transcriptional response strength (Connectivity-Map
reading), and keeping it is the correct modeling choice. The concern's implicit
"if it disappears under normalization it was never learned structure" is a false
dichotomy. See `output/.../null_adjust_test/FINDINGS_multitissue.md`.

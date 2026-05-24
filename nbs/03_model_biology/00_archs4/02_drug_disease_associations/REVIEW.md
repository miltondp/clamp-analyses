# Code Review: `02_drug_disease_associations`

**Branch:** `drug_disease_associations-mp`
**Reviewed:** 2026-05-24
**Scope:** All 14 notebooks under `nbs/03_model_biology/00_archs4/02_drug_disease_associations/`. The supporting Python module `drug_disease_utils` is **not** reviewed (see issue H1 below).

## Pipeline overview

```
00 spredixcan→archs4 ┐                            ┌─ 07 module pred (archs4) ┐
02 spredixcan→gtex   ├─→ projections (LVs)        ├─ 08 module pred (gtex)   ├─→ 10 perf → 11 plots
04 spredixcan→recount2┘                           └─ 09 module pred (recount2)┘     │
01/03/05 lincs→{archs4,gtex,recount2}                                               │
                                                  06 single-gene baseline ──────────┘
                                                  12 LV-group computation → 13 barplots
```

All scoring (gene-based and module-based) uses `score = -1 * drug^T disease` ("`predict_dotprod_neg`"); negative dot product reflects the standard "drug reverses disease signature" convention. The gold standard is **PharmacotherapyDB** (998 pairs, 755 positive / 243 negative, 87 DOIDs).

## Verdict

The pipeline is **internally consistent on the science**: projection direction is correct (`Y ≈ Z·B`, projection returns `LVs × traits`), the scoring metric is consistent across the four prediction notebooks, and the gold-standard universe is filtered consistently. The same pattern is replicated cleanly across ARCHS4 / GTEx / recount2.

The notable problems are **reproducibility-level** (a missing module that breaks end-to-end re-execution, mislabeled output paths, hardcoded paths bypassing `config.R`) and **statistical rigor** (cross-method AUROC differences are reported without any test or paired bootstrap). One genuine logic bug: NB 07 is missing the LV-index alignment assertion that NB 08 and 09 have.

---

## HIGH severity

### H1. `drug_disease_utils` is missing from the repo

NB 06, 07, 08, 09, 12 all do:

```python
sys.path.insert(0, str(here('libs')))
from drug_disease_utils import map_traits_to_doid, _zero_nontop_genes, predict_dotprod_neg
```

But:
- `libs/` does **not** exist in the working tree.
- `libs/drug_disease_utils.py` is **not** tracked in any branch in the repo (`git ls-files`, `git log --all -- libs/drug_disease_utils.py` both empty).
- A *different*-named file `libs/drug_disease.py` (240 lines) exists only on `upstream/bp_coverage` (commit `77d701c`), authored by `msubirana`.

Consequence: the entire prediction half of the pipeline (06–13) cannot be re-executed from a fresh clone. It also means the most safety-critical code — the scoring formula, the top-N gene masking (`_zero_nontop_genes`), the trait→DOID mapping — has **not** been reviewed because it is not in this repo.

**Action:** locate `drug_disease_utils.py` on a collaborator's machine and commit it to `libs/`. Add `libs/` to git tracking (it is not in `.gitignore`). Pin the module's behavior with at least one notebook-level assertion (e.g., a known fixture-based score).

### H2. Output paths label all three datasets as `00_archs4`

For projections that run on **GTEx** and **recount2** models, the notebook output directories still live under `output/03_model_biology/00_archs4/02_drug_disease_associations/...`:

- `03_lincs_projection_gtex.ipynb` → writes to `…/00_archs4/02_drug_disease_associations/03_lincs_projection_gtex/`
- `05_lincs_projection_recount2.ipynb` → writes to `…/00_archs4/02_drug_disease_associations/05_lincs_projection_recount2/`
- Same misnaming applies to `02_spredixcan_projection_gtex.ipynb` and `04_spredixcan_projection_recount2.ipynb`.

This is misleading and a foot-gun: a future grep for "gtex projection artifacts" under `02_gtex/` returns nothing. The notebook *source* lives under `00_archs4/`, but the *outputs* are dataset-specific and should be sorted by dataset.

**Action:** either (a) move the notebooks themselves under per-dataset subdirs, or (b) decouple the output dir from the source dir and key it on dataset (`config$<DATASET>$DATASET_NAME`).

### H3. NB 07 is missing the LV-index alignment assertion

NB 08 and 09 both have:

```python
assert tissue_proj.index.equals(lincs_projection.index)
```

NB 07 (ARCHS4) does not. If LINCS-archs4 and S-PrediXcan-archs4 projections ever fall out of LV-order (e.g., a re-run with different `MULTICLAMP` settings), 07 will silently compute dot products against mismatched LVs and produce **wrong scores** without erroring.

**Action:** add the assertion to NB 07 cell 20. One-line fix.

### H4. No statistical test on cross-method AUROC differences

NB 10 reports aggregate AUROC: gene-based 0.583, ARCHS4 0.625, GTEx 0.602, recount2 0.612. NB 12 adds bootstrap CIs (300 reps) for the LV-group analysis, but **none** of 10/11/13 perform a paired test across methods, and no multiple-testing correction is applied across the 4-way method comparison or the 8-cell LV-group × threshold grid.

The 4-pp AUROC gap between ARCHS4-module and gene-based could be within bootstrap CI overlap. **The headline claim of the analysis is unsupported by a significance test.**

**Action:** add paired bootstrap (resample (drug, disease) pairs, recompute AUROC for each method on the same resample, report 95% CI on the *difference*). Apply BH correction across the 4×3 method-pair grid.

---

## MEDIUM severity

### M1. Hardcoded data paths bypass `config.R`

`config.R` (54 lines total) defines `config$ARCHS4/GTEx/recount2` blocks but has **no** entries for the downstream data this directory consumes:

| Data | Hardcoded in | Suggested config key |
|---|---|---|
| S-PrediXcan TWAS results (49 tissues) | `00_spredixcan_projection_archs4.ipynb` | `config$PHENOMEXCAN$SPREDIXCAN_DIR` |
| LINCS L1000 signatures | `01_lincs_projection_archs4.ipynb`, 03, 05 | `config$LINCS$DATA_FILE` |
| PharmacotherapyDB gold standard | 06, 07, 08, 09 | `config$GOLD_STANDARD$PHARMACOTHERAPYDB` |
| CLAMP model `.rds` artifacts | every projection NB, NB 12 | `config$ARCHS4$MODEL_FILES$<key>` |

NB 12 has the most brittle one — a six-segment hardcoded path to one specific model under `output/01_model_building/04_archs4/06_bp_coverage_rshall/.../CLAMPfull_hall.rds`. If model selection changes, this notebook silently uses the wrong model.

**Action:** extend `config.R` with the new blocks; every notebook here should reference `config$...` not raw paths.

### M2. Asymmetric duplicate-gene handling between archs4 and gtex/recount2

NB 00 drops duplicate gene rows with `data[~data.index.duplicated(keep='first')]` *before* writing the raw output. NB 02 and 04 read that already-deduplicated raw file. Result: the dedup logic for GTEx and recount2 silently inherits a choice made for ARCHS4. If that choice is ever changed to `keep='last'`, downstream results for *all three datasets* drift in lockstep without a corresponding notebook change.

**Action:** deduplicate per-dataset, in each notebook, not by chain. Or, make the dedup step explicit in 02/04 (idempotent) with a comment.

### M3. Silent NaN row drop in NB 00

NB 00 uses `data.dropna(how='any')` against the raw S-PrediXcan stack. Number of dropped rows is not logged. If a future S-PrediXcan release has more NaNs (e.g., new tissues with sparser coverage), gene count drops silently.

**Action:** log `n_dropped` and `n_remaining`; ideally `assert n_dropped < threshold`.

### M4. NB 13 barplot helper sorts ascending so "biggest is at top"

`2df8fad0` in `13_drug_disease_barplots.ipynb` reverses the standard descending convention. Combined with M5 (no significance), the visual order can be misread as "winner is highest bar" without the eye realizing the bar order is inverted relative to most barplots.

**Action:** either use the standard descending sort, or annotate the axis explicitly with rank order.

### M5. Counterintuitive LV-group finding is reported without caveat

NB 12 reports that the "Neither" LV group (no pathway or trait significance) outperforms "Pathway + Trait" (AUROC 0.630 vs 0.616) under the strict threshold. This is potentially interesting science, but it should not be presented as fact without:
- a sample-size disclosure (how many LVs in each group?)
- a permutation/bootstrap test
- a discussion of whether "Neither" includes the bulk of LVs and benefits from sheer dimensionality.

**Action:** add group sizes to the table; report a permutation-test p-value for "Neither vs Pathway+Trait"; if the difference is within CI, label it null.

---

## LOW severity

### L1. CLAMP package and `org.Hs.eg.db` versions not pinned

The projection notebooks rely on `clusterProfiler.bitr()` for SYMBOL→ENSEMBL mapping. Bioconductor annotation DBs change releases; the same input could yield different Ensembl IDs on a future env rebuild. Similarly `CLAMP.projectCLAMP` is called without recording the CLAMP package commit / version. Reproducibility is contingent on conda env state, which is not snapshotted.

**Action:** add a Session-Info cell at the bottom of each notebook (`sessionInfo()` in R, `pip freeze` summary in Python). Commit `envs/clamp-analyses.lock.yaml` (a frozen lockfile from `conda env export`).

### L2. No tissue-file count assertion in 06–09

NB 06–09 build `spredixcan_file_list` via `glob()` and loop over it. There's no `assert len(spredixcan_file_list) == 49`. NB 10 catches downstream coverage gaps, but earlier failure detection costs one line.

**Action:** add `assert len(spredixcan_file_list) == 49, f"Expected 49 tissues, got {len(...)}"`.

### L3. NB 06 logs per-tissue intersection size; 07–09 do not

NB 06 prints `(N common genes)` per tissue. NB 07–09 do not log LV-overlap stats. Cheap parity fix; helps spot upstream drift.

### L4. NB 10 tissue-name regex is over-engineered

Cell `4f069665` in NB 10 tries 5 suffix patterns to extract a tissue name. Likely accreted from successive S-PrediXcan releases. Worth refactoring to a single canonical regex + `assert match is not None`.

### L5. NB 13 loops over `THRESH_LABELS` twice in cell `33553b07`

Minor inefficiency; harmless.

### L6. NB 12 random baseline uses `np.random.seed(42)` declared once at the top

The seed is set once but the random baseline is recomputed in several cells. As long as cells execute in order this is fine; out-of-order execution will produce a different baseline. Tighten by setting the seed immediately before each random-draw cell, or by using a `np.random.default_rng(42)` instance per cell.

---

## What's working well

- **Projection direction is correct** across all 6 projection notebooks. The R call `CLAMP.projectCLAMP(clamp_sub, newdata=r_mat)` consistently returns an `LVs × traits` matrix and the pre/post shapes are consistent with `Y ≈ Z·B`.
- **Scoring metric is consistent** across 06/07/08/09 — same `predict_dotprod_neg` import, same `use_abs=True`, same gold-standard filter, same HDF5 output schema (`metadata` + `prediction` keys).
- **NB 10 has a coverage assertion** (49 tissues × 5 thresholds × 685 pairs = 167,825 rows per method) that would catch most upstream coverage drops.
- **The pipeline is *parallel*-clean across datasets**: 00↔02↔04, 01↔03↔05, 07↔08↔09. Discrepancies are limited to the items in M2 and H3.
- **Reasonable gene-mapping hygiene**: notebooks filter to 1:1 unambiguous SYMBOL↔ENSEMBL mappings (`~dup_symbols & ~dup_ensembl`) consistently.

---

## Recommended fix order

1. **H1** — recover and commit `libs/drug_disease_utils.py`. Without this, the entire downstream pipeline is dead code from any reviewer's perspective.
2. **H3** — one-line assertion fix in NB 07.
3. **H4** — add paired bootstrap CI on cross-method AUROC differences; this is the chart the result hinges on.
4. **H2** — fix output dir mislabeling for GTEx and recount2 (move outputs to `02_gtex/` and `03_recount2/` parents).
5. **M1** — extend `config.R` and migrate hardcoded paths.
6. **M5** — gate the "Neither > Pathway+Trait" finding behind a permutation test before it goes into any figure caption.
7. **M2, M3, L1–L6** — incremental hygiene.

---

## Open questions for the author

- Where does `drug_disease_utils.py` live? (H1)
- Was the cross-method AUROC comparison intentionally left without a significance test, or was it deferred?
- Is the "ARCHS4 module-based" model the intended winner of the headline figure (NB 13)? If so, the gap to single-gene baseline (0.625 vs 0.583) is small and should be defended statistically.
- For NB 12, what's the rationale for the binary AUC/FDR cutoffs (0.7/0.05 strict, 0.6/0.1 loose)? Is there a sensitivity analysis?

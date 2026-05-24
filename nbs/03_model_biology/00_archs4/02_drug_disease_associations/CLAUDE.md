# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Read the **project-root** `CLAUDE.md` first for repo-wide conventions (envs, `config.R`, notebook numbering rule). This file only covers what is specific to this directory.

## What this directory is

A drug→disease association pipeline that projects two external signature sources (S-PrediXcan TWAS results across 49 tissues; LINCS L1000 perturbation signatures) into CLAMP latent space, then scores candidate (drug, disease) pairs against the **PharmacotherapyDB** gold standard (998 pairs: 755 positive / 243 negative, 87 DOIDs).

Despite living under `00_archs4/`, the pipeline evaluates **three** CLAMP models in parallel: ARCHS4, GTEx, recount2. The directory name is the source-tree location only; outputs are per-dataset (see issue H2 below).

All 14 notebooks use the `clamp-analyses` env — they are Python, not R (R is only used inside cells that wrap `CLAMP.projectCLAMP` via `rpy2`).

## Pipeline shape

```
00 spredixcan→archs4 ┐                            ┌─ 07 module pred (archs4) ┐
02 spredixcan→gtex   ├─→ projections (LVs)        ├─ 08 module pred (gtex)   ├─→ 10 perf → 11 plots
04 spredixcan→recount2┘                           └─ 09 module pred (recount2)┘     │
01/03/05 lincs→{archs4,gtex,recount2}                                               │
                                                  06 single-gene baseline ──────────┘
                                                  12 LV-group computation → 13 barplots
```

Three structural parallels are load-bearing — preserve them when editing:
- **Projection trio**: `00↔02↔04` (S-PrediXcan), `01↔03↔05` (LINCS). Same logic, three datasets.
- **Prediction trio**: `07↔08↔09`. Same scoring code, three projection sources.
- **Gold-standard universe** must be filtered identically across `06/07/08/09`, or `10` cross-method AUROC is comparing apples to oranges.

## Scoring convention

All four prediction notebooks (06–09) use `predict_dotprod_neg` with `use_abs=True`:

```python
score = -1 * (drug_vec @ disease_vec)   # negative dot product
```

The minus sign encodes the standard "drug *reverses* disease signature → high score" convention. Don't flip the sign without auditing all four notebooks and `10/11/13`.

Output schema for prediction HDF5s is the same across 06–09: `metadata` key + `prediction` key. Keep this — `10` reads all four into a single comparison frame.

## Known issues (read before extending)

These are documented in full in `REVIEW.md` (commit `dd0e59b`). The ones that affect day-to-day work:

- **`libs/drug_disease_utils.py` is missing from the repo** (REVIEW H1). NB `06/07/08/09/12` all do `from drug_disease_utils import ...` against a `libs/` dir that does **not** exist in the working tree and is not tracked in git. A fresh clone cannot re-execute the prediction half of the pipeline. If you need to run these notebooks, the file must be recovered from a collaborator before anything downstream of projection will work.
- **Output paths mislabel GTEx and recount2 as `00_archs4`** (REVIEW H2). `02/03/04/05` write to `output/03_model_biology/00_archs4/02_drug_disease_associations/...` even though they operate on GTEx/recount2 models. Don't grep for outputs by dataset name — grep by notebook number.
- **NB 07 is missing the LV-index alignment assertion** that NB 08/09 have (REVIEW H3). If you touch any of `07/08/09`, keep `assert tissue_proj.index.equals(lincs_projection.index)` symmetric across the three.
- **Hardcoded paths bypass `config.R`** (REVIEW M1) for S-PrediXcan results, LINCS signatures, PharmacotherapyDB, and the CLAMP `.rds` model files. NB 12 in particular hardcodes a six-segment path to one specific `CLAMPfull_hall.rds` — if model selection changes, NB 12 silently uses the wrong model.

## Conventions specific to this directory

- Gene mapping is **SYMBOL→ENSEMBL via `clusterProfiler.bitr`** (not biomaRt). Filter to 1:1 unambiguous mappings (`~dup_symbols & ~dup_ensembl`) — done consistently across projection NBs; keep it that way.
- NB 00 deduplicates raw S-PrediXcan input with `keep='first'`; NBs `02/04` inherit that file. If you change the dedup policy, change it in all three (REVIEW M2).
- Statistical tests on cross-method AUROC differences (NB 10/11/13) are **not** currently in the pipeline. If you publish a figure that compares ARCHS4 vs GTEx vs recount2 vs gene-baseline AUROC, add a paired bootstrap first (REVIEW H4).

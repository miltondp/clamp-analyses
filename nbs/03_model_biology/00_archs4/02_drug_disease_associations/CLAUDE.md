# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Read the **project-root** `CLAUDE.md` first for repo-wide conventions (envs, `config.R`, notebook numbering rule). This file only covers what is specific to this directory.

## Running commands

If you need to run python or R commands, use the `clamp-analyses` conda environment (e.g. `conda run -n clamp-analyses <cmd>`). Its setup and the pinned CLAMP install are documented in the project-root `README.md`.

## What this directory is

A drug→disease association pipeline that projects two external signature sources (S-PrediXcan TWAS results across 49 tissues; LINCS L1000 perturbation signatures) into CLAMP latent space, then scores candidate (drug, disease) pairs against the **PharmacotherapyDB** gold standard (998 pairs: 755 positive / 243 negative, 87 DOIDs).

Despite living under `00_archs4/`, the pipeline evaluates **three** CLAMP models in parallel: ARCHS4, GTEx, recount2. The directory name is the source-tree location only; outputs are written per-notebook (grep outputs by notebook number, not dataset name).

All 14 notebooks use the `clamp-analyses` env — they are Python, not R (R is only used inside cells that wrap `CLAMP.projectCLAMP` via `rpy2`).

## Pipeline shape

```
00 spredixcan→archs4  ┐                           ┌─ 07 module pred (archs4) ┐
02 spredixcan→gtex    ├─→ projections (LVs)       ├─ 08 module pred (gtex)   ├─→ 10 perf → 11 plots
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

## Open issues to focus on (read before extending)

Full write-ups in `REVIEW.md` (commit `dd0e59b`). In priority order:

1. **No significance test on cross-method AUROC differences** (REVIEW H4) — **top priority.** NB 10/11/13 report aggregate AUROC (gene-based ~0.583, ARCHS4 ~0.625, GTEx ~0.602, recount2 ~0.612) with **no paired test**. The ~4-pp gap that the analysis hinges on could be within bootstrap-CI overlap, so any "model X beats Y" claim is currently unsupported. Before publishing any cross-method figure, add a **paired bootstrap** (resample (drug, disease) pairs, recompute each method's AUROC on the *same* resample, report a 95% CI on the *difference*) and BH-correct across the method-pair grid.
2. **Hardcoded paths bypass `config.R`** (REVIEW M1) for S-PrediXcan results, LINCS signatures, PharmacotherapyDB, and the CLAMP `.rds` model files. NB 12 in particular hardcodes a six-segment path to one specific `CLAMPfull_hall.rds` — if model selection changes, NB 12 silently uses the wrong model.

## Conventions specific to this directory

- Gene mapping is **SYMBOL→ENSEMBL via `clusterProfiler.bitr`** (not biomaRt). Filter to 1:1 unambiguous mappings (`~dup_symbols & ~dup_ensembl`) — done consistently across projection NBs; keep it that way.

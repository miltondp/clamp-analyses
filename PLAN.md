# PLAN: evaluating and improving consensus LVs (metaLVs / MGMs)

This is a living plan on branch `consensus-lvs`. Work through it one step at a time, across sessions.

## How to resume

1. Read this file. Find the first unchecked box in **Status**, and read the **Log**.
2. Read **Key facts** before touching rules, configs or `cluster_metalv.py`.
3. Run `snakemake --profile workflow/profiles/local --dry-run <target>` before any real run. The dry run must list only `cev_*` jobs: no CLAMP fits, no `technical_noise_*`, no `crispercas_*`.
4. When a step is done, tick its box and add a dated line to the **Log**.

## Status

**Blocked on:** the rsync from pivlab01 (Phase 0, step 3). Setup can continue without it, but Phase 1 cannot.

### Phase 0: Setup and data staging (local)
- [ ] 0.1 Create the missing envs (`snakemake`, `gpu-kmeans`), register the kernel, and run `scripts/check_setup.py` (it also runs both `--check` pin scripts).
- [ ] 0.2 Symlink `data/clamp_archs4` to the `ai-comp-exp` fits, and verify the md5s against `MD5SUM.txt`.
- [ ] 0.3 Rsync the noise outputs and the legacy `00_mgm` outputs from pivlab01 (paths below).
- [ ] 0.4 Export Z.csv, summary.csv and manifest.json from the 10 `.rds` files (`cev_export_run`).
- [ ] 0.5 Build the GO BP gmt for ORA if it is missing (`generate_archs4_go_bp_gmt`).

### Phase 1: Reproduce the existing consensus sets (stop if this fails)
- [ ] 1.1 Scaffold the new area: config, `.smk`, Snakefile include and target, staging rules.
- [ ] 1.2 Check that the exported `Z_all.npz` equals the legacy `Z_all.npz` (`cev_zcache_check`).
- [ ] 1.3 Check that `cl_c050_m4_r015` equals `Z_4reps.csv` exactly.
- [ ] 1.4 Check that `cl_c050_m6_r025` equals `Z_6reps.csv`. If not, try the old 6-rep members or flag variants.
- [ ] 1.5 Notebook `00_reproduction`.

### Phase 2: Evaluation framework, applied to the existing sets
- [ ] 2.1 `cev.py` metrics and unit tests (toy data).
- [ ] 2.2 Computational metrics on Z_4reps, Z_6reps and the hierarchical grid (`cev_evaluate`).
- [ ] 2.3 Biological metrics: ORA, member pathway coherence, technical leakage.
- [ ] 2.4 Robustness: split-half and LORO for the hierarchical sets, plus both nulls.
- [ ] 2.5 Semi-synthetic benchmark.
- [ ] 2.6 Notebooks `02_computational` and `03_biological`.

### Phase 3: Ensemble (consensus-index) clustering
- [ ] 3.1 k-means ensemble on the GPU (`km_full`), k-means only.
- [ ] 3.2 Co-association, EAC average linkage, and agreement-vs-k with k selection.
- [ ] 3.3 Build ensemble sets with replication filter, trim and the five combination functions.
- [ ] 3.4 Run the same evaluation on the ensemble sets, including split-half, LORO and synthetic.
- [ ] 3.5 Variants: add hierarchical cuts to the ensemble, PCA-256, and (later) spectral consensus.
- [ ] 3.6 Notebook `01_ensemble`.

### Phase 4: Comparison and recommendation
- [ ] 4.1 Notebook `04_summary`: trade-off curves, size-matched comparisons, recommended set.
- [ ] 4.2 (Later) Connect the chosen set to the `consensus.smk` downstream steps (projection, LINCS, GLS).

### Phase 5: GTEx extension (if 5-vs-5 split-half on ARCHS4 is too noisy)
- [ ] 5.1 Stage the GTEx inputs: TPM gct, FBM, canonical S124 and the existing 10 runs.
- [ ] 5.2 Add a per-compendium `compendia.gtex.extra_seeds: [1008..1017]` in `consensus.yaml` (not global, so recount2's `CON_SEED_PATTERN` is unchanged), and make `con_seed_index` continue past `len(seeds)`.
- [ ] 5.3 Fit the 10 extra runs and their noise screens with the existing `consensus_seed_full` / `consensus_model_noise` rules (48 GB, 8 threads each, two at a time). No `cev_export_run` for GTEx: consensus.smk already writes Z.csv, summary.csv and manifest.json.
- [ ] 5.4 Add `compendia.gtex` to `consensus_eval.yaml`, deriving its 20 models from `con_models("gtex")` plus the extras and its paths from `con_model_dir` / `con_noise_dir`, then run the evaluation with a 10-vs-10 split-half.

## Log

- 2026-11-26: Explored the repo, the PhenoPLIER ensemble code and the ai-comp-exp data, then wrote this plan. Nothing executed yet.
- 2026-11-26: Simplify pass on the plan: derive models/paths/top_n/ORA settings from existing config; one `{subset}` wildcard replaces the split-half/LORO rule variants; cached linkages, uint8 co-association, deduplicated ORA, `Z.npy` storage, `min_frac` replication, manifest.json export, GTEx via the existing consensus.smk chain.

---

## Context

### Current method

`scripts/archs4/crispercas/cluster_metalv.py`, with rules in `workflow/rules/archs4_crispercas.smk:43-130` (ARCHS4) and `workflow/rules/consensus.smk` (GTEx, recount2, CLAMPbase). Config is in `workflow/config/consensus.yaml` and in `archs4.yaml`'s `crispercas:` section; the two copies must be kept in sync by hand.

1. **`cache`**: stack each run's Z.csv (18,423 genes × about 1,728 LVs per run, 17,280 LVs in total) into `Z_all.npz`.
2. **`noise` / `flags`**: flag an LV as technical if either of these holds:
   - `gene_flag`: Pearson correlation of loadings with log10 cDNA length or GC%, with |r| ≥ 0.2, BH FDR < 0.2 and |robust z| ≥ 4. The thresholds are hard-coded in `TH` in `scripts/archs4/technical_noise/run_analysis.R`.
   - `sample_flag`: 5-fold CV R² of the LV's B row on the technical covariates ≥ `noise_r2`, currently 0.15.
3. **`cluster`**: drop flagged LVs and L2-normalize the rest. Run complete linkage on 1−cos and cut at 1−`min_cos` (0.5), so every pair inside a cluster has cos ≥ 0.5. Keep clusters spanning ≥ `min_models` distinct runs (4). The consensus is the gene-wise median of the members' raw loadings. It is not renormalized, and a run that contributes several members counts several times.
4. Evaluation that already exists:
   - Within-cluster cosine (`cluster_summary.csv`).
   - ORA coverage (`run_ora.R` on the top 1% of genes, i.e. 185).
   - supp6: share of run modules that land in an MGM, CLAMPfull vs CLAMPbase.
   - Downstream: LINCS and PhenoPLIER GLS.
   - There is no top-gene Jaccard, no reproducibility analysis and no null.

### Existing consensus sets to compare against

These are in `/home/miltondp/projects/pivlab/ai-comp-exp/data/clamp_archs4/`.

| File | Contents |
|---|---|
| `CLAMPfull_canonical_seed{123,456,789,1001..1007}.rds` | The 10 fits, about 7.6 GB each. `MD5SUM.txt` is alongside. |
| `consensus/Z_4reps.csv` | 18,423 × 1,058. `min_models` 4, `noise_r2` 0.15 (current config). |
| `consensus/Z_6reps.csv` | 18,423 × 618. `min_models` 6 with the **old** `noise_r2` 0.25. Byte-identical to Marc's published set, whose `manifest.json` is in `clamp-gtex-thyroid`. |

### Goal

1. Build a way to measure, computationally and biologically, how well any method combines LVs, and apply it identically to every set.
2. Try consensus-index (ensemble) clustering as in PhenoPLIER: a k-means ensemble over many k and seeds, then a co-association matrix, then a consensus function, then choosing k by ARI agreement with the ensemble.
3. Compare the new approach against Z_4reps and Z_6reps.
4. If 10 ARCHS4 runs aren't enough, add GTEx runs.

### Decisions made (2026-11-25)

- Noise flags are rsynced from pivlab01, not recomputed. Recomputing would need ARCHS4 preprocessing at about 200 GB RAM.
- Everything runs locally, including any extra GTEx fits.
- The work goes in a new Snakemake area, `consensus_eval` (rule prefix `cev_`), following the CLAUDE.md conventions from the start.

---

## Key facts (found during exploration)

- **Z_6reps is not a subset of Z_4reps.**
  - Its column names match the first 618 of Z_4reps, but its values differ: 87 of 618 columns have no Z_4reps match at cos ≥ 0.99, and 36 have none at ≥ 0.9.
  - The cause is commit `3085ab6`, which changed the config from `min_models` 6 / `noise_r2` 0.25 to 4 / 0.15.
  - Reproducing Z_6reps therefore needs flags at 0.25 and `min_models` 6.
- **Assertion trap.** `cluster_metalv.py noise_main` asserts that the S123 flags at `--threshold` equal `lv_technical_classification.csv`, but that file was built with `R2_CUT <- 0.25` (`run_analysis.R:40`). At 0.15 the assertion can fail. Use the `flags` subcommand (`--model NAME=dir`), which has no assertion.
- **Do not edit `cluster_metalv.py`.** Nine rules list it as an input without `ancient()` (`archs4_crispercas.smk:79,95,109`; `consensus.smk:272,291,309,658,673,691`), so any edit reruns the existing consensus chain (ORA, GLS, LINCS) for every compendium. `rerun-triggers: mtime` doesn't help, because an edit changes the file's mtime. Call its subcommands, or load it via importlib like `tests/test_crispr_l1000.py` does. Possible follow-up: wrap those 9 inputs in `ancient()` (precedent: `archs4_crispercas.smk:50`), then factor `cluster_main` into importable functions.
- **`technical_noise_z` takes `Z.csv` without `ancient()`.** Never write a Z.csv into `output/98_final_models/clampfull/canonical/archs4/`: that triggers a rerun which needs the raw h5.
- **Model completeness gating** (`archs4_canonical.smk:12-34`).
  - A model counts as complete only when all of these exist and are non-empty: rds, B.csv, Z.csv, summary.csv, manifest.json, complete.
  - A partly filled folder raises `ValueError` while the job graph is being built.
  - The new area therefore stages its own files under its own root and never fills those folders.
- **What the `.rds` contains.** A list with `B, Z, U, C, L1, L2, Z2, heldOutGenes, Uauc, Up, summary, priorMat, withPrior, call`. `fit_clampfull.R:179-215` writes `write.csv(Z)` with gene rownames, `write.csv(B)`, and `write.csv(summary, row.names=FALSE)` with `LV_index` renamed to `LV`. Exporting the same way reproduces those files exactly. B is about 1,728 × 580k, about 18 GB as CSV. It isn't needed, so skip it.
- **Legacy model order:** `S123, S1001..S1007, S456, S789` (`A4_METALV_SEEDS` is sorted with key `(s<1000, s)`). The order decides linkage tie-breaks and metaLV numbering, so it is required for exact reproduction.
- **ORA** queries the top `ceil(18423 × 0.01)` = **185** positive genes. `run_ora.R` takes `--model-manifest`, like the existing metaLV ORA rules (n_samples is metadata only). The GO BP gmt isn't available locally yet; `generate_archs4_go_bp_gmt` builds it.
- **Local profile limits.** `large_fit_slot=0` blocks only `fit_canonical_clampfull` (archs4). ARCHS4 stability fits request `mem_mb=400000`, more than the local 96000, so they can't be scheduled locally (good). GTEx seed fits request 48 GB and *can* run, so dry-run first.
- **Environment state (2026-11-25).**
  - Conda is miniforge3 at `/home/miltondp/software/miniforge3`. `setup.sh` uses whatever conda is on PATH.
  - `clamp-analyses` and `phenoplier-cli-neo` exist; `snakemake` and `gpu-kmeans` are missing.
  - Don't rerun `setup.sh`: it would reinstall pip, R packages and CLAMP into the existing `clamp-analyses`.
  - Hardware: 32 cores, 125 GB RAM, RTX 4090 with 24 GB.
- **Local data and outputs.**
  - The ARCHS4 h5 exists at `../clamp-technical-artifacts/base/input/human_gene_v2.5.h5`. It isn't needed for this plan.
  - Local `output/` uses an old layout and holds none of the consensus or noise outputs.
- **Absolute symlinks.** `crispercas_seed_noise` writes `gene_properties.csv` and `sample_metadata.tsv.gz` in `seeds/seed*/` as absolute links, which will dangle after the rsync. The new area reads only the `lv_*_canonical.csv` files.
- **Hard-coded run counts.** `cluster_metalv.py:127` prints "/10 models" (text only). In new code, express replication as `min_frac` of the runs in the compendium or run subset (min_models = ceil(min_frac × n_runs)), and keep an absolute `min_models` only for the two frozen reproduction sets. This covers split-half (5 runs), LORO (9) and GTEx (20) with one rule. Likewise, set the k grid relative to the median number of LVs per run.
- **PhenoPLIER code to port** (`/home/miltondp/projects/phenoplier/phenoplier/libs/clustering/`):
  - What to port: `ensembles/utils.py` (`generate_ensemble`, `get_ensemble_distance_matrix`, `run_method_and_compute_agreement`), `ensembles/eac.py`, `ensembles/spectral.py` + `methods.DeltaSpectralClustering`, `ensembles/measures.py` (`aari`), and `utils.compare_arrays` (NaN-aware and weighted).
  - Dependencies are only numpy, pandas, sklearn and tqdm.
  - Required fixes:
    - sklearn ≥ 1.4 renamed `affinity=` to `metric=`.
    - sklearn ≥ 1.6 renamed `force_all_finite=` to `ensure_all_finite=`.
    - Replace the slow Python-callable co-association with a vectorized version.
    - Drop the `IPython.display` import.
  - PhenoPLIER settings, for reference:
    - Ensemble: k = 2..60 × 5 seeds per algorithm, about 4,400 partitions.
    - Consensus functions: EAC (single, complete or average linkage) and SCC (δ in 0.2–0.5).
    - Selection: ARI was chosen because AMI and NMI rise with k. They picked SCC and the top `ari_median` per k, and kept every k with `ari_median` ≥ the 75th percentile.

---

## Design ideas (for discussion)

1. **The main yardstick is reproducibility across independent runs.**
   - If a method finds real modules, consensus sets built from *disjoint* subsets of runs should agree: 5 vs 5 for ARCHS4, 10 vs 10 for GTEx with 20 runs.
   - Match the two sets with the Hungarian algorithm on cosine. Report the fraction matched at cos ≥ 0.9 and the median matched cosine.
   - Leave-one-run-out (LORO) gives a robustness score for each consensus LV.
   - This is the main reason to add GTEx runs.
2. **A semi-synthetic benchmark gives ground truth for the combining step.**
   - Take known modules (clean S123 LVs or Z_6reps) and generate pseudo-runs. In each pseudo-run, a module may be missing, perturbed to a target cosine, split or merged, and spurious LVs are added.
   - Score each method by ARI, module recovery and false consensus LVs.
3. **Nulls calibrate the thresholds.**
   - Shuffling gene labels within each LV, per run, calibrates `min_cos`.
   - Shuffling run labels while keeping clusters calibrates `min_models`: how many clusters span ≥ m runs by chance.
4. **Computational metrics, per set.**
   - Cohesion: within-cluster cosine (mean and min), top-`top_n` (185) Jaccard between members and between each member and the consensus.
   - Icasso Iq: mean similarity within a cluster minus mean similarity to other clusters.
   - Cosine silhouette.
   - Run purity: n_models / n_lvs. Ideally a cluster has one LV per run.
   - Redundancy: the largest off-diagonal cosine and top-gene Jaccard between consensus LVs.
   - Coverage and recall: the share of each run's non-technical LVs represented at cos ≥ t, and the share of "reproducible" LVs captured (best hit at cos ≥ t in ≥ m other runs).
5. **Biological metrics, per set.**
   - ORA with bp, canonical, cellmarker and azimuth: % of LVs with ≥ 1 term, % of pathways recovered, specificity (terms per LV) and redundancy (terms hit by several LVs).
   - Member pathway coherence: Jaccard of each member's significant pathways, from the per-run `summary.csv` (AUC > 0.7, FDR < 0.05).
   - Technical leakage: the gene_flag test applied to the consensus Z.
   - Later: LINCS drug–disease and PhenoPLIER traits.
6. **Read results as trade-offs, not single numbers.** Stricter sets look more cohesive but cover less. Plot cohesion against coverage and reproducibility, like a precision–recall curve, and always compare at matched set sizes too.
7. **Test the combination step separately from clustering.** Candidates, each rescaled to the median member norm:
   - `median`: the current method.
   - `run_median`: average within a run first, then take the median across runs.
   - `mean_l2`: mean of L2-normalized members.
   - `medoid`.
   - `rank1`: rank-1 SVD, with the sign set to positive and clipped at 0.
8. **Pitfalls specific to k-means.**
   - k-means forces every LV into a cluster, so singletons inflate n_models. Trim members far from the cluster medoid (`min_member_cos`, in the spirit of cNMF's outlier filter).
   - Choosing k by ARI against the ensemble leans toward the middle of the k grid. Always show the agreement-vs-k curves and confirm the choice with split-half and synthetic ARI.
   - The expected k is about #LVs / #runs, roughly 1,500 plus singletons, so the grid is about 250–3,000.
9. **Later alternatives:**
   - Spectral clustering as the consensus function (what PhenoPLIER chose).
   - A graph of reciprocal best hits across runs, partitioned with Leiden.
   - A constraint of at most one LV per run in each consensus LV.

---

## Implementation

### Phase 0 details

1. Create only the missing envs:
   ```bash
   conda create -n snakemake  --file envs/snakemake.lock  --override-channels -c conda-forge -y
   conda create -n gpu-kmeans --file envs/gpu-kmeans.lock --override-channels -c conda-forge -y
   conda run -n gpu-kmeans python -m ipykernel install --sys-prefix --name python3
   conda activate snakemake && python scripts/check_setup.py   # also runs both pin --check scripts
   ```
2. Link the fits. `data/` is gitignored, which keeps the personal path out of tracked config:
   ```bash
   ln -s /home/miltondp/projects/pivlab/ai-comp-exp/data/clamp_archs4 data/clamp_archs4
   ```
3. Rsync from pivlab01 with `rsync -a`, which keeps mtimes. The source repo is `/home/marc/Documents/pivlab/clamp-analyses/`; copy into the same relative paths:
   - `output/03_model_biology/02_archs4/06_technical_noise/`, including `gene_properties.csv`, `lv_gene_properties_canonical.csv`, `lv_sample_metadata_canonical.csv`, `lv_technical_classification.csv`, and `seeds/seed*/lv_{gene_properties,sample_metadata}_canonical.csv`.
   - `output/03_model_biology/02_archs4/04_crispercas_ge4r015med/00_mgm/`: `Z_all.npz`, `noise/all_flags.csv`, `cluster_summary.csv`, `cluster_members.csv`, `all_lv_clusters.csv`, `metaZ.csv`.
   - If they still exist, the old 6-rep outputs: `output/03_model_biology/02_archs4/04_crispercas/00_mgm/` and/or `output/98_final_models/consensus_models/archs4/`. They hold the true member lists for Z_6reps.
4. Export with `cev_export_run` (R, `clamp-analyses`, about 32 GB RAM, two at a time). It writes `runs/S{seed}/{Z.csv, summary.csv, manifest.json}`. `manifest.json` uses the `fit_clampfull.R:217-226` keys (`samples` = ncol(B), `latent_variables`, `genes`) plus L2 and md5, so `cev_ora` can pass `--model-manifest` exactly like the existing metaLV ORA rules. No B.csv. Check first how the rds files were made: if they came from `write_clamp_model` (`scripts/archs4/common.R:122-145`) and not from `fit_clampfull.R`, summary.csv is written with row names and a re-export won't match byte for byte.
5. Build the GO BP gmt with the existing rule `generate_archs4_go_bp_gmt`.

### New files

```
scripts/consensus_eval/
  README.md            staging steps, set naming, how to run
  export_run.R         --rds --model --out-dir  → Z.csv, summary.csv, manifest.json
  cev.py               library (numpy/scipy/pandas/sklearn; py3.11 + 3.13; no cuml at top level)
  run.py               CLI subcommands: linkage, hier, eac, build, evaluate, compare, null,
                       synthetic, coherence, leakage, aggregate   (flags come from cluster_metalv.py)
  kmeans_ensemble.py   cuML KMeans ensemble (gpu-kmeans env), optional PCA / gene bagging; reuse the
                       seed scheme base_seed + k*1000 + rep and provenance() from
                       scripts/archs4/projections/benchmark.py:21-31,169-186
  tests/test_cev.py    unittest on toy data (load via importlib like tests/test_crispr_l1000.py)
workflow/config/consensus_eval.yaml
workflow/rules/consensus_eval.smk
nbs/03_model_biology/02_archs4/06_consensus_eval/
  00_reproduction.ipynb  01_ensemble.ipynb  02_computational.ipynb  03_biological.ipynb  04_summary.ipynb
```

Edits to existing files:

- `workflow/Snakefile`: add `configfile: "workflow/config/consensus_eval.yaml"` and `include: "rules/consensus_eval.smk"` **after** `consensus.smk`, because the new rules reuse the `con_*` helpers from `consensus.smk` (`con_models`, `con_noise_dir`, `con_noise_files`, `con_model_dir`) and globals from `archs4_crispercas.smk` (`A4_METALV`, `A4_METALV_TN`, `A4_METALV_ORA`, `a4_metalv_ora_args`, `A4_CRISPERCAS_CFG`). Add a public `rule consensus_evaluation` at the bottom.
- Don't add it to `all` or `archs4_biology` yet; doing so would make `all` depend on the external `.rds` files.

What to reuse:

- **`cluster_metalv.py`**: call its `cache` and `flags` subcommands (`consensus_flags`, `consensus.smk:269-285`, is the template). Don't edit it, and don't use `noise` (see Key facts).
- **Config globals and helpers** instead of new literals: model order `A4_METALV_MODELS`, noise dirs `con_noise_dir("archs4", m)`, `A4_METALV_TN`, `A4_METALV`, top genes `A4_CRISPERCAS_CFG["top_n"]` (185), ORA `--top-pct {A4_CFG[ora][top_pct]}`, ORA databases `A4_METALV_ORA`.
- **`scripts/archs4/coverage/run_ora.R`**: pass `--model-manifest`, and get the database arguments from `a4_metalv_ora_args(db)`. Use `--lv-subset` for ORA deduplication and size-matched subsets instead of writing subset Z files. Aggregate the results the way `consensus_ora_aggregate` (`consensus.smk:480`) does.
- **`run_analysis.R --stage gene_level`**: use it for leakage, with `gene_properties.csv` symlinked. Its output already has r, `_FDR` and `_robust_z` columns for every test, so the variant without the robust-z rule is a column filter (|r| ≥ 0.2 and FDR < 0.2, the `TH` values at `run_analysis.R:39`), not a separate Python implementation.
- **PhenoPLIER clustering code**: port it into `cev.py` with the fixes listed in Key facts.
- **Hungarian matching**: `scipy.optimize.linear_sum_assignment`. The repo has no matcher besides the greedy one in `benchmark.py:49-62`.

`cev.py` contents:

- **Loading and distance:**
  - `load_lvs(npz, flags)` returns `LVTable(X, Xn, models, lvs, genes)`.
  - `condensed_cosine_distance(Xn)`: an exact copy of `cluster_metalv.py:109-117`, so results stay identical. It's inline in `cluster_main`, so it can't be imported without editing that file. The reproduction test guards the copy.
- **Clustering:**
  - `hier_linkage(Xn, method)`, computed once per (noise_r2, linkage, run subset) and saved (n−1 × 4, tiny). `hier_partition(L, min_cos)` is then just `fcluster`, and every min_cos / min_models set reuses the same dendrogram.
  - `coassoc(partitions)`: every partition labels every LV, so the denominator is the scalar P (number of partitions). Accumulate co-membership counts in a dense uint8 matrix (P < 256; about 225 MB at n = 15k), and write 1 − count/P straight into the condensed float64 array in row blocks. Avoid per-partition sparse `H @ H.T`, whose nnz is Σ size² and blows up when one cluster is huge.
  - `eac(D, ks, method="average")`: one linkage, then `fcluster(maxclust=k)` for each k.
  - `eac_threshold`.
  - `compare_arrays`, and `agreement(ensemble, part)` returning ARI mean, median and SD (AMI optionally, on a subsample).
- **Set building:**
  - `replicated(labels, models, m)`.
  - `trim(labels, Xn, min_member_cos)`.
  - Combination functions: `median`, `run_median`, `mean_l2`, `medoid`, `rank1`.
- **Metrics:** `top_genes(Z, n)` (n from `crispercas.top_n`), `cluster_metrics`, `iq`, `silhouette`, `redundancy`, `coverage`, `reproducible_lvs`, `hungarian_match(Za, Zb)`, and `assign_members(Z, Xn, t)` for sets where only Z is available.
- **Robustness:** `run_subsets(models, n_splits, seed)` (full, half{i}{a,b}, loro_{model}), `compare_sets` (Hungarian match for split-half, per-LV recovery for LORO), `gene_perm_null`, `label_perm_null`, `synthetic(...)`. Subset sets are built by the normal set rules, not inside these functions.
- **Storage:** sets store `Z.npy` (float32) plus `lvs.csv`. Write a CSV only where R needs it (`run_ora.R`, `run_analysis.R`).

### Config: `workflow/config/consensus_eval.yaml` (sketch)

Only what is genuinely new goes here. The rules derive the rest:
- model order: `A4_METALV_MODELS`
- noise dirs: `con_noise_dir("archs4", m)`
- gene properties: `{A4_METALV_TN}/gene_properties.csv`
- legacy Z cache: `{A4_METALV}/Z_all.npz`
- top genes: `crispercas.top_n`; ORA top % and databases: `archs4.ora.top_pct` and `A4_METALV_ORA`.

This keeps the plan from adding a third hand-synced copy of the seeds and paths.

```yaml
consensus_eval:
  root: output/03_model_biology/02_archs4/07_consensus_eval   # {root}/{compendium}/...
  notebooks: nbs/03_model_biology/02_archs4/06_consensus_eval
  compendia:
    archs4:
      rds: data/clamp_archs4/CLAMPfull_canonical_seed{seed}.rds
      md5: data/clamp_archs4/MD5SUM.txt
  reproduce:   # frozen legacy sets: absolute min_models, checked against the external Z
    cl_c050_m4_r015: {linkage: complete, min_cos: 0.5, min_models: 4, noise_r2: 0.15,
                      external: data/clamp_archs4/consensus/Z_4reps.csv}
    cl_c050_m6_r025: {linkage: complete, min_cos: 0.5, min_models: 6, noise_r2: 0.25,
                      external: data/clamp_archs4/consensus/Z_6reps.csv}
  grids:       # the .smk expands each grid into set names and the {set} wildcard constraint
    hier:     {linkage: [complete, average], min_cos: [0.4, 0.5, 0.6, 0.7],
               min_frac: [0.3, 0.4, 0.5, 0.6], noise_r2: [0.15], combine: [median]}
    ensemble: {ensemble: [km_full, km_pca256], k: [auto], min_frac: [0.4], noise_r2: [0.15],
               min_member_cos: [0.3], combine: [median, run_median, mean_l2, medoid, rank1]}
  ensembles:   # pure k-means; hierarchical cuts are merged in cev_eac from the CPU linkages
    km_full:   {pca: null, k: [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000],
                seeds: 5, base_seed: 20261125, max_iter: 300, feature_fraction: 1.0}
    km_pca256: {pca: 256, k: [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000],
                seeds: 5, base_seed: 20261125, max_iter: 300, feature_fraction: 1.0}
  hier_cuts: {linkage: [complete, average], min_cos: [0.3, 0.4, 0.5, 0.6, 0.7]}   # added to ensembles in eac
  consensus: {method: eac_average, k: [600, 800, 1000, ..., 3000], select: ari_median}
  eval: {cos_thresholds: [0.5, 0.7, 0.9], match_cos: 0.9, reproducible_cos: 0.5,
         split_half: {n_splits: 10, seed: 1}, loro: true, null: {n_perm: 3, seed: 2},
         synthetic: {n_modules: 600, p_present: [0.5, 0.9], member_cos: [0.6, 0.95],
                     p_split: 0.1, p_merge: 0.05, n_spurious: 400, seeds: [1, 2, 3]}}
  resources: {export: {mem_mb: 32000, runtime: 120}, ensemble: {mem_mb: 16000, gpu: 1, runtime: 600},
              set: {mem_mb: 12000, runtime: 60}, eval: {mem_mb: 12000, runtime: 120}}
```

Set names come from the grid values, for example `hier_complete_c050_f04_r015_median`, so a name can't contradict its parameters. Use only the resource keys that every profile already defines (`mem_mb`, `runtime`, `gpu`). A hierarchical build peaks at about 5 GB (`crispercas_metalv_cluster` requests 16 GB), so 12 GB leaves room. Raise it only where a measurement shows the need. ORA rules use `ora_slots` and `ora_mem_mb`, like `archs4_coverage.smk:365`.

### Rules: `workflow/rules/consensus_eval.smk` (prefix `cev_`)

`{subset}` is `full`, `half{i}{a,b}` or `loro_{model}`, from `run_subsets()`. Split-half and LORO sets go through exactly the same rules as the full sets.

```
cev_export_run (R)            ancient(rds) -> runs/{model}/Z.csv, summary.csv, manifest.json
cev_zcache                    runs/*/Z.csv -> Z_all.npz            (cluster_metalv.py cache, A4_METALV_MODELS order)
cev_zcache_check              Z_all.npz vs ancient(legacy_zcache) -> checks/zcache.ok
cev_flags                     ancient(con_noise_files) -> noise/flags_r{r2}.csv  (cluster_metalv.py flags)
cev_linkage                   Z_all + flags -> linkage/r{r2}/{linkage}/{subset}.npy  (one dendrogram, shared by all cuts)
cev_hier_set                  linkage + Z_all -> sets/{set}/{subset}/{Z.npy, lvs.csv, members.csv, meta.json}
cev_reproduction              reproduce sets vs ancient(external Z) -> report/reproduction.csv
                              (fails unless max|ΔZ| < 1e-6 and names equal)
cev_ensemble     [gpu-kmeans, gpu=1]  Z_all + flags -> ensembles/{ens}/{subset}/partitions.npz   (k-means only)
cev_eac                       partitions + hier cuts from cev_linkage
                              -> ensembles/{ens}/{subset}/{agreement_by_k.csv, consensus_partitions.npz, selected_k.txt}
cev_ensemble_set              -> sets/{set}/{subset}/...   (replication filter -> trim -> combine; CPU)
cev_evaluate                  sets/{set}/full -> sets/{set}/eval/{metrics.csv, per_cluster.csv, coverage.csv}
cev_compare                   sets/{set}/{half*,loro_*} -> sets/{set}/eval/{split_half.csv, loro.csv}   (one CPU rule)
cev_null                      -> null/{gene_perm,label_perm}_r{r2}.csv
cev_synthetic                 -> synthetic/{method}.csv
cev_ora_lvs                   union of unique consensus LVs (hashed by top-gene set) across sets
                              -> ora/{db}/   (run_ora.R --lv-subset, once per db)
cev_coherence                 -> sets/{set}/eval/coherence.csv
cev_leakage                   -> sets/{set}/leakage/
cev_aggregate                 -> report/{metrics_all.csv, ora_all.csv}   (per-set ORA summaries by lookup)
cev_report_00..04 (notebook:) -> report/*.complete
consensus_eval_archs4 (aggregate)  ->  Snakefile: rule consensus_evaluation
```

Rule conventions:

- Wrap every external, rsynced or legacy input in `ancient()`.
- List the new scripts as rule inputs.
- Only `cev_ensemble` runs in `gpu-kmeans`. Everything else, including EAC, runs in `clamp-analyses`, so every set is built with the same scipy.
- Keep `cev.py` free of top-level imports of cupy or cuml.
- The external Z files feed only `cev_reproduction`. They become sets of their own (evaluated separately) only if reproduction fails (step 1.4).
- ORA results per LV don't depend on the set (BH runs within each LV). Nested grid sets share most of their LVs, so ORA runs once on their union.

Split-half and LORO procedure:

- Split-half draws seeded halves of the runs; `min_models = ceil(min_frac × n_runs_in_subset)`.
- `cev_compare` matches the two halves' sets with the Hungarian algorithm.
- LORO reports each consensus LV's recovery rate at cos ≥ 0.9 across the leave-one-run-out sets.
- Ensembles and linkages are keyed by `{subset}`, not by set. Sets that share an ensemble differ only in combine or trim, so they reuse the same GPU work.

### Pitfalls

- **Memory and compute.**
  - For about 15k LVs, the uint8 co-association counts take about 225 MB and the condensed float64 distance about 0.9 GB. Avoid dense `1 − C` and `squareform` temporaries.
  - One cosine-distance-plus-linkage build takes about 12 s at n = 15k (8 s distance, 3 s linkage, measured), which is why linkages are cached per (r2, linkage, subset).
  - k-means on 18k dimensions takes about 1 s per iteration at k = 3000 on CPU, so use cuML on the GPU; PCA-256 is about 70× cheaper.
  - AMI with thousands of clusters is slow, so use ARI as the main criterion.
- **Rescaling.** `mean_l2`, `medoid` and `rank1` change the loading scale. ORA uses ranks only, but projection uses L2 from the reference, so rescale to the median member norm.
- **Size.** Set sizes differ (618 vs 1,058 vs the ensemble sets), and ORA % and coverage depend on size. Report size-matched results: top-N by n_models, plus random subsets (via `run_ora.R --lv-subset` / lookups into the shared ORA results).
- **GTEx notebooks.** Decide before Phase 5 whether GTEx gets sibling notebooks under `nbs/03_model_biology/01_gtex/` or executed-notebook names that include `{compendium}`.
- **GTEx seeds.** Don't append the extra seeds to `consensus.seeds`: ARCHS4 shares that list, and `con_models("archs4")` would point at missing models.

---

## Verification

- `python -m unittest discover -s scripts/consensus_eval/tests` (in `clamp-analyses`). The tests cover:
  - `coassoc` against PhenoPLIER's slow version.
  - `eac` recovering a planted partition.
  - `hier` + `median` reproducing `cluster_metalv.py cluster` exactly.
  - The weighting in `run_median`.
  - `rank1` being non-negative.
  - Iq, purity and coverage on hand-computed cases.
  - Hungarian matching.
  - Ground truth from the synthetic generator.
  - `min_frac` → min_models resolution for full, half and LORO subsets.
- `snakemake --profile workflow/profiles/local --dry-run consensus_evaluation` lists only `cev_*` jobs.
- `cev_zcache_check` and `cev_reproduction` pass.
- Sanity checks:
  - A stricter `min_cos` raises cohesion and lowers coverage.
  - The gene-permutation null gives close to zero replicated clusters.
  - On synthetic data, a perfect method has ARI ≈ 1.
  - The metrics separate Z_4reps from Z_6reps in sensible directions; check this before scoring the ensemble sets.
- Notebooks 00–04 render, and their `.executed.ipynb` files get reviewed.

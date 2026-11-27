# PLAN: evaluating and improving consensus LVs (metaLVs / MGMs)

This is a living plan on branch `consensus-lvs`. Work through it one step at a time, across sessions.

## How to resume

1. Read this file. Find the first unchecked box in **Status**, and read the **Log**.
2. Read **Key facts** before touching rules, configs or `cluster_metalv.py`.
3. Run `snakemake --profile workflow/profiles/local --dry-run <target>` before any real run. The dry run must list only `cev_*` jobs: no CLAMP fits, no `technical_noise_*`, no `crispercas_*`.
4. When a step is done, tick its box and add a dated line to the **Log**.
5. Working style:
   - Present a plan for the step before implementing it. Closely tied steps can be grouped, with a checkpoint between them.
   - After a run, report the results and review them with the user. Before relying on a metric, check that it is informative (correlations across sets, nulls).
   - Commit only when the user asks.

## Status

**Blocked on:** nothing. Every legacy input is now staged in `data/legacy_archs4/` (see Key facts). The CLAMPfull per-run noise dirs were never found, but none of the steps need them any more.

### Phase 0: Setup and data staging (local)
- [x] 0.1 Create the missing envs (`snakemake`, `gpu-kmeans`), register the kernel, and run `scripts/check_setup.py` (it also runs both `--check` pin scripts).
- [x] 0.2 Symlink `data/clamp_archs4` to the `ai-comp-exp` fits (md5 check skipped at the user's request; see Log).
- [x] 0.3 ~~Rsync the noise outputs and the legacy `00_mgm` outputs from pivlab01.~~ They aren't there. Instead, staged the old 4-rep and 6-rep consensus folders, the canonical S123 Z/summary and `gene_properties.csv` under `data/legacy_archs4/`. The 6-rep `noise_flags.csv` is the flag source (see Key facts).
- [x] 0.4 Export Z.csv, summary.csv and manifest.json from the 10 `.rds` files (`cev_export_run`). Check that the S123 export equals `data/legacy_archs4/canonical_S123/{Z,summary}.csv` byte for byte.
- [x] 0.5 Build the GO BP gmt for ORA if it is missing (`generate_archs4_go_bp_gmt`).

### Phase 1: Reproduce the existing consensus sets (stop if this fails)
- [x] 1.1 Scaffold the new area: config, `.smk`, Snakefile include and target, staging rules (`cev_zcache`, `cev_flags`, `cev_reproduce_set`, `cev_reproduction`; CLI `scripts/consensus_eval/run.py`).
- [x] 1.2 Check that the exported `Z_all.npz` has the same (model, LV) keys as the legacy `noise_flags.csv`, plus 18,423 genes with no NaNs (`cev_zcache_check`). There's no legacy `Z_all.npz` to compare against.
- [x] 1.3 Check that `cl_c050_m4_r015` equals `Z_4reps.csv` exactly, and that its members equal `data/legacy_archs4/consensus_4reps/mgm_mapping.csv`.
- [x] 1.4 Check that `cl_c050_m6_r025` equals `Z_6reps.csv`, and that its members equal `data/legacy_archs4/consensus_6reps/mgm_mapping.csv`. If they differ, the member lists show where.
- [x] 1.5 Notebook `00_reproduction`, kept short: an internal check, needed only because the legacy `Z_all.npz` and noise dirs were never found. Show the two `report/reproduction/*.csv` tables and one plot of clusters by number of runs spanned (`reproduce/{set}/cluster_model_count_distribution.csv`). Nothing in Phase 2 depends on it. It gets removed in 4.3.

### Phase 2: Evaluation framework, applied to the existing sets
- [x] 2.1 `cev.py` metrics and unit tests (toy data).
- [x] 2.2 Computational metrics on Z_4reps, Z_6reps and the hierarchical grid (`cev_evaluate`).
- [x] 2.2b Revise the computational metrics after the 2026-11-26 review (see Key facts, "Metric review"):
  - A new "reproducible LV" definition: reciprocal best hit, cos ≥ 0.5, and top-gene Jaccard above the gene-permutation null, in ≥ `ceil(0.4 × n_runs) − 1` other runs.
  - `overlap_to_consensus` replaces `jac_to_consensus`.
  - Drop `iq` and `silhouette`.
  - Write a headline table.
  - Details are in Implementation, "Step 2.2b".
- [x] 2.2c Hierarchical clustering on top-gene overlap (`topj_*` sets): distance = 1 − Jaccard of each LV's top-185 positive genes, instead of cosine on raw loadings. Same grid machinery and the same evaluation. Grid `min_jac` [0.05, 0.1, 0.15, 0.2, 0.3] × both linkages × `min_frac` [0.3–0.6] (40 sets). Also added `recall_multi_c0.9` to the headline. Results are in the Log and Key facts.
- [x] 2.2d Wrap `CEV_SCRIPT` / `CEV_LIB` (`run.py`, `cev.py`) in `ancient()`, so code edits no longer rerun the cev area (166 jobs today, about 2k after 2.4). Results are in the Log and Key facts.
- [x] 2.4 Robustness: split-half and LORO for the hierarchical sets (both `hier_*` and `topj_*`), plus both nulls. Why first, and the decisions taken (2026-11-27): Implementation, "Step 2.4". Numbering is kept, so older references to "2.4" still mean split-half.
- [ ] 2.3 **(next)** Biological metrics: ORA, member pathway coherence, technical leakage, on the 2.4 shortlist only (user, 2026-11-27):
  - `hier_complete_c070_f04`, `hier_complete_c060_f05`, `hier_average_c070_f04`, `hier_average_c060_f05`;
  - `topj_complete_j020_f05`, `topj_average_j030_f04`, `topj_average_j020_f04` (all `_r015_median`);
  - plus Z_4reps and Z_6reps (their twins) as references.
  - Biology has to break the tie: on reproducibility the two families are level for multi-gene modules (see Log).
- [ ] 2.5 Semi-synthetic benchmark.
- [ ] 2.6 Notebooks `02_computational` and `03_biological`.

### Phase 3: Ensemble (consensus-index) clustering
- [ ] 3.1 k-means ensemble on the GPU (`km_full`), k-means only.
- [ ] 3.2 Co-association, EAC average linkage, and agreement-vs-k with k selection.
- [ ] 3.3 Build ensemble sets with replication filter, trim and the five combination functions.
- [ ] 3.4 Run the same evaluation on the ensemble sets, including split-half, LORO and synthetic.
- [ ] 3.5 Variants: add hierarchical cuts to the ensemble, PCA-256, and (later) spectral consensus. After 2.2c: k-means on raw loadings is cosine-like, so it will likely inherit the single-gene problem. Consider adding the `topj` cuts to the ensemble, or k-means on a transformed input (binary top-gene sets, or sqrt/rank loadings), and score it with `frac_single_gene` and `recall_multi_c0.9`.
- [ ] 3.6 Notebook `01_ensemble`.

### Phase 4: Comparison and recommendation
- [ ] 4.1 Notebook `04_summary`: trade-off curves, size-matched comparisons, recommended set.
- [ ] 4.2 (Later) Connect the chosen set to the `consensus.smk` downstream steps (projection, LINCS, GLS).
- [ ] 4.3 Cleanup once everything else is done: remove the `00_reproduction` notebook (and its `.executed.ipynb`) and its `cev_report_00` rule, and drop it from the aggregate. Before that, move any number later notebooks cite (flag counts, set sizes) into `02_computational` or `04_summary`. Keep the `cev_reproduction` / `cev_zcache_check` rules: they are cheap, and they guard the inputs of every set.

### Phase 5: GTEx extension (if 5-vs-5 split-half on ARCHS4 is too noisy)
- [ ] 5.1 Stage the GTEx inputs: TPM gct, FBM, canonical S124 and the existing 10 runs.
- [ ] 5.2 Add a per-compendium `compendia.gtex.extra_seeds: [1008..1017]` in `consensus.yaml` (not global, so recount2's `CON_SEED_PATTERN` is unchanged), and make `con_seed_index` continue past `len(seeds)`.
- [ ] 5.3 Fit the 10 extra runs and their noise screens with the existing `consensus_seed_full` / `consensus_model_noise` rules (48 GB, 8 threads each, two at a time). No `cev_export_run` for GTEx: consensus.smk already writes Z.csv, summary.csv and manifest.json.
- [ ] 5.4 Add `compendia.gtex` to `consensus_eval.yaml`, deriving its 20 models from `con_models("gtex")` plus the extras and its paths from `con_model_dir` / `con_noise_dir`, then run the evaluation with a 10-vs-10 split-half.

## Log

- 2026-11-26: Explored the repo, the PhenoPLIER ensemble code and the ai-comp-exp data, then wrote this plan. Nothing executed yet.
- 2026-11-26: Simplify pass on the plan: derive models/paths/top_n/ORA settings from existing config; one `{subset}` wildcard replaces the split-half/LORO rule variants; cached linkages, uint8 co-association, deduplicated ORA, `Z.npy` storage, `min_frac` replication, manifest.json export, GTEx via the existing consensus.smk chain.
- 2026-11-26: Step 0.1, partly done. Created `snakemake` (8.30.0, SLURM plugin 2.8.0) and `gpu-kmeans` (cuml 25.12, cupy 13.6, CUDA 12.9) from the locks and registered the gpu-kmeans `python3` kernel. A cuML KMeans smoke test on the 4090 passed. `check_setup.py` passes snakemake, R notebook packages, GPU packages and workflow parsing, but fails three `clamp-analyses` checks (CLAMP revision, comparison packages, pip imports); see Key facts. Waiting on a decision about the `clamp-analyses` env.
- 2026-11-26: Step 0.1 done. Per the user's choice, removed and recreated `clamp-analyses` from `envs/clamp-analyses.lock`, then ran setup.sh's remaining clamp-analyses steps. `check_setup.py` now exits 0 with all checks PASS.
- 2026-11-26: Step 0.2 done. Linked `data/clamp_archs4` to the ai-comp-exp fits (10 rds files plus `consensus/Z_{4,6}reps.csv`). md5 verification was skipped at the user's request. `MD5SUM.txt` lists all 10 rds files. Note: `seed123.rds` dates from 2026-09-17, while the other nine are from 2026-11-23.
- 2026-11-26: Step 0.3 probe (read-only, via ssh to pivlab01). The planned source paths don't exist. `/home/marc/Documents/pivlab/clamp-analyses/output/03_model_biology/` holds only `01_gtex`. There's no `02_archs4` (so no `06_technical_noise` and no `04_crispercas*/00_mgm`), and the folder isn't a git checkout. The `pivlab` group gets traverse-only ACLs, so directories can be stat'ed but not listed. The old 6-rep manifest's `reference_model` points at `/home/msubirana/Documents/pivlab/clamp-analyses/...`, which doesn't exist on pivlab01 and is probably on another machine. What *is* on pivlab01 is `output/98_final_models/consensus_models/archs4/`, the old 6-rep consensus model:
  - `Z.csv` has the same md5 as the local `Z_6reps.csv` (`4865841…`).
  - `mgm_mapping.csv` and `cluster_summary.csv` hold the true 6-rep member lists (618 metaLVs).
  - `noise_flags.csv` has columns `model, LV, gene_flag, technical_joint, sample_flag, any_flag` for all 17,280 LVs (1,728 × 10). Checked: `sample_flag == (technical_joint >= 0.25)` and `any_flag == gene_flag | sample_flag` for every row, with no NAs. So flags at any `noise_r2` can be derived from this file: at 0.25, 873 gene, 550 sample and 1,324 any; at 0.15, 1,495 sample and 2,117 any (15,163 kept).
  - Also there: `ora/{archs4_cell_lines,azimuth,bp,canonical,cellmarker}/`, `traits/`, `manifest.json`, `lincs-projection.*`, `phenoplier_*`, and `B.csv` (7.5 GB, not needed).

- 2026-11-26: Step 0.3 closed, as the user chose: fetched the 6-rep folder without `B.csv` to `data/legacy_archs4/consensus_6reps/` (231 MB, 53 files; `Z.csv` md5 checked). The flag source is now the legacy `noise_flags.csv`. Revised 1.2 and 1.4, the `cev_flags` and `cev_zcache_check` rules, and the config sketch to match.
- 2026-11-26: Searched pico (dbmila) at the user's suggestion and found the 4-rep consensus folder (`Z.csv` = `Z_4reps.csv`, true members), the canonical S123 Z/summary and the gene-level `gene_properties.csv`. Staged all three under `data/legacy_archs4/` (see Key facts). Added member checks to 1.3, an export byte-check to 0.4, and config keys. Leakage is no longer blocked.
- 2026-11-26: Step 0.4 done. Added `scripts/consensus_eval/export_run.R`, a minimal `workflow/config/consensus_eval.yaml`, `workflow/rules/consensus_eval.smk` (`cev_export_run`, `cev_export_check`, `consensus_eval_archs4`), and the Snakefile configfile/include plus the public `consensus_evaluation` target (pulled forward from 1.1). The dry run listed only the 10 `cev_export_run` jobs, `cev_export_check` and the target. Ran in about 1.5 min with 6 jobs in parallel. Output is `{root}/archs4/runs/S*/{Z.csv, summary.csv, manifest.json}` (966 MB). `cev_export_check` passed: S123 `Z.csv` and `summary.csv` are byte-identical to `data/legacy_archs4/canonical_S123/`, even though the local `seed123.rds` (re-saved Sep 17) is 9,479 bytes larger than the original fit's rds. The manifest md5 is read from `MD5SUM.txt`, not recomputed: the user had already verified the rds md5s, so don't re-check them.
- 2026-11-26: Step 0.5 done. `data/pathways/go_bp.Hs.symbols.gmt` was missing, so ran `generate_archs4_go_bp_gmt` (the dry run listed only that job; it took about 25 s). It has 15,413 GO BP terms, and its sha256 matches the rule's pin (`33e559d…`).
- 2026-11-26: Steps 1.1–1.4 done. As the user chose, the two reproduction sets are built by the **unchanged** `cluster_metalv.py` (`cache` + `cluster`), so Phase 1 tests only our inputs; cev.py's builder is checked against these outputs in 2.1. Added `run.py` (`flags`, `zcache-check`, `reproduction`), the `reproduce:` config block (each entry carries its `external` Z and legacy `members` file) and rules `cev_zcache`, `cev_zcache_check`, `cev_flags`, `cev_reproduce_set` (→ `{root}/archs4/reproduce/{set}/`) and `cev_reproduction` (→ `report/reproduction/{set}.{csv,log}`). The dry run listed only those 9 jobs; the whole run took about 2.5 min. Results: the cache has 17,280 LVs × 18,423 genes, with keys equal to the legacy flags. Flags: 2,117 at r2 0.15 (15,163 kept) and 1,324 at 0.25 (15,956 kept). `cl_c050_m4_r015` gives 7,860 clusters and keeps 1,058; `cl_c050_m6_r025` gives 8,025 and keeps 618. **Both `metaZ.csv` files are byte-identical to `Z_4reps.csv` / `Z_6reps.csv`** (max|ΔZ| = 0), and every metaLV's member set, `n_models` and `n_lvs` equal the legacy `mgm_mapping.csv`. This confirms the `A4_METALV_MODELS` order and the re-thresholded legacy flags.
- 2026-11-26: Step 1.5 done. Added `nbs/03_model_biology/02_archs4/06_consensus_eval/00_reproduction.ipynb` (R) and its `.executed.ipynb`, run by rule `cev_report_00` (→ `report/00_reproduction.complete`, wired into `consensus_eval_archs4`). Re-added the `notebooks:` config key, and declared `cluster_model_count_distribution.csv` as a `cev_reproduce_set` output. The notebook shows both reproduction tables (all pass) and a bar plot of clusters by runs spanned, with totals of 7,860 / 1,058 and 8,025 / 618. Step 4.3 removes the notebook, its `.executed.ipynb`, `cev_report_00` and its aggregate line. Machine note: a stale fontconfig cache made R/Cairo plots render "sans" with `EBGaramond-Initials`; `fc-cache -f` at user and system level fixed it.
- 2026-11-26: Steps 2.1 + 2.2 done (grouped). Added `scripts/consensus_eval/cev.py` with only what 2.2 uses: pool loading, the copied cosine distance, hier linkage/partition, `resolve_min_models` (rounds first: 0.3 × 10 → 3), `replicated` / `build_set` (`median` only), and metrics (`top_genes` as `run_ora.R`: top 185 *positive* loadings, since about 90% of loadings are 0 and some LVs have < 185 positive genes; `cluster_metrics`, `iq_silhouette`, `redundancy`, `pool_best_hits`, `reproducible`, `coverage`, `summarize`). Added `tests/test_cev.py` (8 tests, all pass). One of them shows the builder equals `cluster_metalv.cluster_main` on toy data at three settings, loaded via importlib, not edited. Added `run.py` subcommands `linkage`, `hier-set`, `set-check`, `pool-hits`, `evaluate` and `aggregate`, the config `grids.hier` / `eval` / resources, and rules `cev_linkage` (3 jobs), `cev_hier_set` (34 sets = 2 reproduce + 32 grid), `cev_set_check`, `cev_pool_hits`, `cev_evaluate` and `cev_aggregate` (→ `report/metrics_all.csv`), all wired into `consensus_eval_archs4`. The whole target took 3 min 14 s, peaking at 4.5 GB; `sets/` is 2.2 GB. **Both `cev_set_check`s pass with float32 Z bit-identical** to `reproduce/{set}/metaZ.csv` and the same names and members. A float64 absolute tolerance of 1e-6 failed at 7.6e-6 on loadings up to about 350 (relative 6e-8), so `set-check` now requires exact float32 equality. `hier_complete_c050_f04_r015_median` has the same metrics as `cl_c050_m4_r015`. The evaluation pool is fixed at r2 0.15 for every set (15,163 LVs, 8,727 reproducible at cos ≥ 0.5 in ≥ 4 runs). Sanity: within a linkage, raising min_cos raises cohesion and lowers coverage, and raising min_frac shrinks the set and lowers coverage. Purity *falls* slightly with min_frac (complete c040: 0.979 → 0.969), because sets spanning more runs pick up more same-run duplicates. Z_6reps vs Z_4reps: 618 vs 1,058 LVs, median mean_cos 0.797 vs 0.759, silhouette 0.71 vs 0.64, coverage at 0.5 0.37 vs 0.54, recall at 0.5 0.63 vs 0.90. Average linkage is much less redundant than complete (frac of LVs with another at cos ≥ 0.5: 0.04–0.18 vs 0.13–0.29), but its members are less tight (median min_cos 0.36–0.72). Note: every rule lists `run.py` as an input, so editing it reruns the whole cev chain, about 3 min, cheap for now. Deferred from the `cev.py` list below: `coassoc`/`eac`/`agreement`/`trim` and the other combine functions (Phase 3), `hungarian_match`/`run_subsets`/nulls (2.4), `synthetic` (2.5). `assign_members` was dropped because every set has members.
- 2026-11-26: Reviewed the 2.2 metrics with the user, using read-only scratch scripts (not kept; 2.2b re-implements their logic). Findings are in Key facts, "Metric review" and "Single-gene-dominated LVs". Decisions:
  - Headline = size (`n_lvs`), fidelity (`overlap_to_consensus`), separation (`frac_redundant_c0.5`) and `recall_c0.9`, with split-half (2.4) as the main score on top.
  - Drop `iq` and `silhouette`.
  - Adopt the new reproducible-LV definition (7,981 LVs vs 8,727).
  - Added step 2.2b (metric revision) and step 2.2c (clustering on top-gene overlap, which the user wants prioritized, before 2.3).
- 2026-11-26: Step 2.2b done. Changes:
  - `cev.pool_best_hits` now also returns the best hit's pool row (`arg`) and its top-gene Jaccard. New functions: `reciprocal`, `jaccard_null` (`rng.permuted` within each top-gene row of run b), `overlap`, `top1_share`, `jaccard_rows` and `quantiles`. `hit_stages` builds the cos / +reciprocal / +Jaccard qualifying-hit masks and `reproducible(ok, min_runs)` counts runs over one mask. The single-gene cutoff is `eval.single_gene` in the config (`--single-gene`).
  - Dropped `iq_silhouette`, `indicator`, `jac_to_consensus`, `iq` and `silhouette`.
  - `cev_pool_hits` now computes the reproducible mask once and writes `best_hits.npz` (`hits, arg, jac, repro, …`; `recip` is `cev.reciprocal(arg, models, runs)`) plus `null.json` (null quantiles, observed Jaccard quantiles, counts at each stage). `cev_evaluate` reads `repro` from that file. `cev_aggregate` also writes `report/headline.csv`; its leading set-definition columns are `CEV_SET_PARAMS` in the smk (derived from the `CEV_SETS` keys).
  - Config: `eval.reproducible` gained `jac_null_q 0.999, null_pairs 3, null_seed 3`, and there is a new `eval.headline`.
  - Ran 10 unit tests (all pass), then the whole chain (84 cev jobs, 1 min 53 s, 4.6 GB peak).
  - Results:
    - Null p99.9 = 0.0306 (max 0.036).
    - The reproducible cascade is **8,727 → 8,071 → 7,981**, exactly as in the exploration.
    - 1.4% of reciprocal cos ≥ 0.5 hits have Jaccard within the null.
    - `cev_set_check` and `cev_reproduction` still pass.
  - Headline:
    - Z_4reps: 1,058 LVs, overlap 0.554, redundant 0.187, recall_c0.9 0.428.
    - Z_6reps: 618, 0.612, 0.100, 0.324.
    - At about 1,040 LVs, average c070_f03 vs complete c050_f04: overlap 0.597 vs 0.554, redundant 0.182 vs 0.187, recall 0.451 vs 0.428.
  - Informativeness (Spearman ρ over the 32 grid sets): `median_overlap_to_consensus` is still the cohesion axis (ρ 0.99 with mean_cos and cos_to_consensus). It becomes independent only for the 2.2c sets. New diagnostic `frac_single_gene` (share of consensus LVs with top1_share > 0.5) is 0.44–0.63, higher than the pool's 28%, and has ρ 0.94 with n_lvs. See Key facts.
- 2026-11-27: Step 2.2c done, run as planned (`topj_*` sets). Changes:
  - `cev.condensed_jaccard_distance` (sparse intersections; distance 1 for empty rows) and `hier_linkage_jac`; `hier_partition(L, min_sim)`.
  - `run.py linkage --dist {cos,jac} --top-n` (npz stores `dist`; `hier-set` asserts it), `grids.topj` (40 sets), `cev_topj_name()`, a `dist` key on every `CEV_SETS` entry, and linkages at `linkage/r{r2}/{dist}_{linkage}/` (old dirs deleted).
  - 2 new tests (Jaccard vs `pdist` over several blocks with empty rows; planted-module recovery), 13 in all, all pass.
  - Run: 166 cev jobs in 4 min 54 s, about 5 GB peak. Both reproductions and set checks still pass bit-identically. Pool top-185 sets: 27 LVs have fewer than 185 positive genes (min 71), none are empty.
  - Set sizes: 198–1,918, nothing degenerate, so no grid cuts were dropped. Only complete j030_f06 (198), average j030_f06 (274) and complete j005_f03 (1,918) fall outside 300–1,700.
  - Size-matched against the hier set of the same linkage, median topj − hier differences:
    - `frac_single_gene` −0.14, lower in every pair (−0.03 to −0.26). topj spans 0.13–0.54, hier 0.44–0.63. **The main question: yes.**
    - `frac_redundant_c0.5` +0.01: lower at loose cuts (j0.05–0.10), higher at strict ones (j0.2–0.3).
    - `recall_c0.9` −0.07, lower in every pair.
    - `median_cos_to_consensus` −0.05, `median_mean_cos` −0.09 (lower at loose cuts, about equal at j ≥ 0.15).
  - The recall drop is mostly single-gene LVs at about 450–550 LVs, but not at Z_4reps' size. So, as the user chose, added `recall_multi_c0.9` (recall over the 4,735 reproducible pool LVs with top1_share ≤ 0.5) to `cev.coverage` / `summarize` and the headline:
    | Set | n_lvs | frac_single_gene | recall_c0.9 | recall_multi_c0.9 |
    |---|---|---|---|---|
    | Z_4reps (= hier_complete_c050_f04) | 1,058 | 0.57 | 0.428 | 0.359 |
    | topj_complete_j010_f04 | 1,077 | 0.45 | 0.343 | 0.290 |
    | hier_average_c070_f03 | 1,040 | 0.59 | 0.451 | 0.382 |
    | topj_average_j015_f03 | 1,091 | 0.50 | 0.332 | 0.254 |
    | hier_complete_c070_f05 | 453 | 0.49 | 0.306 | 0.293 |
    | topj_complete_j020_f05 | 464 | 0.28 | 0.257 | 0.289 |
    | hier_average_c060_f06 | 560 | 0.50 | 0.317 | 0.297 |
    | topj_average_j020_f05 | 526 | 0.31 | 0.254 | 0.266 |
  - Z_4reps vs the size-matched topj sets: 62% (complete j005_f05) and 69% (average j015_f03) of Z_4reps LVs have a topj counterpart at cos ≥ 0.9. The unmatched ones are no more single-gene than the matched ones (top1_share median 0.62 vs 0.66). The topj LVs with no Z_4reps counterpart are clearly multi-gene (median top1_share 0.21–0.28): modules that cosine clustering missed.
  - Informativeness within topj: overlap and cosine cohesion still move together (ρ 0.98). `frac_redundant_c0.5` now follows cohesion (ρ 0.6) rather than size (−0.04), and `recall_c0.9` depends less on size (0.52 vs 0.87 for hier).
  - Decisions (user, 2026-11-27): tick 2.2c, keep `recall_multi_c0.9` in the headline, and **run 2.4 (split-half/LORO) before 2.3**. The rationale is in Implementation, "Step 2.4". Not run yet.
- 2026-11-27: Step 2.2d done. Wrapped `CEV_SCRIPT` / `CEV_LIB` in `ancient()` at their definitions (one edit; all 13 uses are `{input.script}` / `{input.lib}`, so shell commands are unchanged). Checks: the dry run said "Nothing to be done" before the edit and again after `touch scripts/consensus_eval/{run,cev}.py`. `--forcerun cev_evaluate` lists only 74 `cev_evaluate` + `cev_aggregate` + the target (76 jobs), with no linkage, set or pool-hits builds. Rendered commands still show `scripts/consensus_eval/run.py`.
- 2026-11-27: Step 2.4 done. Changes:
  - `cev.load_pool(runs=, gene_perm_seed=)`; `hungarian_match`, `split_half_scores`, `loro_recovery`, `label_perm_null`, `unit_cols` (an all-zero consensus LV gets cosine 0, so it never matches).
  - `run.py hier-sets` (batched; shares `build_hier_set` with `hier-set`), `compare`, `label-null`; `aggregate --robust` writes `report/robustness.csv` and adds `sh_share_c0.9_{all,multi}` and `loro_recovery_c0.9_all` to the headline.
  - smk: `cev_run_subsets` (`report/subsets.json`), `cev_subset_sets`, `cev_null_label`, `cev_compare`. A grid entry may now be a list of grids, which is how the Z_6reps twin is added. `CEV_CHECK_REF` set-checks the twins.
  - Config: `nulls:` (not `null:`, which YAML reads as None).
  - 18 tests, all pass.
  - Run: 395 jobs in 15.5 min (4.8 GB peak). The 4 `topj_average_j005` compares failed: in `half4a` that cut chains 205 LVs into one cluster whose median is all zeros. After the `unit_cols` fix, all 73 compares were rerun (`--forcerun cev_compare`, 2.5 min).
  - All 4 set checks pass: both twins are bit-identical to Z_4reps / Z_6reps.
  - **Disk: `subsets/` is 133 GB.**
  - Results (`report/robustness.csv`):
    - **Precision is enough on ARCHS4.** Split-to-split SD of `sh_share_c0.9_all` has median 0.010 (max 0.016), against a range of 0.11–0.54 across sets and about 0.1 between families. The 10 splits overlap, so the SD understates sampling error, but the gaps are about 10 SD.
    - **Cosine wins split-half at matched size, in every pair.** topj − hier medians: `sh_share_c0.9_all` −0.115, multi −0.144, single −0.076. But the matched top-gene Jaccard favors topj in 60% of pairs (median +0.04; +0.06 to +0.18 at j ≥ 0.2), and LORO is about tied (−0.02; topj is ahead at j0.2–0.3). The split-half score matches on cosine, so it leans toward the cosine family, just as overlap leans toward topj. The tie is only half broken.
    - **The single-gene worry did not hold.** In hier sets at strict cuts, multi-gene LVs reproduce *better* than single-gene ones (c070: 0.58–0.59 vs 0.49–0.50), and hier's lead is largest in the multi stratum.
    - **Gene-permutation null:** 0 consensus LVs for every topj set. Cosine sets get chance clusters at f0.3 (c040 102–125, c050 48–64, c060 21–26, c070 10–11) and at c040–c050 f0.4 (2–7). Shuffled single-gene LVs collide on a hub gene. So min_frac ≥ 0.4 for cosine sets.
    - **Run-label null is weak:** the shuffled-label count is 48–95% of the observed count at the set's m. It shows that run span is mostly a size filter, and it does not calibrate min_frac. Don't headline it.
    - **Prune topj j0.05 and j0.10:** split-half 0.11–0.33, LORO 0.72–0.86, plus the degenerate chain cluster.
    - **Z_4reps is dominated.** hier_complete_c070_f03 (1,023) and hier_average_c070_f03 (1,040) score split-half 0.44–0.45 vs 0.357, multi 0.43–0.44 vs 0.30, and LORO 0.91 vs 0.89, though at f0.3 they carry 10–11 null clusters. Z_6reps' twin (618): 0.496 / 0.488 / 0.839; hier_complete_c060_f05 (598): 0.507 / 0.516 / 0.865 with 0 null clusters.
    - LORO rises as min_frac falls (looser thresholds survive losing a run), so LORO alone can't choose min_frac.
  - Follow-up, as the user chose: a **Jaccard-matched** split-half (Hungarian on top-gene Jaccard, share at Jaccard ≥ `eval.match_jac` 0.25: `jshare_*`) and LORO (`loro_recovery_j0.25_*`). Only `cev_compare` was rerun (4.8 min). Its split SD is 0.009 (range 0.23–0.72).
    - Size-matched topj (j ≥ 0.15) − hier (f ≥ 0.4), 24 pairs:

      | Metric | Median | topj wins |
      |---|---|---|
      | `sh_share_c0.9_all` | −0.107 | 0% |
      | `sh_jshare_j0.25_all` | +0.081 | 100% |
      | `sh_jshare_j0.25_multi` | +0.011 | 75% |
      | `loro_recovery_c0.9_all` | +0.029 | 79% |
      | `loro_recovery_j0.25_all` | +0.046 | 96% |

    - So **each family wins split-half in its own metric**, and for multi-gene LVs the Jaccard-matched gap is about 0. topj's overall Jaccard lead comes from hier's single-gene LVs, whose top-185 tails don't reproduce. LORO slightly favors topj at strict cuts. Reproducibility does not separate the families, so biology (2.3) decides.
    - Shortlist, as `n_lvs`, then sh cos / sh jac / LORO cos (multi-gene in brackets):

      | Set | n_lvs | sh cos | sh jac | LORO cos |
      |---|---|---|---|---|
      | hier_complete_c060_f05 | 598 | 0.51 (0.52) | 0.61 (0.76) | 0.87 |
      | hier_average_c060_f05 | 747 | 0.50 (0.51) | 0.62 (0.78) | 0.87 |
      | topj_complete_j020_f05 | 464 | 0.43 (0.42) | 0.70 (0.78) | 0.86 |
      | topj_average_j030_f04 | 451 | 0.41 (0.39) | 0.69 (0.77) | 0.91 |
      | Z_4reps | 1,058 | 0.36 (0.30) | 0.49 (0.67) | 0.89 |
      | Z_6reps | 618 | 0.50 (0.49) | 0.62 (0.79) | 0.84 |

  - Decisions (user, 2026-11-27): the 2.3 shortlist is in Status. **Phase 5 (GTEx) is deferred:** 5-vs-5 ranks methods with gaps of about 10 SD, so it is needed only later, as a generalization check. Keep the 133 GB of subset sets until Phase 3 is done.
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

- Noise flags are not recomputed: that would need ARCHS4 preprocessing at about 200 GB RAM. (Revised 2026-11-26: the per-run noise dirs aren't on pivlab01, so flags come from the legacy `noise_flags.csv`; see Key facts.)
- Everything runs locally, including any extra GTEx fits.
- The work goes in a new Snakemake area, `consensus_eval` (rule prefix `cev_`), following the CLAUDE.md conventions from the start.

---

## Key facts (found during exploration)

- **Z_6reps is not a subset of Z_4reps.**
  - Its column names match the first 618 of Z_4reps, but its values differ: 87 of 618 columns have no Z_4reps match at cos ≥ 0.99, and 36 have none at ≥ 0.9.
  - The cause is commit `3085ab6`, which changed the config from `min_models` 6 / `noise_r2` 0.25 to 4 / 0.15.
  - Reproducing Z_6reps therefore needs flags at 0.25 and `min_models` 6.
- **Flag source (2026-11-26).** `data/legacy_archs4/consensus_6reps/noise_flags.csv` (fetched from pivlab01, from the old 6-rep consensus model) has exactly `cluster_metalv.py`'s `COLS` (`model, LV, gene_flag, technical_joint, sample_flag, any_flag`) for all 17,280 LVs. It was written at threshold 0.25. `technical_joint` is the 5-fold CV R² with no NAs, so flags at any r2 are `sample_flag = technical_joint >= r2` and `any_flag = gene_flag | sample_flag`, which is exactly `flags()` in `cluster_metalv.py:36-45`. `cev_flags` re-thresholds this file instead of calling `cluster_metalv.py flags`. The result plugs straight into `cluster_metalv.py cluster --flags`. Counts: 2,117 flagged at r2 0.15 (15,163 kept) and 1,324 at 0.25. `gene_flag` (873) is fixed, because it comes from `TH` in `run_analysis.R`.
- **Other legacy files in `data/legacy_archs4/consensus_6reps/`:** `Z.csv` (= `Z_6reps.csv`), `mgm_mapping.csv` / `cluster_summary.csv` (true 6-rep member lists), `manifest.json`, `ora/{bp,canonical,cellmarker,azimuth,archs4_cell_lines}/`, `traits/`, `lincs-projection.*`, `phenoplier_*`. Fetched from pico (`dbmila`, `/pividori_lab/marc_projects/`) on 2026-11-26:
  - `data/legacy_archs4/consensus_4reps/` comes from `clamp-analyses/output/98_final_models/consensus/archs4/`, minus `B.csv` (12.9 GB). `Z.csv` md5 = `Z_4reps.csv` (`44ce257…`). `mgm_mapping.csv` has the true 4-rep members (1,058 metaLVs, 6,591 members). Also present: `cluster_summary.csv`, `manifest.json` (min_models 4), `CLAMPfull_consensus.rds`, `ora/` (7 dbs, incl. reactome and gtex_tissues), `traits/`, `lincs-projection.*`, `projection.json`, `samples.tsv`, `annotations.csv.gz`.
  - `data/legacy_archs4/canonical_S123/{Z,summary}.csv, manifest.json` come from `clamp-analyses/output/98_final_models/clampfull/canonical/archs4/`. They're the reference for the byte-level export check in 0.4.
  - `data/legacy_archs4/noise_inputs/gene_properties.csv` comes from `tmp_clampbase_noise/inputs/`. It's gene-level (cDNA length, GC% per gene, 18,423 genes, mapped from the ARCHS4 h5), so it doesn't depend on the model and is fine for leakage (2.3).
  - Checked: no 4-rep member is flagged under the derived r2 = 0.15 flags, and no 6-rep member at 0.25.
  - **Don't use** `tmp_clampbase_noise/archs4_noise/`, `archs4_flags_r2_0.15.csv` or `archs4_zcache.npz`. They're the **CLAMPbase** screen: technical_joint correlates with CLAMPfull at r = 0.05, and gene_flag agrees on only 93% of LVs.
  - Not reachable: `clamp-analyses/output/98_final_models/consensus_models/` (owner-only `drwx------`). The CLAMPfull ARCHS4 `06_technical_noise/` and `04_crispercas*/00_mgm/` (incl. `Z_all.npz`) weren't found on either server. Marc's original checkout was at `/home/msubirana/Documents/pivlab/clamp-analyses` on an unknown machine.
  - Access: dbmila needs a password, so open a shared connection first (`ssh -fN -o ControlMaster=yes -o ControlPath=~/.ssh/cm-dbmila -o ControlPersist=8h dbmila`), then use `-o ControlPath=~/.ssh/cm-dbmila`. pivlab01 is the sakaiza server.
- **Do not edit `cluster_metalv.py`.** Nine rules list it as an input without `ancient()` (`archs4_crispercas.smk:79,95,109`; `consensus.smk:272,291,309,658,673,691`), so any edit reruns the existing consensus chain (ORA, GLS, LINCS) for every compendium. `rerun-triggers: mtime` doesn't help, because an edit changes the file's mtime. Call its subcommands, or load it via importlib like `tests/test_crispr_l1000.py` does. Possible follow-up: wrap those 9 inputs in `ancient()` (precedent: `archs4_crispercas.smk:50`), then factor `cluster_main` into importable functions.
- **`technical_noise_z` takes `Z.csv` without `ancient()`.** Never write a Z.csv into `output/98_final_models/clampfull/canonical/archs4/`: that triggers a rerun which needs the raw h5.
- **Model completeness gating** (`archs4_canonical.smk:12-34`).
  - A model counts as complete only when all of these exist and are non-empty: rds, B.csv, Z.csv, summary.csv, manifest.json, complete.
  - A partly filled folder raises `ValueError` while the job graph is being built.
  - The new area therefore stages its own files under its own root and never fills those folders.
- **What the `.rds` contains.** A list with `B, Z, U, C, L1, L2, Z2, heldOutGenes, Uauc, Up, summary, priorMatCV, priorMat, withPrior` (+ `call` in S123 only). Z is already an 18,423 × 1,728 data.frame with gene rownames, B a 1,728 × 605,614 data.frame, and summary has `LV`; L2 ≈ 456.94 in every run (one run differs at the 11th digit). Loading takes about 40 s and 12 GB RSS. `fit_clampfull.R:179-215` writes `write.csv(Z)` with gene rownames, `write.csv(B)`, and `write.csv(summary, row.names=FALSE)` with `LV_index` renamed to `LV`. Exporting the same way reproduces those files exactly. B is about 1,728 × 580k, about 18 GB as CSV. It isn't needed, so skip it.
- **Legacy model order:** `S123, S1001..S1007, S456, S789` (`A4_METALV_SEEDS` is sorted with key `(s<1000, s)`). The order decides linkage tie-breaks and metaLV numbering, so it is required for exact reproduction.
- **ORA** queries the top `ceil(18423 × 0.01)` = **185** positive genes. `run_ora.R` takes `--model-manifest`, like the existing metaLV ORA rules (n_samples is metadata only). The GO BP gmt is at `data/pathways/go_bp.Hs.symbols.gmt` (built in 0.5).
- **Local profile limits.** `large_fit_slot=0` blocks only `fit_canonical_clampfull` (archs4). ARCHS4 stability fits request `mem_mb=400000`, more than the local 96000, so they can't be scheduled locally (good). GTEx seed fits request 48 GB and *can* run, so dry-run first.
- **Environment state (2026-11-25).**
  - Conda is miniforge3 at `/home/miltondp/software/miniforge3`. `setup.sh` uses whatever conda is on PATH.
  - As of 2026-11-26 all three repo envs are built from their locks and `scripts/check_setup.py` passes every check. `clamp-analyses` was rebuilt from scratch (old env had CLAMP `818e13b`, no pip deps or comparison packages), then had setup.sh's pip, kernel, `install_analysis_packages.R` and `install_clamp.R` steps run (CLAMP now at the pinned `748559a`). `phenoplier-cli-neo` is untouched.
  - The user-level `ir-clamp` kernel (`~/.local/share/jupyter/kernels/ir-clamp`) points into the `clamp-analyses` prefix, so it now runs the rebuilt env.
  - Note: `conda create --file <explicit lock> --download-only` ignores `--download-only` and installs the env anyway.
  - Don't rerun `setup.sh`: it would reinstall pip, R packages and CLAMP into the existing `clamp-analyses`.
  - Hardware: 32 cores, 125 GB RAM, RTX 4090 with 24 GB.
- **Local data and outputs.**
  - The ARCHS4 h5 exists at `../clamp-technical-artifacts/base/input/human_gene_v2.5.h5`. It isn't needed for this plan.
  - Local `output/` uses an old layout and holds none of the consensus or noise outputs.
- **Hard-coded run counts.** `cluster_metalv.py:127` prints "/10 models" (text only). In new code, express replication as `min_frac` of the runs in the compendium or run subset (min_models = ceil(min_frac × n_runs)), and keep an absolute `min_models` only for the two frozen reproduction sets. This covers split-half (5 runs), LORO (9) and GTEx (20) with one rule. Likewise, set the k grid relative to the median number of LVs per run.
- **Metric review (2026-11-26, the 34 sets of step 2.2).** Spearman ρ across sets shows the 2.2 metrics collapse into three independent axes:
  - **Cohesion:** `mean_cos`, `min_cos`, `cos_to_consensus`, `jac_members`, `jac_to_consensus`, `iq` and `silhouette` correlate at 0.97–1.00.
  - **Size:** `coverage_c*` and `recall_c0.5` correlate at 0.99 with `n_lvs`.
  - **Redundancy:** `frac_redundant_c0.5` has ρ 0.06 with cohesion and 0.55 with size.
  - `recall_c0.9` (ρ 0.88 with size) is the only coverage measure that also asks whether the consensus looks like the run LVs.
  - Separation candidates tried:
    - Member-level ones collapse into cohesion or size: silhouette and Iq (ρ ≥ 0.98 with cohesion), assignment margin (own consensus − best other; ρ 0.93 with cohesion, −0.90 with size), misassigned members (0–1.3% in every set, no spread), and the pool-LV best-minus-second margin (ρ −0.99 with size).
    - Between-consensus separation does carry information: `frac_redundant_c0.5`, and the share of covered pool LVs matching ≥ 2 consensus LVs (ρ 0.96 with `frac_redundant_c0.5`, so the same thing). It separates complete from average linkage (13–29% vs 2–18%).
  - Cohesion metrics are partly circular for cosine-built sets: complete linkage guarantees `min_cos` ≥ the cutoff.
  - Purity barely varies: 8% of Z_4reps clusters have a duplicate run.
  - After 2.2b (32 grid sets): `median_overlap_to_consensus` has ρ 0.99 with `median_mean_cos` and `cos_to_consensus`, and −0.71 with `n_lvs`. For cosine-built sets it is the cohesion axis under another name. It is headlined because it stays non-circular for the 2.2c (Jaccard-built) sets, where cosine becomes the independent check.
- **Null calibration and reproducible LVs (r2 0.15 pool, 15,163 LVs).** The null shuffles gene order within every LV of one run, which keeps each LV's loading distribution and sparsity; it was run on 3 run pairs, seed 0.
  - Null best-hit cosine: median 0.103, p99 0.545, p99.9 0.841, max 0.93. So cos 0.5 is about the 1%-by-chance level for one run pair, and requiring 3 other runs makes chance matches negligible.
  - Observed best-hit cosine: median 0.454.
  - Top-gene Jaccard (top 185) null best hit: p99.9 = max = 0.031, against an observed best-hit median of 0.114 (p25 0.076, p75 0.221). Much better calibrated than cosine.
  - Only 57% of best hits are reciprocal: 71% at cos ≥ 0.3, 88% at ≥ 0.5, 98% at ≥ 0.7.
  - Replication score (cos of the 3rd-best other-run hit) is spread evenly with no break: p10 0.29, p50 0.56, p90 0.90.
  - Reproducible counts with m = 4 runs including its own:

    | Definition | LVs |
    |---|---|
    | best hit cos ≥ 0.5 (2.2) | 8,727 |
    | + reciprocal | 8,071 |
    | + Jaccard > null p99.9 (**adopted**) | 7,981 |
    | reciprocal + Jaccard > null, no cos cutoff | 11,212 (too lenient: 0.031 ≈ 11 shared genes) |

  - Among reciprocal cos ≥ 0.5 hits, only 1.4% have Jaccard within the null, so the Jaccard guard costs little here. It may matter more for 5-run halves or GTEx.
- **Single-gene-dominated LVs.** 28% of pool LVs carry > 50% of their squared norm on one gene (median share 0.24, p90 0.79). CLAMP loadings are non-negative and about 90% zeros.
  - Cosine between such LVs is driven by that gene: the null best hits at cos ≥ 0.5 come from LVs with a median top-gene share of 0.73.
  - 39% of the reproducible LVs (3,387 of 8,727) are single-gene dominated.
  - Consensus sets are *more* single-gene dominated than the pool (2.2b `frac_single_gene`):
    - 0.44–0.63 across the grid, against 28% of pool LVs;
    - Z_4reps 0.57, Z_6reps 0.46;
    - ρ 0.94 with n_lvs, so larger and looser sets admit more of them.
  - Jaccard of reciprocal best hits by cosine bin (2.2b):
    | Cosine bin | Median Jaccard | p10 |
    |---|---|---|
    | [0.3, 0.5) | 0.082 | 0.031 |
    | [0.5, 0.7) | 0.146 | 0.048 |
    | [0.7, 0.9) | 0.289 | 0.118 |
    | ≥ 0.9 | 0.418 | 0.242 |

    Even near-identical-cosine pairs share less than half of their top-185 genes. Over all reciprocal hits, Jaccard quantiles are p25 0.079, p50 0.167, p75 0.321; over those with cos ≥ 0.5, p25 0.138, p50 0.25, p75 0.40.
  - So the legacy method (cosine on raw loadings) partly clusters by top gene. This motivates 2.2c, and it is why the fidelity metric uses top-gene overlap rather than cosine.
- **Top-gene (topj) vs cosine sets (2.2c).** Clustering on 1 − top-185 Jaccard cuts single-gene consensus LVs a lot at matched size (median −0.14 `frac_single_gene`), at no gain in redundancy. Strict cuts (j0.15–0.2, f0.5–0.6, about 450–650 LVs) keep multi-gene recall and cosine cohesion close to the hier sets. At about 1,000 LVs, topj loses multi-gene recall too (0.25–0.29 vs 0.36–0.38). `recall_c0.9` rewards recovering single-gene LVs, so read it with `recall_multi_c0.9` (headline). No 2.2 metric breaks the tie between the families; split-half (2.4) is meant to.
- **Robustness (2.4).** Split-half on ARCHS4 is precise (split SD ≈ 0.01) but metric-dependent: cosine matching favors hier, top-gene matching favors topj, and for multi-gene LVs they tie. Use both `sh_share_c0.9_*` and `sh_jshare_j0.25_*` and compare at matched size. In 5-run halves min_frac collapses (0.3/0.4 → 2, 0.5/0.6 → 3). The gene-permutation null rules out cosine sets at f0.3 (10–125 chance LVs), and topj sets are null-clean. The run-label null is uninformative. A loose average top-gene cut (j0.05) can chain hundreds of LVs into one cluster whose median is all zeros (`unit_cols` handles it).
- **Builder exactness.** Set Z is float32. Compare it to CSVs in float32 (exact equality), not with an absolute float64 tolerance: loadings reach about 350, so the float32→float64 rounding is about 1e-5.
- **Code inputs are `ancient()` (2.2d).** `CEV_SCRIPT` (`run.py`) and `CEV_LIB` (`cev.py`) are wrapped once at their definitions in `consensus_eval.smk`, so new rules inherit it. Editing either file no longer reruns anything. After a real logic change, rerun the affected rule with `--forcerun <rule>`; its downstream jobs follow because their inputs become newer. `export_run.R` and `CEV_METALV_SCRIPT` stay plain inputs on purpose: they change rarely (see "Do not edit `cluster_metalv.py`").
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
4. **Computational metrics, per set** (revised after the 2026-11-26 review; see Key facts):
   - **Headline**, one metric per independent axis:
     - size: `n_lvs`;
     - fidelity: `overlap_to_consensus`, the share of a member's top-185 genes that are in the consensus's top-185, median over members;
     - separation: `frac_redundant_c0.5`, the share of consensus LVs with another consensus LV at cos ≥ 0.5;
     - `recall_c0.9`, the share of reproducible LVs with a consensus LV at cos ≥ 0.9.
   - **Diagnostics** (in `per_cluster.csv` / `metrics.csv`, not headlined): member cosine (mean, min), `cos_to_consensus`, member top-gene Jaccard, purity (n_models / n_lvs), max cos / Jaccard to another consensus LV, and coverage of each run at cos ≥ t.
   - **Dropped:** Icasso Iq and silhouette (identical to cohesion here).
   - **Reproducible LV:** has a *reciprocal* best hit with cos ≥ 0.5 and top-gene Jaccard above the gene-permutation null's p99.9, in ≥ `ceil(min_frac × n_runs) − 1` other runs (min_frac 0.4). With 10 runs that is 3 other runs; in 5-run halves it is 1.
   - Always compare sets at matched size, because coverage and recall rise with `n_lvs`.
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
   - (Promoted to step 2.2c.) Cluster on top-gene overlap instead of cosine on raw loadings. Possible follow-ups: rank-weighted overlap (e.g. RBO) or cosine on transformed loadings (sqrt / rank).
   - Spectral clustering as the consensus function (what PhenoPLIER chose).
   - A graph of reciprocal best hits across runs, partitioned with Leiden.
   - A constraint of at most one LV per run in each consensus LV.

---

## Implementation

### Phase 0 details

Steps 0.1–0.3 are done; the Log and Key facts record how (`data/clamp_archs4` is a symlink because `data/` is gitignored, which keeps the personal path out of tracked config).

4. (Done.) Export with `cev_export_run` (R, `clamp-analyses`; readRDS takes about 40 s and peaks at 12 GB, so the rule requests 16 GB). It writes `runs/S{seed}/{Z.csv, summary.csv, manifest.json}`. `manifest.json` uses the `fit_clampfull.R:217-226` keys (`samples` = ncol(B), `latent_variables`, `genes`) plus L2 and md5, so `cev_ora` can pass `--model-manifest` exactly like the existing metaLV ORA rules. No B.csv. All 10 rds files were saved by `fit_clampfull.R` after its conversions (Z is a data.frame with gene rownames, summary already has `LV`), so the export is plain `write.csv` and matches byte for byte.
5. Build the GO BP gmt with the existing rule `generate_archs4_go_bp_gmt`.

### New files

```
scripts/consensus_eval/
  README.md            staging steps, set naming, how to run
  export_run.R         --rds --model --out-dir  → Z.csv, summary.csv, manifest.json
  cev.py               library (numpy/scipy/pandas/sklearn; py3.11 + 3.13; no cuml at top level)
  run.py               CLI subcommands: linkage, hier, eac, build, evaluate, compare, null,
                       synthetic, coherence, leakage, aggregate   (flags come from cev_flags)
  kmeans_ensemble.py   cuML KMeans ensemble (gpu-kmeans env), optional PCA / gene bagging; reuse the
                       seed scheme base_seed + k*1000 + rep and provenance() from
                       scripts/archs4/projections/benchmark.py:21-31,169-186
  tests/test_cev.py    unittest on toy data (load via importlib like tests/test_crispr_l1000.py)
workflow/config/consensus_eval.yaml
workflow/rules/consensus_eval.smk
nbs/03_model_biology/02_archs4/06_consensus_eval/
  00_reproduction.ipynb (temporary, removed in 4.3)  01_ensemble.ipynb  02_computational.ipynb  03_biological.ipynb  04_summary.ipynb
```

Edits to existing files:

- `workflow/Snakefile`: add `configfile: "workflow/config/consensus_eval.yaml"` and `include: "rules/consensus_eval.smk"` **after** `consensus.smk`, because the new rules reuse the `con_*` helpers from `consensus.smk` (`con_models`, `con_model_dir`, and `con_noise_dir` for GTEx only) and globals from `archs4_crispercas.smk` (`A4_METALV_MODELS`, `A4_METALV_ORA`, `a4_metalv_ora_args`, `A4_CRISPERCAS_CFG`). For ARCHS4, `con_noise_dir` / `A4_METALV_TN` point at the noise dirs that were never found, so don't use them. Add a public `rule consensus_evaluation` at the bottom.
- Don't add it to `all` or `archs4_biology` yet; doing so would make `all` depend on the external `.rds` files.

What to reuse:

- **`cluster_metalv.py`**: call its `cache` subcommand, and feed `cev_flags` output to `cluster --flags`. Don't edit it (see Key facts). Flags come from re-thresholding the legacy file, not from its `flags` or `noise` subcommands.
- **Config globals and helpers** instead of new literals: model order `A4_METALV_MODELS`, top genes `A4_CRISPERCAS_CFG["top_n"]` (185), ORA `--top-pct {A4_CFG[ora][top_pct]}`, ORA databases `A4_METALV_ORA`.
- **`scripts/archs4/coverage/run_ora.R`**: pass `--model-manifest`, and get the database arguments from `a4_metalv_ora_args(db)`. Use `--lv-subset` for ORA deduplication and size-matched subsets instead of writing subset Z files. Aggregate the results the way `consensus_ora_aggregate` (`consensus.smk:480`) does.
- **`run_analysis.R --stage gene_level`**: use it for leakage, with `gene_properties.csv` symlinked. Its output already has r, `_FDR` and `_robust_z` columns for every test, so the variant without the robust-z rule is a column filter (|r| ≥ 0.2 and FDR < 0.2, the `TH` values at `run_analysis.R:39`), not a separate Python implementation.
- **PhenoPLIER clustering code**: port it into `cev.py` with the fixes listed in Key facts.
- **Hungarian matching**: `scipy.optimize.linear_sum_assignment`. The repo has no matcher besides the greedy one in `benchmark.py:49-62`.

`cev.py` contents (2.1 implemented loading/distance, hier, `replicated`, `median` and the metrics except `hungarian_match`; the rest lands with the step that first uses it):

- **Loading and distance:**
  - `load_lvs(npz, flags)` returns `LVTable(X, models, lvs, genes)`; `Xn` is a lazy `cached_property`.
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
- **Metrics:** `top_genes(Z, n)` (n from `crispercas.top_n`), `cluster_metrics`, `redundancy`, `pool_best_hits`, `reproducible`, `coverage`, `summarize`, and `hungarian_match(Za, Zb)` (2.4). `iq_silhouette` is removed in 2.2b, and `assign_members` was dropped.
- **Robustness:** `run_subsets(models, n_splits, seed)` (full, half{i}{a,b}, loro_{model}), `compare_sets` (Hungarian match for split-half, per-LV recovery for LORO), `gene_perm_null`, `label_perm_null`, `synthetic(...)`. Subset sets are built by the normal set rules, not inside these functions.
- **Storage:** sets store `Z.npy` (float32) plus `lvs.csv`. Write a CSV only where R needs it (`run_ora.R`, `run_analysis.R`).

### Step 2.2b: metric revision (done; see Log)

- **`cev.pool_best_hits`** also returns the best-hit index for each (LV, other run), so reciprocity can be checked: LV a's best hit b in run B has a as its best hit in run A. It also returns the top-gene Jaccard of each best-hit pair.
- **`cev.jaccard_null(B, runs, pairs, seed)`** shuffles the gene order within every top-gene row of run b and takes the best-hit Jaccard of each run-a LV, over `null_pairs` run pairs chosen by `null_seed`. It returns the `jac_null_q` quantile.
- **Outputs:** `cev_pool_hits` writes `hits`, `arg`, `jac`, `runs` and pool keys to `best_hits.npz`, plus `null.json` (threshold, quantiles, pairs, seed).
- **`cev.reproducible(hit_stages(...)["cos_reciprocal_jaccard"], m)`** = ≥ m − 1 other runs where the hit is reciprocal, cos ≥ `cos`, and Jaccard > the null threshold. m = `ceil(min_frac × n_runs)`.
- **Per-cluster metrics:**
  - Add `overlap_to_consensus`: |top(member) ∩ top(consensus)| / |top(member)|, the mean over members per consensus LV. Members with fewer than 185 positive genes use their own size.
  - Drop `jac_to_consensus`, `iq` and `silhouette`, and remove `iq_silhouette` and its test.
  - Keep `cos_to_consensus` as a diagnostic: it becomes the non-circular check for the 2.2c sets.
- **Headline table:** `cev_aggregate` also writes `report/headline.csv`, with the set parameters plus the `eval.headline` columns.
- **Tests:** reciprocity on a hand-built case where a broad LV is the non-reciprocal best hit of two specific ones; the Jaccard null on a toy (shuffled rows give a small threshold); and `overlap_to_consensus` on hand-made sets.
- **Run:** `cev_pool_hits`, `cev_evaluate` and `cev_aggregate`; the rest of the chain also reruns because `run.py` changed (about 3 min).
- **Check:** `n_reproducible` ≈ 7,981. The exploration used null seed 0 (Jaccard null p99.9 0.031); the configured `null_seed` may move the count slightly. Record the new headline table in the Log.

### Step 2.2c: hierarchical clustering on top-gene overlap (done; see Log)

The plan below was written 2026-11-27, after 2.2b, and was run as written.

**Idea.** Use the same hierarchical machinery as the `hier_*` sets, with one change to the distance: 1 − Jaccard of the LVs' top-185 positive-gene sets (`cev.top_genes`, the `run_ora.R` selection), instead of 1 − cosine of raw loadings. The consensus is still the median of raw loadings, and the evaluation rules are unchanged.

**Grid** (`grids.topj` in `consensus_eval.yaml`, 40 sets):
```yaml
topj: {linkage: [complete, average], min_jac: [0.05, 0.1, 0.15, 0.2, 0.3],
       min_frac: [0.3, 0.4, 0.5, 0.6], noise_r2: [0.15], combine: [median]}
```
- Why these values (2.2b, Key facts "Single-gene-dominated LVs"):
  - The null maximum is 0.036, so 0.05 is the loosest cut that is still above chance.
  - Reciprocal cos ≥ 0.5 hits have Jaccard p25 0.14, p50 0.25 and p75 0.40.
  - Even cos ≥ 0.9 pairs have a median Jaccard of only 0.42. Complete linkage requires *every* pair to clear the cut, so complete at 0.2–0.3 will probably give very small sets. Average linkage is more forgiving.
- After the run, keep the cuts whose set sizes fall in the cosine grid's range (about 300–1,700 LVs). Drop the cuts that give tiny sets from the grid, and record that in the Log. Don't tune beyond that.
- Name: `topj_{linkage}_j{min_jac*100:03d}_f{min_frac*10:02d}_r{r2*100:03d}_{combine}`, e.g. `topj_average_j010_f04_r015_median`. Add `cev_topj_name()` next to `cev_hier_name()` with the same assert on min_frac, and expand `grids.topj` into `CEV_SETS` in the same loop.

**Code changes:**
- `cev.py`:
  - `condensed_jaccard_distance(B)`: same layout and 256-row blocks as `condensed_cosine_distance`.
    - Intersections come from a float32 matmul of the boolean rows; float32 is exact for counts up to 185.
    - union = |a| + |b| − intersection.
    - The distance is 1 − inter/union, and **1 when union = 0** (an LV with no positive loadings). Check how many pool LVs have fewer than 185 positive genes, or none, and log it.
    - Use a sparse `B @ B[block].T` if the dense version is slow. B has about 185 nonzeros per row out of 18,423.
  - `hier_linkage(X, method)` is unchanged. Add `hier_linkage_jac(B, method)` and keep the cosine path untouched, so the two reproduce sets stay bit-identical.
  - `hier_partition(L, min_sim)`: rename the argument only. It is still `fcluster(t=1 − min_sim)`.
- `run.py`:
  - `linkage` gains `--dist {cos,jac}` and `--top-n`. For `jac`, it builds `B = top_genes(pool.X, top_n)` and stores `dist` in the npz.
  - `hier-set` takes the cut from `p["min_cos"]` or `p["min_jac"]` according to `p.get("dist", "cos")`. It asserts that the linkage npz's `dist` matches.
  - `headline.csv` picks up `dist` and `min_jac` automatically (`CEV_SET_PARAMS` is derived from the `CEV_SETS` keys).
- `consensus_eval.smk`:
  - Every set definition in `CEV_SETS` gets a `dist` key: `cos` for the reproduce and `hier` sets, `jac` for `topj`. Existing set names don't change.
  - The linkage path becomes `linkage/r{r2}/{dist}_{linkage}/{subset}.npz`, with a `dist` wildcard constrained to `cos|jac`. `cev_linkage` passes `--dist {wildcards.dist} --top-n {A4_CRISPERCAS_CFG[top_n]}`, and `cev_hier_set` builds the path from `CEV_SETS[wc.set]["dist"]`.
  - The cosine linkages move to `cos_*`, so the whole chain reruns (about 2 min). Delete the old `linkage/r*/{complete,average}/` dirs after the run.
- No rule edits are needed for evaluation or aggregation: `cev_evaluate` and `cev_aggregate` expand over `CEV_SETS`, so the topj sets join `metrics_all.csv` and `headline.csv` automatically.

**Tests** (`scripts/consensus_eval/tests/test_cev.py`):
- `condensed_jaccard_distance` equals `pdist(B, "jaccard")` on boolean toys. The toys need several blocks (set `BLOCK` small or use n > 256) and include one all-False row, where we expect 1. Note that scipy gives 0 there, so compare that row separately.
- A planted-module toy: the jac linkage recovers the modules at a cut between the within-module and between-module Jaccard.

**Run:**
- Dry run `consensus_evaluation`. Expect only `cev_*` jobs: 5 linkages (3 cos at their new paths, 2 jac), 74 `cev_hier_set`, 74 `cev_evaluate`, plus the checks and the aggregate.
- Do the real run. Record the time and peak memory of a jac linkage in Pitfalls, next to the cosine figure of about 12 s.

**Review with the user (and record in the Log):**
- Set sizes over the grid, and which cuts land in 300–1,700.
- Compare at matched `n_lvs` against the `hier_*` sets: pair each topj set with the hier set of the same linkage closest in size, or plot each headline metric against `n_lvs` by family. Look at:
  - `frac_single_gene`: is it lower than the cosine sets' 0.44–0.63? This is the main question.
  - `frac_redundant_c0.5`: is it lower at equal size?
  - `recall_c0.9`: does it change? It is cosine-based, so it is the independent check here.
  - `median_cos_to_consensus` and `median_mean_cos`: the non-circular cohesion check for Jaccard-built sets. `median_overlap_to_consensus` is now partly circular.
- Informativeness: Spearman ρ of the headline metrics within the topj family, and whether overlap and cosine cohesion come apart there (they had ρ 0.99 for the cosine sets).
- How many of Z_4reps' consensus LVs have a close topj counterpart (cos ≥ 0.9) in the size-matched topj set, and what the unmatched ones look like (`top1_share`).
- Then tick 2.2c, add a Log line, and update Key facts. Split-half (2.4) and biology (2.3) then apply to the topj sets like every other set.

### Step 2.4: split-half and LORO (next; runs before 2.3)

Moved ahead of 2.3 on 2026-11-27, after the 2.2c review. Not run yet. Present a plan for it first (working style); the procedure is already sketched in Design ideas 1 and 3, "Rules" and "Split-half and LORO procedure" below.

**Why it goes first.**
- **It is the agreed main score, and nothing else can settle cosine vs top-gene.** The 2.2 metrics are either size-driven (coverage, recall), circular for one family (cohesion: cosine for `hier_*`, top-gene overlap for `topj_*`), or biased by single-gene LVs (`recall_c0.9` rewards recovering them). After 2.2c the two families trade off (topj: far fewer single-gene LVs; hier: higher recall at Z_4reps' size), and there is no neutral metric to break the tie. Split-half asks directly whether a method's consensus LVs come back when it is run on disjoint runs, and it applies identically to every family, including the Phase 3 ensembles.
- **It shrinks the grid before the expensive step.** 2.3 runs `run_ora.R` on 4 databases plus leakage and coherence. Running it on all 74 sets is wasteful when 2.4 can leave a few candidates per family. ORA is deduplicated by top-gene set, but the grid sets still differ enough that the union is large.
- **It decides Phase 5 early.** If 5-vs-5 halves on ARCHS4 are too noisy to rank methods, the fix is 10 extra GTEx fits (48 GB, 8 threads each). That's the biggest compute decision left, so it should be made before more analysis is built on ARCHS4 alone.
- **It is cheap and needs no new inputs.** It reuses `cev_linkage` / `cev_hier_set` through the `{subset}` wildcard (`CEV_SUBSETS`, today only `full`). A 5-run half has about 7.5k LVs, so each linkage costs about a quarter of the full one.
- **The label-permutation null calibrates `min_frac`.** Shuffling run labels while keeping clusters gives the number of clusters spanning ≥ m runs by chance. That number would justify (or prune) the `min_frac` axis of both grids before 2.3 and Phase 3 inherit it.

**Settle with the user before running.**
- **Scale.** With `split_half.n_splits: 10` there are 20 halves plus 10 LORO subsets, i.e. 30 subsets. That's 150 linkages (5 dist × linkage × r2 combinations) and 2,220 set builds. Each build reloads `Z_all.npz`, so this may take about an hour. Options: start with n_splits 5, build only the `full` set's family (not the reproduce sets), or build every cut of a subset in one job.
- **Reproduce sets.** `cl_c050_m4_r015` has the same metrics as `hier_complete_c050_f04_r015_median`, so leave it out of the subsets (an absolute `min_models` 4 means something else in a 5-run half). `cl_c050_m6_r025` (r2 0.25) has no grid twin. Either skip it, or add a `hier_complete_c050_f06_r025` twin.
- **Score.** Hungarian match of half a vs half b (`scipy.optimize.linear_sum_assignment`). Report the fraction matched at cos ≥ 0.9 and the median matched cosine, **plus the matched top-gene Jaccard**, so that neither family is scored only in its own metric. LORO: each full-set consensus LV's recovery rate at cos ≥ 0.9 across the 10 leave-one-out sets.
- **Stratify by single-gene.** A single-gene LV probably reproduces easily, because the same hub gene recurs in every run, so split-half could favor the cosine sets for the reason 2.2c set out to fix. Report split-half separately for single-gene and multi-gene consensus LVs (`top1_share > eval.single_gene`), in the same way as `recall_multi_c0.9`.
- **Size.** Match rates depend on set size, so compare at matched `n_lvs` as in 2.2c.

**Decisions (user, 2026-11-27).**
- **10 splits, batched builds.** 30 subsets (`half{0..9}{a,b}`, `loro_{model}`) plus a `geneperm0` pseudo-subset. Halves are seeded draws of 5-run combinations that contain the first model, so no split is drawn twice. Subset sets are built by one `cev_subset_sets` job per (r2, dist, linkage, subset), which loads the cache once and writes every cut of that dendrogram to `subsets/r{r2}/{dist}_{linkage}/{subset}/{set}/`. `cev_hier_set` stays `full`-only and unchanged.
- **Reproduce sets** are not rebuilt on subsets. Z_6reps gets a grid twin instead, `hier_complete_c050_f06_r025_median` (a second `grids.hier` entry). `cev_set_check` checks each reproduce set's twin (`hier_complete_c050_f04_r015_median` too) against the `cluster_metalv.py` output.
- **Gene-permutation null: 1 permutation** (`geneperm0`: all runs, gene order shuffled within every pool LV, through the same linkage and set rules; expect about 0 kept clusters). The run-label null (`eval.null`) does the calibration. It shuffles run labels over the full dendrogram's cut and counts clusters that span ≥ m runs.
- **Scores** (`cev_compare` → `sets/{set}/eval/{split_half,loro,robust}.csv`):
  - **Split-half.** Hungarian match on cosine. `sh_share_c0.9` = the share of LVs in both halves whose partner has cos ≥ `match_cos` (= matched pairs / mean size). Also the median paired cosine and top-gene Jaccard. Each of these is reported for all, single-gene and multi-gene LVs (the LV's own `top1_share`), as mean and SD over the 10 splits.
  - **LORO.** Each full-set LV's share of the 10 leave-one-out sets with a best cos ≥ `match_cos`.
- **min_frac collapses in 5-run halves.** `resolve_min_models` gives 2, 2, 3, 3 for min_frac 0.3, 0.4, 0.5, 0.6, so in split-half f03 = f04 and f05 = f06. LORO (9 runs) gives 3, 4, 5, 6, the same as full. So split-half ranks (linkage, cut) at 2 or 3 runs, and LORO is what separates min_frac.

### Config: `workflow/config/consensus_eval.yaml` (sketch)

Only what is genuinely new goes here. The rules derive the rest:
- model order: `A4_METALV_MODELS`
- noise flags: `consensus_eval.compendia.archs4.legacy_flags` (legacy file; `con_noise_dir` isn't available locally)
- gene properties: `consensus_eval.compendia.archs4.gene_properties` (legacy file; used only for leakage)
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
      legacy_flags: data/legacy_archs4/consensus_6reps/noise_flags.csv
      gene_properties: data/legacy_archs4/noise_inputs/gene_properties.csv
      reference_dir: data/legacy_archs4/canonical_{model}   # export check (staged con_ref model)
  reproduce:   # frozen legacy sets: absolute min_models, checked against the external Z
    cl_c050_m4_r015: {linkage: complete, min_cos: 0.5, min_models: 4, noise_r2: 0.15,
                      external: data/clamp_archs4/consensus/Z_4reps.csv,
                      members: data/legacy_archs4/consensus_4reps/mgm_mapping.csv}
    cl_c050_m6_r025: {linkage: complete, min_cos: 0.5, min_models: 6, noise_r2: 0.25,
                      external: data/clamp_archs4/consensus/Z_6reps.csv,
                      members: data/legacy_archs4/consensus_6reps/mgm_mapping.csv}
  grids:       # the .smk expands each grid into set names and the {set} wildcard constraint
    hier:     {linkage: [complete, average], min_cos: [0.4, 0.5, 0.6, 0.7],
               min_frac: [0.3, 0.4, 0.5, 0.6], noise_r2: [0.15], combine: [median]}
    ensemble: {ensemble: [km_full, km_pca256], k: [auto], min_frac: [0.4], noise_r2: [0.15],   # not in the yaml yet
               min_member_cos: [0.3], combine: [median, run_median, mean_l2, medoid, rank1]}
    topj:     {linkage: [complete, average], min_jac: [0.05, 0.1, 0.15, 0.2, 0.3],   # 2.2c, not in the yaml yet
               min_frac: [0.3, 0.4, 0.5, 0.6], noise_r2: [0.15], combine: [median]}
  ensembles:   # pure k-means; hierarchical cuts are merged in cev_eac from the CPU linkages
    km_full:   {pca: null, k: [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000],
                seeds: 5, base_seed: 20261125, max_iter: 300, feature_fraction: 1.0}
    km_pca256: {pca: 256, k: [250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000],
                seeds: 5, base_seed: 20261125, max_iter: 300, feature_fraction: 1.0}
  hier_cuts: {linkage: [complete, average], min_cos: [0.3, 0.4, 0.5, 0.6, 0.7]}   # added to ensembles in eac
  consensus: {method: eac_average, k: [600, 800, 1000, ..., 3000], select: ari_median}
  eval: {pool_r2: 0.15, cos_thresholds: [0.5, 0.7, 0.9], match_cos: 0.9,
         reproducible: {cos: 0.5, min_frac: 0.4, jac_null_q: 0.999, null_pairs: 3, null_seed: 3},
         headline: [n_lvs, median_overlap_to_consensus, frac_redundant_c0.5, recall_c0.9],
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
cev_zcache_check              Z_all.npz keys vs ancient(legacy noise_flags.csv) -> checks/zcache.ok
cev_flags                     ancient(legacy noise_flags.csv) -> noise/flags_r{r2}.csv  (re-threshold technical_joint)
cev_reproduce_set             Z_all + flags -> reproduce/{set}/{metaZ.csv, cluster_summary.csv, ...}
                              (unchanged cluster_metalv.py cluster; the reference that cev_hier_set must match)
cev_linkage                   Z_all + flags -> linkage/r{r2}/{dist}_{linkage}/{subset}.npz  (dist cos|jac from 2.2c; one dendrogram, shared by all cuts)
cev_hier_set                  linkage + Z_all -> sets/{set}/{subset}/{Z.npy, lvs.csv, members.csv, meta.json}
cev_reproduction              reproduce/{set} vs ancient(external Z, legacy members) -> report/reproduction/{set}.csv
                              (fails unless max|ΔZ| < 1e-6, names equal and member sets equal)
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
  - One cosine-distance-plus-linkage build takes about 12 s at n = 15k (8 s distance, 3 s linkage, measured alone), which is why linkages are cached per (r2, dist, linkage, subset). With 5 linkages running at once (2.2c run) cosine took about 67 s and top-gene Jaccard 30 s (average) / 44 s (complete), including `top_genes`; the jac linkage stays within the 16 GB linkage request.
  - k-means on 18k dimensions takes about 1 s per iteration at k = 3000 on CPU, so use cuML on the GPU; PCA-256 is about 70× cheaper.
  - AMI with thousands of clusters is slow, so use ARI as the main criterion.
  - Split-half / LORO (2.4): 155 subset linkages + 155 batched subset builds + 73 compares took about 18 min locally (4.8 GB peak). The subset sets take 133 GB, because halves keep clusters at 2–3 of 5 runs. Only the compare outputs are needed downstream, so the subset `Z.npy` files can be deleted (Snakemake would rebuild them if asked).
- **Rescaling.** `mean_l2`, `medoid` and `rank1` change the loading scale. ORA uses ranks only, but projection uses L2 from the reference, so rescale to the median member norm.
- **Size.** Set sizes differ (618 vs 1,058 vs the ensemble sets), and ORA % and coverage depend on size. Report size-matched results: top-N by n_models, plus random subsets (via `run_ora.R --lv-subset` / lookups into the shared ORA results).
- **GTEx notebooks.** Decide before Phase 5 whether GTEx gets sibling notebooks under `nbs/03_model_biology/01_gtex/` or executed-notebook names that include `{compendium}`.
- **GTEx seeds.** Don't append the extra seeds to `consensus.seeds`: ARCHS4 shares that list, and `con_models("archs4")` would point at missing models.

---

## Verification

- `python -m unittest discover -s scripts/consensus_eval/tests` (in `clamp-analyses`). The tests cover:
  - `coassoc` against PhenoPLIER's slow version.
  - `eac` recovering a planted partition.
  - `hier` + `median` reproducing `cluster_metalv.py cluster` exactly (and `cev_hier_set` for the two reproduce sets equal to `reproduce/{set}/metaZ.csv` and `cluster_summary.csv` members).
  - The weighting in `run_median`.
  - `rank1` being non-negative.
  - Purity, redundancy, coverage, reciprocity and `overlap_to_consensus` on hand-computed cases (Iq/silhouette were dropped in 2.2b).
  - The Jaccard distance against `scipy.spatial.distance.pdist(metric="jaccard")` on boolean toys (2.2c).
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

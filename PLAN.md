# PLAN: Perturb-seq regulator→target benchmark for OneK1K LVs

> **How to resume:** open this file and find the first unchecked step. Do only that step in a
> session. When it's done, tick its box, record the outcome under "Log", and stop. Every step
> starts with a dry-run or a review: never launch heavy jobs without confirming the job list first.

## Why
The OneK1K (`PBMC_1k1k`) benchmark (Fig. 2A) scores each method's LVs by how well the donor scores
(B) track cell-type fractions. That only rewards **composition**, which is a B-side property.
Perturb-seq gives an independent, causal readout on the **gene side** (Z): are an LV's top genes a
co-regulated program, i.e. enriched for the targets of one CRISPRi-perturbed regulator?

It is also non-circular for the prior-based methods. CLAMPfull, PLIER and MOFA-FLEX train on
GO BP, so GO-based ORA would favour them; Perturb-seq targets don't. Nothing in the repo tests
pseudobulk LVs against regulatory or perturbation data yet. The closest existing work is the ARCHS4
L1000 and CRISPR-screen analyses (`archs4_crispercas.smk`, `archs4_drug_diseases.smk`).

## Decisions (agreed 2026-11-25)
- Perturb-seq: **K562 genome-wide (K562GW) only**. It is the only screen with p-values.
- Data: **OneK1K only**. That means the main pseudobulk (mean CPM, 10 methods) plus donor-bulk
  (summed counts, CLAMPfull only).
- Metrics: **regulon ORA + perturbation-response capture**, each against a permutation null.
- OneK1K fits are **not** in the local `output/`. They will be **fitted locally** in the new
  worktree.
- Work happens in a new worktree off **`origin/main`**: `../clamp-analyses-perturbseq`, branch
  `perturbseq-benchmark`.
- Fresh worktrees have empty `output/` and `data/`; they are gitignored.

## Key facts (from exploration, so we don't redo it)

### Perturb-seq source
Folder: `/home/miltondp/projects/pivlab/other/lv-qtl-exploration/data/ota2025`.
- **Release:** Ota et al., Nature 2025, doi:10.1038/s41586-025-09866-3. Zenodo record 17174974 v3,
  file `m_ota_causal_geneRegulation_paper_v3.zip` (7.7 GB, md5
  `148ecd79c5601783ae0db7726f0c356a`). Unzipped, it is about 13 GB under `extracted/`. A local
  data dictionary is in `data/ota2025/CLAUDE.md`.
- **What it contains:** a reprocessing of the Replogle 2022 CRISPRi **cell-line** screens. K562 is
  an erythroleukemia line, not a PBMC. That is a caveat to state in the analysis.
- **Files we need:**
  - `extracted/data/Perturbseq/geneLevel/K562GW/limma_logFC_sum.txt` (1.6 GB) and
    `limma_P_sum.txt` (1.5 GB).
  - Both are tab-separated: first column `GENE` holds **Ensembl IDs, unversioned** (8,248 targets);
    the other columns are **9,866 perturbations named by symbol**.
  - The values are limma-trend log2FC and raw P versus non-targeting controls. There are no
    adjusted p-values, so we BH-correct per column ourselves.
- **ID map:** `extracted/data/gencode_v41_gname_gid_ALL_sorted_onlyID` (Ensembl→symbol).
- **Reusable code:** `lv-qtl-exploration/approaches/05-nf1-programs/src/lvqtl/nf1.py`, with
  `gene_level()` (pyarrow read plus parquet cache), `bh_columns()`, `de_indicator(q=0.05)`,
  `lv_top_indicator()`. In that approach the median perturbation has only 3 DE genes at
  q < 0.05.
- **Prior result there:** its `07_module_mapping.py` found the plain hypergeometric test
  miscalibrated without a peer-knockdown null. That is why both metrics here use nulls.

### OneK1K in this repo
- **Rules:** `workflow/rules/pseudobulk.smk`. Fits run in `clamp_pseudobulk`, `plier_pseudobulk`,
  `pca_nmf_ica_pseudobulk`, `flashier_pseudobulk`, `mofa_flex_pseudobulk`, `gss_pseudobulk` and
  `cogaps_pseudobulk`, aggregated by `full_models_pseudobulk`. The biology aggregate is
  `biology_pseudobulk` (around line 695).
- **Donor bulk:** `workflow/rules/donor_bulk.smk` (`clampfull_donor_bulk`).
- **Config:** `workflow/config/pseudobulk.yaml` (seed 123). k = 130 LVs; GSSig has about k/4 RAVs.
- **Z outputs:** CSV, genes × LVs, **HGNC symbols** (older author-export names), first column is
  the gene.
  - Main pseudobulk: `output/01_model_building/00_pseudobulk/PBMC_1k1k/models/<M>/Z.csv`, where
    `<M>` is one of CLAMPfull, CLAMPbase, PLIER, PCA, NMF, ICA, flashier, MOFA_FLEX_PRIOR, GSS,
    CoGAPS.
  - Donor bulk: `output/01_model_building/00_pseudobulk/donor_bulk/PBMC_1k1k/models/CLAMPfull/Z.csv`.
  - The two strategies have different gene universes: 13,488 genes after filtering versus up to
    35,518 from CELLxGENE.
- **Cell-type benchmark:** `nbs/03_model_biology/00_pseudobulk/00_benchmark.ipynb`, writing
  `output/03_model_biology/00_pseudobulk/00_benchmark/{benchmark_long.csv,
  clampfull_lv_assignments.csv}`.
- **Reusable ORA:** `scripts/archs4/coverage/run_ora.R`.
  - Arguments: `--z --out-dir --dataset --fraction --seed --model --model-manifest --database
    --database-label --database-type gmt --database-path --top-pct --min-size --max-size
    --pvalue-cutoff --qvalue-cutoff --fdr`.
  - It uses **positive loadings only**: top `ceil(top_pct × universe)`, hypergeometric, BH within
    each LV. It writes `enrichment.csv.gz`, `summary.csv` (regulon coverage) and `manifest.json`.
  - The call pattern to copy is in `workflow/rules/consensus.smk:452-477`.
- **Input manifests:** `workflow/config/reference_inputs.json`, downloaded by
  `scripts/inputs/reference_inputs.py` via `workflow/rules/reference_inputs.smk`.

### Module-evaluation metrics (Saelens 2018)
- **Paper:** Saelens, Cannoodt & Saeys, Nat Commun 2018, doi:10.1038/s41467-018-03424-4. It
  benchmarks 42 module-detection methods against known regulatory networks; decomposition (ICA)
  methods came out best.
- **Local copy, untracked:** `lib/module-evaluation/` holds `paper.pdf`, the supplement and the
  code (`code/lib/modulescomparison.py`).
- **Scores:**
  - For E. coli, yeast and synthetic data: recovery/relevance (best-match Jaccard in each
    direction) and recall/precision (extended BCubed).
  - For human data, whose regulatory gold standard is incomplete like ours: **aucodds**
    (`modbindevalscorer`, `cal_aucodds`).
    - For each module × regulator pair: a 2×2 table over the genes, the odds ratio
      (tp·tn)/(fp·fn) and the Fisher right-tail p.
    - Holm–Šidák correction within each module, across regulators.
    - Only pairs with q < 0.05 keep their OR; the rest are set to 0. The code uses 0.05; the paper
      text says 0.1.
    - Per regulator: the maximum OR over modules.
    - aucodds: the fraction of regulators with log10 OR ≥ c, integrated (trapz) over
      c ∈ linspace(0, 3, 100) and divided by 3.
  - **F-aucodds** (`cal_faucodds`): the harmonic mean of aucodds over regulators and aucodds over
    modules (each module's maximum OR). The second half penalizes hub and uninformative modules.
- Modules smaller than 5 genes are dropped.
- Every score is reported as a fold over the same score for permuted known modules (500
  permutations). Our Z-row permutation is the equivalent.
- **How decomposition methods get crisp modules:** genes above 3 SD per component
  (`_ica_zscore`), fdrtool q per component (optionally split by sign), or a percentage cutoff per
  tail (`code/lib/methods/clustering.py`).
- **Decision (2026-11-27):**
  - Adopt aucodds and F-aucodds only. Recovery/relevance is capped when the sizes differ (a
    75-gene LV module vs GATA1's 550 targets), and BCubed needs its Cython ported.
  - Define LV modules two ways, as a sensitivity check: the top 1% per tail (the ORA query) and
    loadings beyond 3 SD per tail.

## Method spec
- **Universe:** Z genes ∩ K562GW genes, mapping Ensembl→symbol via the gencode v41 file. Log the
  overlap, and add HGNC alias rescue only if fewer than about 85% of genes match.
- **Regulons:**
  - BH-adjust within each perturbation column; targets have q < 0.05.
  - Drop the perturbed gene itself, then keep regulons of size 10–1000 after intersecting with the
    universe. The cap was raised from 500 in Step 4 so that GATA1 (550) and MYB (735) stay in.
  - Collapse near-duplicates (Jaccard > 0.5, which catches the ribosomal and essential knockdowns)
    into representatives, and keep the cluster map.
- **Metric 1, regulon ORA:** write the regulons as a GMT and run `run_ora.R` unchanged.
  - Run on Z and on −Z (`<M>` / `<M>_neg`). That covers sign-arbitrary methods (PCA, ICA,
    flashier, GSSig, MOFA-FLEX); nonnegative methods simply have an empty negative tail.
  - Report the % of LVs enriched for ≥ 1 regulon and the % of regulons recovered by ≥ 1 LV.
- **Metric 2, perturbation-response capture:**
  - For each perturbation with ≥ 20 DE genes, compute max |Pearson r| between its logFC vector and
    any LV's loadings over the universe.
  - Report the median of that max-|r| and the % of perturbations above the null 95th percentile.
- **Null (both metrics):** 100 gene-label permutations of Z rows, using the same permutation for
  all LVs, with a seed from config. This preserves each method's k and loading distribution.
  Report observed − null.
- **Stratification:**
  - Split the LVs into those assigned to a cell type (CLAMPfull, from
    `clampfull_lv_assignments.csv`) and the rest. This separates composition from regulation.
  - `clampfull_lv_assignments.csv` comes from `benchmark_pseudobulk`, which needs all 6 datasets.
    Use a OneK1K-only assignment rule instead (see Step 7).
  - Report how much the regulons overlap the GO BP prior (Jaccard).

## Steps (one per session)

- [x] **Step 1: worktree setup.**
  - Worktree `/home/miltondp/projects/clamp/clamp-analyses-perturbseq`, branch
    `perturbseq-benchmark`, created from `origin/main` (not from `consensus-lvs`).
  - `CLAUDE.md` was copied in from the main checkout. This `PLAN.md` lives only in this worktree;
    the main checkout's `PLAN.md` belongs to another project.
- [x] **Step 2: fit OneK1K models locally.**
  - Target list (dry-run first):
    `snakemake --profile workflow/profiles/local output/01_model_building/00_pseudobulk/PBMC_1k1k/models/{CLAMPfull,CLAMPbase,PLIER,PCA,NMF,ICA,flashier,MOFA_FLEX_PRIOR,GSS,CoGAPS}/Z.csv output/01_model_building/00_pseudobulk/donor_bulk/PBMC_1k1k/models/CLAMPfull/Z.csv`
    (16 jobs).
  - `benchmark_long.csv` was **dropped** from the targets: `benchmark_pseudobulk` expands over all
    6 datasets (about 45 GB of raw data, Heart's manual download, days of fits). OneK1K mean r was
    instead checked with a scratch R script that reuses `score_per_ct()` from
    `scripts/pseudobulk/common.R`.
  - Done: every Z.csv exists and all 10 methods reproduce the recorded OneK1K mean r (see Log).
- [x] **Step 3: Perturb-seq input.**
  - The asset went into a **new per-area manifest, `workflow/config/perturbseq_inputs.json`**,
    not `reference_inputs.json`. That manifest's fetcher only handles small single files, and its
    `reference_inputs` target lists every asset, so it would have pulled the 7.7 GB zip.
  - Script `scripts/inputs/perturbseq_inputs.py`, tested in `tests/test_perturbseq_inputs.py`:
    - Downloads the zip, checking md5 and bytes (it reuses `obtain_archive()`).
    - Extracts only the manifest's members and verifies each one's sha256.
    - `--source-root` symlinks an already-extracted copy after verifying it; `--check` only
      verifies.
    - Existing files are never replaced.
  - Rule `ota2025_k562gw` lives in `workflow/rules/perturbseq.smk`. Parameters for Steps 4–7 are
    in `workflow/config/perturbseq.yaml`. The Snakefile includes both.
  - For development, run with `--config perturbseq_source_root=<lv-qtl-exploration>/data/ota2025/extracted`.
- [x] **Step 4: regulons.**
  - Write `scripts/perturbseq/build_regulons.py` (port the BH and loader logic from `nf1.py`) plus
    unit tests in `scripts/perturbseq/tests/`, using a tiny synthetic matrix.
  - Add rule `perturbseq_regulons`, writing to
    `output/03_model_biology/00_pseudobulk/07_perturbseq/regulons/{main,donor_bulk}/` (GMT,
    `logfc.tsv.gz` subset, cluster map, sizes, `report.json`). It is tsv.gz rather than parquet
    because `clamp-analyses` has no pyarrow.
  - Review the regulon size distribution and the ID overlap with the user.
- [x] **Step 5: regulon ORA.**
  - Write `scripts/perturbseq/signed_z.py` and add rules `perturbseq_signed_z` and
    `perturbseq_ora`, which call `run_ora.R` per method, tail and strategy.
  - Add a permutation-null ORA on shuffled Z.
- **Step 6: gene-side response and module metrics.**
  - [x] **6.1 Response capture.**
    - Write `scripts/perturbseq/response_capture.py` (observed values plus permutation null, with
      `--seed`), with tests.
    - Add rule `perturbseq_response`.
  - [x] **6.2 Remove the shared knockdown response.** No peer-knockdown null (decided 2026-11-27).
    - `response_capture.py --remove-pcs` projects the top-n logFC components out before
      correlating. Config `response.remove_pcs`; currently 0 (raw) and 1.
    - New hub statistic `top_lv_share`.
  - [x] **6.3 Review Saelens et al. 2018** and choose metrics to reuse. See "Module-evaluation
    metrics (Saelens 2018)" under Key facts.
  - [ ] **6.4 aucodds / F-aucodds.**
    - New `scripts/perturbseq/module_scores.py` (with tests) and rule `perturbseq_module_scores`
      per strategy and model, plus an aggregate `module_scores_summary.csv`.
    - Reuse `ora_null.Ora`: its sparse regulon membership, annotated universe (N) and `query()`.
    - LV modules, per LV and tail, two definitions:
      - (a) the top `ceil(top_pct × universe)` positive loadings, i.e. the ORA query;
      - (b) loadings > mean + 3 SD (pos) or < mean − 3 SD (neg) of the LV.
      - Both drop modules smaller than 5 genes.
    - 2×2 tables within the annotated universe. The Fisher right tail is the hypergeometric sf, as
      in `Ora`. Holm–Šidák runs within each module across regulons, and only pairs with q < 0.05
      keep their OR. Use OR = inf when fp·fn = 0 and tp > 0.
    - Report aucodds over regulons, aucodds over modules and F-aucodds, per definition.
    - Null: the same 100 Z-row permutations (seed 123). Report observed, null mean,
      observed − null and fold (observed / null mean, the paper's normalization).
    - Config: `perturbseq.modules: {zscore_sd: 3, min_size: 5, q: 0.05, log10_or_max: 3}`.
    - Check against the paper's code: on a toy example, aucodds must equal `cal_aucodds` /
      `cal_faucodds` (a pure-numpy port; no need to build their Cython).
- [ ] **Step 7: report.**
  - Notebook `nbs/03_model_biology/00_pseudobulk/07_perturbseq.ipynb`, run by rule
    `perturbseq_report`. It shows per-method observed − null for regulon ORA, response capture
    (raw and with the shared response removed) and aucodds / F-aucodds, next to the cell-type
    mean r, plus the top CLAMPfull LV↔regulon hits with their cell-type assignment.
  - Add aggregate `perturbseq_pseudobulk` and hook it into `biology_pseudobulk`.
  - The stratification needs OneK1K CLAMPfull LV↔cell-type assignments. Add a small OneK1K-only
    rule that reuses `assign_with_margins()` from `scripts/pseudobulk/common.R` on
    `PBMC_1k1k/models/CLAMPfull/B.csv` and `PBMC_1k1k/pseudobulk/truthFrac_v0.csv`. Don't depend
    on `benchmark_pseudobulk`, which expands over all 6 datasets.
  - Sanity checks:
    - The null hit rate is about the FDR level.
    - Known K562 programs (ribosome/mito, GATA1/TAL1) map to interpretable LVs.
- [ ] **Step 8: discuss results and decide on extensions.** Candidates:
  - Jurkat for metric 2.
  - RPE1 and HepG2 as non-blood negative controls.
  - PBMC_Perez2022 replication.
  - Regulator co-membership: is the regulator in its LV's top genes?
  - Primary CD4+ T-cell Perturb-seq.
  - Whether the result goes into the manuscript.

## Log
- 2026-11-25/26: Explored both repos and agreed the decisions above. Created the worktree from
  `origin/main`, copied `CLAUDE.md` and wrote this plan. No analysis run. Next: Step 2.
- 2026-11-26: **Step 2 done.**
  - Envs (created by another session) pass `scripts/check_setup.py`.
  - All 16 jobs succeeded in 2 h 15 min of wall time. The slowest were CoGAPS (1 h 10 min),
    MOFA-FLEX (1 h 04 min), PLIER (44 min), flashier (23 min), CLAMP (8 min) and donor-bulk
    CLAMPfull (4.6 min).
  - Z shapes (genes × LVs), all with symbol row names:
    - Main pseudobulk: 13,488 genes. LVs: 130 for most methods, CoGAPS 106, GSS 32,
      MOFA-FLEX 662.
    - Donor bulk CLAMPfull: 13,733 × 130.
  - The built `truthFrac_v0.csv` is identical to `references/pseudobulk/PBMC_1k1k/`.
  - OneK1K mean r (9 cell types), observed vs recorded:

    | Method | Observed | Recorded |
    |---|---|---|
    | CLAMPfull | 0.7533 | 0.7533 |
    | CLAMPbase | 0.7374 | 0.7374 |
    | PLIER | 0.7170 | 0.7164 |
    | CoGAPS | 0.7020 | 0.7020 |
    | Flashier | 0.6628 | 0.6628 |
    | NMF | 0.5935 | 0.5934 |
    | PCA | 0.3900 | 0.3900 |
    | ICA | 0.3162 | 0.3134 |
    | GSSig | 0.2728 | 0.2728 |
    | MOFA-FLEX | 0.2099 | 0.2099 |

  - Every method is within 0.003 of its recorded value.
  - Next: Step 3.
- 2026-11-26: **Step 3 done.**
  - The local zip matches Zenodo: 7,730,808,463 B, md5 `148ecd79…`, licence CC-BY-4.0. Its
    sha256 is `eff4ce5c…b739db`.
  - Member sha256 values are recorded in `perturbseq_inputs.json`.
  - I staged the inputs by linking the verified files from `lv-qtl-exploration` with
    `--source-root`. The 7.7 GB download and extraction path was **not** run end-to-end; only the
    synthetic-zip unit tests cover it.
  - `data/perturbseq/ota2025/K562GW/limma_logFC_sum.txt` is 9,866 perturbations × 8,248 genes.
    `--check` passes, and a re-run dry-run reports nothing to do.
  - `reference_inputs` dry-runs as before, with no Ota job. All 41 tests pass.
  - Next: Step 4.
- 2026-11-26: **Step 4 done.**
  - `scripts/perturbseq/build_regulons.py` (10 unit tests in `scripts/perturbseq/tests/`) and rule
    `perturbseq_regulons`, with wildcard `{strategy}` ∈ {main, donor_bulk}. The strategies are in
    `perturbseq.yaml`. Each job takes about 1 min.
  - BH runs over all 8,248 measured targets in each column, before the universe restriction. There
    are no NaN p-values. 8 targets are unmapped and 2 rows share a symbol; those are dropped.
  - Universe, with no alias rescue needed:
    - main: 7,486 genes (90.9% of mapped targets, 55.5% of the 13,488 Z genes). All 10 Z files have
      the same genes; the rule checks this.
    - donor_bulk: 7,763 genes (94.2% of mapped targets).
  - Knockdown efficacy: 7,573 perturbations measure their own transcript, and 75.9% of those have it
    DE. The median own log2FC is −0.16.
  - DE counts across the whole screen are very skewed: median 2, 90th percentile 208, 99th
    percentile 1,297, max 2,949 (MED19).
  - **The size cap went from 500 to 1000 (agreed with the user).** With 500, GATA1 (550 in the
    universe), MYB (735) and 413 other knockdowns were dropped.
  - main, with the 1000 cap:
    - 2,699 regulons in [10, 1000]. 3,575 are below 10 and 147 are above 1000 (MED, ribosome and
      translation knockdowns).
    - Regulon size: median 82, IQR 27–242.
  - The Jaccard > 0.5 collapse barely does anything: 2,699 → 2,685 representatives, 7 clusters with
    more than one member. The largest are NOL7 (5) and TSR1 (4), both ribosome biogenesis. Big
    essential knockdowns are not near-duplicates of each other.
  - donor_bulk: 2,787 regulons → 2,767 representatives.
  - Sanity regulators, main:
    - Kept: GATA1 (550), TAL1 (212), LMO2 (295), LDB1 (448), ZFPM1 (18), SPI1 (35), MYC (50),
      RUNX1 (12).
    - Too small: KLF1 (6), NFE2 (0), TP53 (1).
  - Response matrix for Step 6: 2,353 perturbations with ≥ 20 DE genes in the universe (2,417 for
    donor bulk).
  - **Note for Step 6:** `logfc.tsv.gz` keeps each perturbation's own row. Decide there whether to
    mask it before computing correlations.
  - Next: Step 5.
- 2026-11-27: **Step 5 done.**
  - `scripts/perturbseq/signed_z.py` restricts Z to the regulon universe and writes Z (`pos`), −Z
    (`neg`) and a `samples` manifest. That's needed because `run_ora.R` uses every Z row as its
    universe.
  - `scripts/perturbseq/ora_null.py` reimplements clusterProfiler's enricher, vectorized over LVs
    (decision 2026-11-27). The background is the annotated universe; only regulons with k ≥ 1 are
    tested, and BH runs over those, within each LV. Unpermuted, it must reproduce `run_ora.R`'s
    hits and p-values (rtol 1e-6) or the job fails. **All 11 models pass on both tails.**
  - The null is 100 label permutations from seed 123, shared across LVs, tails and models.
  - Tests: 10 new, in `scripts/perturbseq/tests/`, including an end-to-end comparison against the
    real `run_ora.R`. 20 pass.
  - Rules: `perturbseq_signed_z`, `perturbseq_ora` (`run_ora.R`), `perturbseq_ora_null` and
    `perturbseq_ora_summary`, which writes `07_perturbseq/ora_summary.csv`. There are 45 jobs. The
    slowest `run_ora.R` job took 10 min; the nulls take under 1.5 min each.
  - Two snags:
    - `projections.smk` constrains `{model}` globally, so the new rules set their own `model`
      constraint.
    - `run_ora.R` writes a **truncated header-only gzip** when no LV has hits, because data.table's
      `fwrite(compress="gzip")` leaves it unterminated. `ora_null.py` therefore uses `summary.csv`'s
      `result_rows == 0` instead of reading the file, and Step 7 readers must do the same.
      `run_ora.R` itself is unchanged, since ARCHS4 shares it.
  - **The negative tail is empty for CLAMPfull, CLAMPbase, PLIER, NMF and CoGAPS** (nonnegative Z).
    It has hits only for ICA, PCA, flashier, MOFA-FLEX and GSS.
  - The null is well calibrated: about 2.5–3% of LVs hit ≥ 1 regulon per tail (about 5% for either
    tail). Every method beats its null at perm p = 1/101, the floor.
  - Main pseudobulk, % of LVs enriched in either tail, observed vs null mean (Δ):

    | Method | Observed | Null mean | Δ |
    |---|---|---|---|
    | ICA | 59.2 | 5.2 | 54.0 |
    | PLIER | 38.5 | 3.1 | 35.4 |
    | CLAMPbase | 31.5 | 3.0 | 28.6 |
    | GSSig | 31.3 | 5.5 | 25.8 |
    | CoGAPS | 28.3 | 2.8 | 25.5 |
    | CLAMPfull | 27.7 | 2.6 | 25.1 |
    | NMF | 26.2 | 2.7 | 23.5 |
    | MOFA-FLEX | 26.1 | 5.2 | 20.9 |
    | PCA | 23.8 | 5.2 | 18.7 |
    | flashier | 23.1 | 5.2 | 17.9 |

    Donor-bulk CLAMPfull: 35.4, null 2.5, Δ 32.9. That's higher than main CLAMPfull.
  - **The "% of regulons recovered" metric is saturated** at 74–97% for every method, against a null
    of 1–10%. It doesn't discriminate between methods.
    - The cause is a few promiscuous LVs. CLAMPfull LV37 hits 1,915 of 2,685 regulons and LV119
      hits 1,857. They likely capture K562's shared knockdown response (erythroid, proliferation).
    - LV37 is also the top LV for GATA1 (p.adj 1.6e-78), TAL1, LMO2, LDB1 and MYB, which is
      interpretable but not specific.
    - **For Step 7:** add a specificity view, such as hits per LV, a regulon's best LV, or a hub-LV
      exclusion, before reading regulon recovery.
  - Next: Step 6.
- 2026-11-27: **Step 6 done.**
  - `scripts/perturbseq/response_capture.py` (10 new tests; 30 pass in `scripts/perturbseq/tests/`)
    and rules `perturbseq_response` (11 jobs, about 40 s in total) and `perturbseq_response_summary`
    (`07_perturbseq/response_summary.csv`). `ora_null.summarize()` now takes a statistic list.
  - Each perturbation's own gene is left out of its correlations, exactly (2,003 of 2,353 main
    perturbations; 2,122 of 2,417 donor-bulk). The null uses the same Z row permutations as the ORA
    (seed 123, 100 perms). A perturbation is captured when its perm p ≤ 0.05. The null capture rate
    is 4.0%, not 5%, because the observed value sits in each null member's reference set.
  - Median max |r| over LVs, observed vs null mean (Δ):

    | Method | k | Observed | Null | Δ |
    |---|---|---|---|---|
    | PLIER | 130 | 0.183 | 0.032 | 0.151 |
    | CLAMPfull | 130 | 0.169 | 0.035 | 0.134 |
    | CLAMPbase | 130 | 0.158 | 0.032 | 0.126 |
    | NMF | 130 | 0.155 | 0.031 | 0.124 |
    | flashier | 130 | 0.149 | 0.033 | 0.117 |
    | ICA | 130 | 0.128 | 0.031 | 0.097 |
    | CoGAPS | 106 | 0.128 | 0.031 | 0.097 |
    | PCA | 130 | 0.127 | 0.032 | 0.095 |
    | MOFA-FLEX | 662 | 0.127 | 0.043 | 0.084 |
    | GSSig | 32 | 0.085 | 0.027 | 0.058 |

    Donor-bulk CLAMPfull: 0.180, null 0.035, Δ 0.144. As in the ORA, that's higher than main
    CLAMPfull.
  - The ranking differs from the ORA's. ICA leads the ORA (Δ 54) but is mid-pack here, while PLIER
    and CLAMP lead here.
  - **"% captured" is saturated** at 92–100% for every method (null 4%). The gene-label null only
    tests "any correlation at all", which r ≈ 0.03 over 7.5k genes makes easy.
  - **A common knockdown response dominates.** The first PC of the logFC matrix explains 22.7% of
    its variance. In most methods one hub LV is the best LV for 50–62% of perturbations: CLAMPfull
    LV37 (60%), CLAMPbase LV8, PLIER LV128, PCA PC1, NMF LV16, flashier LV1. Those hubs correlate
    0.29–0.54 with logFC PC1. ICA (IC66, 17%; 86 distinct best LVs) and MOFA-FLEX (26%) spread
    more.
  - Sanity regulators in CLAMPfull: GATA1 (r −0.46), LDB1 (−0.42), MYB, TAL1, LMO2 and MYC all
    pick LV37. SPI1 is barely above null (0.058 vs p95 0.046). That is the same LV37 story as the
    ORA: interpretable, but not specific.
  - **Open decision for Step 7:** remove the shared response before reading this metric. Options:
    regress the top logFC PCs out of the responses, or use a peer-perturbation null, which is what
    `lv-qtl-exploration`'s `07_module_mapping.py` found necessary. Otherwise report a hub-excluded
    view. Median max |r| is usable as is; "% captured" is not.
  - Next: Step 7.
- 2026-11-27: **Step 6 split into substeps; 6.2 and 6.3 done.**
  - Step 6.1 (response capture) was committed as is. The user chose to remove the top shared
    component instead of using a peer-knockdown null.
  - 6.2:
    - `residualize()` centers each perturbation with a mean that leaves its own gene out, zeroes
      the own-gene entry, and projects the top-n left singular vectors out. A unit test caught an
      own-gene leak through the column mean.
    - `top_lv_share` is the % of perturbations whose best LV is the most common one.
    - 34 tests pass. The forced rerun of 11 jobs took 1.6 min.
    - `remove_pcs = 0` reproduces Step 6.1 exactly (difference 0).
  - **Removing PC1 (22.8% of the logFC variance) does not break up the hubs:**

    | Method | Median max \|r\| Δ, raw → PC1 removed | top_lv_share, raw → PC1 removed |
    |---|---|---|
    | PLIER | 0.151 → 0.123 | 60.6 → 57.8 |
    | CLAMPfull | 0.134 → 0.109 | 60.4 → 63.6 |
    | CLAMPbase | 0.126 → 0.096 | 62.1 → 53.6 |
    | NMF | 0.124 → 0.094 | 51.6 → 37.1 |
    | flashier | 0.117 → 0.094 | 49.9 → 41.4 |
    | ICA | 0.097 → 0.079 | 17.4 → 24.9 |
    | CoGAPS | 0.097 → 0.074 | 27.8 → 20.3 |
    | PCA | 0.095 → 0.071 | 56.8 → 46.5 |
    | MOFA-FLEX | 0.084 → 0.068 | 25.9 → 20.1 |
    | GSSig | 0.058 → 0.045 | 40.3 → 29.0 |

    Donor-bulk CLAMPfull: Δ 0.144 → 0.100; top_lv_share 54.4 → 49.6. The method ranking does not
    change.
  - "% captured" stays saturated (87–100%); the null is still 4.0%.
  - CLAMPfull with PC1 removed: GATA1 and LDB1 still pick LV37, but with weaker r (−0.19 and
    −0.15, from −0.46 and −0.42). TAL1, LMO2, MYB and MYC move to other LVs (LV1, LV7, LV108,
    LV65).
  - **The shared response is at least two-dimensional.** A scratch run (5 perms, not in `output/`)
    removed 0/1/2/3/5/10 PCs, which take out 0/23/33/38/43/49% of the variance. Top-LV share:

    | Method | 0 | 1 | 2 | 3 | 5 | 10 |
    |---|---|---|---|---|---|---|
    | CLAMPfull | 60 | 64 | 23 | 19 | 18 | 10 |
    | PLIER | 61 | 58 | 24 | 30 | 33 | 18 |
    | ICA | 17 | 25 | 10 | 9 | 12 | 9 |

    Median max |r| drops steadily: CLAMPfull 0.17 → 0.14 → 0.09 → 0.08 → 0.07 → 0.06. With PCs 1–2
    removed, CLAMPfull's best LVs go from 34 distinct to 44, and PLIER's from 30 to 49.
  - **Open decision:** add `remove_pcs: 2` (and maybe 3) to the config? One PC is not enough to
    remove the shared response.
  - Next: decide `remove_pcs`, then Step 6.4.

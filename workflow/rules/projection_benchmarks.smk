import json
import os
from pathlib import Path

PBENCH = config['projection_benchmark']
PBENCH_ROOT = config.get('projection_benchmark_root', PBENCH['root'])
PBENCH_PY = config.get('projection_benchmark_python', 'python')
PBENCH_R = config.get('projection_benchmark_rscript', 'Rscript')
PBENCH_SCRIPT = 'scripts/archs4/projections/benchmark.py'
PBENCH_DATASETS = PBENCH['datasets']
PBENCH_SEEDS = PBENCH['null_seeds']
PBENCH_METHODS = PBENCH['methods']
PBENCH_SOURCE = config['projections']['archs4_model']
PBENCH_MODEL = os.path.join(PBENCH_SOURCE['root'], PBENCH_SOURCE['rds'])
PBENCH_PB_PROD = config['paths']['production']
PBENCH_GT_PROD = config['gtex']['paths']['production']
PBENCH_GT_BIO = config['gtex']['paths']['biology']
PBENCH_COHORT = config['gtex']['kmeans_clustering']
PBENCH_SCORE_CONDITIONS = PBENCH_METHODS + ['projection'] + [f'null_{s}' for s in PBENCH_SEEDS]
PBENCH_CLUSTER_KEYS = [(k, rep) for k in PBENCH['ks'] for rep in ([0] if k == PBENCH['n_lvs'] else range(PBENCH['n_subsets']))]


def pbench_truth(ds):
    cfg = config['datasets'][ds]
    directory = f'{PBENCH_PB_PROD}/{ds}/pseudobulk' if cfg.get('build_pseudobulk') else cfg['pseudobulk_dir']
    return f'{directory}/truthFrac_v0.csv'


def pbench_expression(ds):
    return f'{PBENCH_GT_PROD}/gtex_fbm_filt.rds' if ds == 'gtex' else f'{PBENCH_PB_PROD}/{ds}/preprocessing/norm.csv'


def pbench_matrix(w):
    if w.condition in PBENCH_METHODS:
        return ancient(f'{PBENCH_PB_PROD}/{w.dataset}/models/{config["methods"][w.condition]}')
    seed = 0 if w.condition == 'projection' else int(w.condition.removeprefix('null_'))
    name = 'real' if seed == 0 else str(seed)
    return f'{PBENCH_ROOT}/projection/{w.dataset}/{name}.rds'


# ============================================================
# Step 1: prepare real and permuted projections
# ============================================================

rule projection_benchmark_model:
    input:
        model=ancient(PBENCH_MODEL),
        script='scripts/archs4/projections/benchmark_project.R',
        common='scripts/archs4/common.R',
    output:
        model=f'{PBENCH_ROOT}/model.rds',
        genes=f'{PBENCH_ROOT}/model_genes.tsv',
        manifest=f'{PBENCH_ROOT}/model.rds.json',
    conda: PBENCH['r_environment']
    resources: mem_mb=40000
    shell:
        'OPENBLAS_NUM_THREADS=4 {PBENCH_R} {input.script:q} --stage prepare --model {input.model:q} '
        '--output {output.model:q} --genes {output.genes:q}'


rule projection_benchmark_project:
    input:
        model=rules.projection_benchmark_model.output.model,
        expression=lambda w: ancient(pbench_expression(w.dataset)),
        backing=lambda w: ancient(f'{PBENCH_GT_PROD}/FBMgtex_preproc_filtered.bk') if w.dataset == 'gtex' else [],
        genes=lambda w: ancient(f'{PBENCH_GT_PROD}/gtex_genes.rds') if w.dataset == 'gtex' else [],
        samples=lambda w: ancient(f'{PBENCH_GT_PROD}/gtex_samples.rds') if w.dataset == 'gtex' else [],
        script='scripts/archs4/projections/benchmark_project.R',
        common='scripts/archs4/common.R',
    output:
        matrix=f'{PBENCH_ROOT}/projection/{{dataset}}/{{seed}}.rds',
        manifest=f'{PBENCH_ROOT}/projection/{{dataset}}/{{seed}}.rds.json',
    params:
        seed=lambda w: 0 if w.seed == 'real' else int(w.seed),
        extra=lambda w: (f' --genes {PBENCH_GT_PROD}/gtex_genes.rds --samples {PBENCH_GT_PROD}/gtex_samples.rds' if w.dataset == 'gtex' else ''),
    wildcard_constraints:
        dataset='|'.join(PBENCH_DATASETS + ['gtex']), seed='real|'+'|'.join(map(str,PBENCH_SEEDS))
    conda: PBENCH['r_environment']
    threads: 4
    resources: mem_mb=20000
    shell:
        'OPENBLAS_NUM_THREADS={threads} {PBENCH_R} {input.script:q} --stage project --model {input.model:q} '
        '--input {input.expression:q} --seed {params.seed} --output {output.matrix:q}{params.extra}'


# ============================================================
# Step 2: score pseudobulk recovery
# ============================================================

rule projection_benchmark_score:
    input:
        matrix=pbench_matrix,
        truth=lambda w: ancient(pbench_truth(w.dataset)),
        script=PBENCH_SCRIPT,
    output:
        targets=f'{PBENCH_ROOT}/scores/{{dataset}}/{{condition}}/targets.tsv',
        null=f'{PBENCH_ROOT}/scores/{{dataset}}/{{condition}}/targets.tsv.null.npz',
        mean_null=f'{PBENCH_ROOT}/scores/{{dataset}}/{{condition}}/targets.tsv.mean_null.tsv.gz',
        manifest=f'{PBENCH_ROOT}/scores/{{dataset}}/{{condition}}/targets.tsv.json',
    params:
        method=lambda w: w.condition if w.condition in PBENCH_METHODS else ('ARCHS4 projection' if w.condition == 'projection' else 'ARCHS4 projection null'),
        matching=lambda w: 'greedy' if w.condition in PBENCH_METHODS else 'max',
        exclusions=lambda w: ' '.join(__import__('shlex').quote(s) for s in config['cell_type_analysis']['excluded_targets'].get(w.dataset,[])),
        seed=lambda w: PBENCH['permutation_seed'] + PBENCH_DATASETS.index(w.dataset),
        B=PBENCH['permutations'],
    wildcard_constraints: dataset='|'.join(PBENCH_DATASETS), condition='|'.join(PBENCH_SCORE_CONDITIONS)
    conda: PBENCH['python_environment']
    threads: 4
    resources: mem_mb=4000
    shell:
        'OPENBLAS_NUM_THREADS={threads} {PBENCH_PY} {input.script:q} score --matrix {input.matrix:q} '
        '--truth {input.truth:q} --dataset {wildcards.dataset} --method {params.method:q} --matching {params.matching} '
        '--seed {params.seed} --permutations {params.B} --exclude {params.exclusions} --output {output.targets:q}'


# ============================================================
# Step 3: evaluate GTEx clustering
# ============================================================

rule projection_benchmark_embedding:
    input:
        matrix=lambda w: f'{PBENCH_ROOT}/projection/gtex/{w.seed}.rds',
        metadata=ancient(config['gtex']['metadata']),
        script=PBENCH_SCRIPT,
    output:
        embedding=f'{PBENCH_ROOT}/gtex/{{seed}}.npy',
        metadata=f'{PBENCH_ROOT}/gtex/{{seed}}.npy.metadata.json',
        manifest=f'{PBENCH_ROOT}/gtex/{{seed}}.npy.json',
    params: minimum=PBENCH_COHORT['ari_min_samples'], clusters=PBENCH_COHORT['n_clusters']
    wildcard_constraints: seed='real|'+'|'.join(map(str,PBENCH_SEEDS))
    conda: PBENCH['python_environment']
    resources: mem_mb=8000
    shell:
        'OPENBLAS_NUM_THREADS=1 {PBENCH_PY} {input.script:q} embedding --matrix {input.matrix:q} --metadata {input.metadata:q} '
        '--min-samples {params.minimum} --n-clusters {params.clusters} --output {output.embedding:q}'


rule projection_benchmark_subset:
    input: metadata=f'{PBENCH_ROOT}/gtex/real.npy.metadata.json', script=PBENCH_SCRIPT
    output: subset=f'{PBENCH_ROOT}/gtex/subsets/{{k}}/{{rep}}.json'
    params: seed=PBENCH['subset_seed']
    wildcard_constraints: k='|'.join(map(str,PBENCH['ks'])), rep=r'\d+'
    conda: PBENCH['python_environment']
    shell:
        '{PBENCH_PY} {input.script:q} subset --metadata {input.metadata:q} --k {wildcards.k} --rep {wildcards.rep} '
        '--seed {params.seed} --output {output.subset:q}'


rule projection_benchmark_cluster:
    input:
        embedding=f'{PBENCH_ROOT}/gtex/real.npy',
        metadata=f'{PBENCH_ROOT}/gtex/real.npy.metadata.json',
        subset=rules.projection_benchmark_subset.output.subset,
        script=PBENCH_SCRIPT,
    output:
        scores=f'{PBENCH_ROOT}/gtex/clusters/{{k}}/{{rep}}.csv',
        manifest=f'{PBENCH_ROOT}/gtex/clusters/{{k}}/{{rep}}.csv.json',
    params: clusters=PBENCH_COHORT['n_clusters'], reps=PBENCH_COHORT['n_reps_per_k'], seed=PBENCH_COHORT['base_seed']
    conda: PBENCH['python_environment']
    resources: gpu=1, mem_mb=4000
    shell:
        'OPENBLAS_NUM_THREADS=1 {PBENCH_PY} {input.script:q} cluster --embedding {input.embedding:q} --metadata {input.metadata:q} '
        '--subset {input.subset:q} --n-clusters {params.clusters} --repetitions {params.reps} --base-seed {params.seed} --output {output.scores:q}'


rule projection_benchmark_null_cluster:
    input:
        embedding=f'{PBENCH_ROOT}/gtex/{{seed}}.npy',
        metadata=f'{PBENCH_ROOT}/gtex/{{seed}}.npy.metadata.json',
        script=PBENCH_SCRIPT,
    output:
        scores=f'{PBENCH_ROOT}/gtex/null_clusters/{{seed}}.csv',
        manifest=f'{PBENCH_ROOT}/gtex/null_clusters/{{seed}}.csv.json',
    params: clusters=PBENCH_COHORT['n_clusters'], reps=PBENCH_COHORT['n_reps_per_k'], seed=PBENCH_COHORT['base_seed']
    wildcard_constraints: seed='|'.join(map(str,PBENCH_SEEDS))
    conda: PBENCH['python_environment']
    resources: gpu=1, mem_mb=4000
    shell:
        'OPENBLAS_NUM_THREADS=1 {PBENCH_PY} {input.script:q} cluster --embedding {input.embedding:q} --metadata {input.metadata:q} '
        '--n-clusters {params.clusters} --repetitions {params.reps} --base-seed {params.seed} --output {output.scores:q}'


# ============================================================
# Step 4: combine benchmark results
# ============================================================

rule projection_benchmark_settings:
    output: f'{PBENCH_ROOT}/settings.json'
    params: settings=dict(PBENCH, cluster_repetitions=PBENCH_COHORT['n_reps_per_k'], clustering=PBENCH_COHORT)
    run:
        Path(output[0]).parent.mkdir(parents=True,exist_ok=True)
        Path(output[0]).write_text(json.dumps(params.settings,indent=2))


rule projection_benchmark_tables:
    input:
        scores=expand(f'{PBENCH_ROOT}/scores/{{dataset}}/{{condition}}/targets.tsv',dataset=PBENCH_DATASETS,condition=PBENCH_SCORE_CONDITIONS),
        clustering=[f'{PBENCH_ROOT}/gtex/clusters/{k}/{rep}.csv' for k,rep in PBENCH_CLUSTER_KEYS],
        null_clusters=expand(f'{PBENCH_ROOT}/gtex/null_clusters/{{seed}}.csv',seed=PBENCH_SEEDS),
        settings=rules.projection_benchmark_settings.output,
        local=ancient(f'{PBENCH_GT_BIO}/00_kmeans_clustering/ari_data.csv'),
        script=PBENCH_SCRIPT,
    output:
        pseudobulk=f'{PBENCH_ROOT}/pseudobulk_adjusted.tsv',
        gtex=f'{PBENCH_ROOT}/gtex_full_ari.tsv',
        repetitions=f'{PBENCH_ROOT}/gtex_sensitivity_repetitions.tsv',
        subsets=f'{PBENCH_ROOT}/gtex_sensitivity_subsets.tsv',
        summary=f'{PBENCH_ROOT}/gtex_sensitivity_summary.tsv',
        manifest=f'{PBENCH_ROOT}/pseudobulk_adjusted.tsv.json',
    conda: PBENCH['python_environment']
    shell:
        '{PBENCH_PY} {input.script:q} aggregate --root {PBENCH_ROOT:q} --settings {input.settings:q} --local-ari {input.local:q}'


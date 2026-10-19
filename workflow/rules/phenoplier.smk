import os


PHENOPLIER_CFG = config["phenoplier"]
PHENOPLIER_TARGET = config.get("phenoplier_target", PHENOPLIER_CFG["target"])
PHENOPLIER_CLUSTER = PHENOPLIER_CFG["clusters"][PHENOPLIER_TARGET]
PHENOPLIER_WORKERS = PHENOPLIER_CFG["workers"]
PHENOPLIER_GLS_RES = PHENOPLIER_CFG["resources"]["gls"]
PHENOPLIER_STORE_RES = PHENOPLIER_CFG["resources"]["store"]
PHENOPLIER_REPORT = PHENOPLIER_CFG["report_root"]
PHENOPLIER_TRAIT_FILTER = PHENOPLIER_CFG["trait_filter"]
PHENOPLIER_EXPECTED_PHENOTYPES = int(
    PHENOPLIER_CFG["expected_phenotypes"][PHENOPLIER_TRAIT_FILTER]
)


PHENOPLIER_WORKSPACE_ARG = (
    f"--workspace {PHENOPLIER_CFG['workspace']} " if PHENOPLIER_CFG.get("workspace") else ""
)


A4_TRAIT_DIRS = A4_CFG["traits"]
PHENOPLIER_COV_BASE_DIR = A4_TRAIT_DIRS["coverage"]["clampbase_dir"]
PHENOPLIER_SAT_BASE_DIR = A4_TRAIT_DIRS["saturation"]["clampbase_dir"]
PHENOPLIER_FIN_BASE_DIR = A4_TRAIT_DIRS["finals"]["clampbase_dir"]


PHENOPLIER_SAT_KS = [int(k) for k in PHENOPLIER_CFG["saturation_k_values"]]
PHENOPLIER_SAT_K_PATTERN = "|".join(map(str, PHENOPLIER_SAT_KS))
PHENOPLIER_SAT_CELLS = [
    (fraction, k, seed)
    for (fraction, k, seed) in A4_SAT_CELLS
    if k in PHENOPLIER_SAT_KS
]


A4_FINAL_DATASETS = list(A4_COV_CFG["final_datasets"])
A4_FINAL_DATASET_PATTERN = "|".join(A4_FINAL_DATASETS)
A4_FINAL_LINKED = ["archs4"]
A4_FINAL_LINKED_PATTERN = "|".join(A4_FINAL_LINKED)
A4_FINAL_RUN_DATASETS = [d for d in A4_FINAL_DATASETS if d not in A4_FINAL_LINKED]
A4_FINAL_RUN_PATTERN = "|".join(A4_FINAL_RUN_DATASETS)


def a4_final_base_seed(dataset):
    return int(A4_COV_CFG["clampbase_seed"][dataset])


def a4_final_base_rds(dataset):
    if dataset in A4_FINAL_LINKED:
        return f"{a4_cov_cell(100, a4_final_base_seed(dataset))}/CLAMPbase.rds"
    return A4_COV_CFG["clampbase_source"][dataset]["model"]


def phenoplier_final_name(dataset, model):

    if dataset in A4_FINAL_LINKED:
        return f"cov_rs100_seed{a4_final_base_seed(dataset)}_{model}"
    return f"final_{dataset}_{model}"


def phenoplier_store_dir(summary_dir):

    return os.path.join(os.path.dirname(summary_dir), "stores")


PHENOPLIER_COV_BASE = [
    f"{PHENOPLIER_COV_BASE_DIR}/cov_rs{fraction}_seed{seed}.tsv.gz"
    for fraction in A4_COV_LEVELS
    for seed in a4_cov_seeds_for(fraction)
]
PHENOPLIER_SAT_BASE = [
    f"{PHENOPLIER_SAT_BASE_DIR}/sat_rs{fraction}_k{k}_seed{seed}.tsv.gz"
    for (fraction, k, seed) in PHENOPLIER_SAT_CELLS
]
PHENOPLIER_FIN_BASE = [
    f"{PHENOPLIER_FIN_BASE_DIR}/final_{dataset}.tsv.gz" for dataset in A4_FINAL_DATASETS
]
PHENOPLIER_CANONICAL = [
    f"{A4_CAN_FINAL_ROOT}/{dataset}/traits/canon_{dataset}.tsv.gz"
    for dataset in A4_CAN_COMPENDIA
]
PHENOPLIER_FIN_STORES = [
    f"{phenoplier_store_dir(PHENOPLIER_FIN_BASE_DIR)}/final_{dataset}_CLAMPbase.h5"
    for dataset in A4_FINAL_DATASETS
]
PHENOPLIER_COMP_BASE_DIR = A4_TRAIT_DIRS["coverage"]["comparators_clampbase_dir"]
PHENOPLIER_COMP_BASE = [
    f"{PHENOPLIER_COMP_BASE_DIR}/gtex_rs100_seed{rep}.tsv.gz" for rep in A4_COV_GTEX_REPLICATES
]


PHENOPLIER_GLS_SHELL = "".join([
    "bash {input.script} --rds {params.rds} --name {params.name} ",
    "--summary-out {output.summary} --n-jobs {threads} ",
    PHENOPLIER_WORKSPACE_ARG,
    f"--conda-env {PHENOPLIER_CFG['conda_env']} ",
    f"--require-version {PHENOPLIER_CFG['version']} ",
    f"--namespace {PHENOPLIER_CFG['namespace']} ",
    f"--cohort {PHENOPLIER_CFG['cohort']} ",
    f"--lv-percentile {PHENOPLIER_CFG['lv_percentile']} ",
    f"--trait-filter {PHENOPLIER_TRAIT_FILTER} ",
    f"--expected-phenotypes {PHENOPLIER_EXPECTED_PHENOTYPES} ",
    f"--executor {PHENOPLIER_CLUSTER['executor']} ",
    f"--cluster '{PHENOPLIER_CLUSTER.get('cluster', '')}' ",
    f"--workers-step6 {PHENOPLIER_WORKERS['step6']} ",
    f"--blas-step6 {PHENOPLIER_WORKERS['blas_step6']} ",
    f"--workers-step7 {PHENOPLIER_WORKERS['step7']} ",
    f"--blas-step7 {PHENOPLIER_WORKERS['blas_step7']} ",
    "> {log} 2>&1",
])

PHENOPLIER_STORE_SHELL = "".join([
    "bash {input.script} --rds {params.rds} --name {params.name} ",
    "--store-out {output.store} ",
    PHENOPLIER_WORKSPACE_ARG,
    f"--conda-env {PHENOPLIER_CFG['conda_env']} ",
    f"--require-version {PHENOPLIER_CFG['version']} ",
    f"--cohort {PHENOPLIER_CFG['cohort']} ",
    "> {log} 2>&1",
])


# ============================================================
# Step 1: test module-trait associations
# ============================================================

rule gls_bp_coverage_base:
    input:
        rds=lambda wc: f"{a4_cov_cell(wc.fraction, wc.seed)}/CLAMPbase.rds",
        script="scripts/phenoplier/run_gls.sh",
        verify="scripts/phenoplier/verify_summary.py",
    output:
        summary=f"{PHENOPLIER_COV_BASE_DIR}/cov_rs{{fraction}}_seed{{seed}}.tsv.gz",
    log:
        f"{PHENOPLIER_COV_BASE_DIR}/cov_rs{{fraction}}_seed{{seed}}.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: f"cov_rs{wc.fraction}_seed{wc.seed}_CLAMPbase",
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        fraction=A4_COV_LEVEL_PATTERN,
        seed=A4_COV_SEED_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_bp_saturation_base:
    input:
        rds=lambda wc: ancient(a4_sat_base_rds(wc.fraction, wc.k, wc.seed)),
        script="scripts/phenoplier/run_gls.sh",
        verify="scripts/phenoplier/verify_summary.py",
    output:
        summary=f"{PHENOPLIER_SAT_BASE_DIR}/sat_rs{{fraction}}_k{{k}}_seed{{seed}}.tsv.gz",
    log:
        f"{PHENOPLIER_SAT_BASE_DIR}/sat_rs{{fraction}}_k{{k}}_seed{{seed}}.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: f"sat_rs{wc.fraction}_k{wc.k}_seed{wc.seed}_CLAMPbase",
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=PHENOPLIER_SAT_K_PATTERN,
        seed=A4_SAT_SEED_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_gtex_coverage_base:
    input:
        rds=lambda wc: f"{a4_cov_gtex_rep_dir(wc.replicate)}/CLAMPbase.rds",
        script="scripts/phenoplier/run_gls.sh",
        verify="scripts/phenoplier/verify_summary.py",
    output:
        summary=f"{PHENOPLIER_COMP_BASE_DIR}/gtex_rs100_seed{{replicate}}.tsv.gz",
    log:
        f"{PHENOPLIER_COMP_BASE_DIR}/gtex_rs100_seed{{replicate}}.log",
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: f"gtex_rs100_seed{wc.replicate}_CLAMPbase",
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        replicate="|".join(map(str, A4_COV_GTEX_REPLICATES)),
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_final_base:
    input:
        rds=lambda wc: a4_final_base_rds(wc.dataset),
        script="scripts/phenoplier/run_gls.sh",
        verify="scripts/phenoplier/verify_summary.py",
    output:
        summary=f"{PHENOPLIER_FIN_BASE_DIR}/final_{{dataset}}.tsv.gz",
    log:
        f"{PHENOPLIER_FIN_BASE_DIR}/final_{{dataset}}.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: phenoplier_final_name(wc.dataset, "CLAMPbase"),
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        dataset=A4_FINAL_RUN_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule link_final_base:
    input:
        summary=lambda wc: (
            f"{PHENOPLIER_COV_BASE_DIR}/cov_rs100_seed{a4_final_base_seed(wc.dataset)}.tsv.gz"
        ),
    output:
        summary=f"{PHENOPLIER_FIN_BASE_DIR}/final_{{dataset}}.tsv.gz",
    wildcard_constraints:
        dataset=A4_FINAL_LINKED_PATTERN,
    shell:
        "cp -f {input.summary} {output.summary}"


rule gls_canonical:
    input:
        rds=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/{A4_CAN_MODEL_NAME}.rds",
        script="scripts/phenoplier/run_gls.sh",
        verify="scripts/phenoplier/verify_summary.py",
    output:
        summary=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/traits/canon_{{dataset}}.tsv.gz",
    log:
        f"{A4_CAN_FINAL_ROOT}/{{dataset}}/traits/canon_{{dataset}}.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: f"canon_{wc.dataset}_{A4_CAN_MODEL_NAME}",
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        dataset=A4_CAN_COMPENDIUM_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


# ============================================================
# Step 2: build final-model study stores
# ============================================================

rule store_final_base:
    input:
        summary=f"{PHENOPLIER_FIN_BASE_DIR}/final_{{dataset}}.tsv.gz",
        model=lambda wc: a4_final_base_rds(wc.dataset),
        script="scripts/phenoplier/store_build.sh",
    output:
        store=f"{phenoplier_store_dir(PHENOPLIER_FIN_BASE_DIR)}/final_{{dataset}}_CLAMPbase.h5",
    log:
        f"{phenoplier_store_dir(PHENOPLIER_FIN_BASE_DIR)}/final_{{dataset}}_CLAMPbase.log"
    params:
        rds=lambda wc: a4_final_base_rds(wc.dataset),
        name=lambda wc: phenoplier_final_name(wc.dataset, "CLAMPbase"),
    resources:
        mem_mb=int(PHENOPLIER_STORE_RES["mem_mb"]),
        runtime=int(PHENOPLIER_STORE_RES["runtime"]),
    wildcard_constraints:
        dataset=A4_FINAL_DATASET_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_STORE_SHELL


# ============================================================
# Step 3: combine trait results
# ============================================================

rule aggregate_phenoplier_traits:
    input:
        clampbase=PHENOPLIER_COV_BASE + PHENOPLIER_SAT_BASE + PHENOPLIER_FIN_BASE,
        canonical=PHENOPLIER_CANONICAL,
        script="scripts/phenoplier/aggregate_traits.py",
        wrapper="scripts/phenoplier/aggregate_traits.sh",
    output:
        long=f"{PHENOPLIER_REPORT}/phenoplier_traits_long.csv.gz",
    log:
        f"{PHENOPLIER_REPORT}/aggregate.log"
    params:
        conda_env=PHENOPLIER_CFG["conda_env"],
    resources:
        mem_mb=8000,
        runtime=720,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        "bash {input.wrapper} {params.conda_env} --out {output.long} "
        "--clampbase {input.clampbase} "
        "--canonical {input.canonical} > {log} 2>&1"


# ============================================================
# Step 4: collect completed analyses
# ============================================================

rule phenoplier_coverage:
    input:
        PHENOPLIER_COV_BASE + PHENOPLIER_COMP_BASE,


rule phenoplier_saturation:
    input:
        PHENOPLIER_SAT_BASE,


rule phenoplier_finals:
    input:
        PHENOPLIER_FIN_BASE,


rule phenoplier_canonical:
    input:
        PHENOPLIER_CANONICAL,


rule phenoplier_final_stores:
    input:
        PHENOPLIER_FIN_STORES,


rule phenoplier_traits:
    input:
        PHENOPLIER_COV_BASE
        + PHENOPLIER_SAT_BASE
        + PHENOPLIER_FIN_BASE
        + PHENOPLIER_CANONICAL,


rule phenoplier_long_table:
    input:
        rules.aggregate_phenoplier_traits.output.long,

import os
import csv
import re
import shlex


PHENOPLIER_CFG = config["phenoplier"]
PHENOPLIER_CONDA_ENV = os.path.join(REPO_ROOT, PHENOPLIER_CFG["conda_env_file"])
# Some shell managers prepend shims even during Conda activation.
PHENOPLIER_ENV_SHELL = 'export PATH="${{CONDA_PREFIX:?Use --use-conda}}/bin:$PATH"; '
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
    f"--workspace {shlex.quote(PHENOPLIER_CFG['workspace'])} " if PHENOPLIER_CFG.get("workspace") else ""
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


PHENOPLIER_COV_SEEDS = A4_PW_INDICES
PHENOPLIER_COV_SEED_PATTERN = "|".join(map(str, PHENOPLIER_COV_SEEDS))
PHENOPLIER_GTEX_REPLICATES = sorted(A4_COV_GTEX_REPLICATES) + [
    r for r in A4_PW_GTEX_NEW if r > max(A4_COV_GTEX_REPLICATES)
]
PHENOPLIER_Z_CELLS = [
    a4_cov_cell(fraction, seed)
    for fraction in A4_COV_LEVELS
    for seed in PHENOPLIER_COV_SEEDS
    if seed not in A4_COV_SEEDS
] + [a4_cov_gtex_rep_dir(r) for r in PHENOPLIER_GTEX_REPLICATES if r not in A4_COV_GTEX_REPLICATES]


def phenoplier_cov_model(fraction, seed):
    cell = a4_cov_cell(fraction, seed)
    if int(seed) in A4_COV_SEEDS:
        return f"{cell}/CLAMPbase.rds"
    return f"{cell}/CLAMPbase/phenoplier_Z.pkl"


def phenoplier_cov_name(fraction, seed):
    fraction, seed = int(fraction), int(seed)
    if fraction == 100 and seed > 1:
        label = seed if seed in A4_COV_SEEDS else A4_PW_FULL_SEEDS[seed - 1]
        return f"archs4_rs100_seed{label}_CLAMPbase"
    return f"cov_rs{fraction}_seed{seed}_CLAMPbase"


def phenoplier_gtex_model(replicate):
    if int(replicate) in A4_COV_GTEX_REPLICATES:
        return f"{a4_cov_gtex_rep_dir(replicate)}/CLAMPbase.rds"
    return f"{a4_cov_gtex_rep_dir(replicate)}/CLAMPbase/phenoplier_Z.pkl"


A4_FINAL_DATASETS = list(A4_COV_CFG["final_datasets"])
PHENOPLIER_TRAIT_FINALS = list(PHENOPLIER_CFG["trait_finals"])
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
    for seed in PHENOPLIER_COV_SEEDS
]
PHENOPLIER_SAT_BASE = [
    f"{PHENOPLIER_SAT_BASE_DIR}/sat_rs{fraction}_k{k}_seed{seed}.tsv.gz"
    for (fraction, k, seed) in PHENOPLIER_SAT_CELLS
]
PHENOPLIER_FIN_BASE = [
    f"{PHENOPLIER_FIN_BASE_DIR}/final_{dataset}.tsv.gz" for dataset in PHENOPLIER_TRAIT_FINALS
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
    f"{PHENOPLIER_COMP_BASE_DIR}/gtex_rs100_seed{rep}.tsv.gz" for rep in PHENOPLIER_GTEX_REPLICATES
]


PHENOPLIER_GLS_SHELL = PHENOPLIER_ENV_SHELL + "".join([
    "bash {input.script:q} --rds {params.rds:q} --name {params.name:q} ",
    "--summary-out {output.summary:q} --n-jobs {threads} --mem-mb {resources.mem_mb} ",
    PHENOPLIER_WORKSPACE_ARG,
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
    "> {log:q} 2>&1",
])

PHENOPLIER_STORE_SHELL = PHENOPLIER_ENV_SHELL + "".join([
    "bash {input.script:q} --rds {params.rds:q} --name {params.name:q} ",
    "--store-out {output.store:q} ",
    f"--namespace {PHENOPLIER_CFG['namespace']} ",
    f"--lv-percentile {PHENOPLIER_CFG['lv_percentile']} ",
    f"--trait-filter {PHENOPLIER_TRAIT_FILTER} ",
    PHENOPLIER_WORKSPACE_ARG,
    f"--require-version {PHENOPLIER_CFG['version']} ",
    f"--cohort {PHENOPLIER_CFG['cohort']} ",
    "> {log:q} 2>&1",
])


# ============================================================
# Step 1: test module-trait associations
# ============================================================

rule phenoplier_z_input:
    input:
        z=ancient("{cell}/CLAMPbase/Z.csv"),
    output:
        temp("{cell}/CLAMPbase/phenoplier_Z.pkl"),
    wildcard_constraints:
        cell="|".join(re.escape(cell) for cell in PHENOPLIER_Z_CELLS),
    conda: "clamp-analyses"
    shell:
        "python -c \"import pandas as pd; pd.read_csv('{input.z}', index_col=0).to_pickle('{output}')\""


rule gls_bp_coverage_base:
    input:
        rds=lambda wc: phenoplier_cov_model(wc.fraction, wc.seed),
        script="scripts/phenoplier/run_gls.sh",
        adapter="scripts/phenoplier/workflow.py",
        verify="scripts/phenoplier/verify_summary.py",
    output:
        summary=f"{PHENOPLIER_COV_BASE_DIR}/cov_rs{{fraction}}_seed{{seed}}.tsv.gz",
    log:
        f"{PHENOPLIER_COV_BASE_DIR}/cov_rs{{fraction}}_seed{{seed}}.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: phenoplier_cov_name(wc.fraction, wc.seed),
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        fraction=A4_COV_LEVEL_PATTERN,
        seed=PHENOPLIER_COV_SEED_PATTERN,
    conda: PHENOPLIER_CONDA_ENV
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_bp_saturation_base:
    input:
        rds=lambda wc: a4_sat_base_rds(wc.fraction, wc.k, wc.seed),
        script="scripts/phenoplier/run_gls.sh",
        adapter="scripts/phenoplier/workflow.py",
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
    conda: PHENOPLIER_CONDA_ENV
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_gtex_coverage_base:
    input:
        rds=lambda wc: phenoplier_gtex_model(wc.replicate),
        script="scripts/phenoplier/run_gls.sh",
        adapter="scripts/phenoplier/workflow.py",
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
        replicate="|".join(map(str, PHENOPLIER_GTEX_REPLICATES)),
    conda: PHENOPLIER_CONDA_ENV
    shell:
        PHENOPLIER_GLS_SHELL


rule gls_final_base:
    input:
        rds=lambda wc: a4_final_base_rds(wc.dataset),
        script="scripts/phenoplier/run_gls.sh",
        adapter="scripts/phenoplier/workflow.py",
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
    conda: PHENOPLIER_CONDA_ENV
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
        adapter="scripts/phenoplier/workflow.py",
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
    conda: PHENOPLIER_CONDA_ENV
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
        adapter="scripts/phenoplier/workflow.py",
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
    conda: PHENOPLIER_CONDA_ENV
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
        check_env="scripts/phenoplier/check_env.py",
    output:
        long=f"{PHENOPLIER_REPORT}/phenoplier_traits_long.csv.gz",
    log:
        f"{PHENOPLIER_REPORT}/aggregate.log"
    params:
        version=PHENOPLIER_CFG["version"],
    resources:
        mem_mb=8000,
        runtime=720,
    conda: PHENOPLIER_CONDA_ENV
    shell:
        PHENOPLIER_ENV_SHELL +
        "python {input.check_env:q} --require-version {params.version:q} && "
        "python {input.script:q} --out {output.long:q} "
        "--clampbase {input.clampbase} "
        "--canonical {input.canonical} > {log:q} 2>&1"


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


# October 2026 campaign: 192 delivered units, checked on pico on 2026-11-02.
PHENOPLIER_CAMPAIGN = PHENOPLIER_CFG["campaign"]
with open(PHENOPLIER_CAMPAIGN["manifest"]) as handle:
    _campaign_rows = list(csv.DictReader(handle, delimiter="\t"))
PHENOPLIER_MODELS = {row["model"]: row for row in _campaign_rows}
if len(PHENOPLIER_MODELS) != len(_campaign_rows):
    raise ValueError("Duplicate model names in the PhenoPLIER campaign manifest")
PHENOPLIER_CAMPAIGN_MODE = config.get("phenoplier_campaign_mode", PHENOPLIER_CAMPAIGN["mode"])
if PHENOPLIER_CAMPAIGN_MODE not in {"recorded", "recompute"}:
    raise ValueError("phenoplier_campaign_mode must be recorded or recompute")
PHENOPLIER_RESULTS_ROOT = config.get("phenoplier_results_root", PHENOPLIER_CAMPAIGN["results_root"])
PHENOPLIER_CAMPAIGN_ROOT = PHENOPLIER_CAMPAIGN["output_root"]
if PHENOPLIER_CAMPAIGN_MODE == "recompute":
    PHENOPLIER_CAMPAIGN_ROOT += "-recomputed"
PHENOPLIER_MODEL_PATTERN = "|".join(re.escape(name) for name in PHENOPLIER_MODELS)


def phenoplier_campaign_source(wc):
    row = PHENOPLIER_MODELS[wc.model]
    root = PHENOPLIER_CAMPAIGN[
        "historical_root" if row["batch"] == "historical" else "source_root"
    ]
    return os.path.normpath(os.path.join(root, row["source"]))


def phenoplier_campaign_args(wc):
    row = PHENOPLIER_MODELS[wc.model]
    args = [
        "--name", wc.model, "--namespace", PHENOPLIER_CFG["namespace"] + "/oct2026",
        "--require-version", PHENOPLIER_CFG["version"],
        "--cohort", PHENOPLIER_CFG["cohort"],
        "--lv-percentile", PHENOPLIER_CFG["lv_percentile"],
        "--fit-trait-filter", row["fit_trait_filter"],
        "--trait-filter", "phenomexcan-curated", "--expected-phenotypes", 2824,
        "--expected-genes", row["genes"], "--expected-lvs", row["lvs"],
        "--executor", PHENOPLIER_CLUSTER["executor"],
        "--workers-step6", PHENOPLIER_WORKERS["step6"],
        "--workers-step7", PHENOPLIER_WORKERS["step7"],
    ]
    if PHENOPLIER_CFG.get("workspace"):
        args += ["--workspace", PHENOPLIER_CFG["workspace"]]
    if PHENOPLIER_CLUSTER.get("cluster"):
        args += ["--cluster", PHENOPLIER_CLUSTER["cluster"]]
    return shlex.join(str(arg) for arg in args)


if PHENOPLIER_CAMPAIGN_MODE == "recorded":
    rule collect_campaign:
        input:
            summary=f"{PHENOPLIER_RESULTS_ROOT}/{{model}}/summary/gls-summary-phenomexcan_rapid_gwas-phenomexcan-curated.tsv.gz",
            pickle=f"{PHENOPLIER_RESULTS_ROOT}/{{model}}/summary/gls-summary-phenomexcan_rapid_gwas-phenomexcan-curated.pkl.gz",
            store=f"{PHENOPLIER_RESULTS_ROOT}/{{model}}/study-phenomexcan-curated.h5",
            excluded=f"{PHENOPLIER_RESULTS_ROOT}/{{model}}/summary/trait-filter-excluded-phenomexcan-curated.tsv",
            manifest=PHENOPLIER_CAMPAIGN["manifest"],
            script="scripts/phenoplier/collect_campaign.py",
            check_env="scripts/phenoplier/check_env.py",
            verify="scripts/phenoplier/verify_summary.py",
        output:
            summary=f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/curated.tsv.gz",
            pickle=f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/curated.pkl.gz",
            store=f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/study-phenomexcan-curated.h5",
            excluded=f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/curated_excluded.tsv",
            provenance=f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/source.json",
        params:
            version=PHENOPLIER_CFG["version"],
            source=lambda wc: os.path.join(PHENOPLIER_RESULTS_ROOT, wc.model),
            dest=lambda wc: os.path.join(PHENOPLIER_CAMPAIGN_ROOT, wc.model),
            lvs=lambda wc: PHENOPLIER_MODELS[wc.model]["lvs"],
        log:
            f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/collect.log",
        resources:
            mem_mb=32000,
            runtime=60,
        wildcard_constraints:
            model=PHENOPLIER_MODEL_PATTERN,
        conda: PHENOPLIER_CONDA_ENV
        shell:
            PHENOPLIER_ENV_SHELL +
            "python {input.check_env:q} --require-version {params.version:q} > {log:q} 2>&1 && "
            "python {input.verify:q} {input.pickle:q} 2824 {params.lvs} >> {log:q} 2>&1 && "
            "python {input.script:q} --source {params.source:q} --dest {params.dest:q} "
            "--manifest {input.manifest:q} --model {wildcards.model:q} >> {log:q} 2>&1"

else:
    rule gls_campaign:
        input:
            model=phenoplier_campaign_source,
            manifest=PHENOPLIER_CAMPAIGN["manifest"],
            script="scripts/phenoplier/run_gls.sh",
            adapter="scripts/phenoplier/workflow.py",
            verify="scripts/phenoplier/verify_summary.py",
        output:
            summary=f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/curated.tsv.gz",
            pickle=f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/curated.pkl.gz",
        params:
            args=phenoplier_campaign_args,
        log:
            f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/gls.log",
        threads: int(PHENOPLIER_GLS_RES["threads"])
        resources:
            mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
            runtime=int(PHENOPLIER_GLS_RES["runtime"]),
        wildcard_constraints:
            model=PHENOPLIER_MODEL_PATTERN,
        conda: PHENOPLIER_CONDA_ENV
        shell:
            "bash {input.script:q} --rds {input.model:q} --summary-out {output.summary:q} "
            "--n-jobs {threads} --mem-mb {resources.mem_mb} {params.args} > {log:q} 2>&1"


    rule store_campaign:
        input:
            summary=rules.gls_campaign.output.summary,
            pickle=rules.gls_campaign.output.pickle,
            model=phenoplier_campaign_source,
            script="scripts/phenoplier/store_build.sh",
            adapter="scripts/phenoplier/workflow.py",
        output:
            store=f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/study-phenomexcan-curated.h5",
        params:
            args=phenoplier_campaign_args,
        log:
            f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/store.log",
        resources:
            mem_mb=int(PHENOPLIER_STORE_RES["mem_mb"]),
            runtime=int(PHENOPLIER_STORE_RES["runtime"]),
        wildcard_constraints:
            model=PHENOPLIER_MODEL_PATTERN,
        conda: PHENOPLIER_CONDA_ENV
        shell:
            "bash {input.script:q} --rds {input.model:q} --store-out {output.store:q} "
            "{params.args} > {log:q} 2>&1"

rule phenoplier_campaign:
    input:
        expand(f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/study-phenomexcan-curated.h5",
               model=list(PHENOPLIER_MODELS)),


rule phenoplier_campaign_new:
    input:
        expand(f"{PHENOPLIER_CAMPAIGN_ROOT}/{{model}}/study-phenomexcan-curated.h5",
               model=[name for name, row in PHENOPLIER_MODELS.items() if row["batch"] == "new"]),


# Install/check the isolated CLI without model inputs, reference data or GLS.
rule phenoplier_environment:
    input:
        script="scripts/phenoplier/check_env.py",
    output:
        report="output/03_model_biology/03_phenoplier/environment.json",
    params:
        version=PHENOPLIER_CFG["version"],
    conda: PHENOPLIER_CONDA_ENV
    shell:
        PHENOPLIER_ENV_SHELL +
        "python {input.script:q} --require-version {params.version:q} --output {output.report:q}"

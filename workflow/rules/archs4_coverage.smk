import json
import re

A4_COV_CFG = A4_CFG["coverage"]
A4_COV_BIO = f"{A4_BIO}/00_coverage"
A4_COV_NB = f"{A4_BIO_NB}/00_coverage"
A4_COV_ORA_ROOT = A4_COV_CFG["ora_root"]
A4_COV_INPUT_ROOT = A4_COV_CFG["input_root"]
A4_COV_LEVELS = [int(level) for level in A4_COV_CFG["levels"]]
A4_COV_SEEDS = [int(seed) for seed in A4_COV_CFG["seeds"]]
A4_COV_REFERENCE_SEED = 1
A4_COV_RNG_SEEDS = [int(seed) for seed in A4_COV_CFG["rng_seeds"]]
A4_COV_DATABASES = list(A4_COV_CFG["databases"])
A4_COV_COMPARATORS = list(A4_COV_CFG["comparators"])
A4_COV_LEVEL_PATTERN = "|".join(map(str, A4_COV_LEVELS))
A4_COV_SEED_PATTERN = "|".join(map(str, A4_COV_SEEDS))
A4_COV_DATABASE_PATTERN = "|".join(A4_COV_DATABASES)
A4_COV_COMPARATOR_PATTERN = "|".join(A4_COV_COMPARATORS)
A4_COV_RETRY_FACTORS = [float(value) for value in A4_COV_CFG["routing"]["retry_mem_factors"]]
A4_COV_RETRY_CAP = int(A4_COV_CFG["routing"]["retry_mem_cap_mb"])

A4_COV_FULL_REPLICATES = {int(k): int(v) for k, v in A4_COV_CFG["full_data_replicates"].items()}
A4_COV_GTEX_REPLICATES = {
    int(k): int(v) for k, v in A4_COV_CFG["comparators"]["gtex"]["replicate_seeds"].items()
}
A4_COV_GTEX_REP_ROOT = A4_COV_CFG["comparators"]["gtex"]["replicate_root"]


def a4_cov_gtex_rep_dir(replicate):
    return f"{A4_COV_GTEX_REP_ROOT}/rs100_seed_{replicate}"


def a4_cov_cell(fraction, seed):
    fraction = str(fraction)
    level_dir = A4_COV_CFG["levels"][fraction]
    return (
        f"{A4_COV_INPUT_ROOT}/{level_dir}/"
        f"study_coverage_rs{fraction}_seed_{seed}"
    )


def a4_cov_ora_dir(dataset, fraction, seed, model, database):
    return f"{A4_COV_ORA_ROOT}/{dataset}/rs{fraction}/seed{seed}/{model}/{database}"


def a4_cov_rng_seed(seed):
    return A4_PW_SUB_SEEDS[int(seed) - 1]


def a4_cov_initial_mem(fraction):
    return int(A4_COV_CFG["resources"][str(fraction)]["mem_mb"])


def a4_cov_retry_mem(fraction, attempt):
    factor = A4_COV_RETRY_FACTORS[min(int(attempt), len(A4_COV_RETRY_FACTORS)) - 1]
    return min(int(round(a4_cov_initial_mem(fraction) * factor)), A4_COV_RETRY_CAP)


def a4_cov_db(wildcards):
    return A4_CFG["ora"]["databases"][wildcards.database]


def a4_cov_base_z(wildcards):
    seed = int(wildcards.seed)
    if wildcards.dataset == "archs4":
        return f"{a4_cov_cell(wildcards.fraction, wildcards.seed)}/CLAMPbase/Z.csv"
    if wildcards.dataset == "gtex" and (seed in A4_COV_GTEX_REPLICATES or seed in A4_PW_GTEX_NEW):
        return f"{a4_cov_gtex_rep_dir(seed)}/CLAMPbase/Z.csv"
    return A4_COV_CFG["comparators"][wildcards.dataset]["base_z"]


def a4_cov_sample_arg(wildcards, full):
    if wildcards.dataset == "archs4":
        return f"--subsample-info {a4_cov_cell(wildcards.fraction, wildcards.seed)}/subsample_info.rds"
    samples = A4_COV_CFG["comparators"][wildcards.dataset]["samples_rds"]
    return f"--samples-rds {samples}"


def a4_cov_seeds_for(fraction):
    return A4_COV_SEEDS


def a4_cov_fbm(fraction, seed):
    if int(fraction) == 100:
        return f"{A4_PREP}/fbm_filtered.bk"
    return f"{a4_cov_cell(fraction, seed)}/fbm_subsampled.bk"

def a4_cov_fbm_dep(fraction, seed):
    if int(fraction) == 100:
        return f"{A4_PREP}/fbm_filtered.bk"
    return f"{a4_cov_cell(fraction, seed)}/fbm_subsampled.bk"

A4_COV_CELLS = {
    a4_cov_cell(fraction, seed): (fraction, seed)
    for fraction in A4_COV_LEVELS
    for seed in a4_cov_seeds_for(fraction)
}
A4_PW_CFG = A4_COV_CFG["pathway"]
A4_PW_SUB_SEEDS = [int(s) for s in A4_PW_CFG["subsample_seeds"]]
A4_PW_FULL_SEEDS = [int(s) for s in A4_PW_CFG["archs4_full_seeds"]]
A4_PW_GTEX_SEEDS = [int(s) for s in A4_PW_CFG["gtex_seeds"]]
A4_PW_INDICES = list(range(1, len(A4_PW_SUB_SEEDS) + 1))
A4_PW_NEW = [i for i in A4_PW_INDICES if i not in A4_COV_SEEDS]
A4_PW_GTEX_NEW = [i for i in A4_PW_INDICES if i not in A4_COV_GTEX_REPLICATES]
A4_PW_SUB_LEVELS = [f for f in A4_COV_LEVELS if f < 100]
A4_PW_CELLS = {
    a4_cov_cell(fraction, index): (fraction, index)
    for fraction in A4_PW_SUB_LEVELS
    for index in A4_PW_NEW
}
A4_PW_FULL_CELLS = {a4_cov_cell(100, index): index for index in A4_PW_NEW}
A4_COV_ALL_CELLS = {**A4_COV_CELLS, **A4_PW_CELLS, **{cell: (100, index) for cell, index in A4_PW_FULL_CELLS.items()}}
A4_COV_CELL_PATTERN = "|".join(re.escape(cell) for cell in A4_COV_ALL_CELLS)


A4_COV_FULL_REP_CELLS = {a4_cov_cell(100, replicate): replicate for replicate in A4_COV_FULL_REPLICATES}
A4_COV_FULL_REP_CELL_PATTERN = "|".join(re.escape(cell) for cell in A4_COV_FULL_REP_CELLS)
A4_COV_SUBSAMPLE_CELL_PATTERN = "|".join(
    re.escape(cell) for cell in A4_COV_CELLS if cell not in A4_COV_FULL_REP_CELLS
)
A4_COV_SVD_CELL_PATTERN = A4_COV_SUBSAMPLE_CELL_PATTERN + "|" + "|".join(re.escape(cell) for cell in A4_PW_CELLS)


def a4_cov_fraction_of(cell):
    return A4_COV_ALL_CELLS[str(cell)][0]


def a4_cov_seed_of(cell):
    return A4_COV_ALL_CELLS[str(cell)][1]


A4_COV_CELL_TARGETS = [f"{cell}/CLAMPbase.rds" for cell in A4_COV_CELLS]
A4_COV_GTEX_REP_TARGETS = [
    f"{a4_cov_gtex_rep_dir(replicate)}/CLAMPbase/Z.csv" for replicate in A4_COV_GTEX_REPLICATES
]

A4_COV_ARCH_BASE_ORA = [
    a4_cov_ora_dir("archs4", fraction, seed, "CLAMPbase", database)
    for fraction in A4_COV_LEVELS
    for seed in a4_cov_seeds_for(fraction)
    for database in A4_COV_DATABASES
]
A4_COV_COMPARATOR_BASE_ORA = [
    a4_cov_ora_dir(dataset, 100, 0, "CLAMPbase", database)
    for dataset in A4_COV_COMPARATORS
    for database in A4_COV_DATABASES
] + [
    a4_cov_ora_dir("gtex", 100, seed, "CLAMPbase", database)
    for seed in A4_COV_GTEX_REPLICATES
    for database in A4_COV_DATABASES
]
A4_COV_BASE_ORA = A4_COV_ARCH_BASE_ORA + A4_COV_COMPARATOR_BASE_ORA


# ============================================================
# Step 1: prepare study subsamples
# ============================================================

rule subsample_bp_coverage_cell:
    input:
        h5=A4_CFG["raw_h5"],
        metadata=f"{A4_PREP}/metadata_filtered.rds",
        fbm=f"{A4_PREP}/fbm_filtered.bk",
        script="scripts/archs4/coverage/subsample_cells.R",
        common="scripts/archs4/common.R",
    output:
        info=f"{{cell}}/subsample_info.rds",
        metadata=f"{{cell}}/metadata_cell.rds",
        samples=f"{{cell}}/samples.rds",
        fbm=f"{{cell}}/fbm_subsampled.bk",
    log:
        f"{{cell}}/subsample.log"
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.info),
        fraction=lambda wc: a4_cov_fraction_of(wc.cell),
        seed=lambda wc: a4_cov_rng_seed(a4_cov_seed_of(wc.cell)),
        run_index=lambda wc: a4_cov_seed_of(wc.cell),
        sc_max=A4_CFG["preprocess"]["single_cell_probability_max"],
        block_size=A4_CFG["raw_data"]["block_size"],
    resources:
        mem_mb=lambda wc: a4_cov_initial_mem(a4_cov_fraction_of(wc.cell)),
        runtime=5760,
    wildcard_constraints:
        cell=A4_COV_CELL_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --h5 {input.h5} --metadata-filtered {input.metadata} "
        "--fbm-filtered {input.fbm} --out-dir {params.out_dir} "
        "--fraction {params.fraction} --seed {params.seed} --run-index {params.run_index} "
        "--single-cell-probability-max {params.sc_max} --block-size {params.block_size} "
        "> {log} 2>&1"


rule svd_bp_coverage_cell:
    input:
        metadata=rules.subsample_bp_coverage_cell.output.metadata,
        fbm=lambda wc: a4_cov_fbm_dep(a4_cov_fraction_of(wc.cell), a4_cov_seed_of(wc.cell)),
        script="scripts/archs4/svd.R",
        common="scripts/archs4/common.R",
    output:
        svd=f"{{cell}}/svd.rds",
        k=f"{{cell}}/CLAMP_K.rds",
    log:
        f"{{cell}}/svd.log"
    params:
        fbm_prefix=lambda wc: a4_cov_fbm(a4_cov_fraction_of(wc.cell), a4_cov_seed_of(wc.cell))[: -len(".bk")],
        n_cores=A4_CFG["svd"]["n_cores"],
        k_multiplier=A4_CFG["clamp"]["k_multiplier"],
        seed=A4_CFG["svd"]["seed"],
    threads: A4_CFG["svd"]["n_cores"]
    resources:
        mem_mb=lambda wc: a4_cov_initial_mem(a4_cov_fraction_of(wc.cell)),
        runtime=6000,
    wildcard_constraints:
        cell=A4_COV_SVD_CELL_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --metadata {input.metadata} --fbm {params.fbm_prefix} "
        "--svd {output.svd} --k {output.k} --n-cores {params.n_cores} "
        "--k-multiplier {params.k_multiplier} --seed {params.seed} "
        "> {log} 2>&1"


# ============================================================
# Step 2: fit CLAMPbase models
# ============================================================

rule clampbase_bp_coverage_cell:
    input:
        metadata=rules.subsample_bp_coverage_cell.output.metadata,
        samples=rules.subsample_bp_coverage_cell.output.samples,
        fbm=lambda wc: a4_cov_fbm_dep(a4_cov_fraction_of(wc.cell), a4_cov_seed_of(wc.cell)),
        svd=rules.svd_bp_coverage_cell.output.svd,
        k=rules.svd_bp_coverage_cell.output.k,
        script=ancient("scripts/archs4/clamp.R"),
        common="scripts/archs4/common.R",
    output:
        rds=f"{{cell}}/CLAMPbase.rds",
        b=f"{{cell}}/CLAMPbase/B.csv",
        z=f"{{cell}}/CLAMPbase/Z.csv",
    log:
        f"{{cell}}/clampbase.log"
    params:
        fbm_prefix=lambda wc: a4_cov_fbm(a4_cov_fraction_of(wc.cell), a4_cov_seed_of(wc.cell))[: -len(".bk")],
        out_dir=lambda wildcards, output: os.path.dirname(output.rds),
        seed=A4_CFG["clamp"]["seed"],
    resources:
        mem_mb=lambda wc: a4_cov_retry_mem(a4_cov_fraction_of(wc.cell), 1),
        runtime=12000,
    wildcard_constraints:
        cell=A4_COV_SUBSAMPLE_CELL_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --stage base --metadata {input.metadata} "
        "--samples {input.samples} --fbm {params.fbm_prefix} --svd {input.svd} "
        "--k {input.k} --out-dir {params.out_dir} --seed {params.seed} "
        "> {log} 2>&1"


rule clampbase_full_replicate_cell:
    input:
        metadata=f"{A4_PREP}/metadata_filtered.rds",
        fbm=f"{A4_PREP}/fbm_filtered.bk",
        subsample=lambda wc: ancient(f"{wc.cell}/subsample_info.rds"),
        svd_script="scripts/archs4/svd.R",
        base_script=ancient("scripts/archs4/saturation/fit_clampbase.R"),
        common="scripts/archs4/common.R",
    output:
        svd=f"{{cell}}/svd.rds",
        k=f"{{cell}}/CLAMP_K.rds",
        rds=f"{{cell}}/CLAMPbase.rds",
        b=f"{{cell}}/CLAMPbase/B.csv",
        z=f"{{cell}}/CLAMPbase/Z.csv",
        manifest=f"{{cell}}/manifest.json",
    log:
        f"{{cell}}/clampbase_replicate.log"
    params:
        fbm_prefix=lambda wildcards, input: input.fbm[: -len(".bk")],
        out_dir=lambda wildcards, output: os.path.dirname(output.rds),
        replicate=lambda wc: A4_COV_FULL_REP_CELLS[wc.cell],
        rng_seed=lambda wc: A4_COV_FULL_REPLICATES[A4_COV_FULL_REP_CELLS[wc.cell]],
        svd_k=A4_COV_CFG["full_data_svd_k"],
        n_cores=A4_CFG["svd"]["n_cores"],
    threads: 16
    resources:
        mem_mb=lambda wc: a4_cov_initial_mem(100),
        runtime=12000,
    wildcard_constraints:
        cell=A4_COV_FULL_REP_CELL_PATTERN,
    conda: "clamp-analyses"
    shell:
        "(Rscript {input.svd_script} --metadata {input.metadata} --fbm {params.fbm_prefix} "
        "--svd {output.svd} --k {output.k} --svd-k {params.svd_k} --n-cores {params.n_cores} "
        "--seed {params.rng_seed} && touch {params.out_dir}/svd.complete && "
        "Rscript {input.base_script} --metadata {input.metadata} --subsample-info {input.subsample} "
        "--fbm-backing {input.fbm} --svd {output.svd} --k {output.k} --out-dir {params.out_dir} "
        "--seed {params.rng_seed} --seed-index {params.replicate} --fraction 100 "
        "--requested-cpus {threads} --requested-mem-mb {resources.mem_mb}) > {log} 2>&1"


rule clampbase_gtex_replicate:
    input:
        fbm=ancient(A4_COV_CFG["comparators"]["gtex"]["fbm_rds"]),
        backing=ancient(A4_COV_CFG["comparators"]["gtex"]["fbm_backing"]),
        genes=ancient(A4_COV_CFG["comparators"]["gtex"]["genes_rds"]),
        samples=ancient(A4_COV_CFG["comparators"]["gtex"]["samples_rds"]),
        k=ancient(A4_COV_CFG["comparators"]["gtex"]["k"]),
        runner="scripts/gtex/run_clampbase_seed.sh",
        script="scripts/gtex/fit_clampbase_seed.R",
    output:
        svd=f"{A4_COV_GTEX_REP_ROOT}/rs100_seed_{{replicate}}/svd.rds",
        rds=f"{A4_COV_GTEX_REP_ROOT}/rs100_seed_{{replicate}}/CLAMPbase.rds",
        b=f"{A4_COV_GTEX_REP_ROOT}/rs100_seed_{{replicate}}/CLAMPbase/B.csv",
        z=f"{A4_COV_GTEX_REP_ROOT}/rs100_seed_{{replicate}}/CLAMPbase/Z.csv",
        complete=f"{A4_COV_GTEX_REP_ROOT}/rs100_seed_{{replicate}}/complete",
    log:
        f"{A4_COV_GTEX_REP_ROOT}/logs/rs100_seed_{{replicate}}.log"
    params:
        rng_seed=lambda wc: A4_COV_GTEX_REPLICATES[int(wc.replicate)],
        out_dir=lambda wildcards, output: os.path.dirname(output.rds),
    threads: 8
    resources:
        mem_mb=96000,
        runtime=2880,
    wildcard_constraints:
        replicate="|".join(map(str, A4_COV_GTEX_REPLICATES)),
    conda: "clamp-analyses"
    shell:
        "CLAMP_THREADS={threads} bash {input.runner} {params.rng_seed} {wildcards.replicate} "
        "{params.out_dir} > {log} 2>&1"


# ============================================================
# Step 3: score pathway recovery and report results
# ============================================================

rule ora_bp_base_coverage:
    input:
        z=a4_cov_base_z,
        samples=lambda wc: (f"{a4_cov_cell(wc.fraction, wc.seed)}/subsample_info.rds"
                            if wc.dataset == "archs4" else A4_COV_CFG["comparators"][wc.dataset]["samples_rds"]),
        database=lambda wc: a4_cov_db(wc)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(
            f"{A4_COV_ORA_ROOT}/{{dataset}}/rs{{fraction}}/seed{{seed}}/CLAMPbase/{{database}}"
        ),
    log:
        f"{A4_COV_ORA_ROOT}/{{dataset}}/rs{{fraction}}/seed{{seed}}/CLAMPbase/{{database}}.log"
    params:
        sample_arg=lambda wc: a4_cov_sample_arg(wc, False),
        database_type=lambda wc: a4_cov_db(wc)["type"],
        database_label=lambda wc: a4_cov_db(wc)["label"],
        exclude=lambda wc: (
            f"--exclude-term-regex '{a4_cov_db(wc)['exclude_term_regex']}'"
            if "exclude_term_regex" in a4_cov_db(wc) else ""
        ),
        sheet=lambda wc: f"--sheet {a4_cov_db(wc)['sheet']}" if "sheet" in a4_cov_db(wc) else "",
        columns=lambda wc: (
            f"--term-column {a4_cov_db(wc)['term_column']} --gene-column {a4_cov_db(wc)['gene_column']}"
            if "term_column" in a4_cov_db(wc) else ""
        ),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        ora_mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
        ora_slots=1,
    wildcard_constraints:
        dataset="archs4|" + A4_COV_COMPARATOR_PATTERN,
        fraction=A4_COV_LEVEL_PATTERN,
        seed="0|" + "|".join(map(str, A4_PW_INDICES)),
        database=A4_COV_DATABASE_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset {wildcards.dataset} --fraction {wildcards.fraction} --seed {wildcards.seed} "
        "--model CLAMPbase {params.sample_arg} --database {wildcards.database} "
        "--database-label '{params.database_label}' --database-type {params.database_type} "
        "--database-path {input.database} {params.exclude} {params.sheet} {params.columns} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CFG[ora][min_size]} "
        "--max-size {A4_CFG[ora][max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


rule aggregate_bp_coverage:
    input:
        A4_COV_BASE_ORA,
        pathway=lambda wc: A4_PW_ORA,
        script="scripts/archs4/coverage/aggregate_coverage.R",
    output:
        coverage_long=f"{A4_COV_BIO}/coverage_long.csv",
        cross_dataset=f"{A4_COV_BIO}/cross_dataset_coverage.csv",
        panel_ready=f"{A4_COV_BIO}/coverage_panel_ready.csv",
    log:
        f"{A4_COV_BIO}/aggregate.log"
    params:
        databases=",".join(A4_COV_DATABASES),
        labels=lambda wc: A4_PW_DATABASE_LABELS,
    resources:
        mem_mb=16000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --ora-root {A4_COV_ORA_ROOT} "
        "--coverage-out {output.coverage_long} --cross-out {output.cross_dataset} "
        "--panel-out {output.panel_ready} --databases {params.databases} "
        "--database-labels '{params.labels}' "
        "> {log} 2>&1"


rule coverage_report_archs4:
    input:
        coverage_long=rules.aggregate_bp_coverage.output.coverage_long,
        cross_dataset=rules.aggregate_bp_coverage.output.cross_dataset,
        panel_ready=rules.aggregate_bp_coverage.output.panel_ready,
        notebook=f"{A4_COV_NB}/00_coverage.ipynb",
    output:
        complete=touch(f"{A4_COV_BIO}/notebook.complete"),
    log:
        notebook=f"{A4_COV_NB}/00_coverage.executed.ipynb",
    params:
        coverage_dir=A4_COV_BIO,
    resources:
        mem_mb=16000,
        runtime=1440,
    conda: "clamp-analyses"
    notebook:
        f"{A4_COV_NB}/00_coverage.ipynb"


rule coverage_bp_cells:
    input:
        A4_COV_CELL_TARGETS,
        A4_COV_GTEX_REP_TARGETS,


rule coverage_bp_ora:
    input:
        A4_COV_BASE_ORA,


rule coverage_bp_archs4_base_ora:
    input:
        A4_COV_ARCH_BASE_ORA,


rule coverage_bp_comparator_base_ora:
    input:
        A4_COV_COMPARATOR_BASE_ORA,


rule coverage_bp:
    input:
        rules.coverage_report_archs4.output.complete,


rule archs4_coverage:
    input:
        rules.coverage_report_archs4.output.complete,


A4_PW_DATABASE_LABELS = ",".join(
    f"{db}={A4_CFG['ora']['databases'][db]['label']}" for db in A4_COV_DATABASES
)
A4_PW_ORA = [
    a4_cov_ora_dir("archs4", fraction, index, "CLAMPbase", database)
    for fraction in A4_COV_LEVELS
    for index in A4_PW_NEW
    for database in A4_COV_DATABASES
] + [
    a4_cov_ora_dir("gtex", 100, index, "CLAMPbase", database)
    for index in A4_PW_GTEX_NEW
    for database in A4_COV_DATABASES
]


def a4_pw_inputs(target, inputs):
    return {} if os.path.exists(target) else inputs


# ============================================================
# Step 4: fit additional pathway replicates
# ============================================================

rule clampbase_pathway_cell:
    input:
        unpack(lambda wc: a4_pw_inputs(f"{wc.cell}/CLAMPbase/Z.csv", {
            "metadata": f"{A4_PREP}/metadata_filtered.rds",
            "fbm": a4_cov_fbm(*A4_PW_CELLS[wc.cell]),
            "subsample": f"{wc.cell}/subsample_info.rds",
            "svd": f"{wc.cell}/svd.rds",
            "k": f"{wc.cell}/CLAMP_K.rds",
            "script": "scripts/archs4/saturation/fit_clampbase.R",
        })),
    output:
        z=f"{{cell}}/CLAMPbase/Z.csv",
    log:
        f"{{cell}}/clampbase.log"
    params:
        fraction=lambda wc: A4_PW_CELLS[wc.cell][0],
        index=lambda wc: A4_PW_CELLS[wc.cell][1],
        rng_seed=lambda wc: a4_cov_rng_seed(A4_PW_CELLS[wc.cell][1]),
        fbm=lambda wc: a4_cov_fbm(*A4_PW_CELLS[wc.cell]),
    threads: 16
    resources:
        mem_mb=lambda wc: a4_cov_initial_mem(A4_PW_CELLS[wc.cell][0]),
        runtime=lambda wc: int(A4_COV_CFG["resources"][str(A4_PW_CELLS[wc.cell][0])]["runtime"]),
    wildcard_constraints:
        cell="|".join(re.escape(cell) for cell in A4_PW_CELLS),
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/archs4/saturation/fit_clampbase.R --metadata {A4_PREP}/metadata_filtered.rds "
        "--subsample-info {wildcards.cell}/subsample_info.rds --fbm-backing {params.fbm} "
        "--svd {wildcards.cell}/svd.rds --k {wildcards.cell}/CLAMP_K.rds --out-dir {wildcards.cell} "
        "--seed {params.rng_seed} --seed-index {params.index} --fraction {params.fraction} "
        "--initialization svd --z-only true --requested-cpus {threads} > {log} 2>&1"


rule clampbase_gtex_pathway_replicate:
    input:
        unpack(lambda wc: a4_pw_inputs(f"{a4_cov_gtex_rep_dir(wc.replicate)}/CLAMPbase/Z.csv", {
            "fbm": A4_COV_CFG["comparators"]["gtex"]["fbm_rds"],
            "backing": A4_COV_CFG["comparators"]["gtex"]["fbm_backing"],
            "genes": A4_COV_CFG["comparators"]["gtex"]["genes_rds"],
            "samples": A4_COV_CFG["comparators"]["gtex"]["samples_rds"],
            "k": A4_COV_CFG["comparators"]["gtex"]["k"],
            "runner": "scripts/gtex/run_clampbase_seed.sh",
            "script": "scripts/gtex/fit_clampbase_seed.R",
        })),
    output:
        z=f"{A4_COV_GTEX_REP_ROOT}/rs100_seed_{{replicate}}/CLAMPbase/Z.csv",
    log:
        f"{A4_COV_GTEX_REP_ROOT}/logs/rs100_seed_{{replicate}}.log"
    params:
        rng_seed=lambda wc: A4_PW_GTEX_SEEDS[int(wc.replicate) - 1],
        out_dir=lambda wc: a4_cov_gtex_rep_dir(wc.replicate),
    threads: 8
    resources:
        mem_mb=96000,
        runtime=2880,
    wildcard_constraints:
        replicate="|".join(map(str, A4_PW_GTEX_NEW)),
    conda: "clamp-analyses"
    shell:
        "CLAMP_THREADS={threads} bash scripts/gtex/run_clampbase_seed.sh {params.rng_seed} "
        "{wildcards.replicate} {params.out_dir} > {log} 2>&1"


rule clampbase_full_pathway_cell:
    input:
        unpack(lambda wc: a4_pw_inputs(f"{wc.cell}/CLAMPbase/Z.csv", {
            "metadata": f"{A4_PREP}/metadata_filtered.rds",
            "fbm": f"{A4_PREP}/fbm_filtered.bk",
            "subsample": f"{wc.cell}/subsample_info.rds",
            "metadata": f"{A4_PREP}/metadata_filtered.rds",
            "fbm": f"{A4_PREP}/fbm_filtered.bk",
            "subsample": f"{wc.cell}/subsample_info.rds",
            "svd_script": "scripts/archs4/svd.R",
            "base_script": "scripts/archs4/saturation/fit_clampbase.R",
        })),
    output:
        z=f"{{cell}}/CLAMPbase/Z.csv",
    log:
        f"{{cell}}/clampbase_replicate.log"
    params:
        index=lambda wc: A4_PW_FULL_CELLS[wc.cell],
        rng_seed=lambda wc: A4_PW_FULL_SEEDS[A4_PW_FULL_CELLS[wc.cell] - 1],
        svd_k=A4_COV_CFG["full_data_svd_k"],
        n_cores=A4_CFG["svd"]["n_cores"],
    threads: 16
    resources:
        mem_mb=lambda wc: a4_cov_initial_mem(100),
        runtime=12000,
    wildcard_constraints:
        cell="|".join(re.escape(c) for c in A4_PW_FULL_CELLS),
    conda: "clamp-analyses"
    shell:
        "(Rscript scripts/archs4/svd.R --metadata {A4_PREP}/metadata_filtered.rds --fbm {A4_PREP}/fbm_filtered "
        "--svd {wildcards.cell}/svd.rds --k {wildcards.cell}/CLAMP_K.rds --svd-k {params.svd_k} "
        "--n-cores {params.n_cores} --seed {params.rng_seed} && "
        "Rscript scripts/archs4/saturation/fit_clampbase.R --metadata {A4_PREP}/metadata_filtered.rds "
        "--subsample-info {wildcards.cell}/subsample_info.rds --fbm-backing {A4_PREP}/fbm_filtered.bk "
        "--svd {wildcards.cell}/svd.rds --k {wildcards.cell}/CLAMP_K.rds --out-dir {wildcards.cell} "
        "--seed {params.rng_seed} --seed-index {params.index} --fraction 100 --initialization random "
        "--z-only true --requested-cpus {threads}) > {log} 2>&1"


rule coverage_pathway_ora:
    input:
        A4_PW_ORA,

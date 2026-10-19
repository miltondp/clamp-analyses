A4_SAT_CFG = A4_CFG["saturation"]
A4_SAT_BIO = f"{A4_BIO}/01_saturation"
A4_SAT_NB = f"{A4_BIO_NB}/01_saturation"
A4_SAT_ORA_ROOT = A4_SAT_CFG["ora_root"]
A4_SAT_INPUT_ROOT = A4_SAT_CFG["input_root"]
A4_SAT_BASE_ROOT = A4_SAT_CFG["base_root"]
A4_SAT_FRACTIONS = [int(f) for f in A4_SAT_CFG["fractions"]]
A4_SAT_KS = [int(k) for k in A4_SAT_CFG["k_values"]]
A4_SAT_SEEDS = [int(s) for s in A4_SAT_CFG["seeds"]]
A4_SAT_RNG_SEEDS = [int(s) for s in A4_SAT_CFG["rng_seeds"]]
A4_SAT_DATABASES = ["canonical", "cellmarker"]
A4_SAT_FRACTION_PATTERN = "|".join(map(str, A4_SAT_FRACTIONS))
A4_SAT_K_PATTERN = "|".join(map(str, A4_SAT_KS))
A4_SAT_SEED_PATTERN = "|".join(map(str, A4_SAT_SEEDS))
A4_SAT_DATABASE_PATTERN = "|".join(A4_SAT_DATABASES)
A4_SAT_RETRY_FACTORS = [float(v) for v in A4_SAT_CFG["routing"]["retry_mem_factors"]]
A4_SAT_RETRY_CAP = int(A4_SAT_CFG["routing"]["retry_mem_cap_mb"])
A4_SAT_THREADS = int(A4_SAT_CFG["routing"]["threads_per_fit"])


def a4_sat_input_cell(fraction, seed):

    level = A4_SAT_CFG["input_levels"][str(fraction)]
    leaf = A4_SAT_CFG["input_leaf"].format(fraction=fraction, seed=seed)
    return f"{A4_SAT_INPUT_ROOT}/{level}/{leaf}"


def a4_sat_base_cell(fraction, k, seed):
    leaf = A4_SAT_CFG["base_leaf"].format(fraction=fraction, k=k, seed=seed)
    return f"{A4_SAT_BASE_ROOT}/{leaf}"


def a4_sat_base_rds(fraction, k, seed):
    return f"{a4_sat_base_cell(fraction, k, seed)}/CLAMPbase.rds"


def a4_sat_base_z(fraction, k, seed):
    return f"{a4_sat_base_cell(fraction, k, seed)}/CLAMPbase/Z.csv"


def a4_sat_ora_dir(fraction, k, seed, model, database):
    return f"{A4_SAT_ORA_ROOT}/rs{fraction}/k{k}/seed{seed}/{model}/{database}"


def a4_sat_rng_seed(seed):
    return A4_SAT_RNG_SEEDS[int(seed) - 1]


def a4_sat_resource(fraction, k, key):
    return int(A4_SAT_CFG["resources"][str(fraction)][str(k)][key])


def a4_sat_retry_mem(fraction, k, attempt):
    factor = A4_SAT_RETRY_FACTORS[min(int(attempt), len(A4_SAT_RETRY_FACTORS)) - 1]
    return min(int(round(a4_sat_resource(fraction, k, "mem_mb") * factor)), A4_SAT_RETRY_CAP)


def a4_sat_db(wildcards):
    return A4_CFG["ora"]["databases"][wildcards.database]


A4_SAT_EXCLUDED = {
    (int(f), int(k), int(s))
    for f, k, s in A4_SAT_CFG.get("exclude_cells", [])
}

A4_SAT_CELLS = [
    (f, k, s)
    for f in A4_SAT_FRACTIONS
    for k in A4_SAT_KS
    for s in A4_SAT_SEEDS
    if (f, k, s) not in A4_SAT_EXCLUDED
]


def a4_sat_base_present():

    return [c for c in A4_SAT_CELLS if os.path.exists(a4_sat_base_z(*c))]


def a4_sat_base_missing():
    present = set(a4_sat_base_present())
    return [c for c in A4_SAT_CELLS if c not in present]


A4_SAT_BASE_ORA_PRESENT = [
    a4_sat_ora_dir(f, k, s, "CLAMPbase", db)
    for (f, k, s) in a4_sat_base_present()
    for db in A4_SAT_DATABASES
]
A4_SAT_BASE_ORA_ALL = [
    a4_sat_ora_dir(f, k, s, "CLAMPbase", db)
    for (f, k, s) in A4_SAT_CELLS
    for db in A4_SAT_DATABASES
]


# ============================================================
# Step 1: fit and score fixed-rank models
# ============================================================

rule fit_bp_saturation_base:
    input:
        svd=lambda wc: ancient(f"{a4_sat_input_cell(wc.fraction, wc.seed)}/svd.rds"),
        subsample=lambda wc: ancient(f"{a4_sat_input_cell(wc.fraction, wc.seed)}/subsample_info.rds"),
        metadata=ancient(f"{A4_PROD}/00_preprocess/metadata_filtered.rds"),
        fbm=lambda wc: f"{a4_sat_input_cell(wc.fraction, wc.seed)}/fbm_subsampled.bk",
        script=ancient("scripts/archs4/saturation/fit_clampbase.R"),
        runner=ancient("scripts/archs4/coverage/run_with_metrics.py"),
    output:
        base_rds=f"{A4_SAT_BASE_ROOT}/{A4_SAT_CFG['base_leaf']}/CLAMPbase.rds",
        base_z=f"{A4_SAT_BASE_ROOT}/{A4_SAT_CFG['base_leaf']}/CLAMPbase/Z.csv",
    log:
        f"{A4_SAT_BASE_ROOT}/{A4_SAT_CFG['base_leaf']}/clampbase.log"
    params:
        fbm=lambda wc: f"{a4_sat_input_cell(wc.fraction, wc.seed)}/fbm_subsampled.bk",
        out_dir=lambda wc: a4_sat_base_cell(wc.fraction, wc.k, wc.seed),
        metrics_dir=lambda wc: f"{a4_sat_base_cell(wc.fraction, wc.k, wc.seed)}/attempts",
        rng_seed=lambda wc: a4_sat_rng_seed(wc.seed),
        initial_mem=lambda wc: a4_sat_resource(wc.fraction, wc.k, "mem_mb"),
    threads: A4_SAT_THREADS
    resources:
        mem_mb=lambda wc, attempt: a4_sat_retry_mem(wc.fraction, wc.k, attempt),
        fit_mem_mb=lambda wc, attempt: a4_sat_retry_mem(wc.fraction, wc.k, attempt),
        runtime=lambda wc: a4_sat_resource(wc.fraction, wc.k, "runtime"),
        fit_slots=1,
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=A4_SAT_K_PATTERN,
        seed=A4_SAT_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.runner} --metrics-dir {params.metrics_dir} "
        "--requested-mem-mb {resources.mem_mb} --initial-mem-mb {params.initial_mem} "
        "--requested-cpus {threads} --requested-runtime-min {resources.runtime} "
        "--seed {params.rng_seed} -- "
        "Rscript {input.script} --fraction {wildcards.fraction} --clamp-k {wildcards.k} "
        "--seed-index {wildcards.seed} --seed {params.rng_seed} "
        "--metadata {input.metadata} --subsample-info {input.subsample} "
        "--fbm-backing {params.fbm} --svd {input.svd} "
        "--requested-mem-mb {resources.mem_mb} --requested-cpus {threads} "
        "--requested-runtime-min {resources.runtime} "
        "--out-dir {params.out_dir} > {log} 2>&1"


rule ora_bp_saturation_base:
    input:
        z=lambda wc: ancient(a4_sat_base_z(wc.fraction, wc.k, wc.seed)),
        subsample=lambda wc: f"{a4_sat_input_cell(wc.fraction, wc.seed)}/subsample_info.rds",
        database=lambda wc: a4_sat_db(wc)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(
            f"{A4_SAT_ORA_ROOT}/rs{{fraction}}/k{{k}}/seed{{seed}}/CLAMPbase/{{database}}"
        ),
    log:
        f"{A4_SAT_ORA_ROOT}/rs{{fraction}}/k{{k}}/seed{{seed}}/CLAMPbase/{{database}}.log"
    params:
        subsample=lambda wc: f"{a4_sat_input_cell(wc.fraction, wc.seed)}/subsample_info.rds",
        database_type=lambda wc: a4_sat_db(wc)["type"],
        database_label=lambda wc: a4_sat_db(wc)["label"],
        exclude=lambda wc: (
            f"--exclude-term-regex '{a4_sat_db(wc)['exclude_term_regex']}'"
            if "exclude_term_regex" in a4_sat_db(wc) else ""
        ),
        sheet=lambda wc: f"--sheet {a4_sat_db(wc)['sheet']}" if "sheet" in a4_sat_db(wc) else "",
        columns=lambda wc: (
            f"--term-column {a4_sat_db(wc)['term_column']} --gene-column {a4_sat_db(wc)['gene_column']}"
            if "term_column" in a4_sat_db(wc) else ""
        ),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        ora_mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
        ora_slots=1,
    wildcard_constraints:
        fraction=A4_SAT_FRACTION_PATTERN,
        k=A4_SAT_K_PATTERN,
        seed="|".join(map(str, A4_PW_INDICES)),
        database=A4_SAT_DATABASE_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset archs4 --fraction {wildcards.fraction} --seed {wildcards.seed} "
        "--model CLAMPbase --subsample-info {params.subsample} --database {wildcards.database} "
        "--database-label '{params.database_label}' --database-type {params.database_type} "
        "--database-path {input.database} {params.exclude} {params.sheet} {params.columns} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CFG[ora][min_size]} "
        "--max-size {A4_CFG[ora][max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


A4_SAT_FULL_DATA_ORA_ROOT = f"{A4_COV_ORA_ROOT}/archs4/rs100"
A4_SAT_FULL_DATA_ORA = [
    a4_cov_ora_dir("archs4", 100, seed, model, database)
    for seed in a4_cov_seeds_for(100)
    for model in ("CLAMPbase",)
    for database in A4_SAT_DATABASES
]
A4_SAT_PW_CELLS = [
    (f, k, i) for f in A4_SAT_FRACTIONS for k in A4_SAT_KS for i in A4_PW_NEW
] + sorted(A4_SAT_EXCLUDED)
A4_SAT_PW_ORA = [
    a4_sat_ora_dir(f, k, i, "CLAMPbase", db) for (f, k, i) in A4_SAT_PW_CELLS for db in A4_SAT_DATABASES
]
A4_SAT_PW_FULL_ORA = [
    a4_cov_ora_dir("archs4", 100, i, "CLAMPbase", db) for i in A4_PW_NEW for db in A4_SAT_DATABASES
]
A4_SAT_PW_LEAVES = {
    A4_SAT_CFG["base_leaf"].format(fraction=f, k=k, seed=i): (f, k, i) for (f, k, i) in A4_SAT_PW_CELLS
}
A4_SAT_EXTENDED_SVD = {int(f): [int(k) for k in ks] for f, ks in A4_PW_CFG["extended_svd"].items()}


def a4_sat_pw_svd(fraction, k, index):
    cell = a4_sat_input_cell(fraction, index)
    if int(k) in A4_SAT_EXTENDED_SVD.get(int(fraction), []):
        return f"{cell}/extra_svd/k{k}.rds"
    return f"{cell}/svd.rds"


# ============================================================
# Step 2: fit additional pathway replicates
# ============================================================

rule svd_saturation_extended:
    input:
        metadata=f"{{cell}}/metadata_cell.rds",
        subsample=f"{{cell}}/subsample_info.rds",
        fbm=f"{{cell}}/fbm_subsampled.bk",
        script="scripts/archs4/svd.R",
        common="scripts/archs4/common.R",
    output:
        svd=f"{{cell}}/extra_svd/k{{k}}.rds",
        k=f"{{cell}}/extra_svd/k{{k}}.K.rds",
    log:
        f"{{cell}}/extra_svd/k{{k}}.log"
    params:
        fbm_prefix=lambda wc: f"{wc.cell}/fbm_subsampled",
        n_cores=A4_CFG["svd"]["n_cores"],
        seed=A4_CFG["svd"]["seed"],
    threads: A4_CFG["svd"]["n_cores"]
    resources:
        mem_mb=lambda wc: a4_cov_initial_mem(a4_cov_fraction_of(wc.cell)),
        runtime=6000,
    wildcard_constraints:
        cell="|".join(re.escape(a4_sat_input_cell(f, i)) for f in A4_SAT_EXTENDED_SVD for i in A4_PW_INDICES),
        k="|".join(str(k) for ks in A4_SAT_EXTENDED_SVD.values() for k in ks),
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --metadata {input.metadata} --fbm {params.fbm_prefix} "
        "--svd {output.svd} --k {output.k} --svd-k {wildcards.k} --n-cores {params.n_cores} "
        "--seed {params.seed} > {log} 2>&1"


rule fit_saturation_pathway_base:
    input:
        unpack(lambda wc: a4_pw_inputs(f"{A4_SAT_BASE_ROOT}/{wc.pathway_leaf}/CLAMPbase/Z.csv", {
            "svd": a4_sat_pw_svd(*A4_SAT_PW_LEAVES[wc.pathway_leaf]),
            "subsample": f"{a4_sat_input_cell(A4_SAT_PW_LEAVES[wc.pathway_leaf][0], A4_SAT_PW_LEAVES[wc.pathway_leaf][2])}/subsample_info.rds",
            "metadata": f"{A4_PREP}/metadata_filtered.rds",
            "fbm": f"{a4_sat_input_cell(A4_SAT_PW_LEAVES[wc.pathway_leaf][0], A4_SAT_PW_LEAVES[wc.pathway_leaf][2])}/fbm_subsampled.bk",
            "script": "scripts/archs4/saturation/fit_clampbase.R",
        })),
    output:
        z=f"{A4_SAT_BASE_ROOT}/{{pathway_leaf}}/CLAMPbase/Z.csv",
    log:
        f"{A4_SAT_BASE_ROOT}/{{pathway_leaf}}/clampbase.log"
    params:
        cell=lambda wc: A4_SAT_PW_LEAVES[wc.pathway_leaf],
        input_cell=lambda wc: a4_sat_input_cell(A4_SAT_PW_LEAVES[wc.pathway_leaf][0], A4_SAT_PW_LEAVES[wc.pathway_leaf][2]),
        svd=lambda wc: a4_sat_pw_svd(*A4_SAT_PW_LEAVES[wc.pathway_leaf]),
        rng_seed=lambda wc: a4_cov_rng_seed(A4_SAT_PW_LEAVES[wc.pathway_leaf][2]),
        out_dir=lambda wc: f"{A4_SAT_BASE_ROOT}/{wc.pathway_leaf}",
    threads: A4_SAT_THREADS
    resources:
        mem_mb=lambda wc: a4_sat_resource(*A4_SAT_PW_LEAVES[wc.pathway_leaf][:2], "mem_mb"),
        runtime=lambda wc: a4_sat_resource(*A4_SAT_PW_LEAVES[wc.pathway_leaf][:2], "runtime"),
        fit_slots=1,
    wildcard_constraints:
        pathway_leaf="|".join(re.escape(leaf) for leaf in A4_SAT_PW_LEAVES),
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/archs4/saturation/fit_clampbase.R --fraction {params.cell[0]} --clamp-k {params.cell[1]} "
        "--seed-index {params.cell[2]} --seed {params.rng_seed} --initialization random --z-only true "
        "--metadata {A4_PREP}/metadata_filtered.rds --subsample-info {params.input_cell}/subsample_info.rds "
        "--fbm-backing {params.input_cell}/fbm_subsampled.bk --svd {params.svd} "
        "--requested-cpus {threads} --out-dir {params.out_dir} > {log} 2>&1"


ruleorder: fit_saturation_pathway_base > fit_bp_saturation_base


# ============================================================
# Step 3: summarize recovery and render the report
# ============================================================

rule aggregate_bp_saturation:
    input:
        A4_SAT_BASE_ORA_ALL,
        full_data=[ancient(d) for d in A4_SAT_FULL_DATA_ORA],
        pathway=lambda wc: A4_SAT_PW_ORA + A4_SAT_PW_FULL_ORA,
        script="scripts/archs4/saturation/aggregate_saturation.R",
    output:
        saturation_long=f"{A4_SAT_BIO}/saturation_long.csv",
        panel_ready=f"{A4_SAT_BIO}/saturation_panel_ready.csv",
    log:
        f"{A4_SAT_BIO}/aggregate.log"
    resources:
        mem_mb=16000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --ora-root {A4_SAT_ORA_ROOT} "
        "--full-ora-root {A4_SAT_FULL_DATA_ORA_ROOT} "
        "--saturation-out {output.saturation_long} --panel-out {output.panel_ready} "
        "--database-labels '{A4_PW_DATABASE_LABELS}' > {log} 2>&1"


rule saturation_report_archs4:
    input:
        saturation_long=rules.aggregate_bp_saturation.output.saturation_long,
        panel_ready=rules.aggregate_bp_saturation.output.panel_ready,
        notebook=f"{A4_SAT_NB}/00_saturation.ipynb",
    output:
        complete=touch(f"{A4_SAT_BIO}/notebook.complete"),
    log:
        notebook=f"{A4_SAT_NB}/00_saturation.executed.ipynb",
    params:
        saturation_dir=A4_SAT_BIO,
        ora_root=A4_SAT_ORA_ROOT,
        full_ora_root=A4_SAT_FULL_DATA_ORA_ROOT,
    resources:
        mem_mb=16000,
        runtime=60,
    conda: "clamp-analyses"
    notebook:
        f"{A4_SAT_NB}/00_saturation.ipynb"


rule saturation_bp_base_ora:
    input:
        A4_SAT_BASE_ORA_PRESENT,


rule saturation_bp_base_models:
    input:
        [a4_sat_base_rds(*c) for c in a4_sat_base_missing()],


rule saturation_bp_ora:
    input:
        A4_SAT_BASE_ORA_ALL,


rule saturation_pathway_ora:
    input:
        A4_SAT_PW_ORA,


rule saturation_bp:
    input:
        rules.saturation_report_archs4.output.complete,


rule archs4_saturation:
    input:
        rules.saturation_report_archs4.output.complete,

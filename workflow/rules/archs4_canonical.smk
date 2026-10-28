A4_CAN_CFG = A4_CFG["canonical"]
A4_CAN_MODEL_NAME = A4_CAN_CFG["model_name"]
A4_CAN_FINAL_ROOT = A4_CAN_CFG["final_root"]
A4_CAN_COMPENDIA = list(A4_CAN_CFG["compendia"])
A4_CAN_DATABASES = list(A4_CAN_CFG["ora"]["databases"])
A4_CAN_COMPENDIUM_PATTERN = "|".join(A4_CAN_COMPENDIA)
A4_CAN_DATABASE_PATTERN = "|".join(A4_CAN_DATABASES)


A4_CAN_PRIOR = A4_CFG["drug_diseases"]["prior"]
A4_CAN_PRIOR_MANIFEST = os.path.join(os.path.dirname(A4_CAN_PRIOR["output_gmt"]), "canonical_prior_manifest.json")
A4_CAN_FILES = [f"{A4_CAN_MODEL_NAME}.rds", "B.csv", "Z.csv", "summary.csv", "manifest.json", "complete"]


def a4_model_complete(out_dir, files):
    return all(os.path.isfile(os.path.join(out_dir, name)) and os.path.getsize(os.path.join(out_dir, name)) > 0 for name in files)


def a4_model_inputs(out_dir, files, inputs):
    present = [os.path.isfile(os.path.join(out_dir, name)) and os.path.getsize(os.path.join(out_dir, name)) > 0 for name in files]
    if all(present):
        return {}
    if any(present):
        missing = [name for name, exists in zip(files, present) if not exists]
        raise ValueError(f"Incomplete model bundle at {out_dir}; restore or move it before rebuilding. Missing: {', '.join(missing)}")
    if os.path.isdir(out_dir) and os.listdir(out_dir):
        raise ValueError(f"Refusing to fit into nonempty directory: {out_dir}; move it before rebuilding")
    return inputs


A4_CAN_FIT_PATTERN = "|".join(
    dataset for dataset in A4_CAN_COMPENDIA
    if not a4_model_complete(f"{A4_CAN_FINAL_ROOT}/{dataset}", A4_CAN_FILES)
) or r"(?!)"


def a4_can_fit_inputs(wc):
    inputs = dict(
        prior=A4_CAN_PRIOR["output_gmt"],
        prior_manifest=A4_CAN_PRIOR_MANIFEST,
        script="scripts/archs4/coverage/fit_clampfull.R",
        common="scripts/archs4/common.R",
    )
    if wc.dataset == "archs4":
        cell = a4_cov_cell(100, 1)
        inputs.update(metadata=f"{A4_PREP}/metadata_filtered.rds", raw_metadata=f"{A4_PREP}/metadata.rds",
                      subsample=f"{cell}/subsample_info.rds", fbm=f"{A4_PREP}/fbm_filtered.bk",
                      svd=f"{cell}/svd.rds", k=f"{cell}/CLAMP_K.rds", base=f"{cell}/CLAMPbase.rds")
    else:
        c = A4_COV_CFG["comparators"][wc.dataset]
        inputs.update(fbm=c["fbm_rds"], backing=c["fbm_backing"], genes=c["genes_rds"],
                      samples=c["samples_rds"], svd=c["svd"], k=c["k"], base=c["base"])
    return a4_model_inputs(f"{A4_CAN_FINAL_ROOT}/{wc.dataset}", A4_CAN_FILES, inputs)


def a4_can_expression_args(wc):
    if wc.dataset == "archs4":
        return (f"--metadata {A4_PREP}/metadata_filtered.rds "
                f"--subsample-info {a4_cov_cell(100, 1)}/subsample_info.rds "
                f"--fbm-backing {A4_PREP}/fbm_filtered.bk")
    c = A4_COV_CFG["comparators"][wc.dataset]
    return f"--fbm-rds {c['fbm_rds']} --genes-rds {c['genes_rds']} --samples-rds {c['samples_rds']}"


def a4_can_training_file(wc, name):
    if wc.dataset == "archs4":
        return f"{a4_cov_cell(100, 1)}/{dict(svd='svd.rds', k='CLAMP_K.rds', base='CLAMPbase.rds')[name]}"
    return A4_COV_CFG["comparators"][wc.dataset][name]


# ============================================================
# Step 1: build the canonical prior
# ============================================================

rule build_canonical_prior:
    input:
        c2_all=A4_CAN_PRIOR["c2_all"],
        reactome=A4_CAN_PRIOR["reactome"],
        cellmarker=A4_CAN_PRIOR["cellmarker"],
        script="scripts/archs4/drug_diseases/build_canonical_prior.py",
    output:
        gmt=A4_CAN_PRIOR["output_gmt"],
        manifest=A4_CAN_PRIOR_MANIFEST,
    conda: "clamp-analyses"
    shell:
        "python {input.script:q} --c2-all {input.c2_all:q} --reactome {input.reactome:q} "
        "--cellmarker {input.cellmarker:q} --sheet {A4_CAN_PRIOR[sheet]:q} "
        "--term-column {A4_CAN_PRIOR[term_column]:q} --gene-column {A4_CAN_PRIOR[gene_column]:q} "
        "--output-gmt {output.gmt:q} --manifest {output.manifest:q}"


# ============================================================
# Step 2: fit canonical models
# ============================================================

rule fit_canonical_clampfull:
    input:
        unpack(a4_can_fit_inputs),
    output:
        model=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/{A4_CAN_MODEL_NAME}.rds",
        b=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/B.csv",
        z=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/Z.csv",
        summary=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/summary.csv",
        manifest=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/manifest.json",
        complete=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/complete",
    log:
        f"{A4_PROD}/canonical_logs/{{dataset}}.log",
    params:
        out_dir=lambda wc: f"{A4_CAN_FINAL_ROOT}/{wc.dataset}",
        expression=a4_can_expression_args,
        svd=lambda wc: a4_can_training_file(wc, "svd"),
        k=lambda wc: a4_can_training_file(wc, "k"),
        base=lambda wc: a4_can_training_file(wc, "base"),
        seed=lambda wc: A4_CAN_CFG["fits"][wc.dataset]["seed"],
        max_iter=lambda wc: A4_CAN_CFG["fits"][wc.dataset]["max_iter"],
    threads: 16
    resources:
        mem_mb=lambda wc: A4_CAN_CFG["fits"][wc.dataset]["mem_mb"],
        runtime=lambda wc: A4_CAN_CFG["fits"][wc.dataset]["runtime"],
        large_fit_slot=lambda wc: 1 if wc.dataset == "archs4" else 0,
    wildcard_constraints:
        dataset=A4_CAN_FIT_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript scripts/archs4/coverage/fit_clampfull.R --dataset {wildcards.dataset} "
        "{params.expression} --prior-gmt {A4_CAN_PRIOR[output_gmt]:q} "
        "--svd {params.svd:q} --k {params.k:q} --base-model {params.base:q} "
        "--out-dir {params.out_dir:q} --model-name {A4_CAN_MODEL_NAME} --seed {params.seed} "
        "--seed-index 1 --fraction 100 --max-iter {params.max_iter} --multiplier 100 "
        "--requested-cpus {threads} --requested-mem-mb {resources.mem_mb} "
        "--requested-runtime-min {resources.runtime} > {log:q} 2>&1"


rule canonical_models:
    input:
        expand(rules.fit_canonical_clampfull.output, dataset=A4_CAN_COMPENDIA)


def a4_can_database(wildcards):
    return A4_CFG["ora"]["databases"][wildcards.database]


# ============================================================
# Step 3: evaluate pathway recovery
# ============================================================

rule ora_canonical:
    input:
        model=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/{A4_CAN_MODEL_NAME}.rds",
        z=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/Z.csv",
        manifest=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/manifest.json",
        database=lambda wc: a4_can_database(wc)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora/{{database}}"),
    log:
        f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora/{{database}}.log"
    params:
        database_type=lambda wc: a4_can_database(wc)["type"],
        database_label=lambda wc: a4_can_database(wc)["label"],
        exclude=lambda wc: (
            f"--exclude-term-regex '{a4_can_database(wc)['exclude_term_regex']}'"
            if "exclude_term_regex" in a4_can_database(wc) else ""
        ),
        sheet=lambda wc: f"--sheet {a4_can_database(wc)['sheet']}" if "sheet" in a4_can_database(wc) else "",
        columns=lambda wc: (
            f"--term-column {a4_can_database(wc)['term_column']} "
            f"--gene-column {a4_can_database(wc)['gene_column']}"
            if "term_column" in a4_can_database(wc) else ""
        ),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
    wildcard_constraints:
        dataset=A4_CAN_COMPENDIUM_PATTERN,
        database=A4_CAN_DATABASE_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset {wildcards.dataset} --fraction 100 --seed 1 "
        "--model {A4_CAN_MODEL_NAME} --model-manifest {input.manifest} "
        "--database {wildcards.database} --database-label '{params.database_label}' "
        "--database-type {params.database_type} --database-path {input.database} "
        "{params.exclude} {params.sheet} {params.columns} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {A4_CFG[ora][min_size]} "
        "--max-size {A4_CFG[ora][max_size]} --pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} "
        "--qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} --fdr {A4_CFG[ora][fdr][0]} "
        "> {log} 2>&1"


rule aggregate_canonical_ora:
    wildcard_constraints:
        dataset=A4_CAN_COMPENDIUM_PATTERN,
    input:
        ora=lambda wc: [
            f"{A4_CAN_FINAL_ROOT}/{wc.dataset}/ora/{database}"
            for database in A4_CAN_DATABASES
        ],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        coverage_long=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora.csv",
        cross_dataset=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora_cross.csv",
        panel_ready=f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora_panel.csv",
    log:
        f"{A4_CAN_FINAL_ROOT}/{{dataset}}/ora_aggregate.log"
    params:
        databases=",".join(A4_CAN_DATABASES)
    resources:
        mem_mb=16000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --ora-root {A4_CAN_FINAL_ROOT}/{wildcards.dataset} "
        "--coverage-out {output.coverage_long} --cross-out {output.cross_dataset} "
        "--panel-out {output.panel_ready} --databases {params.databases} "
        "> {log} 2>&1"


rule archs4_canonical_ora:
    input:
        expand(rules.aggregate_canonical_ora.output, dataset=A4_CAN_COMPENDIA)

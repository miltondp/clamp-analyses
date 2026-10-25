CON_CFG = config["consensus"]
CON_ROOT = CON_CFG["root"]
CON_NAME = CON_CFG["model_name"]
CON_SEEDS = sorted((int(s) for s in CON_CFG["seeds"]), key=lambda s: (s < 1000, s))
CON_COMPENDIA = list(CON_CFG["compendia"])
CON_FIT_COMPENDIA = [c for c in CON_COMPENDIA if c != "archs4"]
CON_DATABASES = list(CON_CFG["ora_databases"])
CON_RES = CON_CFG["resources"]

CON_GTEX_FROM_COVERAGE = {seed: rep for rep, seed in A4_COV_GTEX_REPLICATES.items()}
CON_STAB = {c: CON_CFG["compendia"][c]["stability_root"] for c in CON_FIT_COMPENDIA}
CON_STAB_OF = {root: c for c, root in CON_STAB.items()}
CON_STAB_PATTERN = "|".join(re.escape(root) for root in CON_STAB_OF)
CON_COMPENDIUM_PATTERN = "|".join(CON_COMPENDIA)
CON_FIT_PATTERN = "|".join(CON_FIT_COMPENDIA)
CON_SEED_PATTERN = "|".join(map(str, CON_SEEDS))


def con_ref(c):
    return CON_CFG["compendia"][c]["reference"]


def con_models(c):
    return [con_ref(c)] + [f"S{s}" for s in CON_SEEDS]


def con_seed_index(seed):
    return [int(s) for s in CON_CFG["seeds"]].index(int(seed)) + 2


def con_model_dir(c, model):
    if model == con_ref(c):
        return f"{A4_CAN_FINAL_ROOT}/{c}"
    if c == "archs4":
        return a4_seed_model(model[1:])
    return f"{CON_STAB[c]}/seed{model[1:]}/CLAMPfull_canonical"


def con_noise_dir(c, model):
    if c == "archs4":
        return A4_METALV_TN if model == con_ref(c) else f"{A4_METALV_NOISE}/seed{model[1:]}"
    return f"{CON_ROOT}/{c}/noise/{model}"


def con_expression_args(c):
    if c == "archs4":
        return f"--preprocess-dir {A4_PREP}"
    comp = A4_COV_CFG["comparators"][c]
    return (f"--fbm-rds {comp['fbm_rds']} --fbm-backing {comp['fbm_backing']} "
            f"--genes-rds {comp['genes_rds']} --samples-rds {comp['samples_rds']}")


def con_expression_inputs(c):
    if c == "archs4":
        return dict(metadata=ancient(f"{A4_PREP}/metadata_filtered.rds"), fbm=ancient(f"{A4_PREP}/fbm_filtered.bk"))
    comp = A4_COV_CFG["comparators"][c]
    return dict(fbm=ancient(comp["fbm_rds"]), backing=ancient(comp["fbm_backing"]),
                genes=ancient(comp["genes_rds"]), samples=ancient(comp["samples_rds"]))


def con_published_k(c):
    return A4_COV_CFG["comparators"][c].get("k") or f"{R2_DIR}/CLAMP_K_recount2.rds"


def con_covariate_args(c):
    cov = CON_CFG["compendia"][c]["covariates"]
    return (f"--categorical {','.join(cov['categorical'])} --numeric {','.join(cov['numeric'])} "
            f"--study {cov['study']} --min-level-n {cov['min_level_n']}")


# ============================================================
# Step 1: fit independent CLAMPbase replicates
# ============================================================

rule consensus_gtex_base_from_coverage:
    input:
        complete=lambda wc: ancient(f"{a4_cov_gtex_rep_dir(CON_GTEX_FROM_COVERAGE[int(wc.seed)])}/complete"),
    output:
        svd=f"{CON_STAB['gtex']}/seed{{seed}}/svd.rds",
        k=f"{CON_STAB['gtex']}/seed{{seed}}/CLAMP_K.rds",
        rds=f"{CON_STAB['gtex']}/seed{{seed}}/CLAMPbase.rds",
        z=f"{CON_STAB['gtex']}/seed{{seed}}/CLAMPbase/Z.csv",
        b=f"{CON_STAB['gtex']}/seed{{seed}}/CLAMPbase/B.csv",
    params:
        src=lambda wc: a4_cov_gtex_rep_dir(CON_GTEX_FROM_COVERAGE[int(wc.seed)]),
    wildcard_constraints:
        seed="|".join(map(str, CON_GTEX_FROM_COVERAGE)),
    shell:
        "for f in svd.rds CLAMP_K.rds CLAMPbase.rds CLAMPbase/Z.csv CLAMPbase/B.csv; do "
        "ln -sfr {params.src}/$f $(dirname {output.svd})/$f; done"


rule consensus_seed_base:
    input:
        unpack(lambda wc: con_expression_inputs(CON_STAB_OF[wc.root])),
        k=lambda wc: ancient(con_published_k(CON_STAB_OF[wc.root])),
        script=ancient("scripts/gtex/fit_clampbase_seed.R"),
    output:
        svd="{root}/seed{seed}/svd.rds",
        k="{root}/seed{seed}/CLAMP_K.rds",
        rds="{root}/seed{seed}/CLAMPbase.rds",
        z="{root}/seed{seed}/CLAMPbase/Z.csv",
        b="{root}/seed{seed}/CLAMPbase/B.csv",
    log:
        "{root}/logs/seed{seed}.base.log"
    params:
        compendium=lambda wc: CON_STAB_OF[wc.root],
        index=lambda wc: con_seed_index(wc.seed),
    threads: int(CON_RES["base"]["threads"])
    resources:
        mem_mb=int(CON_RES["base"]["mem_mb"]),
        runtime=int(CON_RES["base"]["runtime"]),
    wildcard_constraints:
        root=CON_STAB_PATTERN,
        seed=CON_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "export OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} OMP_NUM_THREADS={threads} && "
        "Rscript {input.script} --dataset {params.compendium} --fbm-rds {input.fbm} --backing {input.backing} "
        "--genes-rds {input.genes} --samples-rds {input.samples} --out-dir {wildcards.root}/seed{wildcards.seed} "
        "--seed {wildcards.seed} --seed-index {params.index} "
        "--expected-k $(Rscript -e 'cat(as.integer(readRDS(commandArgs(TRUE)[1])))' {input.k}) > {log} 2>&1"


ruleorder: consensus_gtex_base_from_coverage > consensus_seed_base


# ============================================================
# Step 2: fit canonical CLAMPfull replicates
# ============================================================

rule consensus_seed_full:
    input:
        unpack(lambda wc: a4_model_inputs(f"{wc.root}/seed{wc.seed}/CLAMPfull_canonical", A4_CAN_FILES, {
            **con_expression_inputs(CON_STAB_OF[wc.root]),
            "prior": ancient(A4_CAN_PRIOR["output_gmt"]),
            "svd": f"{wc.root}/seed{wc.seed}/svd.rds",
            "k": f"{wc.root}/seed{wc.seed}/CLAMP_K.rds",
            "base": f"{wc.root}/seed{wc.seed}/CLAMPbase.rds",
            "script": ancient("scripts/archs4/coverage/fit_clampfull.R"),
            "common": ancient("scripts/archs4/common.R"),
        })),
    output:
        model=f"{{root}}/seed{{seed}}/CLAMPfull_canonical/{A4_CAN_MODEL_NAME}.rds",
        b="{root}/seed{seed}/CLAMPfull_canonical/B.csv",
        z="{root}/seed{seed}/CLAMPfull_canonical/Z.csv",
        summary="{root}/seed{seed}/CLAMPfull_canonical/summary.csv",
        manifest="{root}/seed{seed}/CLAMPfull_canonical/manifest.json",
        complete="{root}/seed{seed}/CLAMPfull_canonical/complete",
    log:
        "{root}/logs/seed{seed}.full.log"
    params:
        compendium=lambda wc: CON_STAB_OF[wc.root],
        expression=lambda wc: con_expression_args(CON_STAB_OF[wc.root]),
        index=lambda wc: con_seed_index(wc.seed),
        max_iter=lambda wc: CON_CFG["compendia"][CON_STAB_OF[wc.root]]["fit"]["max_iter"],
    threads: lambda wc: int(CON_CFG["compendia"][CON_STAB_OF[wc.root]]["fit"]["threads"])
    resources:
        mem_mb=lambda wc: int(CON_CFG["compendia"][CON_STAB_OF[wc.root]]["fit"]["mem_mb"]),
        runtime=lambda wc: int(CON_CFG["compendia"][CON_STAB_OF[wc.root]]["fit"]["runtime"]),
    wildcard_constraints:
        root=CON_STAB_PATTERN,
        seed=CON_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "export OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} OMP_NUM_THREADS={threads} && "
        "Rscript scripts/archs4/coverage/fit_clampfull.R --dataset {params.compendium} {params.expression} "
        "--prior-gmt {A4_CAN_PRIOR[output_gmt]:q} --svd {wildcards.root}/seed{wildcards.seed}/svd.rds "
        "--k {wildcards.root}/seed{wildcards.seed}/CLAMP_K.rds "
        "--base-model {wildcards.root}/seed{wildcards.seed}/CLAMPbase.rds "
        "--out-dir {wildcards.root}/seed{wildcards.seed}/CLAMPfull_canonical --model-name {A4_CAN_MODEL_NAME} "
        "--seed {wildcards.seed} --seed-index {params.index} --fraction 100 --max-iter {params.max_iter} "
        "--multiplier 100 --requested-cpus {threads} --requested-mem-mb {resources.mem_mb} "
        "--requested-runtime-min {resources.runtime} > {log} 2>&1"


# ============================================================
# Step 3: screen technical noise
# ============================================================

rule consensus_gene_properties:
    input:
        z=lambda wc: ancient(f"{con_model_dir(wc.compendium, con_ref(wc.compendium))}/Z.csv"),
        ensembl=expand(f"{A4_TN_ENS}/{{name}}", name=A4_TN_URLS),
        script=A4_TN_SCRIPT,
        common=A4_TN_COMMON,
    output:
        f"{CON_ROOT}/{{compendium}}/noise/gene_properties.csv",
    log:
        f"{CON_ROOT}/{{compendium}}/noise/gene_properties.log"
    params:
        model=lambda wc: con_model_dir(wc.compendium, con_ref(wc.compendium)),
        out=lambda wc, output: os.path.dirname(output[0]),
    resources:
        mem_mb=int(CON_RES["noise"]["mem_mb"]),
        runtime=int(CON_RES["noise"]["runtime"]),
    wildcard_constraints:
        compendium=CON_FIT_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --stage gene_properties --model-dir {params.model} --out-dir {params.out} "
        "--ensembl-dir {A4_TN_ENS} --gene-map gtf > {log} 2>&1"


rule consensus_sample_covariates:
    input:
        samples=lambda wc: ancient(A4_COV_CFG["comparators"][wc.compendium]["samples_rds"]),
        source=lambda wc: ancient(CON_CFG["compendia"][wc.compendium].get("attributes")
                                  or CON_CFG["compendia"][wc.compendium]["rse_dir"]),
        script="scripts/consensus/sample_covariates.R",
    output:
        f"{CON_ROOT}/{{compendium}}/noise/sample_metadata.tsv.gz",
    log:
        f"{CON_ROOT}/{{compendium}}/noise/sample_metadata.log"
    params:
        source=lambda wc, input: (f"--attributes {input.source}" if wc.compendium == "gtex"
                                  else f"--rse-dir {input.source}"),
    threads: int(CON_RES["noise"]["threads"])
    resources:
        mem_mb=int(CON_RES["noise"]["mem_mb"]),
        runtime=int(CON_RES["noise"]["runtime"]),
    wildcard_constraints:
        compendium=CON_FIT_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --dataset {wildcards.compendium} --samples-rds {input.samples} "
        "{params.source} --n-cores {threads} --out {output} > {log} 2>&1"


rule consensus_model_noise:
    input:
        z=lambda wc: ancient(f"{con_model_dir(wc.compendium, wc.model)}/Z.csv"),
        b=lambda wc: ancient(f"{con_model_dir(wc.compendium, wc.model)}/B.csv"),
        annotation=lambda wc: ancient(f"{con_model_dir(wc.compendium, wc.model)}/summary.csv"),
        gene_properties=rules.consensus_gene_properties.output,
        sample_metadata=rules.consensus_sample_covariates.output,
        script=A4_TN_SCRIPT,
    output:
        gene=f"{CON_ROOT}/{{compendium}}/noise/{{model}}/lv_gene_properties_canonical.csv",
        sample=f"{CON_ROOT}/{{compendium}}/noise/{{model}}/lv_sample_metadata_canonical.csv",
    log:
        gene=f"{CON_ROOT}/{{compendium}}/noise/{{model}}/gene_level.log",
        sample=f"{CON_ROOT}/{{compendium}}/noise/{{model}}/sample_level.log",
    params:
        model=lambda wc: con_model_dir(wc.compendium, wc.model),
        out=lambda wc, output: os.path.dirname(output.gene),
        covariates=lambda wc: con_covariate_args(wc.compendium),
    threads: int(CON_RES["noise"]["threads"])
    resources:
        mem_mb=int(CON_RES["noise"]["mem_mb"]),
        runtime=int(CON_RES["noise"]["runtime"]),
    wildcard_constraints:
        compendium=CON_FIT_PATTERN,
        model=r"S\d+",
    conda: "clamp-analyses"
    shell:
        "ln -sfr {input.gene_properties} {params.out}/gene_properties.csv && "
        "ln -sfr {input.sample_metadata} {params.out}/sample_metadata.tsv.gz && "
        "Rscript {input.script} --stage gene_level --model-dir {params.model} --out-dir {params.out} > {log.gene} 2>&1 && "
        "Rscript {input.script} --stage sample_level --model-dir {params.model} --out-dir {params.out} "
        "{params.covariates} --n-cores {threads} > {log.sample} 2>&1"


def con_noise_files(c, model):
    return [f"{con_noise_dir(c, model)}/lv_gene_properties_canonical.csv",
            f"{con_noise_dir(c, model)}/lv_sample_metadata_canonical.csv"]


rule consensus_flags:
    input:
        noise=lambda wc: [f for m in con_models(wc.compendium) for f in con_noise_files(wc.compendium, m)],
        script="scripts/archs4/crispercas/cluster_metalv.py",
    output:
        f"{CON_ROOT}/{{compendium}}/noise/all_flags.csv",
    params:
        models=lambda wc: " ".join(f"{m}={con_noise_dir(wc.compendium, m)}" for m in con_models(wc.compendium)),
    wildcard_constraints:
        compendium=CON_FIT_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.script} flags --model {params.models} --threshold {CON_CFG[noise_r2]} --out {output}"


# ============================================================
# Step 4: cluster gene modules across models
# ============================================================

rule consensus_zcache:
    input:
        z=lambda wc: [ancient(f"{con_model_dir(wc.compendium, m)}/Z.csv") for m in con_models(wc.compendium)],
        script="scripts/archs4/crispercas/cluster_metalv.py",
    output:
        f"{CON_ROOT}/{{compendium}}/Z_all.npz",
    params:
        models=lambda wc: " ".join(con_models(wc.compendium)),
    resources:
        mem_mb=16000,
    wildcard_constraints:
        compendium=CON_FIT_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.script} cache --models {params.models} --z {input.z} --out {output}"


rule consensus_cluster:
    input:
        z=rules.consensus_zcache.output,
        flags=rules.consensus_flags.output,
        script="scripts/archs4/crispercas/cluster_metalv.py",
    output:
        metaz=f"{CON_ROOT}/{{compendium}}/metaZ.csv",
        summary=f"{CON_ROOT}/{{compendium}}/cluster_summary.csv",
        members=f"{CON_ROOT}/{{compendium}}/cluster_members.csv",
        all_lvs=f"{CON_ROOT}/{{compendium}}/all_lv_clusters.csv",
    log:
        f"{CON_ROOT}/{{compendium}}/cluster.log"
    params:
        out=lambda wc, output: os.path.dirname(output.metaz),
        reference=lambda wc: con_ref(wc.compendium),
    resources:
        mem_mb=16000,
    wildcard_constraints:
        compendium=CON_FIT_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.script} cluster --z-cache {input.z} --flags {input.flags} "
        "--min-cos {CON_CFG[min_cos]} --min-models {CON_CFG[min_models]} --reference {params.reference} "
        "--out-dir {params.out} > {log} 2>&1"


# ============================================================
# Step 5: build consensus models
# ============================================================

rule consensus_model:
    input:
        metaz=rules.consensus_cluster.output.metaz,
        summary=rules.consensus_cluster.output.summary,
        reference=lambda wc: ancient(f"{A4_CAN_FINAL_ROOT}/{wc.compendium}/{A4_CAN_MODEL_NAME}.rds"),
        manifest=lambda wc: ancient(f"{A4_CAN_FINAL_ROOT}/{wc.compendium}/manifest.json"),
        script="scripts/consensus/build_consensus_model.R",
    output:
        z=f"{CON_ROOT}/{{compendium}}/Z.csv",
        mapping=f"{CON_ROOT}/{{compendium}}/mgm_mapping.csv",
        manifest=f"{CON_ROOT}/{{compendium}}/manifest.json",
    params:
        out=lambda wc, output: os.path.dirname(output.z),
        models=lambda wc: ",".join(con_models(wc.compendium)),
    resources:
        mem_mb=32000,
    wildcard_constraints:
        compendium=CON_COMPENDIUM_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --dataset {wildcards.compendium} --metaz {input.metaz} --summary {input.summary} "
        "--reference-model {input.reference} --reference-manifest {input.manifest} --models {params.models} "
        "--min-cos {CON_CFG[min_cos]} --min-models {CON_CFG[min_models]} --model-name {CON_NAME} "
        "--out-dir {params.out} && ln -sfr {input.metaz} {output.z}"


CON_MGM_FILES = {
    "metaZ.csv": rules.crispercas_metalv_cluster.output.metaz,
    "cluster_summary.csv": rules.crispercas_metalv_cluster.output.summary,
    "cluster_members.csv": rules.crispercas_metalv_cluster.output.members,
    "all_lv_clusters.csv": rules.crispercas_metalv_cluster.output.all_lvs,
    "noise/all_flags.csv": rules.crispercas_noise_flags.output[0],
}


rule consensus_archs4_from_mgm:
    input:
        lambda wc: CON_MGM_FILES[wc.name],
    output:
        f"{CON_ROOT}/archs4/{{name}}",
    wildcard_constraints:
        name="|".join(re.escape(name) for name in CON_MGM_FILES),
    shell:
        "ln -sfr {input} {output}"


rule consensus_archs4_ora_from_mgm:
    input:
        lambda wc: f"{A4_METALV}/ora/{wc.database}",
    output:
        directory(f"{CON_ROOT}/archs4/ora/{{database}}"),
    wildcard_constraints:
        database="|".join(CON_DATABASES),
    shell:
        "ln -sfnr {input} {output}"


rule consensus_model_links:
    input:
        z=lambda wc: [ancient(f"{con_model_dir(wc.compendium, m)}/Z.csv") for m in con_models(wc.compendium)],
    output:
        touch(f"{CON_ROOT}/{{compendium}}/models/links.complete"),
    params:
        links=lambda wc: " ".join(f"{m}={con_model_dir(wc.compendium, m)}" for m in con_models(wc.compendium)),
        out=lambda wc, output: os.path.dirname(output[0]),
    wildcard_constraints:
        compendium=CON_COMPENDIUM_PATTERN,
    shell:
        "for spec in {params.links}; do ln -sfnr ${{spec#*=}} {params.out}/${{spec%%=*}}; done"


# ============================================================
# Step 6: project samples, enrich pathways, project drugs and traits
# ============================================================

rule consensus_project_samples:
    input:
        unpack(lambda wc: con_expression_inputs(wc.compendium)),
        z=rules.consensus_model.output.z,
        manifest=rules.consensus_model.output.manifest,
        script="scripts/consensus/project_samples.R",
    output:
        f"{CON_ROOT}/{{compendium}}/B.csv",
    log:
        f"{CON_ROOT}/{{compendium}}/projection.log"
    params:
        expression=lambda wc: con_expression_args(wc.compendium),
    threads: int(CON_RES["projection"]["threads"])
    resources:
        mem_mb=int(CON_RES["projection"]["mem_mb"]),
        runtime=int(CON_RES["projection"]["runtime"]),
    wildcard_constraints:
        compendium=CON_COMPENDIUM_PATTERN,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --manifest {input.manifest} {params.expression} --n-cores {threads} "
        "--out {output} > {log} 2>&1"


rule consensus_ora:
    input:
        z=rules.consensus_model.output.z,
        manifest=rules.consensus_model.output.manifest,
        database=lambda wc: a4_metalv_ora_args(wc.database)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(f"{CON_ROOT}/{{compendium}}/ora/{{database}}"),
    log:
        f"{CON_ROOT}/{{compendium}}/ora/{{database}}.log"
    params:
        args=lambda wc: a4_metalv_ora_args(wc.database),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
    wildcard_constraints:
        compendium=CON_FIT_PATTERN,
        database="|".join(CON_DATABASES),
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset {wildcards.compendium} --fraction 100 --seed 1 "
        "--model {CON_NAME} --model-manifest {input.manifest} "
        "--database {wildcards.database} --database-label '{params.args[label]}' "
        "--database-type {params.args[type]} --database-path {input.database} {params.args[extra]} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {params.args[min_size]} --max-size {params.args[max_size]} "
        "--pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} --qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} "
        "--fdr {A4_CFG[ora][fdr][0]} > {log} 2>&1"


rule consensus_ora_aggregate:
    input:
        ora=lambda wc: [f"{CON_ROOT}/{wc.compendium}/ora/{db}" for db in CON_DATABASES],
    output:
        f"{CON_ROOT}/{{compendium}}/ora.csv",
    params:
        summaries=lambda wc, input: " ".join(f"{d}/summary.csv" for d in input.ora),
    wildcard_constraints:
        compendium=CON_COMPENDIUM_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python -c \"import sys, pandas as pd; "
        "pd.concat([pd.read_csv(f) for f in sys.argv[2:]]).to_csv(sys.argv[1], index=False)\" "
        "{output} {params.summaries}"


rule consensus_lincs:
    input:
        z=rules.consensus_model.output.z,
        manifest=rules.consensus_model.output.manifest,
        lincs=f"{A4_DD_DATA}/lincs-data.pkl",
        script="scripts/archs4/drug_diseases/project_lincs.py",
    output:
        f"{CON_ROOT}/{{compendium}}/lincs-projection.pkl",
    resources:
        mem_mb=16000,
        runtime=360,
    wildcard_constraints:
        compendium=CON_COMPENDIUM_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.script} --z {input.z} "
        "--l2 $(python -c \"import json; print(repr(json.load(open('{input.manifest}'))['L2']))\") "
        "--lincs-file {input.lincs} --compendium {wildcards.compendium} --output {output}"


rule consensus_lincs_tsv:
    input:
        lincs=rules.consensus_lincs.output,
        drugbank=rules.fetch_drugbank_drug_diseases.output,
        script="scripts/archs4/drug_diseases/export_lincs_projection_tsv.py",
    output:
        f"{CON_ROOT}/{{compendium}}/lincs-projection.tsv",
    resources:
        mem_mb=8000,
        runtime=60,
    wildcard_constraints:
        compendium=CON_COMPENDIUM_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.script} --lincs-projection {input.lincs} --drugbank {input.drugbank} --output {output}"


rule consensus_gls_input:
    input:
        z=rules.consensus_model.output.z,
    output:
        z=temp(f"{CON_ROOT}/{{compendium}}/traits/phenoplier_Z.pkl"),
        mapping=f"{CON_ROOT}/{{compendium}}/traits/lv_mapping.csv",
    wildcard_constraints:
        compendium=CON_COMPENDIUM_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python -c \"import pandas as pd; z = pd.read_csv('{input.z}', index_col=0); "
        "m = pd.DataFrame({{'phenoplier_lv': [f'LV{{i + 1}}' for i in range(z.shape[1])], 'metaLV': z.columns}}); "
        "z.columns = m.phenoplier_lv; z.to_pickle('{output.z}'); m.to_csv('{output.mapping}', index=False)\""


rule consensus_gls:
    input:
        rds=rules.consensus_gls_input.output.z,
        script="scripts/phenoplier/run_gls.sh",
        verify="scripts/phenoplier/verify_summary.py",
    output:
        summary=f"{CON_ROOT}/{{compendium}}/traits/gls_summary.tsv.gz",
    log:
        f"{CON_ROOT}/{{compendium}}/traits/gls.log"
    params:
        rds=lambda wc, input: input.rds,
        name=lambda wc: f"consensus_{wc.compendium}_{CON_NAME}",
    threads: int(PHENOPLIER_GLS_RES["threads"])
    resources:
        mem_mb=int(PHENOPLIER_GLS_RES["mem_mb"]),
        runtime=int(PHENOPLIER_GLS_RES["runtime"]),
    wildcard_constraints:
        compendium=CON_COMPENDIUM_PATTERN,
    conda: PHENOPLIER_CFG["conda_env"]
    shell:
        PHENOPLIER_GLS_SHELL


rule consensus_traits:
    input:
        summary=rules.consensus_gls.output.summary,
        mapping=rules.consensus_gls_input.output.mapping,
        script="scripts/phenoplier/metalv/collect_summary.py",
    output:
        f"{CON_ROOT}/{{compendium}}/traits/consensus_{{compendium}}.tsv.gz",
    wildcard_constraints:
        compendium=CON_COMPENDIUM_PATTERN,
    conda: "clamp-analyses"
    shell:
        "python {input.script} --summary {input.summary} --mapping {input.mapping} --out {output}"


CON_CORE = [
    f"{CON_ROOT}/{c}/{name}"
    for c in CON_COMPENDIA
    for name in ["models/links.complete", "B.csv", "ora.csv", "lincs-projection.tsv", "manifest.json"]
]


rule consensus_models:
    input:
        CON_CORE,

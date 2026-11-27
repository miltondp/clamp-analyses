import re

# Evaluation of consensus LVs (metaLVs) and ensemble clustering (rule prefix cev_).
# Included after consensus.smk / archs4_crispercas.smk: reuses A4_METALV_MODELS and con_ref.
# All outputs live under consensus_eval.root, never in the canonical/stability model
# folders (see the completeness gating in archs4_canonical.smk).

CEV_CFG = config["consensus_eval"]
CEV_RES = CEV_CFG["resources"]
CEV_A4_CFG = CEV_CFG["compendia"]["archs4"]
CEV_A4 = f"{CEV_CFG['root']}/archs4"
CEV_A4_RUNS = f"{CEV_A4}/runs"
CEV_A4_REF = con_ref("archs4")
CEV_A4_REF_DIR = CEV_A4_CFG["reference_dir"].format(model=CEV_A4_REF)
CEV_REPRO = CEV_CFG["reproduce"]
CEV_A4_R2 = sorted({str(v["noise_r2"]) for v in CEV_REPRO.values()})
CEV_SCRIPT = "scripts/consensus_eval/run.py"
CEV_METALV_SCRIPT = "scripts/archs4/crispercas/cluster_metalv.py"
assert all(v["linkage"] == "complete" for v in CEV_REPRO.values()), "cluster_metalv.py only does complete linkage"


# ============================================================
# Step 0: export Z / summary / manifest from the saved fits
# ============================================================

rule cev_export_run:
    input:
        rds=ancient(lambda wc: CEV_A4_CFG["rds"].format(seed=wc.model[1:])),
        md5=ancient(CEV_A4_CFG["md5"]),
        script="scripts/consensus_eval/export_run.R",
    output:
        z=f"{CEV_A4_RUNS}/{{model}}/Z.csv",
        summary=f"{CEV_A4_RUNS}/{{model}}/summary.csv",
        manifest=f"{CEV_A4_RUNS}/{{model}}/manifest.json",
    wildcard_constraints:
        model="|".join(A4_METALV_MODELS),
    params:
        out_dir=lambda wc, output: os.path.dirname(output.z),
    resources:
        mem_mb=CEV_RES["export"]["mem_mb"],
        runtime=CEV_RES["export"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --rds {input.rds} --model {wildcards.model} "
        "--md5-file {input.md5} --out-dir {params.out_dir}"


rule cev_export_check:
    """The reference-model export must equal the published model byte for byte."""
    input:
        z=f"{CEV_A4_RUNS}/{CEV_A4_REF}/Z.csv",
        summary=f"{CEV_A4_RUNS}/{CEV_A4_REF}/summary.csv",
        ref_z=ancient(f"{CEV_A4_REF_DIR}/Z.csv"),
        ref_summary=ancient(f"{CEV_A4_REF_DIR}/summary.csv"),
    output:
        f"{CEV_A4}/checks/export_{CEV_A4_REF}.ok",
    shell:
        "cmp {input.z} {input.ref_z} && cmp {input.summary} {input.ref_summary} && "
        "echo '{CEV_A4_REF} Z.csv and summary.csv match {CEV_A4_REF_DIR}' > {output}"


# ============================================================
# Step 1: stage the LV cache and noise flags
# ============================================================

rule cev_zcache:
    """Stack the run Z.csv files in A4_METALV_MODELS order (decides tie-breaks and metaLV numbering)."""
    input:
        z=expand(rules.cev_export_run.output.z, model=A4_METALV_MODELS),
        script=CEV_METALV_SCRIPT,
    output:
        f"{CEV_A4}/Z_all.npz",
    resources:
        mem_mb=CEV_RES["zcache"]["mem_mb"],
        runtime=CEV_RES["zcache"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} cache --models {A4_METALV_MODELS} --z {input.z} --out {output}"


rule cev_zcache_check:
    """The cache must hold exactly the (model, LV) pairs of the legacy flags, all genes, no NaNs."""
    input:
        z=rules.cev_zcache.output,
        flags=ancient(CEV_A4_CFG["legacy_flags"]),
        script=CEV_SCRIPT,
    output:
        f"{CEV_A4}/checks/zcache.ok",
    resources:
        mem_mb=CEV_RES["zcache"]["mem_mb"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} zcache-check --z-cache {input.z} --flags {input.flags} "
        "--n-genes {CEV_A4_CFG[genes]} --out {output}"


rule cev_flags:
    """Re-threshold the legacy technical_joint (5-fold CV R2) at noise_r2; gene_flag is fixed."""
    input:
        legacy=ancient(CEV_A4_CFG["legacy_flags"]),
        script=CEV_SCRIPT,
    output:
        f"{CEV_A4}/noise/flags_r{{r2}}.csv",
    wildcard_constraints:
        r2="|".join(re.escape(r) for r in CEV_A4_R2),
    conda: "clamp-analyses"
    shell:
        "python {input.script} flags --legacy {input.legacy} --r2 {wildcards.r2} --out {output}"


# ============================================================
# Step 2: reproduce the legacy consensus sets with cluster_metalv.py
# ============================================================

rule cev_reproduce_set:
    """Unchanged legacy clustering; the reference that the cev.py builder must match (PLAN 2.1)."""
    input:
        z=rules.cev_zcache.output,
        flags=lambda wc: f"{CEV_A4}/noise/flags_r{CEV_REPRO[wc.set]['noise_r2']}.csv",
        check=rules.cev_zcache_check.output,
        script=CEV_METALV_SCRIPT,
    output:
        metaz=f"{CEV_A4}/reproduce/{{set}}/metaZ.csv",
        summary=f"{CEV_A4}/reproduce/{{set}}/cluster_summary.csv",
        members=f"{CEV_A4}/reproduce/{{set}}/cluster_members.csv",
        all_lvs=f"{CEV_A4}/reproduce/{{set}}/all_lv_clusters.csv",
    log:
        f"{CEV_A4}/reproduce/{{set}}/cluster.log",
    params:
        out=lambda wc, output: os.path.dirname(output.metaz),
        min_cos=lambda wc: CEV_REPRO[wc.set]["min_cos"],
        min_models=lambda wc: CEV_REPRO[wc.set]["min_models"],
    wildcard_constraints:
        set="|".join(CEV_REPRO),
    resources:
        mem_mb=CEV_RES["set"]["mem_mb"],
        runtime=CEV_RES["set"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} cluster --z-cache {input.z} --flags {input.flags} "
        "--min-cos {params.min_cos} --min-models {params.min_models} --reference {CEV_A4_REF} "
        "--out-dir {params.out} > {log} 2>&1"


rule cev_reproduction:
    """Fails unless names match, max|dZ| < 1e-6 and every metaLV has the legacy members."""
    input:
        metaz=rules.cev_reproduce_set.output.metaz,
        summary=rules.cev_reproduce_set.output.summary,
        external=lambda wc: ancient(CEV_REPRO[wc.set]["external"]),
        members=lambda wc: ancient(CEV_REPRO[wc.set]["members"]),
        script=CEV_SCRIPT,
    output:
        f"{CEV_A4}/report/reproduction/{{set}}.csv",
    log:
        f"{CEV_A4}/report/reproduction/{{set}}.log",
    wildcard_constraints:
        set="|".join(CEV_REPRO),
    resources:
        mem_mb=CEV_RES["set"]["mem_mb"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} reproduction --set {wildcards.set} --metaz {input.metaz} "
        "--summary {input.summary} --external {input.external} --members {input.members} "
        "--out {output} > {log} 2>&1"


# ============================================================
# Aggregate
# ============================================================

rule consensus_eval_archs4:
    input:
        expand(rules.cev_export_run.output, model=A4_METALV_MODELS),
        rules.cev_export_check.output,
        rules.cev_zcache_check.output,
        expand(rules.cev_reproduction.output, set=CEV_REPRO),

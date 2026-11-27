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
# Aggregate
# ============================================================

rule consensus_eval_archs4:
    input:
        expand(rules.cev_export_run.output, model=A4_METALV_MODELS),
        rules.cev_export_check.output,

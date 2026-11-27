import itertools
import json
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
CEV_NB = os.path.join(REPO_ROOT, CEV_CFG["notebooks"])
CEV_A4_REF = con_ref("archs4")
CEV_A4_REF_DIR = CEV_A4_CFG["reference_dir"].format(model=CEV_A4_REF)
CEV_REPRO = CEV_CFG["reproduce"]
CEV_EVAL = CEV_CFG["eval"]
CEV_SCRIPT = "scripts/consensus_eval/run.py"
CEV_LIB = "scripts/consensus_eval/cev.py"
CEV_SUBSETS = ["full"]   # split-half / LORO subsets arrive in PLAN 2.4
CEV_METALV_SCRIPT = "scripts/archs4/crispercas/cluster_metalv.py"
assert all(v["linkage"] == "complete" for v in CEV_REPRO.values()), "cluster_metalv.py only does complete linkage"


def _cev_grid_name(prefix, p, sim_key, sim_tag):
    name = (f"{prefix}_{p['linkage']}_{sim_tag}{round(p[sim_key] * 100):03d}_f{round(p['min_frac'] * 10):02d}"
            f"_r{round(p['noise_r2'] * 100):03d}_{p['combine']}")
    assert round(p["min_frac"] * 10) == p["min_frac"] * 10, f"min_frac {p['min_frac']} does not fit the set name"
    return name


def cev_hier_name(p):
    return _cev_grid_name("hier", p, "min_cos", "c")


def cev_topj_name(p):
    return _cev_grid_name("topj", p, "min_jac", "j")


# Every hierarchical set: the frozen reproduce sets (absolute min_models) plus the grids (min_frac).
# dist picks the linkage: cos (1 - cosine of raw loadings) or jac (1 - top-gene Jaccard).
CEV_SETS = {s: {"dist": "cos", "linkage": v["linkage"], "min_cos": v["min_cos"], "min_models": v["min_models"],
                "noise_r2": v["noise_r2"], "combine": "median"} for s, v in CEV_REPRO.items()}
for _grid_key, _dist, _namer in [("hier", "cos", cev_hier_name), ("topj", "jac", cev_topj_name)]:
    _grid = CEV_CFG["grids"][_grid_key]
    for _vals in itertools.product(*_grid.values()):
        _p = {"dist": _dist, **dict(zip(_grid, _vals))}
        _name = _namer(_p)
        assert _name not in CEV_SETS, f"duplicate set name {_name}"
        CEV_SETS[_name] = _p
# Set-definition columns of metrics.csv (meta.json keys), in first-seen order: they lead headline.csv.
CEV_SET_PARAMS = ["set", "subset", *dict.fromkeys(k for v in CEV_SETS.values() for k in v), "min_models_resolved"]
CEV_A4_R2 = sorted({str(v["noise_r2"]) for v in CEV_SETS.values()} | {str(CEV_EVAL["pool_r2"])})


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
        mem_mb=CEV_RES["check"]["mem_mb"],
        runtime=CEV_RES["check"]["runtime"],
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
        "python {input.script} flags --legacy {input.legacy} --legacy-r2 {CEV_A4_CFG[legacy_flags_r2]} --r2 {wildcards.r2} --out {output}"


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
        counts=f"{CEV_A4}/reproduce/{{set}}/cluster_model_count_distribution.csv",
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
    """Fails unless names match, max|dZ| < reproduce_tol and every metaLV has the legacy members."""
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
        mem_mb=CEV_RES["check"]["mem_mb"],
        runtime=CEV_RES["check"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} reproduction --set {wildcards.set} --metaz {input.metaz} "
        "--summary {input.summary} --external {input.external} --members {input.members} "
        "--tol {CEV_CFG[reproduce_tol]} "
        "--out {output} > {log} 2>&1"


# ============================================================
# Step 3: hierarchical sets (cev.py builder) and computational metrics
# ============================================================

rule cev_linkage:
    """One dendrogram per (noise_r2, dist, linkage, subset), shared by every similarity / min_frac cut."""
    input:
        z=rules.cev_zcache.output,
        flags=f"{CEV_A4}/noise/flags_r{{r2}}.csv",
        check=rules.cev_zcache_check.output,
        script=CEV_SCRIPT,
        lib=CEV_LIB,
    output:
        f"{CEV_A4}/linkage/r{{r2}}/{{dist}}_{{linkage}}/{{subset}}.npz",
    log:
        f"{CEV_A4}/linkage/r{{r2}}/{{dist}}_{{linkage}}/{{subset}}.log",
    wildcard_constraints:
        r2="|".join(re.escape(r) for r in CEV_A4_R2),
        dist="cos|jac",
        linkage="|".join(sorted({v["linkage"] for v in CEV_SETS.values()})),
        subset="|".join(CEV_SUBSETS),
    threads: CEV_RES["linkage"]["threads"]
    resources:
        mem_mb=CEV_RES["linkage"]["mem_mb"],
        runtime=CEV_RES["linkage"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} python {input.script} linkage "
        "--z-cache {input.z} --flags {input.flags} --method {wildcards.linkage} --dist {wildcards.dist} "
        "--top-n {A4_CRISPERCAS_CFG[top_n]} --out {output} > {log} 2>&1"


rule cev_hier_set:
    input:
        z=rules.cev_zcache.output,
        flags=lambda wc: f"{CEV_A4}/noise/flags_r{CEV_SETS[wc.set]['noise_r2']}.csv",
        linkage=lambda wc: (f"{CEV_A4}/linkage/r{CEV_SETS[wc.set]['noise_r2']}/"
                            f"{CEV_SETS[wc.set]['dist']}_{CEV_SETS[wc.set]['linkage']}/{wc.subset}.npz"),
        script=CEV_SCRIPT,
        lib=CEV_LIB,
    output:
        z=f"{CEV_A4}/sets/{{set}}/{{subset}}/Z.npy",
        lvs=f"{CEV_A4}/sets/{{set}}/{{subset}}/lvs.csv",
        members=f"{CEV_A4}/sets/{{set}}/{{subset}}/members.csv",
        meta=f"{CEV_A4}/sets/{{set}}/{{subset}}/meta.json",
    log:
        f"{CEV_A4}/sets/{{set}}/{{subset}}/build.log",
    params:
        out=lambda wc, output: os.path.dirname(output.z),
        json=lambda wc: json.dumps(CEV_SETS[wc.set]),
    wildcard_constraints:
        set="|".join(CEV_SETS),
        subset="|".join(CEV_SUBSETS),
    resources:
        mem_mb=CEV_RES["set"]["mem_mb"],
        runtime=CEV_RES["set"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} hier-set --set {wildcards.set} --subset {wildcards.subset} "
        "--z-cache {input.z} --flags {input.flags} --linkage {input.linkage} "
        "--params '{params.json}' --out-dir {params.out} > {log} 2>&1"


rule cev_set_check:
    """cev.py builder on real data == unchanged cluster_metalv.py output (the reproduce sets)."""
    input:
        set=lambda wc: expand(rules.cev_hier_set.output, set=wc.set, subset="full"),
        reproduce=rules.cev_reproduce_set.output,
        script=CEV_SCRIPT,
    output:
        f"{CEV_A4}/checks/set_{{set}}.ok",
    params:
        set_dir=lambda wc, input: os.path.dirname(input.set[0]),
        reproduce_dir=lambda wc, input: os.path.dirname(input.reproduce[0]),
    wildcard_constraints:
        set="|".join(CEV_REPRO),
    resources:
        mem_mb=CEV_RES["check"]["mem_mb"],
        runtime=CEV_RES["check"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} set-check --set-dir {params.set_dir} --reproduce-dir {params.reproduce_dir} "
        "--out {output}"


rule cev_pool_hits:
    """Best hit of every evaluation-pool LV in each other run, its Jaccard null and the reproducible LVs."""
    input:
        z=rules.cev_zcache.output,
        flags=f"{CEV_A4}/noise/flags_r{{r2}}.csv",
        script=CEV_SCRIPT,
        lib=CEV_LIB,
    output:
        hits=f"{CEV_A4}/pool/r{{r2}}/best_hits.npz",
        null=f"{CEV_A4}/pool/r{{r2}}/null.json",
    params:
        repro=CEV_EVAL["reproducible"],
    log:
        f"{CEV_A4}/pool/r{{r2}}/best_hits.log",
    wildcard_constraints:
        r2="|".join(re.escape(r) for r in CEV_A4_R2),
    threads: CEV_RES["eval"]["threads"]
    resources:
        mem_mb=CEV_RES["eval"]["mem_mb"],
        runtime=CEV_RES["eval"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} python {input.script} pool-hits "
        "--z-cache {input.z} --flags {input.flags} --top-n {A4_CRISPERCAS_CFG[top_n]} "
        "--reproducible-cos {params.repro[cos]} --reproducible-frac {params.repro[min_frac]} "
        "--jac-null-q {params.repro[jac_null_q]} --null-pairs {params.repro[null_pairs]} "
        "--null-seed {params.repro[null_seed]} --out {output.hits} --null-out {output.null} > {log} 2>&1"


rule cev_evaluate:
    input:
        set=lambda wc: expand(rules.cev_hier_set.output, set=wc.set, subset="full"),
        z=rules.cev_zcache.output,
        hits=f"{CEV_A4}/pool/r{CEV_EVAL['pool_r2']}/best_hits.npz",
        script=CEV_SCRIPT,
        lib=CEV_LIB,
    output:
        metrics=f"{CEV_A4}/sets/{{set}}/eval/metrics.csv",
        per_cluster=f"{CEV_A4}/sets/{{set}}/eval/per_cluster.csv",
        coverage=f"{CEV_A4}/sets/{{set}}/eval/coverage.csv",
    log:
        f"{CEV_A4}/sets/{{set}}/eval/evaluate.log",
    params:
        set_dir=lambda wc, input: os.path.dirname(input.set[0]),
        out=lambda wc, output: os.path.dirname(output.metrics),
        thresholds=" ".join(str(t) for t in CEV_EVAL["cos_thresholds"]),
    wildcard_constraints:
        set="|".join(CEV_SETS),
    threads: CEV_RES["eval"]["threads"]
    resources:
        mem_mb=CEV_RES["eval"]["mem_mb"],
        runtime=CEV_RES["eval"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} python {input.script} evaluate "
        "--set-dir {params.set_dir} --z-cache {input.z} --best-hits {input.hits} "
        "--pool-r2 {CEV_EVAL[pool_r2]} --top-n {A4_CRISPERCAS_CFG[top_n]} --cos-thresholds {params.thresholds} "
        "--single-gene {CEV_EVAL[single_gene]} --out-dir {params.out} > {log} 2>&1"


rule cev_aggregate:
    input:
        metrics=expand(rules.cev_evaluate.output.metrics, set=CEV_SETS),
        script=CEV_SCRIPT,
    output:
        all=f"{CEV_A4}/report/metrics_all.csv",
        headline=f"{CEV_A4}/report/headline.csv",
    params:
        set_params=" ".join(CEV_SET_PARAMS),
        headline=" ".join(CEV_EVAL["headline"]),
    conda: "clamp-analyses"
    shell:
        "python {input.script} aggregate --metrics {input.metrics} --params {params.set_params} "
        "--headline {params.headline} "
        "--out {output.all} --headline-out {output.headline}"


rule cev_report_00:
    """Temporary internal check, removed in PLAN 4.3 (legacy cache/noise dirs were never found)."""
    input:
        reports=expand(rules.cev_reproduction.output, set=CEV_REPRO),
        counts=expand(rules.cev_reproduce_set.output.counts, set=CEV_REPRO),
        notebook=f"{CEV_NB}/00_reproduction.ipynb",
    output:
        complete=touch(f"{CEV_A4}/report/00_reproduction.complete"),
    log:
        notebook=f"{CEV_NB}/00_reproduction.executed.ipynb",
    params:
        sets=list(CEV_REPRO),
        min_models=[CEV_REPRO[s]["min_models"] for s in CEV_REPRO],
    resources:
        mem_mb=CEV_RES["check"]["mem_mb"],
        runtime=CEV_RES["check"]["runtime"],
    conda: "clamp-analyses"
    notebook:
        f"{CEV_NB}/00_reproduction.ipynb"


# ============================================================
# Aggregate
# ============================================================

rule consensus_eval_archs4:
    input:
        expand(rules.cev_export_run.output, model=A4_METALV_MODELS),
        rules.cev_export_check.output,
        rules.cev_zcache_check.output,
        expand(rules.cev_reproduction.output, set=CEV_REPRO),
        rules.cev_report_00.output,   # temporary, removed in PLAN 4.3
        expand(rules.cev_set_check.output, set=CEV_REPRO),
        rules.cev_aggregate.output,

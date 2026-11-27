import itertools
import json
import math
import random
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
# ancient(): editing the code must not rerun the area; after a logic change, rerun with --forcerun <rule>.
CEV_SCRIPT = ancient("scripts/consensus_eval/run.py")
CEV_LIB = ancient("scripts/consensus_eval/cev.py")
CEV_METALV_SCRIPT = "scripts/archs4/crispercas/cluster_metalv.py"
assert all(v["linkage"] == "complete" for v in CEV_REPRO.values()), "cluster_metalv.py only does complete linkage"


def cev_grid_name(prefix, p):
    """{prefix}_{linkage}_{dist initial}{min_<dist>}_f{min_frac}_r{noise_r2}_{combine}, e.g. hier_complete_c050_f04_r015_median."""
    name = (f"{prefix}_{p['linkage']}_{p['dist'][0]}{round(p['min_' + p['dist']] * 100):03d}"
            f"_f{round(p['min_frac'] * 10):02d}_r{round(p['noise_r2'] * 100):03d}_{p['combine']}")
    assert round(p["min_frac"] * 10) == p["min_frac"] * 10, f"min_frac {p['min_frac']} does not fit the set name"
    return name


def cev_min_models(min_frac, n_runs):
    """As cev.resolve_min_models: ceil(min_frac x n_runs), rounded first so 0.3 x 10 gives 3."""
    return math.ceil(round(min_frac * n_runs, 9))


# Every hierarchical set: the frozen reproduce sets (absolute min_models) plus the grids (min_frac).
# dist picks the linkage: cos (1 - cosine of raw loadings) or jac (1 - top-gene Jaccard), cut at min_<dist>.
CEV_SETS = {s: {"dist": "cos", "linkage": v["linkage"], "min_cos": v["min_cos"], "min_models": v["min_models"],
                "noise_r2": v["noise_r2"], "combine": "median"} for s, v in CEV_REPRO.items()}
CEV_GRID_SETS = {}   # grid sets only: these get split-half / LORO / null subsets (PLAN 2.4)
for _prefix, _grids in CEV_CFG["grids"].items():
    for _grid in (_grids if isinstance(_grids, list) else [_grids]):
        for _vals in itertools.product(*_grid.values()):
            _p = dict(zip(_grid, _vals))
            _name = cev_grid_name(_prefix, _p)
            assert _name not in CEV_SETS, f"duplicate set name {_name}"
            CEV_SETS[_name] = CEV_GRID_SETS[_name] = _p
# Set-definition columns of metrics.csv (meta.json keys), in first-seen order: they lead headline.csv.
CEV_SET_PARAMS = ["set", "subset", *dict.fromkeys(k for v in CEV_SETS.values() for k in v), "min_models_resolved"]
CEV_A4_R2 = sorted({str(v["noise_r2"]) for v in CEV_SETS.values()} | {str(CEV_EVAL["pool_r2"])})
# Reference (cluster_metalv.py) output each set is checked against: the reproduce sets and their grid twins
# (same cosine cut, linkage, r2 and min_frac x 10 runs = min_models).
CEV_CHECK_REF = {s: s for s in CEV_REPRO}
for _r, _v in CEV_REPRO.items():
    for _name, _p in CEV_GRID_SETS.items():
        if ((_p["dist"], _p["linkage"], _p.get("min_cos"), _p["noise_r2"], _p["combine"])
                == ("cos", _v["linkage"], _v["min_cos"], _v["noise_r2"], "median")
                and cev_min_models(_p["min_frac"], len(A4_METALV_MODELS)) == _v["min_models"]):
            CEV_CHECK_REF[_name] = _r


def cev_run_subsets(models, n_splits, seed, loro):
    """{subset: runs}. Halves: n_splits seeded draws of half the runs that contain models[0] (so half a
    and half b are never swapped copies of another split), plus the complement. Runs keep model order."""
    combos = [c for c in itertools.combinations(models, len(models) // 2) if models[0] in c]
    subsets = {}
    for i, a in enumerate(random.Random(seed).sample(combos, n_splits)):
        subsets[f"half{i}a"] = list(a)
        subsets[f"half{i}b"] = [m for m in models if m not in a]
    if loro:
        subsets.update({f"loro_{m}": [x for x in models if x != m] for m in models})
    return subsets


CEV_SUBSET_RUNS = cev_run_subsets(A4_METALV_MODELS, CEV_EVAL["split_half"]["n_splits"],
                                  CEV_EVAL["split_half"]["seed"], CEV_EVAL["loro"])
CEV_HALVES = [(f"half{i}a", f"half{i}b") for i in range(CEV_EVAL["split_half"]["n_splits"])]
CEV_LORO = [s for s in CEV_SUBSET_RUNS if s.startswith("loro_")]
_gp = CEV_EVAL["nulls"]["gene_perm"]
CEV_GENEPERM = {f"geneperm{i}": _gp["seed"] + i for i in range(_gp["n_perm"])}   # all runs, genes shuffled
CEV_SUBSETS = ["full", *CEV_SUBSET_RUNS, *CEV_GENEPERM]


def cev_pool_args(subset):
    """run.py pool selection for a subset: none for full, --runs for split-half / LORO, --gene-perm-seed for nulls."""
    if subset in CEV_GENEPERM:
        return f"--gene-perm-seed {CEV_GENEPERM[subset]}"
    if subset in CEV_SUBSET_RUNS:
        return "--runs " + " ".join(CEV_SUBSET_RUNS[subset])
    return ""


def cev_dendrogram(p):
    """r{r2}/{dist}_{linkage}: the linkage (and subset-set folder) a set is cut from."""
    return f"r{p['noise_r2']}/{p['dist']}_{p['linkage']}"


def cev_dendrogram_sets(r2, dist, linkage):
    """Grid sets cut from one linkage, as {name: definition}."""
    key = f"r{r2}/{dist}_{linkage}"
    return {n: p for n, p in CEV_GRID_SETS.items() if cev_dendrogram(p) == key}


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
    params:
        pool=lambda wc: cev_pool_args(wc.subset),
    wildcard_constraints:
        r2="|".join(re.escape(r) for r in CEV_A4_R2),
        dist="|".join(sorted({v["dist"] for v in CEV_SETS.values()})),
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
        "--top-n {A4_CRISPERCAS_CFG[top_n]} {params.pool} --out {output} > {log} 2>&1"


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
        subset="full",   # subset sets: cev_subset_sets
    resources:
        mem_mb=CEV_RES["set"]["mem_mb"],
        runtime=CEV_RES["set"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} hier-set --set {wildcards.set} --subset {wildcards.subset} "
        "--z-cache {input.z} --flags {input.flags} --linkage {input.linkage} "
        "--params '{params.json}' --out-dir {params.out} > {log} 2>&1"


rule cev_set_check:
    """cev.py builder on real data == unchanged cluster_metalv.py output (the reproduce sets and their grid twins)."""
    input:
        set=lambda wc: expand(rules.cev_hier_set.output, set=wc.set, subset="full"),
        reproduce=lambda wc: expand(rules.cev_reproduce_set.output, set=CEV_CHECK_REF[wc.set]),
        script=CEV_SCRIPT,
    output:
        f"{CEV_A4}/checks/set_{{set}}.ok",
    params:
        set_dir=lambda wc, input: os.path.dirname(input.set[0]),
        reproduce_dir=lambda wc, input: os.path.dirname(input.reproduce[0]),
    wildcard_constraints:
        set="|".join(CEV_CHECK_REF),
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


# ============================================================
# Step 4: robustness (PLAN 2.4): split-half, LORO and nulls for the grid sets
# ============================================================

rule cev_subsets:
    """Provenance: the runs of every split-half / LORO subset and the gene-permutation seeds."""
    output:
        f"{CEV_A4}/report/subsets.json",
    run:
        with open(output[0], "w") as f:
            json.dump({"runs": CEV_SUBSET_RUNS, "gene_perm_seeds": CEV_GENEPERM}, f, indent=2)


rule cev_subset_sets:
    """Every grid set cut from one subset linkage, from one pool load (same builder as cev_hier_set)."""
    input:
        z=rules.cev_zcache.output,
        flags=f"{CEV_A4}/noise/flags_r{{r2}}.csv",
        linkage=f"{CEV_A4}/linkage/r{{r2}}/{{dist}}_{{linkage}}/{{subset}}.npz",
        script=CEV_SCRIPT,
        lib=CEV_LIB,
    output:
        directory(f"{CEV_A4}/subsets/r{{r2}}/{{dist}}_{{linkage}}/{{subset}}"),
    log:
        f"{CEV_A4}/subsets/r{{r2}}/{{dist}}_{{linkage}}/{{subset}}.log",
    params:
        json=lambda wc: json.dumps(cev_dendrogram_sets(wc.r2, wc.dist, wc.linkage)),
        pool=lambda wc: cev_pool_args(wc.subset),
    wildcard_constraints:
        r2="|".join(re.escape(r) for r in CEV_A4_R2),
        subset="|".join(s for s in CEV_SUBSETS if s != "full"),
    resources:
        mem_mb=CEV_RES["set"]["mem_mb"],
        runtime=CEV_RES["set"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} hier-sets --subset {wildcards.subset} --z-cache {input.z} --flags {input.flags} "
        "{params.pool} --linkage {input.linkage} --params '{params.json}' --out-dir {output} > {log} 2>&1"


rule cev_null_label:
    """Clusters spanning >= m runs at every grid cut of a full linkage, observed vs run labels shuffled."""
    input:
        linkage=f"{CEV_A4}/linkage/r{{r2}}/{{dist}}_{{linkage}}/full.npz",
        script=CEV_SCRIPT,
        lib=CEV_LIB,
    output:
        f"{CEV_A4}/null/label_perm/r{{r2}}/{{dist}}_{{linkage}}.csv",
    log:
        f"{CEV_A4}/null/label_perm/r{{r2}}/{{dist}}_{{linkage}}.log",
    params:
        cuts=lambda wc: " ".join(str(c) for c in sorted({p[f"min_{wc.dist}"] for p in
                                                         cev_dendrogram_sets(wc.r2, wc.dist, wc.linkage).values()})),
        null=CEV_EVAL["nulls"]["label_perm"],
    wildcard_constraints:
        r2="|".join(re.escape(r) for r in CEV_A4_R2),
    resources:
        mem_mb=CEV_RES["check"]["mem_mb"],
        runtime=CEV_RES["check"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} label-null --linkage {input.linkage} --cuts {params.cuts} "
        "--n-perm {params.null[n_perm]} --seed {params.null[seed]} --out {output} > {log} 2>&1"


def cev_subset_set_dirs(wc, subsets):
    return [f"{CEV_A4}/subsets/{cev_dendrogram(CEV_GRID_SETS[wc.set])}/{s}/{wc.set}" for s in subsets]


rule cev_compare:
    """Split-half (Hungarian match of half a vs half b), LORO recovery and null counts of one grid set."""
    input:
        full=lambda wc: expand(rules.cev_hier_set.output, set=wc.set, subset="full"),
        subsets=lambda wc: [f"{CEV_A4}/subsets/{cev_dendrogram(CEV_GRID_SETS[wc.set])}/{s}"
                            for s in [*CEV_SUBSET_RUNS, *CEV_GENEPERM]],
        label_null=lambda wc: f"{CEV_A4}/null/label_perm/{cev_dendrogram(CEV_GRID_SETS[wc.set])}.csv",
        script=CEV_SCRIPT,
        lib=CEV_LIB,
    output:
        split_half=f"{CEV_A4}/sets/{{set}}/eval/split_half.csv",
        loro=f"{CEV_A4}/sets/{{set}}/eval/loro.csv",
        robust=f"{CEV_A4}/sets/{{set}}/eval/robust.csv",
    log:
        f"{CEV_A4}/sets/{{set}}/eval/compare.log",
    params:
        full=lambda wc, input: os.path.dirname(input.full[0]),
        half_a=lambda wc: " ".join(cev_subset_set_dirs(wc, [a for a, _ in CEV_HALVES])),
        half_b=lambda wc: " ".join(cev_subset_set_dirs(wc, [b for _, b in CEV_HALVES])),
        loro=lambda wc: " ".join(cev_subset_set_dirs(wc, CEV_LORO)),
        geneperm=lambda wc: " ".join(cev_subset_set_dirs(wc, CEV_GENEPERM)),
        out=lambda wc, output: os.path.dirname(output.robust),
    wildcard_constraints:
        set="|".join(CEV_GRID_SETS),
    resources:
        mem_mb=CEV_RES["check"]["mem_mb"],
        runtime=CEV_RES["check"]["runtime"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} compare --full-dir {params.full} --half-a {params.half_a} --half-b {params.half_b} "
        "--loro {params.loro} --geneperm {params.geneperm} --label-null {input.label_null} "
        "--top-n {A4_CRISPERCAS_CFG[top_n]} --match-cos {CEV_EVAL[match_cos]} --match-jac {CEV_EVAL[match_jac]} --single-gene {CEV_EVAL[single_gene]} "
        "--out-dir {params.out} > {log} 2>&1"


rule cev_aggregate:
    input:
        metrics=expand(rules.cev_evaluate.output.metrics, set=CEV_SETS),
        robust=expand(rules.cev_compare.output.robust, set=CEV_GRID_SETS),
        script=CEV_SCRIPT,
    output:
        all=f"{CEV_A4}/report/metrics_all.csv",
        headline=f"{CEV_A4}/report/headline.csv",
        robust=f"{CEV_A4}/report/robustness.csv",
    params:
        set_params=" ".join(CEV_SET_PARAMS),
        headline=" ".join(CEV_EVAL["headline"]),
    conda: "clamp-analyses"
    shell:
        "python {input.script} aggregate --metrics {input.metrics} --params {params.set_params} "
        "--headline {params.headline} --robust {input.robust} "
        "--out {output.all} --headline-out {output.headline} --robust-out {output.robust}"


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
        expand(rules.cev_set_check.output, set=CEV_CHECK_REF),
        rules.cev_subsets.output,
        rules.cev_aggregate.output,

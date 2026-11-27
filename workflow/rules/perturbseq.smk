# Perturb-seq regulator->target benchmark of the OneK1K LVs (PLAN.md).
import json
import shlex

PS_CFG = config["perturbseq"]
with open(PS_CFG["inputs_manifest"]) as handle:
    PS_MEMBERS = {key: member["path"] for key, member in json.load(handle)["members"].items()}
PS_SOURCE = config["perturbseq_source_root"]


rule ota2025_k562gw:
    output:
        logfc=PS_MEMBERS["k562gw_logfc"],
        p=PS_MEMBERS["k562gw_p"],
        ids=PS_MEMBERS["gencode_v41_ids"],
    log:
        os.path.join(os.path.dirname(PS_MEMBERS["gencode_v41_ids"]), "ota2025_k562gw.log")
    params:
        manifest=PS_CFG["inputs_manifest"],
        source=("--source-root " + shlex.quote(PS_SOURCE)) if PS_SOURCE else "",
    resources:
        runtime=240,  # 7.7 GB download plus ~3 GB extraction; seconds when linking.
    conda: "snakemake"
    shell:
        "python scripts/inputs/perturbseq_inputs.py --manifest {params.manifest:q} {params.source} > {log:q} 2>&1"


PS_STRATEGIES = PS_CFG["strategies"]
# Each strategy reads the Z.csv files written by its fitting area (pseudobulk.smk, donor_bulk.smk).
PS_Z_ROOTS = {
    "main": config["paths"]["production"],
    "donor_bulk": config["donor_bulk"]["paths"]["production"],
}
PS_REGULONS = os.path.join(PS_CFG["output_root"], "regulons", "{strategy}")


def ps_z_files(wildcards):
    root = os.path.join(PS_Z_ROOTS[wildcards.strategy], PS_CFG["dataset"], "models")
    return [os.path.join(root, model, "Z.csv") for model in PS_STRATEGIES[wildcards.strategy]["models"]]


rule perturbseq_regulons:
    input:
        logfc=rules.ota2025_k562gw.output.logfc,
        p=rules.ota2025_k562gw.output.p,
        ids=rules.ota2025_k562gw.output.ids,
        z=ps_z_files,
    output:
        gmt=os.path.join(PS_REGULONS, "regulons.gmt"),
        clusters=os.path.join(PS_REGULONS, "clusters.csv"),
        sizes=os.path.join(PS_REGULONS, "sizes.csv"),
        logfc=os.path.join(PS_REGULONS, "logfc.tsv.gz"),
        report=os.path.join(PS_REGULONS, "report.json"),
    log:
        os.path.join(PS_REGULONS, "build_regulons.log"),
    wildcard_constraints:
        strategy="|".join(PS_STRATEGIES),
    params:
        out_dir=lambda wc, output: os.path.dirname(output.gmt),
        cfg=PS_CFG["regulons"],
        min_de=PS_CFG["response"]["min_de"],
    resources:
        mem_mb=16000,
        runtime=30,
    conda: "clamp-analyses"
    shell:
        "python scripts/perturbseq/build_regulons.py --logfc {input.logfc:q} --pvalues {input.p:q} "
        "--id-map {input.ids:q} --z {input.z:q} --q {params.cfg[q]} --min-size {params.cfg[min_size]} "
        "--max-size {params.cfg[max_size]} --jaccard {params.cfg[jaccard]} --min-de {params.min_de} "
        "--out-dir {params.out_dir:q} > {log:q} 2>&1"


# Step 5: regulon ORA. run_ora.R gives the observed hits per tail (Z and -Z, restricted to the
# regulon universe); ora_null.py reproduces them and adds the gene-label permutation null.
PS_ORA = os.path.join(PS_CFG["output_root"], "ora", "{strategy}", "{model}")
PS_SIGNED = os.path.join(PS_CFG["output_root"], "signed_z", "{strategy}", "{model}")
PS_TAILS = ["pos", "neg"]
# projections.smk constrains {model} globally; these rules override it with the Perturb-seq models.
PS_MODEL_PATTERN = "|".join(sorted({m for cfg in PS_STRATEGIES.values() for m in cfg["models"]}))


def ps_model_dir(wildcards):
    return os.path.join(PS_Z_ROOTS[wildcards.strategy], PS_CFG["dataset"], "models", wildcards.model)


def ps_strategy_models():
    return [(s, m) for s, cfg in PS_STRATEGIES.items() for m in cfg["models"]]


rule perturbseq_signed_z:
    input:
        z=lambda wc: os.path.join(ps_model_dir(wc), "Z.csv"),
        b=lambda wc: os.path.join(ps_model_dir(wc), "B.csv"),
        universe=rules.perturbseq_regulons.output.logfc,
    output:
        pos=os.path.join(PS_SIGNED, "pos.csv"),
        neg=os.path.join(PS_SIGNED, "neg.csv"),
        manifest=os.path.join(PS_SIGNED, "manifest.json"),
    wildcard_constraints:
        strategy="|".join(PS_STRATEGIES),
        model=PS_MODEL_PATTERN,
    resources:
        mem_mb=4000,
        runtime=10,
    conda: "clamp-analyses"
    shell:
        "python scripts/perturbseq/signed_z.py --z {input.z:q} --b {input.b:q} --universe {input.universe:q} "
        "--out-pos {output.pos:q} --out-neg {output.neg:q} --out-manifest {output.manifest:q}"


rule perturbseq_ora:
    input:
        z=os.path.join(PS_SIGNED, "{tail}.csv"),
        manifest=rules.perturbseq_signed_z.output.manifest,
        gmt=rules.perturbseq_regulons.output.gmt,
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(os.path.join(PS_ORA, "{tail}")),
    log:
        os.path.join(PS_ORA, "{tail}.log"),
    params:
        ora=PS_CFG["ora"],
        reg=PS_CFG["regulons"],
    wildcard_constraints:
        strategy="|".join(PS_STRATEGIES),
        model=PS_MODEL_PATTERN,
        tail="|".join(PS_TAILS),
    resources:
        mem_mb=int(PS_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(PS_CFG["ora"]["resources"]["runtime"]),
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z:q} --out-dir {output.ora_dir:q} "
        "--dataset {PS_CFG[dataset]} --fraction 100 --seed {PS_CFG[seed]} "
        "--model {wildcards.model}_{wildcards.tail} --model-manifest {input.manifest:q} "
        "--database {params.ora[database]} --database-label '{params.ora[database_label]}' "
        "--database-type gmt --database-path {input.gmt:q} "
        "--top-pct {params.ora[top_pct]} --min-size {params.reg[min_size]} --max-size {params.reg[max_size]} "
        "--pvalue-cutoff {params.ora[pvalue_cutoff]} --qvalue-cutoff {params.ora[qvalue_cutoff]} "
        "--fdr {params.ora[fdr]} > {log:q} 2>&1"


rule perturbseq_ora_null:
    input:
        pos=os.path.join(PS_SIGNED, "pos.csv"),
        neg=os.path.join(PS_SIGNED, "neg.csv"),
        ora_pos=os.path.join(PS_ORA, "pos"),
        ora_neg=os.path.join(PS_ORA, "neg"),
        gmt=rules.perturbseq_regulons.output.gmt,
    output:
        null=os.path.join(PS_ORA, "null.csv.gz"),
        lv_hits=os.path.join(PS_ORA, "lv_hits.csv"),
        summary=os.path.join(PS_ORA, "summary.csv"),
    log:
        os.path.join(PS_ORA, "ora_null.log"),
    params:
        out_dir=lambda wc, output: os.path.dirname(output.summary),
        ora=PS_CFG["ora"],
        reg=PS_CFG["regulons"],
    wildcard_constraints:
        strategy="|".join(PS_STRATEGIES),
        model=PS_MODEL_PATTERN,
    resources:
        mem_mb=int(PS_CFG["ora"]["null_resources"]["mem_mb"]),
        runtime=int(PS_CFG["ora"]["null_resources"]["runtime"]),
    conda: "clamp-analyses"
    shell:
        "python scripts/perturbseq/ora_null.py --z-pos {input.pos:q} --z-neg {input.neg:q} --gmt {input.gmt:q} "
        "--ora-pos {input.ora_pos:q} --ora-neg {input.ora_neg:q} "
        "--strategy {wildcards.strategy} --model {wildcards.model} "
        "--top-pct {params.ora[top_pct]} --min-size {params.reg[min_size]} --max-size {params.reg[max_size]} "
        "--fdr {params.ora[fdr]} --n-perm {PS_CFG[n_perm]} --seed {PS_CFG[seed]} "
        "--out-dir {params.out_dir:q} > {log:q} 2>&1"


rule perturbseq_ora_summary:
    input:
        [os.path.join(PS_CFG["output_root"], "ora", s, m, "summary.csv") for s, m in ps_strategy_models()],
    output:
        os.path.join(PS_CFG["output_root"], "ora_summary.csv"),
    conda: "clamp-analyses"
    shell:
        "python -c \"import sys, pandas as pd; "
        "pd.concat([pd.read_csv(f) for f in sys.argv[2:]]).to_csv(sys.argv[1], index=False)\" "
        "{output:q} {input:q}"

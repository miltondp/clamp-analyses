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

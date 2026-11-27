# Perturb-seq regulator->target benchmark of the OneK1K LVs (PLAN.md).
import json
import shlex

PS_CFG = config["perturbseq"]
with open(PS_CFG["inputs_manifest"]) as handle:
    PS_MEMBERS = {key: member["path"] for key, member in json.load(handle)["members"].items()}


def perturbseq_source_flag(wildcards):
    source = config.get("perturbseq_source_root") or ""
    return f"--source-root {shlex.quote(source)}" if source else ""


rule ota2025_k562gw:
    output:
        logfc=PS_MEMBERS["k562gw_logfc"],
        p=PS_MEMBERS["k562gw_p"],
        ids=PS_MEMBERS["gencode_v41_ids"],
    log:
        f"{PS_CFG['data_root']}/ota2025_k562gw.log"
    params:
        manifest=PS_CFG["inputs_manifest"],
        source=perturbseq_source_flag,
    resources:
        runtime=240,  # 7.7 GB download plus ~3 GB extraction; seconds when linking.
    conda: "snakemake"
    shell:
        "python scripts/inputs/perturbseq_inputs.py --manifest {params.manifest:q} {params.source} > {log:q} 2>&1"

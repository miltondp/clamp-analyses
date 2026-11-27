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

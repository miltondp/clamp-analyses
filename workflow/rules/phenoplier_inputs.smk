import shlex

# Public PhenoPLIER inputs, pinned to upstream content. No GLS fits run here.
with open("workflow/config/phenoplier_inputs.json") as handle:
    PHENO_INPUT_ASSETS = json.load(handle)["assets"]
PHENO_INPUT_PATHS = {asset["path"]: key for key, asset in PHENO_INPUT_ASSETS.items()}
PHENO_INPUT_SOURCE = config.get("phenoplier_input_root", os.environ.get("PHENOPLIER_DATA_ROOT", ""))

rule fetch_phenoplier_input:
    input:
        lambda wc: [PHENO_INPUT_ASSETS[PHENO_INPUT_PATHS[wc.phenoplier_file]]["local_source"]]
        if PHENO_INPUT_ASSETS[PHENO_INPUT_PATHS[wc.phenoplier_file]].get("local_source") else [],
    output:
        "{phenoplier_file}",
    wildcard_constraints:
        phenoplier_file="|".join(re.escape(path) for path in PHENO_INPUT_PATHS),
    params:
        key=lambda wc: PHENO_INPUT_PATHS[wc.phenoplier_file],
        source_arg=("--source-root " + shlex.quote(PHENO_INPUT_SOURCE)) if PHENO_INPUT_SOURCE else "",
    conda: "snakemake"
    shell:
        "python scripts/inputs/reference_inputs.py --manifest workflow/config/phenoplier_inputs.json "
        "--only {params.key:q} {params.source_arg}"

rule link_spredixcan_inputs:
    output:
        directory(config["archs4"]["drug_diseases"]["spredixcan_dir"]),
    params:
        source_arg=("--source-root " + shlex.quote(PHENO_INPUT_SOURCE)) if PHENO_INPUT_SOURCE else "",
    resources:
        mem_mb=8000,
    conda: "clamp-analyses"
    shell:
        "python scripts/inputs/link_spredixcan_inputs.py --destination {output:q} {params.source_arg}"

rule drug_disease_inputs:
    input:
        config["archs4"]["drug_diseases"]["drugbank_tsv"],
        rules.link_spredixcan_inputs.output,
        [asset["path"] for asset in PHENO_INPUT_ASSETS.values() if asset["group"] == "drug_disease"],

rule crispr_inputs:
    input:
        [asset["path"] for asset in PHENO_INPUT_ASSETS.values() if asset["group"] in {"crispr", "phenomexcan"}],

rule phenoplier_inputs:
    input:
        rules.link_spredixcan_inputs.output,
        [asset["path"] for asset in PHENO_INPUT_ASSETS.values()],

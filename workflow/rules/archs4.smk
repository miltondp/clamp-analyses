import os

# Shared ARCHS4 workflow context. Final models are published under
# output/98_final_models and consumed by their dedicated analyses.
A4_CFG = config["archs4"]
A4_PROD = A4_CFG["paths"]["production"]
A4_BIO = A4_CFG["paths"]["biology"]
A4_BIO_NB = os.path.join(REPO_ROOT, A4_CFG["paths"]["biology_notebooks"])
A4_PREP = f"{A4_PROD}/00_preprocess"

rule download_archs4:
    output:
        A4_CFG["raw_h5"]
    params:
        url=A4_CFG["raw_url"],
    conda: "clamp-analyses"
    shell:
        "mkdir -p $(dirname {output}) && "
        "curl --fail --location --retry 3 '{params.url}' --output {output}.download && "
        "mv {output}.download {output}"


rule preprocess_archs4:
    input:
        h5=ancient(A4_CFG["raw_h5"]),
        gene_lengths=ancient(A4_CFG["gene_lengths"]),
        script=ancient("scripts/archs4/preprocess.R"),
        common=ancient("scripts/archs4/common.R"),
    output:
        metadata=f"{A4_PREP}/metadata.rds",
        metadata_filtered=f"{A4_PREP}/metadata_filtered.rds",
        all_samples=f"{A4_PREP}/all_samples.rds",
        fbm_filtered=f"{A4_PREP}/fbm_filtered.bk",
    log:
        f"{A4_PREP}/preprocess.log"
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.metadata),
        mean_cutoff=A4_CFG["preprocess"]["mean_cutoff"],
        var_cutoff=A4_CFG["preprocess"]["var_cutoff"],
        sc_max=A4_CFG["preprocess"]["single_cell_probability_max"],
        block_size=A4_CFG["raw_data"]["block_size"],
        n_cores=A4_CFG["raw_data"]["n_cores"],
        seed=A4_CFG["preprocess"]["seed"],
    threads: A4_CFG["raw_data"]["n_cores"]
    resources:
        mem_mb=200000,
        runtime=18000,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --h5 {input.h5} --gene-lengths {input.gene_lengths} "
        "--out-dir {params.out_dir} --mean-cutoff {params.mean_cutoff} "
        "--var-cutoff {params.var_cutoff} "
        "--single-cell-probability-max {params.sc_max} "
        "--block-size {params.block_size} --n-cores {params.n_cores} --seed {params.seed} "
        "> {log} 2>&1"


rule generate_archs4_go_bp_gmt:
    input:
        script="scripts/archs4/generate_go_bp_gmt.R",
    output:
        A4_CFG["coverage"]["prior_gmt"],
    params:
        sha256="33e559df968d0d8c2028c0f735624555dc2edbebe955648c3351dfa31cfd2e2a",
    resources:
        mem_mb=8000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "mkdir -p $(dirname {output:q}) && "
        "tmp=$(mktemp {output:q}.XXXXXX) && "
        "trap 'if test -e \"$tmp\"; then unlink \"$tmp\"; fi' EXIT && "
        "Rscript {input.script:q} \"$tmp\" && "
        "printf '%s  %s\\n' '{params.sha256}' \"$tmp\" | sha256sum --check --status && "
        "mv \"$tmp\" {output:q}"

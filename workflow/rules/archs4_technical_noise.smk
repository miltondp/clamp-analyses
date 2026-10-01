import re

# Technical-noise screen for the published canonical ARCHS4 model.
A4_TN_CFG = A4_CFG["technical_noise"]
A4_TN_ROOT = A4_TN_CFG["compute_root"]
A4_TN_NB = os.path.join(REPO_ROOT, A4_TN_CFG["notebooks"])
A4_TN_ENS = A4_TN_CFG["ensembl_dir"]
A4_TN_URLS = A4_TN_CFG["ensembl_urls"]
A4_TN_MODEL = f"{A4_CAN_FINAL_ROOT}/archs4"
A4_TN_TRAITS = f"{A4_TN_MODEL}/traits/canon_archs4.tsv.gz"
A4_TN_SCRIPT = "scripts/archs4/technical_noise/run_analysis.R"
A4_TN_COMMON = "scripts/archs4/common.R"


rule technical_noise_fetch_ensembl:
    output:
        f"{A4_TN_ENS}/{{name}}",
    log:
        f"{A4_TN_ENS}/{{name}}.log"
    params:
        url=lambda wc: A4_TN_URLS[wc.name],
    wildcard_constraints:
        name="|".join(re.escape(name) for name in A4_TN_URLS),
    shell:
        "curl --fail --location --retry 3 '{params.url}' --output {output} 2> {log}"


rule technical_noise_gene_properties:
    input:
        z=f"{A4_TN_MODEL}/Z.csv",
        h5=A4_CFG["raw_h5"],
        ensembl=expand(f"{A4_TN_ENS}/{{name}}", name=A4_TN_URLS),
        script=A4_TN_SCRIPT,
        common=A4_TN_COMMON,
    output:
        f"{A4_TN_ROOT}/gene_properties.csv",
    log:
        f"{A4_TN_ROOT}/gene_properties.log"
    params:
        model=A4_TN_MODEL,
        out=A4_TN_ROOT,
        ensembl=A4_TN_ENS,
    resources:
        mem_mb=24000,
        runtime=120,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --stage gene_properties --model-dir {params.model} "
        "--out-dir {params.out} --h5 {input.h5} --ensembl-dir {params.ensembl} "
        "> {log} 2>&1"


rule technical_noise_sample_qc:
    input:
        b=f"{A4_TN_MODEL}/B.csv",
        h5=A4_CFG["raw_h5"],
        metadata=f"{A4_PREP}/metadata_filtered.rds",
        script=A4_TN_SCRIPT,
        common=A4_TN_COMMON,
    output:
        metadata=f"{A4_TN_ROOT}/sample_metadata.tsv.gz",
        unreadable=f"{A4_TN_ROOT}/unreadable_samples.csv",
    log:
        f"{A4_TN_ROOT}/sample_qc.log"
    params:
        model=A4_TN_MODEL,
        out=A4_TN_ROOT,
        prep=A4_PREP,
    threads: 8
    resources:
        mem_mb=96000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --stage sample_qc --model-dir {params.model} "
        "--out-dir {params.out} --h5 {input.h5} --preprocess-dir {params.prep} "
        "--n-cores {threads} > {log} 2>&1"


rule technical_noise_z:
    input:
        z=f"{A4_TN_MODEL}/Z.csv",
        annotation=f"{A4_TN_MODEL}/summary.csv",
        properties=rules.technical_noise_gene_properties.output,
        script=A4_TN_SCRIPT,
        common=A4_TN_COMMON,
    output:
        f"{A4_TN_ROOT}/lv_gene_properties_canonical.csv",
    log:
        f"{A4_TN_ROOT}/lv_gene_properties.log"
    params:
        model=A4_TN_MODEL,
        out=A4_TN_ROOT,
    resources:
        mem_mb=16000,
        runtime=180,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --stage gene_level --model-dir {params.model} "
        "--out-dir {params.out} > {log} 2>&1"


rule technical_noise_b:
    input:
        b=f"{A4_TN_MODEL}/B.csv",
        metadata=rules.technical_noise_sample_qc.output.metadata,
        script=A4_TN_SCRIPT,
        common=A4_TN_COMMON,
    output:
        f"{A4_TN_ROOT}/lv_sample_metadata_canonical.csv",
    log:
        f"{A4_TN_ROOT}/lv_sample_metadata.log"
    params:
        model=A4_TN_MODEL,
        out=A4_TN_ROOT,
    threads: 8
    resources:
        mem_mb=96000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --stage sample_level --model-dir {params.model} "
        "--out-dir {params.out} --n-cores {threads} > {log} 2>&1"


rule technical_noise_summary:
    input:
        z=rules.technical_noise_z.output,
        b=rules.technical_noise_b.output,
        annotation=f"{A4_TN_MODEL}/summary.csv",
        traits=A4_TN_TRAITS,
        script=A4_TN_SCRIPT,
        common=A4_TN_COMMON,
    output:
        classification=f"{A4_TN_ROOT}/lv_technical_classification.csv",
        affected=f"{A4_TN_ROOT}/affected_lv_summary.csv",
        groups=f"{A4_TN_ROOT}/affected_lv_group_summary.csv",
    log:
        f"{A4_TN_ROOT}/summary.log"
    params:
        model=A4_TN_MODEL,
        out=A4_TN_ROOT,
    resources:
        mem_mb=16000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --stage summary --model-dir {params.model} "
        "--out-dir {params.out} --traits {input.traits} > {log} 2>&1"


rule technical_noise_report_archs4:
    input:
        groups=rules.technical_noise_summary.output.groups,
        affected=rules.technical_noise_summary.output.affected,
        notebook=f"{A4_TN_NB}/00_technical_noise.ipynb",
    output:
        complete=touch(f"{A4_TN_ROOT}/notebook.complete"),
    log:
        notebook=f"{A4_TN_NB}/00_technical_noise.executed.ipynb",
    resources:
        mem_mb=4000,
        runtime=30,
    conda: "clamp-analyses"
    notebook:
        f"{A4_TN_NB}/00_technical_noise.ipynb"

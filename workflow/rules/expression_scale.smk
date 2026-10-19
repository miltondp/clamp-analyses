ES_ROOT = "output/01_model_building/04_expression_scale"
ES_NB = "nbs/01_model_building/04_expression_scale"


# ============================================================
# Step 1: sample linear-scale expression matrices
# ============================================================

rule expression_scale_sample:
    input:
        gtex_bk=ancient(f"{GTEX_PROD}/FBMgtex.bk"),
        gtex_samples=ancient(f"{GTEX_PROD}/gtex_samples.rds"),
        gtex_genes=ancient(f"{GTEX_PROD}/gtex_genes.rds"),
        recount2_bk=ancient(f"{R2_DIR}/FBMrecount2.bk"),
        recount2_samples=ancient(f"{R2_DIR}/recount2_samples.rds"),
        recount2_genes=ancient(f"{R2_DIR}/recount2_genes.rds"),
        archs4_h5=ancient(A4_CFG["raw_h5"]),
        archs4_gene_lengths=ancient(A4_CFG["gene_lengths"]),
        archs4_metadata=ancient(f"{A4_PREP}/metadata_filtered.rds"),
        pseudobulk=[ancient(f) for f in expand(f"{PROD}/{{dataset}}/preprocessing/cpm_filt.csv", dataset=DATASETS)],
        script="scripts/expression_scale/sample_linear.R",
    output:
        f"{ES_ROOT}/linear_matrices.rds",
    log:
        f"{ES_ROOT}/sample_linear.log",
    params:
        pseudobulk=lambda wc, input: ",".join(input.pseudobulk),
        sc_max=A4_CFG["preprocess"]["single_cell_probability_max"],
    resources:
        mem_mb=16000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --gtex-bk {input.gtex_bk} --gtex-samples {input.gtex_samples} "
        "--gtex-genes {input.gtex_genes} --recount2-bk {input.recount2_bk} "
        "--recount2-samples {input.recount2_samples} --recount2-genes {input.recount2_genes} "
        "--archs4-h5 {input.archs4_h5} --archs4-gene-lengths {input.archs4_gene_lengths} "
        "--archs4-metadata {input.archs4_metadata} --archs4-sc-max {params.sc_max} --pseudobulk {params.pseudobulk} "
        "--n-samples 2000 --seed 123 --out {output} > {log} 2>&1"


# ============================================================
# Step 2: compare expression scales
# ============================================================

rule expression_scale_report:
    input:
        matrices=rules.expression_scale_sample.output,
        notebook=f"{ES_NB}/00_log2_transform.ipynb",
    output:
        summary=f"{ES_ROOT}/summary.csv",
        mean_sd=f"{ES_ROOT}/meanSdPlot.png",
        mean_var=f"{ES_ROOT}/plotMeanVar.png",
        skewness=f"{ES_ROOT}/skewness.png",
        complete=touch(f"{ES_ROOT}/notebook.complete"),
    log:
        notebook=f"{ES_NB}/00_log2_transform.executed.ipynb",
    resources:
        mem_mb=16000,
        runtime=30,
    conda: "clamp-analyses"
    notebook:
        f"{ES_NB}/00_log2_transform.ipynb"

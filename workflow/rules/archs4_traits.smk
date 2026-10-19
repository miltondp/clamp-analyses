A4_TRAIT_CFG = A4_CFG["traits"]


# ============================================================
# Step 1: summarize trait recovery by training size
# ============================================================

rule aggregate_archs4_coverage_traits:
    input:
        clampbase=PHENOPLIER_COV_BASE,
        finals_clampbase=PHENOPLIER_FIN_BASE,
        comparators_clampbase=PHENOPLIER_COMP_BASE,
        script="scripts/archs4/traits/aggregate_coverage_traits.R",
    output:
        traits_long=f"{A4_COV_BIO}/coverage_trait_recovery.csv",
        finals_long=f"{A4_COV_BIO}/final_trait_recovery.csv",
    log:
        f"{A4_COV_BIO}/traits_aggregate.log"
    params:
        clampbase_dir=lambda wc, input: os.path.dirname(input.clampbase[0]),
        finals_clampbase_dir=lambda wc, input: os.path.dirname(input.finals_clampbase[0]),
        comparators_clampbase_dir=A4_TRAIT_CFG["coverage"]["comparators_clampbase_dir"],
        fdr=A4_TRAIT_CFG["fdr"],
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --clampbase-dir {params.clampbase_dir} "
        "--finals-clampbase-dir {params.finals_clampbase_dir} "
        "--comparators-clampbase-dir {params.comparators_clampbase_dir} --fdr {params.fdr} "
        "--traits-out {output.traits_long} --finals-out {output.finals_long} > {log} 2>&1"


rule coverage_traits_report_archs4:
    input:
        traits_long=rules.aggregate_archs4_coverage_traits.output.traits_long,
        finals_long=rules.aggregate_archs4_coverage_traits.output.finals_long,
        notebook=f"{A4_COV_NB}/01_coverage_traits.ipynb",
    output:
        complete=touch(f"{A4_COV_BIO}/traits.complete"),
    log:
        notebook=f"{A4_COV_NB}/01_coverage_traits.executed.ipynb",
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    notebook:
        f"{A4_COV_NB}/01_coverage_traits.ipynb"


# ============================================================
# Step 2: summarize trait recovery at fixed rank
# ============================================================

rule aggregate_archs4_saturation_traits:
    input:
        clampbase=PHENOPLIER_SAT_BASE,
        full_clampbase=[p for p in PHENOPLIER_COV_BASE if "/cov_rs100_" in p],
        script="scripts/archs4/traits/aggregate_saturation_traits.R",
    output:
        traits_long=f"{A4_SAT_BIO}/saturation_trait_recovery_k1728.csv",
    log:
        f"{A4_SAT_BIO}/traits_aggregate.log"
    params:
        clampbase_dir=lambda wc, input: os.path.dirname(input.clampbase[0]),
        full_clampbase_dir=A4_TRAIT_CFG["coverage"]["clampbase_dir"],
        fdr=A4_TRAIT_CFG["fdr"],
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --clampbase-dir {params.clampbase_dir} "
        "--full-clampbase-dir {params.full_clampbase_dir} --fdr {params.fdr} "
        "--traits-out {output.traits_long} > {log} 2>&1"


rule saturation_traits_report_archs4:
    input:
        traits_long=rules.aggregate_archs4_saturation_traits.output.traits_long,
        notebook=f"{A4_SAT_NB}/01_saturation_traits.ipynb",
    output:
        complete=touch(f"{A4_SAT_BIO}/traits.complete"),
    log:
        notebook=f"{A4_SAT_NB}/01_saturation_traits.executed.ipynb",
    resources:
        mem_mb=8000,
        runtime=30,
    conda: "clamp-analyses"
    notebook:
        f"{A4_SAT_NB}/01_saturation_traits.ipynb"


rule archs4_traits:
    input:
        rules.coverage_traits_report_archs4.output.complete,
        rules.saturation_traits_report_archs4.output.complete,

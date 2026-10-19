A4_STAB_CFG = A4_CFG["stability"]
A4_STAB_ROOT = A4_STAB_CFG["root"]
A4_STAB_SEEDS = [int(seed) for seed in A4_STAB_CFG["seeds"]]

A4_STAB_FROM_COVERAGE = {rng: a4_cov_cell(100, rep) for rep, rng in A4_COV_FULL_REPLICATES.items()}
A4_STAB_OWN_SEEDS = [seed for seed in A4_STAB_SEEDS if seed not in A4_STAB_FROM_COVERAGE]

A4_STAB_FIT_SEED_PATTERN = "|".join(
    str(seed) for seed in A4_STAB_SEEDS
    if not a4_model_complete(f"{A4_STAB_ROOT}/seed{seed}/CLAMPfull_canonical", A4_CAN_FILES)
) or r"(?!)"


# ============================================================
# Step 1: reuse or fit independent CLAMPbase replicates
# ============================================================

rule archs4_stability_base_from_coverage:
    input:
        svd=lambda wc: f"{A4_STAB_FROM_COVERAGE[int(wc.seed)]}/svd.rds",
        k=lambda wc: f"{A4_STAB_FROM_COVERAGE[int(wc.seed)]}/CLAMP_K.rds",
        rds=lambda wc: f"{A4_STAB_FROM_COVERAGE[int(wc.seed)]}/CLAMPbase.rds",
        z=lambda wc: f"{A4_STAB_FROM_COVERAGE[int(wc.seed)]}/CLAMPbase/Z.csv",
        b=lambda wc: f"{A4_STAB_FROM_COVERAGE[int(wc.seed)]}/CLAMPbase/B.csv",
    output:
        svd=f"{A4_STAB_ROOT}/seed{{seed}}/svd.rds",
        k=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMP_K.rds",
        rds=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPbase.rds",
        z=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPbase/Z.csv",
        b=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPbase/B.csv",
    wildcard_constraints:
        seed="|".join(map(str, A4_STAB_FROM_COVERAGE)),
    shell:
        "mkdir -p $(dirname {output.z:q}) && "
        "ln -sfr {input.svd:q} {output.svd:q} && "
        "ln -sfr {input.k:q} {output.k:q} && "
        "ln -sfr {input.rds:q} {output.rds:q} && "
        "ln -sfr {input.z:q} {output.z:q} && ln -sfr {input.b:q} {output.b:q}"


rule archs4_stability_base:
    input:
        metadata=ancient(f"{A4_PREP}/metadata_filtered.rds"),
        subsample=ancient(f"{a4_cov_cell(100, 1)}/subsample_info.rds"),
        fbm=ancient(f"{A4_PREP}/fbm_filtered.bk"),
        svd_script="scripts/archs4/svd.R",
        base_script=ancient("scripts/archs4/saturation/fit_clampbase.R"),
        common="scripts/archs4/common.R",
        raw_metadata=f"{A4_PREP}/metadata.rds",
    output:
        svd=f"{A4_STAB_ROOT}/seed{{seed}}/svd.rds",
        k=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMP_K.rds",
        rds=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPbase.rds",
        z=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPbase/Z.csv",
        b=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPbase/B.csv",
    log:
        f"{A4_STAB_ROOT}/logs/seed{{seed}}.base.log",
    params:
        fbm_prefix=lambda wildcards, input: input.fbm[: -len(".bk")],
        out_dir=lambda wildcards, output: os.path.dirname(output.rds),
        index=lambda wc: A4_STAB_SEEDS.index(int(wc.seed)) + 2,
        svd_k=A4_STAB_CFG["svd_k"],
    threads: int(A4_STAB_CFG["resources"]["threads"])
    resources:
        mem_mb=int(A4_STAB_CFG["resources"]["mem_mb"]),
        runtime=int(A4_STAB_CFG["resources"]["runtime"]),
    wildcard_constraints:
        seed="|".join(map(str, A4_STAB_OWN_SEEDS)),
    conda: "clamp-analyses"
    shell:
        "(Rscript {input.svd_script} --metadata {input.metadata} --fbm {params.fbm_prefix} "
        "--svd {output.svd} --k {output.k} --svd-k {params.svd_k} --n-cores {threads} "
        "--seed {wildcards.seed} && touch {params.out_dir}/svd.complete && "
        "Rscript {input.base_script} --metadata {input.metadata} --subsample-info {input.subsample} "
        "--fbm-backing {input.fbm} --svd {output.svd} --k {output.k} --out-dir {params.out_dir} "
        "--seed {wildcards.seed} --seed-index {params.index} --fraction 100 "
        "--requested-cpus {threads} --requested-mem-mb {resources.mem_mb}) > {log} 2>&1"


# ============================================================
# Step 2: fit canonical CLAMPfull replicates
# ============================================================

rule archs4_stability_full:
    input:
        unpack(lambda wc: a4_model_inputs(f"{A4_STAB_ROOT}/seed{wc.seed}/CLAMPfull_canonical", A4_CAN_FILES, {
            "metadata": ancient(f"{A4_PREP}/metadata_filtered.rds"),
            "subsample": ancient(f"{a4_cov_cell(100, 1)}/subsample_info.rds"),
            "fbm": ancient(f"{A4_PREP}/fbm_filtered.bk"),
            "prior": ancient(A4_CAN_PRIOR["output_gmt"]),
            "svd": f"{A4_STAB_ROOT}/seed{wc.seed}/svd.rds",
            "k": f"{A4_STAB_ROOT}/seed{wc.seed}/CLAMP_K.rds",
            "base": f"{A4_STAB_ROOT}/seed{wc.seed}/CLAMPbase.rds",
            "script": "scripts/archs4/coverage/fit_clampfull.R",
            "common": "scripts/archs4/common.R",
            "raw_metadata": f"{A4_PREP}/metadata.rds",
        })),
    output:
        full=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPfull_canonical/complete",
        model=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPfull_canonical/CLAMPfull_canonical.rds",
        z=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPfull_canonical/Z.csv",
        b=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPfull_canonical/B.csv",
        summary=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPfull_canonical/summary.csv",
        manifest=f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPfull_canonical/manifest.json",
    log:
        f"{A4_STAB_ROOT}/logs/seed{{seed}}.full.log",
    params:
        subsample=f"{a4_cov_cell(100, 1)}/subsample_info.rds",
        base_dir=lambda wc: f"{A4_STAB_ROOT}/seed{wc.seed}",
        out_dir=lambda wildcards, output: os.path.dirname(output.full),
        index=lambda wc: A4_STAB_SEEDS.index(int(wc.seed)) + 2,
    threads: int(A4_STAB_CFG["resources"]["threads"])
    resources:
        mem_mb=int(A4_STAB_CFG["resources"]["mem_mb"]),
        runtime=int(A4_STAB_CFG["resources"]["runtime"]),
    wildcard_constraints:
        seed=A4_STAB_FIT_SEED_PATTERN,
    conda: "clamp-analyses"
    shell:
        "(Rscript scripts/archs4/coverage/fit_clampfull.R --dataset archs4 --metadata {A4_PREP}/metadata_filtered.rds "
        "--subsample-info {params.subsample} --fbm-backing {A4_PREP}/fbm_filtered.bk --prior-gmt {A4_CAN_PRIOR[output_gmt]} "
        "--svd {params.base_dir}/svd.rds --base-model {params.base_dir}/CLAMPbase.rds --k {params.base_dir}/CLAMP_K.rds --out-dir {params.out_dir} "
        "--model-name CLAMPfull_canonical --seed {wildcards.seed} --seed-index {params.index} "
        "--fraction 100 --max-iter 5000 --multiplier 100 --requested-cpus {threads} "
        "--requested-mem-mb {resources.mem_mb} ) > {log} 2>&1"


# ============================================================
# Step 3: collect completed models
# ============================================================

rule archs4_stability_complete:
    input:
        f"{A4_STAB_ROOT}/seed{{seed}}/CLAMPfull_canonical/complete",
    output:
        touch(f"{A4_STAB_ROOT}/seed{{seed}}/complete"),
    wildcard_constraints:
        seed="|".join(map(str, A4_STAB_SEEDS)),


rule archs4_stability:
    input:
        expand(f"{A4_STAB_ROOT}/seed{{seed}}/complete", seed=A4_STAB_SEEDS),

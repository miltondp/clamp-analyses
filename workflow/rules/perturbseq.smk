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
# Each strategy reads the fits written by its fitting area (pseudobulk.smk, donor_bulk.smk).
PS_Z_ROOTS = {
    "main": config["paths"]["production"],
    "donor_bulk": config["donor_bulk"]["paths"]["production"],
}
PS_REGULONS = os.path.join(PS_CFG["output_root"], "regulons", "{strategy}")

# Step 7: seed and k refits of the main fits, scored as strategy `variants` (model <M>__<variant>)
# against the main regulons. They are not in PS_STRATEGIES, so the main summaries leave them out.
PS_VAR = PS_CFG["variants"]
PS_VAR_ROOT = PS_VAR["root"]
PS_VARIANT_MODELS = (
    [f"{m}__seed{s}" for s in PS_VAR["seeds"] for m in PS_VAR["seed_models"]]
    + [f"{m}__k{k}" for k in PS_VAR["ks"] for m in PS_VAR["k_models"]]
)


def ps_regulon_strategy(strategy):
    """The strategy whose regulons (and gene universe) a strategy is scored against."""
    return "main" if strategy == "variants" else strategy


def ps_model_dir(strategy, model):
    if strategy == "variants":
        base, variant = model.split("__")
        return os.path.join(PS_VAR_ROOT, PS_CFG["dataset"], variant, "models", base)
    return os.path.join(PS_Z_ROOTS[strategy], PS_CFG["dataset"], "models", model)


def ps_genes_file(strategy):
    """The strategy's preprocessing gene list: the regulon universe source, which every Z must match."""
    return os.path.join(PS_Z_ROOTS[ps_regulon_strategy(strategy)], PS_CFG["dataset"],
                        "preprocessing", "row_stats.csv")


def ps_regulon_file(name):
    """Input function: a perturbseq_regulons output (`gmt`, `logfc`, `sizes`) for the scored strategy."""
    def _file(wildcards):
        return getattr(rules.perturbseq_regulons.output, name).format(
            strategy=ps_regulon_strategy(wildcards.strategy))
    return _file


rule perturbseq_regulons:
    input:
        logfc=rules.ota2025_k562gw.output.logfc,
        p=rules.ota2025_k562gw.output.p,
        ids=rules.ota2025_k562gw.output.ids,
        # Regulons depend on the screen and the preprocessing gene list only, never on a fit.
        universe=lambda wc: ps_genes_file(wc.strategy),
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
        "--id-map {input.ids:q} --universe {input.universe:q} --q {params.cfg[q]} --min-size {params.cfg[min_size]} "
        "--max-size {params.cfg[max_size]} --jaccard {params.cfg[jaccard]} --min-de {params.min_de} "
        "--out-dir {params.out_dir:q} > {log:q} 2>&1"


# Step 5: regulon ORA. run_ora.R gives the observed hits per tail (Z and -Z, restricted to the
# regulon universe); ora_null.py reproduces them and adds the gene-label permutation null.
PS_ORA = os.path.join(PS_CFG["output_root"], "ora", "{strategy}", "{model}")
PS_SIGNED = os.path.join(PS_CFG["output_root"], "signed_z", "{strategy}", "{model}")
PS_TAILS = ["pos", "neg"]
# projections.smk constrains {model} globally; these rules override it with the Perturb-seq models.
PS_WC = dict(
    strategy="|".join([*PS_STRATEGIES, "variants"]),
    model="|".join(sorted({m for cfg in PS_STRATEGIES.values() for m in cfg["models"]}
                          | set(PS_VARIANT_MODELS))),
)


rule perturbseq_signed_z:
    input:
        z=lambda wc: os.path.join(ps_model_dir(wc.strategy, wc.model), "Z.csv"),
        b=lambda wc: os.path.join(ps_model_dir(wc.strategy, wc.model), "B.csv"),
        genes=lambda wc: ps_genes_file(wc.strategy),
        universe=ps_regulon_file("logfc"),
    output:
        pos=os.path.join(PS_SIGNED, "pos.csv"),
        neg=os.path.join(PS_SIGNED, "neg.csv"),
        manifest=os.path.join(PS_SIGNED, "manifest.json"),
    wildcard_constraints:
        **PS_WC,
    resources:
        mem_mb=4000,
        runtime=10,
    conda: "clamp-analyses"
    shell:
        "python scripts/perturbseq/signed_z.py --z {input.z:q} --b {input.b:q} --genes {input.genes:q} "
        "--universe {input.universe:q} "
        "--out-pos {output.pos:q} --out-neg {output.neg:q} --out-manifest {output.manifest:q}"


rule perturbseq_ora:
    input:
        z=os.path.join(PS_SIGNED, "{tail}.csv"),
        manifest=rules.perturbseq_signed_z.output.manifest,
        gmt=ps_regulon_file("gmt"),
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(os.path.join(PS_ORA, "{tail}")),
    log:
        os.path.join(PS_ORA, "{tail}.log"),
    params:
        ora=PS_CFG["ora"],
        reg=PS_CFG["regulons"],
    wildcard_constraints:
        **PS_WC,
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
        z=rules.perturbseq_signed_z.output.pos,
        ora_pos=os.path.join(PS_ORA, "pos"),
        ora_neg=os.path.join(PS_ORA, "neg"),
        gmt=ps_regulon_file("gmt"),
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
        **PS_WC,
    resources:
        mem_mb=int(PS_CFG["ora"]["null_resources"]["mem_mb"]),
        runtime=int(PS_CFG["ora"]["null_resources"]["runtime"]),
    conda: "clamp-analyses"
    shell:
        "python scripts/perturbseq/ora_null.py --z {input.z:q} --gmt {input.gmt:q} "
        "--ora-pos {input.ora_pos:q} --ora-neg {input.ora_neg:q} "
        "--strategy {wildcards.strategy} --model {wildcards.model} "
        "--top-pct {params.ora[top_pct]} --min-size {params.reg[min_size]} --max-size {params.reg[max_size]} "
        "--fdr {params.ora[fdr]} --n-perm {PS_CFG[n_perm]} --seed {PS_CFG[seed]} "
        "--out-dir {params.out_dir:q} > {log:q} 2>&1"


def ps_per_model(pattern):
    """`pattern` expanded over every strategy's models."""
    return [expand(pattern, strategy=s, model=cfg["models"]) for s, cfg in PS_STRATEGIES.items()]


# Shell command that stacks the {input} CSVs into {output}.
PS_CONCAT_CSV = (
    "python -c \"import sys, pandas as pd; "
    "pd.concat([pd.read_csv(f) for f in sys.argv[2:]]).to_csv(sys.argv[1], index=False)\" "
    "{output:q} {input:q}"
)


rule perturbseq_ora_summary:
    input:
        ps_per_model(rules.perturbseq_ora_null.output.summary),
    output:
        os.path.join(PS_CFG["output_root"], "ora_summary.csv"),
    conda: "clamp-analyses"
    shell:
        PS_CONCAT_CSV


# Step 6: perturbation-response capture. max |r| between each perturbation's logFC (raw, and without
# its top shared components) and the LVs, against the same gene-label permutation null as the ORA.
PS_RESPONSE = os.path.join(PS_CFG["output_root"], "response", "{strategy}", "{model}")


rule perturbseq_response:
    input:
        z=rules.perturbseq_signed_z.output.pos,
        logfc=ps_regulon_file("logfc"),
        sizes=ps_regulon_file("sizes"),
    output:
        perturbations=os.path.join(PS_RESPONSE, "perturbations.csv"),
        null=os.path.join(PS_RESPONSE, "null.csv.gz"),
        summary=os.path.join(PS_RESPONSE, "summary.csv"),
    log:
        os.path.join(PS_RESPONSE, "response_capture.log"),
    params:
        out_dir=lambda wc, output: os.path.dirname(output.summary),
        alpha=PS_CFG["response"]["alpha"],
        remove_pcs=PS_CFG["response"]["remove_pcs"],  # a list: the shell formatter space-joins it
    wildcard_constraints:
        **PS_WC,
    resources:
        mem_mb=int(PS_CFG["response"]["resources"]["mem_mb"]),
        runtime=int(PS_CFG["response"]["resources"]["runtime"]),
    conda: "clamp-analyses"
    shell:
        "python scripts/perturbseq/response_capture.py --z {input.z:q} --logfc {input.logfc:q} "
        "--sizes {input.sizes:q} --strategy {wildcards.strategy} --model {wildcards.model} "
        "--remove-pcs {params.remove_pcs} --n-perm {PS_CFG[n_perm]} --seed {PS_CFG[seed]} --alpha {params.alpha} "
        "--out-dir {params.out_dir:q} > {log:q} 2>&1"


rule perturbseq_response_summary:
    input:
        ps_per_model(rules.perturbseq_response.output.summary),
    output:
        os.path.join(PS_CFG["output_root"], "response_summary.csv"),
    conda: "clamp-analyses"
    shell:
        PS_CONCAT_CSV


# Step 6.4: aucodds / F-aucodds (Saelens et al. 2018) of the LV modules against the regulons, with
# the same gene-label permutation null.
PS_MODULES = os.path.join(PS_CFG["output_root"], "modules", "{strategy}", "{model}")


rule perturbseq_module_scores:
    input:
        z=rules.perturbseq_signed_z.output.pos,
        gmt=ps_regulon_file("gmt"),
    output:
        null=os.path.join(PS_MODULES, "null.csv.gz"),
        module_hits=os.path.join(PS_MODULES, "module_hits.csv"),
        regulon_hits=os.path.join(PS_MODULES, "regulon_hits.csv"),
        summary=os.path.join(PS_MODULES, "summary.csv"),
    log:
        os.path.join(PS_MODULES, "module_scores.log"),
    params:
        out_dir=lambda wc, output: os.path.dirname(output.summary),
        mod=PS_CFG["modules"],
        ora=PS_CFG["ora"],
        reg=PS_CFG["regulons"],
    wildcard_constraints:
        **PS_WC,
    resources:
        mem_mb=int(PS_CFG["modules"]["resources"]["mem_mb"]),
        runtime=int(PS_CFG["modules"]["resources"]["runtime"]),
    conda: "clamp-analyses"
    shell:
        "python scripts/perturbseq/module_scores.py --z {input.z:q} --gmt {input.gmt:q} "
        "--strategy {wildcards.strategy} --model {wildcards.model} "
        "--top-pct {params.ora[top_pct]} --min-size {params.reg[min_size]} --max-size {params.reg[max_size]} "
        "--module-min-size {params.mod[min_size]} --zscore-sd {params.mod[zscore_sd]} --q {params.mod[q]} "
        "--log10-or-max {params.mod[log10_or_max]} --n-perm {PS_CFG[n_perm]} --seed {PS_CFG[seed]} "
        "--out-dir {params.out_dir:q} > {log:q} 2>&1"


rule perturbseq_module_scores_summary:
    input:
        ps_per_model(rules.perturbseq_module_scores.output.summary),
    output:
        os.path.join(PS_CFG["output_root"], "module_scores_summary.csv"),
    conda: "clamp-analyses"
    shell:
        PS_CONCAT_CSV


# Step 7.1: seed and k refits of the main OneK1K fits (PLAN.md). Each reuses its production rule
# (shell, conda, resources and the other params inherited) and overrides only the outputs, the
# output directories, --seed and the k file. `use rule ... with:` replaces a whole directive, so
# ps_override() rebuilds it from the production rule.
PS_VAR_MODELS = os.path.join(PS_VAR_ROOT, "{dataset}", "{variant}", "models")
PS_VAR_WC = dict(dataset=PS_CFG["dataset"], variant=r"seed\d+|k\d+")


def ps_override(named, **changes):
    out = dict(named.items())
    out.update(changes)
    return out


def ps_variant_seed(production_seed):
    """Params function: the seed of a seed<s> variant; k variants keep the production seed."""
    def _seed(wildcards):
        v = wildcards.variant
        return int(v[len("seed"):]) if v.startswith("seed") else production_seed
    return _seed


def ps_variant_k(wildcards):
    """The k file: the variant's own for k<k>, the production one for seed<s>."""
    if wildcards.variant.startswith("k"):
        return os.path.join(PS_VAR_ROOT, wildcards.dataset, wildcards.variant, "k.csv")
    return expand(rules.preprocess_pseudobulk.output.k, dataset=wildcards.dataset)[0]


def ps_variant_dir(model):
    return lambda wc: os.path.join(PS_VAR_ROOT, wc.dataset, wc.variant, "models", model)


rule perturbseq_variant_k:
    output:
        os.path.join(PS_VAR_ROOT, "{dataset}", "k{k}", "k.csv"),
    wildcard_constraints:
        dataset=PS_CFG["dataset"],
        k=r"\d+",
    shell:
        "printf 'k\\n%s\\n' {wildcards.k} > {output:q}"


use rule clamp_pseudobulk as perturbseq_variant_clamp with:
    input:
        **ps_override(rules.clamp_pseudobulk.input, k=ps_variant_k),
    output:
        base_B=os.path.join(PS_VAR_MODELS, "CLAMPbase", "B.csv"),
        base_Z=os.path.join(PS_VAR_MODELS, "CLAMPbase", "Z.csv"),
        base_rds=os.path.join(PS_VAR_MODELS, "CLAMPbase", "CLAMPbase.rds"),
        full_B=os.path.join(PS_VAR_MODELS, "CLAMPfull", "B.csv"),
        full_Z=os.path.join(PS_VAR_MODELS, "CLAMPfull", "Z.csv"),
        full_rds=os.path.join(PS_VAR_MODELS, "CLAMPfull", "CLAMPfull.rds"),
        l2=os.path.join(PS_VAR_MODELS, "CLAMPfull", "L2.csv"),
        summary=os.path.join(PS_VAR_MODELS, "CLAMPfull", "summary.csv"),
    params:
        **ps_override(rules.clamp_pseudobulk.params, base_dir=ps_variant_dir("CLAMPbase"),
                      full_dir=ps_variant_dir("CLAMPfull"), seed=ps_variant_seed(config["clamp"]["seed"])),
    wildcard_constraints:
        **PS_VAR_WC,


use rule plier_pseudobulk as perturbseq_variant_plier with:
    input:
        **ps_override(rules.plier_pseudobulk.input, k=ps_variant_k),
    output:
        B=os.path.join(PS_VAR_MODELS, "PLIER", "B.csv"),
        Z=os.path.join(PS_VAR_MODELS, "PLIER", "Z.csv"),
        rds=os.path.join(PS_VAR_MODELS, "PLIER", "PLIER.rds"),
        summary=os.path.join(PS_VAR_MODELS, "PLIER", "summary.csv"),
    params:
        **ps_override(rules.plier_pseudobulk.params, out_dir=ps_variant_dir("PLIER"),
                      seed=ps_variant_seed(config["preprocess"]["seed"])),
    wildcard_constraints:
        **PS_VAR_WC,


use rule flashier_pseudobulk as perturbseq_variant_flashier with:
    input:
        **ps_override(rules.flashier_pseudobulk.input, k=ps_variant_k),
    output:
        B=os.path.join(PS_VAR_MODELS, "flashier", "B.csv"),
        Z=os.path.join(PS_VAR_MODELS, "flashier", "Z.csv"),
        rds=os.path.join(PS_VAR_MODELS, "flashier", "flashier_model.rds"),
    params:
        **ps_override(rules.flashier_pseudobulk.params, out_dir=ps_variant_dir("flashier"),
                      seed=ps_variant_seed(config["preprocess"]["seed"])),
    wildcard_constraints:
        **PS_VAR_WC,


use rule mofa_flex_pseudobulk as perturbseq_variant_mofa_flex with:
    input:
        **ps_override(rules.mofa_flex_pseudobulk.input, k=ps_variant_k),
    output:
        B=os.path.join(PS_VAR_MODELS, "MOFA_FLEX_PRIOR", "B.csv"),
        Z=os.path.join(PS_VAR_MODELS, "MOFA_FLEX_PRIOR", "Z.csv"),
        model=os.path.join(PS_VAR_MODELS, "MOFA_FLEX_PRIOR", "model.pkl"),
    params:
        **ps_override(rules.mofa_flex_pseudobulk.params, out_dir=ps_variant_dir("MOFA_FLEX_PRIOR"),
                      seed=ps_variant_seed(config["mofa_flex"]["seed"])),
    wildcard_constraints:
        **PS_VAR_WC,


use rule cogaps_pseudobulk as perturbseq_variant_cogaps with:
    input:
        **ps_override(rules.cogaps_pseudobulk.input, k=ps_variant_k),
    output:
        B=os.path.join(PS_VAR_MODELS, "CoGAPS", "B.csv"),
        Z=os.path.join(PS_VAR_MODELS, "CoGAPS", "Z.csv"),
        rds=os.path.join(PS_VAR_MODELS, "CoGAPS", "cogaps_model.rds"),
    params:
        **ps_override(rules.cogaps_pseudobulk.params, out_dir=ps_variant_dir("CoGAPS"),
                      seed=ps_variant_seed(config["cogaps"]["seed"])),
    wildcard_constraints:
        **PS_VAR_WC,


# The production rule fits PCA, NMF and ICA in one job, so ICA and NMF variants get their own rule
# with the same script. NMF seed variants use a seeded random start (as in cNMF replicates); NMF k
# variants keep the production nndsvd init.
rule perturbseq_variant_ica_nmf:
    input:
        norm=rules.pca_nmf_ica_pseudobulk.input.norm,
        k=ps_variant_k,
    output:
        B=os.path.join(PS_VAR_MODELS, "{method}", "B.csv"),
        Z=os.path.join(PS_VAR_MODELS, "{method}", "Z.csv"),  # the script also writes <method>_model.pkl
    params:
        out_dir=lambda wc, output: os.path.dirname(output.Z),
        seed=ps_variant_seed(config["preprocess"]["seed"]),
        nmf_init=lambda wc: "random" if wc.method == "NMF" and wc.variant.startswith("seed") else "nndsvd",
    wildcard_constraints:
        **PS_VAR_WC,
        method="ICA|NMF",
    resources:
        mem_mb=48000,
        runtime=720,
    conda: "clamp-analyses"
    shell:
        "python scripts/pseudobulk/pca_nmf_ica.py --norm {input.norm} --k {input.k} "
        "--models-dir {params.out_dir} --method {wildcards.method} --flat-output "
        "--seed {params.seed} --nmf-init {params.nmf_init}"


def ps_variant_fit_files(wildcards=None):
    return [os.path.join(ps_model_dir("variants", m), "Z.csv") for m in PS_VARIANT_MODELS]


rule perturbseq_variant_fits:
    input:
        ps_variant_fit_files(),


# OneK1K cell-type mean r of each fit (score_per_ct(), as in the Fig. 2A benchmark), for production
# and variant fits alike. Main and variants only: donor-bulk samples are donors, not the pseudobulk
# samples in truthFrac_v0.csv.
PS_CELLTYPE = os.path.join(PS_CFG["output_root"], "celltype", "{strategy}", "{model}.csv")


rule perturbseq_celltype_score:
    input:
        b=lambda wc: os.path.join(ps_model_dir(wc.strategy, wc.model), "B.csv"),
        truth=expand(rules.build_pseudobulk.output.truth0, dataset=PS_CFG["dataset"]),
        script="scripts/perturbseq/celltype_score.R",
    output:
        PS_CELLTYPE,
    log:
        PS_CELLTYPE.replace(".csv", ".log"),
    params:
        exclude=",".join(config["cell_type_analysis"]["excluded_targets"].get(PS_CFG["dataset"], [])),
    wildcard_constraints:
        strategy="main|variants",
        model=PS_WC["model"],
    resources:
        mem_mb=4000,
        runtime=10,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --b {input.b:q} --truth {input.truth:q} --dataset {PS_CFG[dataset]} "
        "--strategy {wildcards.strategy} --model {wildcards.model} --exclude '{params.exclude}' "
        "--out {output:q} > {log:q} 2>&1"


rule perturbseq_celltype_summary:
    input:
        expand(PS_CELLTYPE, strategy="main", model=PS_STRATEGIES["main"]["models"]),
    output:
        os.path.join(PS_CFG["output_root"], "celltype_summary.csv"),
    conda: "clamp-analyses"
    shell:
        PS_CONCAT_CSV

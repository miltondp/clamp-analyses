import os


A4_CRISPERCAS_CFG = A4_CFG["crispercas"]
A4_CRISPERCAS_ROOT = A4_CRISPERCAS_CFG["compute_root"]
A4_METALV = A4_CRISPERCAS_CFG["metalv_root"]
A4_METALV_NOISE = A4_CRISPERCAS_CFG["noise_seeds_root"]
A4_METALV_SEEDS = sorted(A4_CFG["stability"]["seeds"], key=lambda s: (s < 1000, s))
A4_METALV_MODELS = ["S123"] + [f"S{s}" for s in A4_METALV_SEEDS]
A4_METALV_TN = A4_CFG["technical_noise"]["compute_root"]
A4_METALV_ORA = ["bp", "canonical", "cellmarker", "archs4_cell_lines", "azimuth"]


def a4_seed_model(seed):
    return f"{A4_CFG['stability']['root']}/seed{seed}/CLAMPfull_canonical"


def a4_metalv_z(model):
    if model == "S123":
        return f"{A4_CAN_FINAL_ROOT}/archs4/Z.csv"
    return f"{a4_seed_model(model[1:])}/Z.csv"


def a4_metalv_ora_args(db):
    if db == "azimuth":
        return dict(path=A4_CRISPERCAS_CFG["azimuth_gmt"], type="gmt", label="Azimuth", extra="",
                    min_size=A4_CRISPERCAS_CFG["azimuth_min_size"], max_size=A4_CRISPERCAS_CFG["azimuth_max_size"])
    if db == "archs4_cell_lines":
        return dict(path=A4_CRISPERCAS_CFG["archs4_cell_lines_gmt"], type="gmt", label="ARCHS4 Cell-lines", extra="",
                    min_size=A4_CFG["ora"]["min_size"], max_size=A4_CRISPERCAS_CFG["archs4_cell_lines_max_size"])
    d = A4_CFG["ora"]["databases"][db]
    extra = f"--exclude-term-regex '{d['exclude_term_regex']}' " if "exclude_term_regex" in d else ""
    if "sheet" in d:
        extra += f"--sheet {d['sheet']} --term-column {d['term_column']} --gene-column {d['gene_column']}"
    return dict(path=d["path"], type=d["type"], label=d["label"], extra=extra,
                min_size=A4_CFG["ora"]["min_size"], max_size=A4_CFG["ora"]["max_size"])


# ============================================================
# Step 1: identify technical noise
# ============================================================

rule crispercas_seed_noise:
    input:
        z=lambda wc: ancient(f"{a4_seed_model(wc.seed)}/Z.csv"),
        b=lambda wc: ancient(f"{a4_seed_model(wc.seed)}/B.csv"),
        annotation=lambda wc: ancient(f"{a4_seed_model(wc.seed)}/summary.csv"),
        gene_properties=ancient(f"{A4_METALV_TN}/gene_properties.csv"),
        sample_metadata=ancient(f"{A4_METALV_TN}/sample_metadata.tsv.gz"),
        script=ancient("scripts/archs4/technical_noise/run_analysis.R"),
    output:
        gene=f"{A4_METALV_NOISE}/seed{{seed}}/lv_gene_properties_canonical.csv",
        sample=f"{A4_METALV_NOISE}/seed{{seed}}/lv_sample_metadata_canonical.csv",
    log:
        gene=f"{A4_METALV_NOISE}/seed{{seed}}/gene_level.log",
        sample=f"{A4_METALV_NOISE}/seed{{seed}}/sample_level.log",
    params:
        model=lambda wc: a4_seed_model(wc.seed),
        out=lambda wc: f"{A4_METALV_NOISE}/seed{wc.seed}",
    threads: 4
    resources:
        mem_mb=96000,
        runtime=1440,
    conda: "clamp-analyses"
    shell:
        "ln -sf $(readlink -f {input.gene_properties}) {params.out}/gene_properties.csv && "
        "ln -sf $(readlink -f {input.sample_metadata}) {params.out}/sample_metadata.tsv.gz && "
        "Rscript {input.script} --stage gene_level --model-dir {params.model} --out-dir {params.out} > {log.gene} 2>&1 && "
        "Rscript {input.script} --stage sample_level --model-dir {params.model} --out-dir {params.out} "
        "--n-cores {threads} > {log.sample} 2>&1"


rule crispercas_noise_flags:
    input:
        seeds=expand(rules.crispercas_seed_noise.output, seed=A4_METALV_SEEDS),
        reference=[ancient(f"{A4_METALV_TN}/{f}") for f in ["lv_gene_properties_canonical.csv",
                                                             "lv_sample_metadata_canonical.csv",
                                                             "lv_technical_classification.csv"]],
        script="scripts/archs4/crispercas/cluster_metalv.py",
    output:
        f"{A4_METALV}/noise/all_flags.csv",
    conda: "clamp-analyses"
    shell:
        "python {input.script} noise --reference-dir {A4_METALV_TN} --seeds-dir {A4_METALV_NOISE} "
        "--seeds {A4_METALV_SEEDS} --threshold {A4_CRISPERCAS_CFG[noise_r2]} --out {output}"


# ============================================================
# Step 2: cluster gene modules across models
# ============================================================

rule crispercas_metalv_zcache:
    input:
        z=[ancient(a4_metalv_z(m)) for m in A4_METALV_MODELS],
        script="scripts/archs4/crispercas/cluster_metalv.py",
    output:
        f"{A4_METALV}/Z_all.npz",
    resources:
        mem_mb=16000,
    conda: "clamp-analyses"
    shell:
        "python {input.script} cache --models {A4_METALV_MODELS} --z {input.z} --out {output}"


rule crispercas_metalv_cluster:
    input:
        z=rules.crispercas_metalv_zcache.output,
        flags=rules.crispercas_noise_flags.output,
        script="scripts/archs4/crispercas/cluster_metalv.py",
    output:
        metaz=f"{A4_METALV}/metaZ.csv",
        summary=f"{A4_METALV}/cluster_summary.csv",
        members=f"{A4_METALV}/cluster_members.csv",
        all_lvs=f"{A4_METALV}/all_lv_clusters.csv",
        counts=f"{A4_METALV}/cluster_model_count_distribution.csv",
    log:
        f"{A4_METALV}/cluster.log"
    resources:
        mem_mb=16000,
    conda: "clamp-analyses"
    shell:
        "python {input.script} cluster --z-cache {input.z} --flags {input.flags} "
        "--min-cos {A4_CRISPERCAS_CFG[min_cos]} --min-models {A4_CRISPERCAS_CFG[min_models]} "
        "--out-dir {A4_METALV} > {log} 2>&1"


# ============================================================
# Step 3: test screen, pathway, and trait associations
# ============================================================

rule crispercas_metalv_fgsea:
    input:
        z=rules.crispercas_metalv_cluster.output.metaz,
        lipid_deg=A4_CRISPERCAS_CFG["lipid_deg"],
        script="scripts/archs4/crispercas/gene_enrichment_crisprcas.R",
    output:
        f"{A4_METALV}/gsea/fgsea_all_lvs.csv",
    log:
        f"{A4_METALV}/gsea/fgsea.log"
    resources:
        mem_mb=16000,
        runtime=60,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --lipid-deg {input.lipid_deg} --out {output} "
        "--n-reps {A4_CRISPERCAS_CFG[n_reps]} > {log} 2>&1"


rule crispercas_metalv_gsea_best:
    input:
        fgsea=rules.crispercas_metalv_fgsea.output,
        script="scripts/archs4/crispercas/select_metalvs.py",
    output:
        f"{A4_METALV}/gsea/gsea_best.csv",
    conda: "clamp-analyses"
    shell:
        "python {input.script} gsea-best --fgsea {input.fgsea} --out {output}"


rule fetch_azimuth_gmt:
    output:
        A4_CRISPERCAS_CFG["azimuth_gmt"],
    log:
        A4_CRISPERCAS_CFG["azimuth_gmt"] + ".log"
    conda: "snakemake"
    shell:
        "python scripts/inputs/reference_inputs.py --only azimuth --destination {output:q} > {log:q} 2>&1"


rule fetch_archs4_cell_lines_gmt:
    output:
        A4_CRISPERCAS_CFG["archs4_cell_lines_gmt"],
    log:
        A4_CRISPERCAS_CFG["archs4_cell_lines_gmt"] + ".log"
    conda: "snakemake"
    shell:
        "python scripts/inputs/reference_inputs.py --only archs4_cell_lines --destination {output:q} > {log:q} 2>&1"


rule crispercas_metalv_ora:
    input:
        z=rules.crispercas_metalv_cluster.output.metaz,
        manifest=f"{A4_CAN_FINAL_ROOT}/archs4/manifest.json",
        database=lambda wc: a4_metalv_ora_args(wc.database)["path"],
        script="scripts/archs4/coverage/run_ora.R",
    output:
        ora_dir=directory(f"{A4_METALV}/ora/{{database}}"),
    log:
        f"{A4_METALV}/ora/{{database}}.log"
    params:
        args=lambda wc: a4_metalv_ora_args(wc.database),
    wildcard_constraints:
        database="|".join(A4_METALV_ORA),
    resources:
        mem_mb=int(A4_CFG["ora"]["resources"]["mem_mb"]),
        runtime=int(A4_CFG["ora"]["resources"]["runtime"]),
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --z {input.z} --out-dir {output.ora_dir} "
        "--dataset archs4 --fraction 100 --seed 1 "
        "--model {A4_CAN_MODEL_NAME} --model-manifest {input.manifest} "
        "--database {wildcards.database} --database-label '{params.args[label]}' "
        "--database-type {params.args[type]} --database-path {input.database} {params.args[extra]} "
        "--top-pct {A4_CFG[ora][top_pct]} --min-size {params.args[min_size]} --max-size {params.args[max_size]} "
        "--pvalue-cutoff {A4_CFG[ora][pvalue_cutoff]} --qvalue-cutoff {A4_CFG[ora][qvalue_cutoff]} "
        "--fdr {A4_CFG[ora][fdr][0]} > {log} 2>&1"


# ============================================================
# Step 4: select and interpret metaLVs
# ============================================================

rule crispercas_symbol_ensembl_mapping:
    input:
        script="scripts/archs4/crispercas/extract_symbol_ensembl.R",
    output:
        mapping=A4_CRISPERCAS_CFG["symbol_ensembl"],
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} {output.mapping}"


rule crispercas_metalv_select:
    input:
        summary=rules.crispercas_metalv_cluster.output.summary,
        metaz=rules.crispercas_metalv_cluster.output.metaz,
        gsea=rules.crispercas_metalv_gsea_best.output,
        gls=A4_CRISPERCAS_CFG["gls_summary"],
        canonical=f"{A4_METALV}/ora/canonical",
        smultixcan=A4_CRISPERCAS_CFG["smultixcan_zscores"],
        mapping=ancient(A4_CRISPERCAS_CFG["symbol_ensembl"]),
        script="scripts/archs4/crispercas/select_metalvs.py",
    output:
        funnel=f"{A4_METALV}/selection/funnel.csv",
        selected=f"{A4_METALV}/selection/selected_metalvs.csv",
        pairs=f"{A4_METALV}/selection/screen_gene_trait_pairs.csv",
        gls=f"{A4_METALV}/selection/gls_metalv.csv",
    conda: "clamp-analyses"
    shell:
        "python {input.script} select --cluster-summary {input.summary} --metaz {input.metaz} --gsea-best {input.gsea} "
        "--gls {input.gls} --canonical-ora {input.canonical}/enrichment.csv.gz --smultixcan-zscores {input.smultixcan} "
        "--symbol-ensembl {input.mapping} --fdr {A4_CRISPERCAS_CFG[fdr]} --top-n {A4_CRISPERCAS_CFG[top_n]} "
        "--lipid-regex '{A4_CRISPERCAS_CFG[lipid_regex]}' --out-dir {A4_METALV}/selection"


rule crispercas_fig5_selection:
    input:
        summary=rules.crispercas_metalv_cluster.output.summary,
        select=rules.crispercas_metalv_select.output.selected,
    output:
        candidates=f"{A4_METALV}/fig5_selection/candidates.csv",
        selected=f"{A4_METALV}/fig5_selection/selected.csv",
    params:
        candidates=",".join(A4_CRISPERCAS_CFG["fig5_candidates"]),
        selected=",".join(A4_CRISPERCAS_CFG["fig5_metalvs"]),
    conda: "clamp-analyses"
    shell:
        "python -c \"import sys, pandas as pd; s = pd.read_csv(sys.argv[1]).set_index('metaLV'); "
        "c, f = sys.argv[3].split(','), sys.argv[4].split(','); assert set(f) <= set(c) <= set(s.index), (c, f); "
        "s.loc[c, ['n_models']].reset_index().to_csv(sys.argv[5], index=False); "
        "s.loc[f, ['n_models']].reset_index().to_csv(sys.argv[6], index=False)\" "
        "{input.summary} {input.select} {params.candidates} {params.selected} {output.candidates} {output.selected}"


rule crispercas_metalv_pathways:
    input:
        selected=rules.crispercas_fig5_selection.output.selected,
        ora=expand(f"{A4_METALV}/ora/{{database}}", database=["bp", "canonical"]),
        bp_gmt=A4_CFG["ora"]["databases"]["bp"]["path"],
        script="scripts/archs4/crispercas/merge_pathways.R",
    output:
        members=f"{A4_METALV}/pathways/canonical_bp_members.csv",
        top=f"{A4_METALV}/pathways/canonical_bp_top.csv",
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --selected {input.selected} --ora-dir {A4_METALV}/ora --bp-gmt {input.bp_gmt} "
        "--max-size {A4_CRISPERCAS_CFG[pathway_max_size]} --jaccard {A4_CRISPERCAS_CFG[pathway_jaccard]} "
        "--members-out {output.members} --top-out {output.top}"


# ============================================================
# Step 5: score drug signatures
# ============================================================

rule crispercas_metalv_l1000:
    input:
        metaz=rules.crispercas_metalv_cluster.output.metaz,
        selected=rules.crispercas_fig5_selection.output.candidates,
        model=ancient(f"{A4_CAN_FINAL_ROOT}/archs4/{A4_CAN_MODEL_NAME}.rds"),
        lincs=ancient(f"{A4_DD_DATA}/lincs-data.pkl"),
        drugbank=ancient(A4_DD_CFG["drugbank_tsv"]),
        mapping=ancient(A4_CRISPERCAS_CFG["symbol_ensembl"]),
        shortlist="workflow/config/fig5_drug_shortlist.csv",
        script="scripts/archs4/crispercas/l1000_metalv.py",
    output:
        scores=f"{A4_METALV}/l1000/l1000_scores.csv",
        ranking=f"{A4_METALV}/l1000/l1000_fig5_ranking.csv",
        unique=f"{A4_METALV}/l1000/l1000_unique_signatures.csv",
        curated=f"{A4_METALV}/l1000/l1000_fig5_curated.csv",
        metadata=f"{A4_METALV}/l1000/l1000_metadata.json",
    threads: 4
    resources:
        mem_mb=16000,
    conda: "clamp-analyses"
    shell:
        "OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} python {input.script} --metaz {input.metaz} --selected {input.selected} "
        "--reference-model {input.model} --lincs {input.lincs} --drugbank {input.drugbank} --symbol-ensembl {input.mapping} "
        "--shortlist {input.shortlist} --fdr {A4_CRISPERCAS_CFG[fdr]} "
        "--n-perm {A4_CRISPERCAS_CFG[l1000_n_perm]} --seed {A4_CRISPERCAS_CFG[l1000_seed]} --out-dir {A4_METALV}/l1000"


# ============================================================
# Step 6: project metaLV activity and prepare Figure 5
# ============================================================

rule crispercas_metalv_projection:
    input:
        metaz=rules.crispercas_metalv_cluster.output.metaz,
        selected=rules.crispercas_fig5_selection.output.candidates,
        model=ancient(f"{A4_CAN_FINAL_ROOT}/archs4/{A4_CAN_MODEL_NAME}.rds"),
        metadata=ancient(f"{A4_PREP}/metadata_filtered.rds"),
        samples=ancient(f"{A4_METALV_TN}/sample_metadata.tsv.gz"),
        h5=ancient(A4_CFG["raw_h5"]),
        script="scripts/archs4/crispercas/project_metalv.R",
    output:
        f"{A4_METALV}/projection/archs4_metalv_B.tsv.gz",
    log:
        f"{A4_METALV}/projection/projection.log"
    threads: 8
    resources:
        mem_mb=32000,
    conda: "clamp-analyses"
    shell:
        "Rscript {input.script} --metaz {input.metaz} --selected {input.selected} --model {input.model} "
        "--preprocess-dir {A4_PREP} --sample-metadata {input.samples} --h5 {input.h5} "
        "--n-cores {threads} --out {output} > {log} 2>&1"


rule crispercas_metalv_composition:
    input:
        b=rules.crispercas_metalv_projection.output,
        selected=rules.crispercas_fig5_selection.output.selected,
        script="scripts/archs4/crispercas/top_sample_composition.py",
    output:
        composition=f"{A4_METALV}/projection/composition.csv",
        enrichment=f"{A4_METALV}/projection/top_sample_enrichment.csv",
    params:
        top=A4_CRISPERCAS_CFG["fig5_top_samples"],
    conda: "clamp-analyses"
    shell:
        "python {input.script} --b {input.b} --selected {input.selected} --top-fraction {params.top[top_fraction]} "
        "--min-pct {params.top[min_pct]} --n-boot {params.top[n_boot]} --seed {params.top[seed]} "
        "--enrichment-out {output.enrichment} --out {output.composition}"


A4_FIG5_DATA = ["modules.csv", "module_genes.csv", "pathway_edges.csv", "module_trait_groups.csv",
                "module_traits_all.csv", "gene_trait_edges.csv", "ring_pvalues.csv", "twas_rcp_matrix.csv",
                "projection_composition.csv", "module_drugs.csv"]


rule crispercas_fig5_data:
    input:
        metaz=rules.crispercas_metalv_cluster.output.metaz,
        summary=rules.crispercas_metalv_cluster.output.summary,
        phenomexcan=A4_CRISPERCAS_CFG["phenomexcan_dir"],
        selected=rules.crispercas_fig5_selection.output.selected,
        pairs=rules.crispercas_metalv_select.output.pairs,
        gls=rules.crispercas_metalv_select.output.gls,
        members=rules.crispercas_metalv_pathways.output.members,
        top=rules.crispercas_metalv_pathways.output.top,
        composition=rules.crispercas_metalv_composition.output.composition,
        l1000=rules.crispercas_metalv_l1000.output.curated,
        script="scripts/archs4/crispercas/fig5_data.py",
    output:
        expand(f"{A4_METALV}/fig5/{{f}}", f=A4_FIG5_DATA),
    conda: "clamp-analyses"
    shell:
        "python {input.script} figure --metaz {input.metaz} --cluster-summary {input.summary} --selected {input.selected} "
        "--pairs {input.pairs} --pathway-members {input.members} --pathway-top {input.top} --gls {input.gls} "
        "--phenotypes-info {A4_CRISPERCAS_CFG[phenomexcan_dir]}/phenotypes_info.tsv "
        "--phenomexcan {A4_CRISPERCAS_CFG[phenomexcan_dir]} --composition {input.composition} --l1000 {input.l1000} "
        "--min-cases {A4_CRISPERCAS_CFG[min_cases]} --fdr {A4_CRISPERCAS_CFG[fdr]} "
        "--top-n {A4_CRISPERCAS_CFG[top_n]} --reference-genes " + ",".join(A4_CRISPERCAS_CFG["fig5_reference_genes"]) + " --out-dir {A4_METALV}/fig5"


rule crispercas_mgm_report:
    input:
        flags=rules.crispercas_noise_flags.output[0],
        summary=rules.crispercas_metalv_cluster.output.summary,
        ora=expand(rules.crispercas_metalv_ora.output.ora_dir, database=A4_METALV_ORA),
        funnel=rules.crispercas_metalv_select.output.funnel,
        pathways=rules.crispercas_metalv_pathways.output.top,
        modules=f"{A4_METALV}/fig5/modules.csv",
        traits=f"{A4_METALV}/fig5/module_trait_groups.csv",
        composition=rules.crispercas_metalv_composition.output.composition,
        ranking=rules.crispercas_metalv_l1000.output.ranking,
        curated=rules.crispercas_metalv_l1000.output.curated,
        unique=rules.crispercas_metalv_l1000.output.unique,
        notebook=f"{A4_BIO_NB}/04_crispercas/00_mgm.ipynb",
    output:
        complete=touch(f"{A4_METALV}/notebook.complete"),
    log:
        notebook=f"{A4_BIO_NB}/04_crispercas/00_mgm.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{A4_BIO_NB}/04_crispercas/00_mgm.ipynb"


rule biology_archs4_crispercas:
    input:
        rules.crispercas_fig5_data.output,
        rules.crispercas_mgm_report.output,
        expand(rules.crispercas_metalv_ora.output.ora_dir, database=A4_METALV_ORA),

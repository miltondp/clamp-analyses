PANELS = config["paths"]["panels"]
PANELS_NB = os.path.join(REPO_ROOT, config["paths"]["panels_notebooks"])


def frozen_gtex(relative_path):
    return ancient(f"{GTEX_BIO}/{relative_path}")


# ============================================================
# Step 1: render pseudobulk and GTEx panels
# ============================================================

rule fig2_panel:
    input:
        benchmark_long=rules.benchmark_pseudobulk.output.long,
        holdout_predictions=rules.grouped_cv_analysis_pseudobulk.output.predictions,
        holdout_thresholded_metrics=rules.grouped_cv_analysis_pseudobulk.output.thresholded_metrics,
        holdout_thresholded_summary=rules.grouped_cv_analysis_pseudobulk.output.thresholded_summary,
        marker=rules.disentangle_pseudobulk.output.module_panel_ready,
        hard=rules.hard_cell_types_pseudobulk.output.group_dot_ready,
        hard_cross_grid=rules.hard_cell_types_pseudobulk.output.cross_grid,
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        hard_pair_loading_ready=rules.hard_pair_loadings_pseudobulk.output.panel_ready,
        hard_pair_loading_stats=rules.hard_pair_loadings_pseudobulk.output.stats,
        subtissue_confusion=frozen_gtex("04_subtissues/subtissue_confusion_matrices.tsv"),
        subtissue_anatomical=frozen_gtex("04_subtissues/anatomical_subtissue_results.tsv"),
        liver_xcell_scatter=frozen_gtex("04_liver_disentangle_xcell_rf_true_labels/liver_xcell_scatter.csv"),
        liver_lv_pathways=frozen_gtex("04_liver_disentangle_xcell_rf_true_labels/liver_lv_pathways.csv"),
        multitissue_panel_ready=frozen_gtex("06_multitissue_xcell_recovery/multitissue_recovery_panel_ready.csv"),
        donor_bulk_predictions=rules.donor_bulk_recovery_report.output.fig2_predictions,
        donor_bulk_statistics=rules.donor_bulk_recovery_report.output.fig2_statistics,
        donor_bulk_perez_annotation=rules.donor_bulk_recovery_report.output.perez_annotation,
        donor_bulk_perez_activity=rules.donor_bulk_recovery_report.output.perez_activity,
        donor_bulk_perez_summary=rules.donor_bulk_recovery_report.output.perez_summary,
        style="scripts/panels/style.R",
        notebook=f"{PANELS_NB}/fig2.ipynb",
    output:
        png=f"{PANELS}/fig2/fig2.png",
        pdf=f"{PANELS}/fig2/fig2.pdf",
        svg=f"{PANELS}/fig2/fig2.svg",
        complete=touch(f"{PANELS}/fig2/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/fig2.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/fig2.ipynb"


rule supp1_panel:
    input:
        benchmark_long=rules.benchmark_pseudobulk.output.long,
        labels="data/pseudobulk/cell_type_labels.csv",
        bootstrap=rules.benchmark_pseudobulk.output.bootstrap,


        heatmap_long=rules.single_cell_recovery_pseudobulk.output.heatmap,
        pseudobulk_recovery=rules.single_cell_recovery_pseudobulk.output.recovery,
        pseudobulk_overall=rules.single_cell_recovery_pseudobulk.output.overall,
        corr_full=rules.benchmark_pseudobulk.output.corr,
        assignments=rules.benchmark_pseudobulk.output.assignments,
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        ari_data=frozen_gtex("00_kmeans_clustering/ari_data.csv"),
        ari_comparisons=frozen_gtex("00_kmeans_clustering/ari_comparisons.csv"),
        gene_fraction_ari_data=frozen_gtex("00_kmeans_clustering/gene_fraction_ari_data.csv"),
        gene_fraction_ari_comparisons=frozen_gtex("00_kmeans_clustering/gene_fraction_ari_comparisons.csv"),
        tissue_subtissue_heatmap=frozen_gtex("02_LV_importance_rf_true_labels_biology/tissue_subtissue_heatmap.csv"),
        z_matrix_subtissue=frozen_gtex("07_global_alignment_rf_true_labels/cumulative20/gtex_global_alignment_summary.csv"),
        marker=rules.disentangle_pseudobulk.output.module_panel_ready,
        hard=rules.hard_cell_types_pseudobulk.output.group_dot_ready,
        orthogonality_per_term=frozen_gtex("07_geneset_orthogonality/orthogonality_per_term.csv"),
        orthogonality_summary=frozen_gtex("07_geneset_orthogonality/orthogonality_summary.csv"),
        style="scripts/panels/style.R",
        notebook=f"{PANELS_NB}/supp1.ipynb",
    output:
        png=f"{PANELS}/supp1/supp1.png",
        pdf=f"{PANELS}/supp1/supp1.pdf",
        svg=f"{PANELS}/supp1/supp1.svg",
        complete=touch(f"{PANELS}/supp1/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp1.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp1.ipynb"


rule supp2_panel:
    input:
        purity=rules.donor_bulk_recovery_report.output.supp2_purity,
        umap_cells=rules.donor_bulk_recovery_report.output.supp2_umap_cells,
        umap_lvs=rules.donor_bulk_recovery_report.output.supp2_umap_lvs,
        donor_bulk_recovery=rules.donor_bulk_analysis_tables.output.recovery,
        donor_bulk_overall=rules.donor_bulk_analysis_tables.output.recovery_overall,
        style="scripts/panels/style.R",
        notebook=f"{PANELS_NB}/supp2.ipynb",
    output:
        png=f"{PANELS}/supp2/supp2.png",
        pdf=f"{PANELS}/supp2/supp2.pdf",
        svg=f"{PANELS}/supp2/supp2.svg",
        complete=touch(f"{PANELS}/supp2/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp2.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp2.ipynb"


rule supp3_panel:
    input:
        benchmark_long=rules.benchmark_pseudobulk.output.long,
        labels="data/pseudobulk/cell_type_labels.csv",
        bootstrap=rules.benchmark_pseudobulk.output.bootstrap,
        heatmap_long=rules.single_cell_recovery_pseudobulk.output.heatmap,
        corr_full=rules.benchmark_pseudobulk.output.corr,
        assignments=rules.benchmark_pseudobulk.output.assignments,
        related_corr=rules.hard_cell_types_pseudobulk.output.correlations,
        ari_data=frozen_gtex("00_kmeans_clustering/ari_data.csv"),
        ari_comparisons=frozen_gtex("00_kmeans_clustering/ari_comparisons.csv"),
        gene_fraction_ari_data=frozen_gtex("00_kmeans_clustering/gene_fraction_ari_data.csv"),
        gene_fraction_ari_comparisons=frozen_gtex("00_kmeans_clustering/gene_fraction_ari_comparisons.csv"),
        gtex_top1_detail=frozen_gtex("07_global_alignment_rf_true_labels/top1/gtex_global_alignment_detail.csv"),
        gtex_ora_per_lv=frozen_gtex("07_global_alignment_rf_true_labels/gtex_tissue_ora_per_lv.csv"),
        marker=rules.disentangle_pseudobulk.output.module_panel_ready,
        hard=rules.hard_cell_types_pseudobulk.output.group_dot_ready,
        orthogonality_per_term=frozen_gtex("07_geneset_orthogonality/orthogonality_per_term.csv"),
        orthogonality_summary=frozen_gtex("07_geneset_orthogonality/orthogonality_summary.csv"),
        style="scripts/panels/style.R",
        notebook=f"{PANELS_NB}/supp3.ipynb",
    output:
        png=f"{PANELS}/supp3/supp3.png",
        pdf=f"{PANELS}/supp3/supp3.pdf",
        svg=f"{PANELS}/supp3/supp3.svg",
        complete=touch(f"{PANELS}/supp3/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp3.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp3.ipynb"


rule supp4_panel:
    input:
        runtime_seed_totals=rules.computational_timing_report_pseudobulk.output.seed_totals,
        runtime_gtex_per_fit=rules.computational_timing_report_gtex.output.per_fit,
        style="scripts/panels/style.R",
        notebook=f"{PANELS_NB}/supp4.ipynb",
    output:
        png=f"{PANELS}/supp4/supp4.png",
        pdf=f"{PANELS}/supp4/supp4.pdf",
        svg=f"{PANELS}/supp4/supp4.svg",
        complete=touch(f"{PANELS}/supp4/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp4.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp4.ipynb"


rule panels_pseudobulk:
    input:
        rules.fig2_panel.output.complete,
        rules.supp1_panel.output.complete,
        rules.supp2_panel.output.complete,
        rules.supp3_panel.output.complete,
        rules.supp4_panel.output.complete,


# ============================================================
# Step 2: render ARCHS4 supplementary panels
# ============================================================

rule supp5_panel:
    input:
        coverage=rules.aggregate_bp_coverage.output.coverage_long,
        saturation=rules.aggregate_bp_saturation.output.saturation_long,
        traits=rules.aggregate_archs4_coverage_traits.output.traits_long,
        finals=rules.aggregate_archs4_coverage_traits.output.finals_long,
        saturation_traits=rules.aggregate_archs4_saturation_traits.output.traits_long,
        notebook=f"{PANELS_NB}/supp5.ipynb",
        style="scripts/panels/style.R",
    output:
        png=f"{PANELS}/supp5/supp5.png",
        pdf=f"{PANELS}/supp5/supp5.pdf",
        svg=f"{PANELS}/supp5/supp5.svg",
        complete=touch(f"{PANELS}/supp5/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp5.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp5.ipynb"


rule supp6_panel:
    input:
        rules.projection_benchmark_figure.output,
    output:
        complete=touch(f"{PANELS}/supp6/notebook.complete"),


rule supp7_panel:
    input:
        aggregate=rules.aggregate_drug_diseases.output,
        drug_diseases=rules.prepare_three_model_comparison_data.output,
        notebook=f"{PANELS_NB}/supp7.ipynb",
        style="scripts/panels/style.R",
    output:
        png=f"{PANELS}/supp7/supp7.png",
        pdf=f"{PANELS}/supp7/supp7.pdf",
        svg=f"{PANELS}/supp7/supp7.svg",
        complete=touch(f"{PANELS}/supp7/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp7.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp7.ipynb"


rule supp8_panel:
    input:
        z=rules.crispercas_metalv_zcache.output,
        members=rules.crispercas_metalv_cluster.output.members,
        gsea=rules.crispercas_metalv_gsea_best.output,
        funnel=rules.crispercas_metalv_select.output.funnel,
        notebook=f"{PANELS_NB}/supp8.ipynb",
    output:
        png=f"{PANELS}/supp8/supp8.png",
        pdf=f"{PANELS}/supp8/supp8.pdf",
        svg=f"{PANELS}/supp8/supp8.svg",
        complete=touch(f"{PANELS}/supp8/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/supp8.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/supp8.ipynb"


# ============================================================
# Step 3: render ARCHS4 main figures
# ============================================================

rule fig3_panel:
    input:
        notebook=f"{PANELS_NB}/fig3.ipynb",
        helpers="scripts/panels/main_figure.R",
        coverage=rules.aggregate_bp_coverage.output.coverage_long,
        saturation=rules.aggregate_bp_saturation.output.saturation_long,
        traits=ancient(f"{A4_COV_BIO}/coverage_trait_recovery.csv"),
        finals=ancient(f"{A4_COV_BIO}/final_trait_recovery.csv"),
    output:
        png=f"{PANELS}/fig3/fig3.png",
        pdf=f"{PANELS}/fig3/fig3.pdf",
        svg=f"{PANELS}/fig3/fig3.svg",
        complete=touch(f"{PANELS}/fig3/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/fig3.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/fig3.ipynb"


rule fig4_panel:
    input:
        projections=rules.aggregate_projections.output.agg,
        drug_diseases=rules.prepare_three_model_comparison_data.output,
        scheme="scripts/panels/assets/projection_scheme.svg",
        mechanism_plots="scripts/archs4/projections/plots.R",
        notebook=f"{PANELS_NB}/fig4.ipynb",
        helpers="scripts/panels/main_figure.R",
    output:
        png=f"{PANELS}/fig4/fig4.png",
        pdf=f"{PANELS}/fig4/fig4.pdf",
        svg=f"{PANELS}/fig4/fig4.svg",
        complete=touch(f"{PANELS}/fig4/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/fig4.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/fig4.ipynb"


rule fig5_panel:
    input:
        data=rules.crispercas_fig5_data.output,
        helpers="scripts/panels/main_figure.R",
        notebook=f"{PANELS_NB}/fig5.ipynb",
    output:
        png=f"{PANELS}/fig5/fig5.png",
        pdf=f"{PANELS}/fig5/fig5.pdf",
        svg=f"{PANELS}/fig5/fig5.svg",
        complete=touch(f"{PANELS}/fig5/notebook.complete"),
    log:
        notebook=f"{PANELS_NB}/fig5.executed.ipynb",
    conda: "clamp-analyses"
    notebook:
        f"{PANELS_NB}/fig5.ipynb"


# ============================================================
# Step 4: build supplementary tables
# ============================================================

rule supplementary_tables:
    input:
        metalv=rules.crispercas_fig5_data.output,
        coverage_traits=PHENOPLIER_COV_BASE,
        projection_sources=rules.projection_benchmark_figure.output.sources,
        phenomexcan=[f"{A4_CRISPERCAS_CFG['phenomexcan_dir']}/{name}" for name in
                    ["phenotypes_info.tsv", "gene_assoc/smultixcan-mashr-pvalues.pkl",
                     "gene_assoc/smultixcan-mashr-zscores.pkl", "genes_metadata/genes_mapping_name_to_id.pkl"]],
        funnel=rules.crispercas_metalv_select.output.funnel,
        flags=rules.crispercas_noise_flags.output,
        ora=expand(rules.crispercas_metalv_ora.output.ora_dir, database=A4_METALV_ORA),
        lipid_deg=A4_CRISPERCAS_CFG["lipid_deg"],
        l1000=rules.crispercas_metalv_l1000.output.ranking,
        script="scripts/tables/build_supplementary_tables.py",
    output:
        expand(f"{PANELS}/tables/supplementary_table_{{n}}.xlsx", n=range(1, 9)),
    conda: "clamp-analyses"
    shell:
        "python {input.script:q} --phenomexcan-dir {A4_CRISPERCAS_CFG[phenomexcan_dir]:q}"


PUBLISHED_FIGURES = ["fig2", "fig3", "fig4", "fig5"] + [f"supp{i}" for i in range(1, 9)]


# ============================================================
# Step 5: publish figures and tables
# ============================================================

rule publish_figures:
    input:
        figures=expand(f"{PANELS}/{{name}}/{{name}}.{{ext}}", name=PUBLISHED_FIGURES, ext=["pdf", "png", "svg"]),
        tables=rules.supplementary_tables.output,
    output:
        figures=expand("figures/{name}/{name}.{ext}", name=PUBLISHED_FIGURES, ext=["pdf", "png", "svg"]),
        tables=expand("figures/tables/supplementary_table_{n}.xlsx", n=range(1, 9)),
    run:
        import shutil
        for src, dst in zip(list(input.figures) + list(input.tables), list(output.figures) + list(output.tables)):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)

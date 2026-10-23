# metaLV GLS run (CRISPR-Cas9 Figure 5)

PhenoPLIER GLS of the GSEA-significant CRISPR metaLVs against the PhenomeXcan
traits. The result is the declared input `crispercas.gls_summary` of
`workflow/rules/archs4_crispercas.smk`; `crispercas_metalv_select` recomputes the
Benjamini-Hochberg FDR over the GSEA-significant metaLVs of the current run.

1. Locally, after `crispercas_metalv_gsea_best`:

   ```bash
   python scripts/phenoplier/metalv/export_metalv_z.py \
     --metaz output/03_model_biology/02_archs4/04_crispercas/00_mgm/metaZ.csv \
     --gsea-best output/03_model_biology/02_archs4/04_crispercas/00_mgm/gsea/gsea_best.csv \
     --out-pkl metalv_z.pkl \
     --out-mapping output/03_model_biology/02_archs4/04_crispercas/00_mgm/phenoplier/lv_mapping.csv
   ```

   PhenoPLIER needs contiguous LV names, so the metaLVs are renamed `LV1..LVn`
   and the mapping is kept.

2. On pico, copy `metalv_z.pkl` and run (one node, local executor, as for the
   final models; step-6/7 pools are sized from the 32-CPU cgroup):

   ```bash
   PHENOPLIER_HOME=/pividori_lab/<user>/metalv_phenoplier/workspace \
     sbatch run_gls_metalv.sbatch metalv_z.pkl metalv
   python ../verify_summary.py <workspace>/projects/metalv/.../gls-summary.tsv.gz
   ```

3. Copy the GLS summary back and name the LVs:

   ```bash
   python scripts/phenoplier/metalv/collect_summary.py --summary gls-summary.tsv.gz \
     --mapping output/03_model_biology/02_archs4/04_crispercas/00_mgm/phenoplier/lv_mapping.csv \
     --out output/03_model_biology/02_archs4/04_crispercas/00_mgm/phenoplier/gls_summary_metalv.tsv.gz
   ```

The manuscript run (PhenoPLIER v0.5.2, biomedical trait filter, 2,366 traits,
`--lv-percentile 0.01`) used the 40 metaLVs significant when clusters were kept
at >= 5/10 seeds; the 27 metaLVs significant at >= 6/10 seeds are a subset with
the same clusters and loadings, and GLS is fitted per metaLV, so their P values
are unchanged.

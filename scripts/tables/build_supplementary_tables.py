#!/usr/bin/env python
# Build the eight manuscript supplementary workbooks.
from pathlib import Path

import argparse
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "output/99_panels/tables"
ARCHS4 = ROOT / "output/03_model_biology/02_archs4"
parser = argparse.ArgumentParser()
parser.add_argument("--phenomexcan-dir", type=Path, default=ROOT / "data/phenomexcan")
PHENOMEXCAN = parser.parse_args().phenomexcan_dir
COVERAGE_SUMMARIES = ROOT / "output/01_model_building/02_archs4/clampbase_coverage/summaries"
SUPP6_SOURCE = ROOT / "output/99_panels/supp6/source_data"
CRISPR = ARCHS4 / "04_crispercas"
FDR = 0.05


def mgm(x):

    if not isinstance(x, str):
        return x
    for a, b in [("meta-modules (metaLVs)", "meta-gene modules (MGMs)"), ("meta-modules", "meta-gene modules"),
                 ("metaLV", "MGM")]:
        x = x.replace(a, b)
    return x


LV_WORD = re.compile(r"(?<![A-Za-z])(?:LV|lv)(?=s?(?![A-Za-z]))")
LV_IDS = re.compile(r"(?:S\d+:)?LV\d+(?:\s*[,;/]\s*(?:S\d+:)?LV\d+)*")


def gm(x):
    return LV_WORD.sub("GM", mgm(x)) if isinstance(x, str) else x


def gm_ids(x):
    if isinstance(x, str) and LV_IDS.fullmatch(x.strip()):
        return LV_WORD.sub("GM", x)
    return mgm(x)


def mgm_frame(df):
    df = df.copy()
    for c in df.columns[df.dtypes.map(lambda t: isinstance(t, pd.CategoricalDtype))]:
        df[c] = df[c].astype(object)
    for c in df.columns[df.dtypes == object]:
        df[c] = df[c].map(gm_ids)
    return df.rename(columns=gm)


README_TEXT = [("CLAMP-full", "CLAMPfull"), ("normalised", "normalized"), ("a MGM", "an MGM"),
               ("Supplementary Figure ", "Supplementary Fig. "), ("Figure ", "Fig. ")]
COLUMN_NAMES = {"S123_GMs": "reference_model_GMs", "GeneRatio": "gene_ratio", "BgRatio": "bg_ratio",
                "FoldEnrichment": "fold_enrichment", "traits_FDR_0.05": "traits_FDR_0_05",
                "genes_SMultiXcan_FDR_0.05": "genes_SMultiXcan_FDR_0_05"}


def tidy_readme(x):
    x = gm(x)
    for a, b in README_TEXT:
        x = x.replace(a, b)
    return x


def tidy_columns(df):
    df = df.rename(columns=COLUMN_NAMES)
    return df.rename(columns=lambda c: c.replace("lipid-increasing", "lipid_increasing")
                     .replace("lipid-decreasing", "lipid_decreasing"))


def write_workbook(path, readme, sheets):

    sheets = {gm(k): tidy_columns(mgm_frame(v)) for k, v in sheets.items()}
    readme = [tidy_readme(r) for r in readme]
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        pd.DataFrame({"README": readme}).to_excel(xl, sheet_name="README", index=False, header=False)
        for name, df in sheets.items():
            df.to_excel(xl, sheet_name=name, index=False)
    print(f"{path.name}: " + ", ".join(f"{k}={len(v):,} rows" for k, v in sheets.items()))


def table_traits():
    traits = set()
    for f in sorted(COVERAGE_SUMMARIES.glob("cov_rs*_seed*.tsv.gz")):
        traits |= set(pd.read_csv(f, sep="\t", usecols=["phenotype"], dtype=str).phenotype)
    assert len(traits) == 2366, len(traits)

    info = pd.read_csv(PHENOMEXCAN / "phenotypes_info.tsv", sep="\t", dtype={"short_code": str})
    info = info[info.short_code.isin(traits)].copy()
    assert len(info) == len(traits)

    df = info.rename(columns={
        "short_code": "trait_id", "full_code": "phenomexcan_code", "description": "description",
        "type": "trait_type", "n": "n_total", "n_cases": "n_cases", "n_controls": "n_controls",
        "source": "source"})
    df = df[["trait_id", "phenomexcan_code", "description", "trait_type", "source",
             "n_total", "n_cases", "n_controls"]].sort_values("trait_id")
    readme = [
        "Supplementary Table 1. Traits evaluated in the genetic-trait recovery analyses.",
        "All 2,366 PhenomeXcan traits tested by the PhenoPLIER generalized least-squares (GLS) regression in the training-size, fixed-rank and full-compendium analyses (Figure 3). Per-model trait recovery counts are provided as Figure 3 source data.",
        "trait_id / phenomexcan_code / description / trait_type / source: PhenomeXcan trait annotation (phenotypes_info.tsv).",
        "n_total, n_cases, n_controls: GWAS sample sizes reported by PhenomeXcan (cases and controls are empty for continuous traits).",
    ]
    return readme, {"traits": df}


STUDIES = {
    "cytokines": ("IFN-alpha-treated EndoC-betaH1 beta cells (GSE133218)", "Colli et al. 2020",
                  "10.1038/s41467-020-16327-0"),
    "monocyte": ("LPS-treated human monocyte-derived macrophages (GSE193336)", "Li et al. 2022",
                 "10.1038/s41598-022-16005-9"),
    "placenta": ("Sorted cells from the first-trimester maternal-fetal interface (E-MTAB-6701)",
                 "Vento-Tormo et al. 2018", "10.1038/s41586-018-0698-6"),
}
NONE = "None reported (transcriptomic evidence only)"
EVIDENCE = {

    "Interferon alpha response": "Proteomics; western blot of STAT1/STAT2 phosphorylation and IRF1; siRNA knockdown of IRF1, STAT1 and STAT2",
    "IRF/STAT response": "Proteomics; western blot of STAT1/STAT2 phosphorylation and IRF1; siRNA knockdown of IRF1, STAT1 and STAT2",
    "MHC-I antigen presentation": "Proteomics; flow cytometry of MHC class I surface protein",
    "ER stress": "Proteomics",
    "Protein modification/degradation": "Proteomics",
    "Alternative splicing": "RT-PCR validation of splicing events (RNA level)",
    "PD-L1/HLA-E immune protection": "Flow cytometry of HLA-E surface protein",
    "Pancreatic beta cells": "Experimental design (beta-cell line)",


    "NAD+ depletion / SIRT1 inhibition": "Plasma metabolomics (nicotinamide, NAD+ salvage) and MALDI-MS imaging of NAD+ in paws of mice with collagen-induced arthritis",
    "Increased glucose transport / glycolytic reprogramming": "Western blot of c-Myc and the lactate/pyruvate transporter MCT1 in joints of mice with collagen-induced arthritis",
    "Reduced TCA cycle / OXPHOS": "Plasma metabolomics (itaconate, isocitrate, malate) in mice with collagen-induced arthritis",
    "Oxidative stress / ROS": "Plasma metabolomics (gamma-glutamyl amino acids) in mice with collagen-induced arthritis",
    "HIF1A / NF-κB inflammatory signaling": NONE,
    "Itaconate / macrophage metabolic rewiring": "Plasma metabolomics (itaconate) in mice with collagen-induced arthritis",
    "Glutathione / redox metabolism": "Plasma metabolomics (gamma-glutamyl cycle metabolites) in mice with collagen-induced arthritis",
    "Neutrophil activation": "Plasma CXCL1 protein (Luminex) in mice with collagen-induced arthritis",
    "Macrophage": "Experimental design (monocyte-derived macrophages)",
    "Neutrophil": NONE,


    "Trophoblast differentiation into EVT / SCT": NONE,
    "EVT invasion / EMT": "Immunohistochemistry of HLA-G+ EVT at the invasion site",
    "Spiral-artery remodelling": "Immunohistochemistry and multiplexed smFISH of perivascular cells",
    "dNK-trophoblast HLA-C / HLA-E / HLA-G interactions": "Flow cytometry of KIR2DL1 on dNK1 cells",
    "dNK chemokine / immunomodulatory signaling": "Flow cytometry of dNK granule proteins",
    "Immune checkpoint / maternal-fetal immune tolerance": "Immunohistochemistry of PDL1 and HLA-G on EVT",
    "Adenosine-mediated immunoregulation": NONE,
    "dNK1 glycolytic metabolic priming": NONE,
    "Trophoblast": "Immunohistochemistry (HLA-G)",
    "Extravillous trophoblast (EVT)": "Immunohistochemistry (HLA-G)",
    "Syncytiotrophoblast (SCT)": NONE,
    "Villous cytotrophoblast (VCT)": NONE,
    "Decidual NK / dNK": "Flow cytometry of three dNK subsets",
    "Decidual stromal cells": "Immunohistochemistry and multiplexed smFISH of decidual layers",
    "Decidual macrophages / maternal myeloid cells": "Flow cytometry sorting",
}


def table_projection_targets():
    rows = []
    for group, f in [("cytokines", "00_cytokines_manual_recovery.csv"),
                     ("monocyte", "01_monocyte_manual_recovery.csv"),
                     ("placenta", "02_placenta_manual_recovery.csv")]:
        d = pd.read_csv(SUPP6_SOURCE / f)
        for target, g in d.groupby("mechanism", sort=False):
            a = g[g.model == "ARCHS4"].iloc[0]
            l = g[g.model == "local"].iloc[0]
            dataset, study, doi = STUDIES[group]
            rows.append({
                "dataset": dataset, "source_study": study, "doi": doi, "target": target,
                "target_type": "Cell type" if a.category == "cell_type" else "Pathway or process",
                "non_transcriptomic_evidence_in_source_study": EVIDENCE[target],
                "ARCHS4_projection_recovered": bool(a.recovered),
                "ARCHS4_LV": a.LV, "ARCHS4_gene_set": a.term, "ARCHS4_ORA_FDR": a.fdr,
                "dataset_trained_recovered": bool(l.recovered),
                "dataset_trained_LV": l.LV, "dataset_trained_gene_set": l.term,
                "dataset_trained_ORA_FDR": l.fdr,
            })
    df = pd.DataFrame(rows)
    assert len(df) == 33, len(df)
    assert df.ARCHS4_projection_recovered.sum() == 27
    assert df.dataset_trained_recovered.sum() == 12
    readme = [
        "Supplementary Table 2. Curated pathway, process and cell-type targets for the external-dataset projection comparison (Figure 4b and Supplementary Figure 6c,d).",
        "Each target is a pathway, cellular process or cell type reported by the source study. A target was recovered when a module's top 1% of gene loadings was enriched (ORA FDR < 0.05) for a gene set matched to it; LV and gene set give the best match, blank when not recovered.",
        "ARCHS4_*: projection into the frozen ARCHS4 CLAMPfull model. dataset_trained_*: CLAMPfull model trained on the evaluated dataset.",
        "non_transcriptomic_evidence_in_source_study: measurements other than RNA-seq that the source study reported for the target, curated from its full text. 'Experimental design' marks cell types defined by the cells profiled.",
        "Totals: 33 targets; 27 recovered by ARCHS4 projection and 12 by the dataset-trained models.",
    ]
    return readme, {"targets": df}


METALV = CRISPR / "00_mgm"
FIG5_LVS = ["metaLV175", "metaLV335", "metaLV588", "metaLV91"]
SEEDS = ["S123", "S1001", "S1002", "S1003", "S1004", "S1005", "S1006", "S1007", "S456", "S789"]
GENE_SET = {"gene_set_increase": "lipid-increasing", "gene_set_decrease": "lipid-decreasing"}
CELL_TYPE = {"metaLV175": "Hepatocyte", "metaLV335": "Macrophage",
             "metaLV588": "Adrenocortical cells, oligodendrocyte", "metaLV91": "Adrenocortical cells, oligodendrocyte"}
ORA_LIBRARIES = {"canonical": "Canonical pathways (MSigDB C2 CP, no Reactome)", "bp": "GO Biological Process",
                 "cellmarker": "CellMarker 2.0", "azimuth": "Azimuth 2023", "archs4_cell_lines": "ARCHS4 cell lines"}
TOP_N = 185


def sci(x):
    return f"{x:.1e}"


def load_ora():
    with open(ROOT / "data/pathways/go_bp.Hs.symbols.gmt") as fh:
        go = dict(line.split("\t", 2)[:2] for line in fh)
    ora = {}
    for db in ORA_LIBRARIES:
        e = pd.read_csv(METALV / f"ora/{db}/enrichment.csv.gz")
        if db == "bp":
            e["Description"] = e.ID.map(go).fillna(e.Description)
        ora[db] = e[e["p.adjust"] < FDR]
    return ora


def table_metalv_selection():
    summ = pd.read_csv(METALV / "cluster_summary.csv")
    funnel = pd.read_csv(METALV / "selection/funnel.csv")
    gsea = pd.read_csv(METALV / "gsea/gsea_best.csv")
    assert len(summ) == len(funnel) == 618
    t = summ.drop(columns="cluster").rename(columns={"n_models": "n_seeds", "n_lvs": "n_member_LVs",
                                                     "min_cos": "min_cosine", "mean_cos": "mean_cosine"})
    t["members"] = t.members.str.replace(";", ", ")
    t["S123_lvs"] = t.S123_lvs.fillna("").str.replace(";", ", ")
    for k, label in GENE_SET.items():
        g = gsea[gsea.pathway == k].set_index("lv")
        t[f"GSEA_{label}_NES"] = t.metaLV.map(g.NES)
        t[f"GSEA_{label}_padj"] = t.metaLV.map(g.padj)
        t[f"GSEA_{label}_leading_edge"] = t.metaLV.map(g.leadingEdge)
    z_all = pd.read_csv(METALV / "metaZ.csv", index_col=0)
    screen = pd.read_csv(ROOT / "data/archs4/crispr_cas9/lipid_DEG.csv").set_index("gene_name")["lipid effect"]

    def screen_top(m):
        z = z_all[m][z_all[m] > 0]
        r = z.rank(ascending=False, method="min")
        r = r[(r <= TOP_N) & r.index.isin(screen.index)].sort_values()
        return ", ".join(f"{g} (rank {int(k)}, {screen[g]})" for g, k in r.items())

    t["CRISPR_screen_genes_top1pct"] = t.metaLV.map(screen_top)
    steps = ["pass_gsea", "pass_trait", "pass_canonical", "pass_screen_gene_top", "pass_missed_by_twas", "pass_lipid_pathway"]
    t = t.merge(funnel[["metaLV", *steps, "n_lipid_pathways"]], on="metaLV")
    t["figure5_module"] = t.metaLV.isin(FIG5_LVS)
    assert set(t[t.pass_lipid_pathway].metaLV) == set(FIG5_LVS)
    counts = pd.DataFrame({"step": ["metaLVs (>= 6/10 seeds)", *[s.replace("pass_", "") for s in steps]],
                           "n_metaLVs": [len(t), *[int(t[s].sum()) for s in steps]]})

    z = z_all[FIG5_LVS]
    top = pd.read_csv(METALV / "pathways/canonical_bp_top.csv")
    pairs = pd.read_csv(METALV / "selection/screen_gene_trait_pairs.csv")
    comp = pd.read_csv(METALV / "projection/composition.csv")
    comp = comp[(comp.set == "top 1%") & (comp.category != "other")].sort_values("pct", ascending=False)
    traits = pd.read_csv(METALV / "fig5/module_traits_all.csv")
    ring = pd.read_csv(METALV / "fig5/ring_pvalues.csv")
    ring = ring[pd.to_numeric(ring.fdr) < FDR]
    ora = load_ora()

    def top_hit(db, m):
        x = ora[db][ora[db].LV == m].nsmallest(1, "p.adjust")
        return f"{x.Description.iloc[0]} padj {sci(x['p.adjust'].iloc[0])}" if len(x) else ""

    rows = []
    for m in FIG5_LVS:
        s = t.set_index("metaLV").loc[m]
        hit = pairs[(pairs.metaLV == m) & pairs.in_top & pairs.missed_by_twas]
        sig = gsea[(gsea.lv == m) & (gsea.padj < FDR)].sort_values("padj")
        p3 = top[top.metaLV == m].head(3)
        tr = traits[traits.metaLV == m].sort_values("pvalue")
        rg = ring[ring.metaLV == m]
        c = comp[comp.metaLV == m]
        rows.append({
            "metaLV": m, "n_seeds": f"{s.n_seeds}/10", "n_member_LVs": s.n_member_LVs, "member_LVs": s.members,
            "min_cosine": round(s.min_cosine, 2), "mean_cosine": round(s.mean_cosine, 2),
            "CRISPR_gene_set": ", ".join(GENE_SET[x] for x in sig.pathway),
            "CRISPR_GSEA_NES": ", ".join(f"{x:.2f}" for x in sig.NES),
            "CRISPR_GSEA_padj": ", ".join(f"{x:.4g}" for x in sig.padj),
            "CRISPR_GSEA_leading_edge": "; ".join(sig.leadingEdge),
            "CRISPR_screen_genes": ", ".join(dict.fromkeys(hit.gene)),
            "CRISPR_screen_genes_top1pct": s.CRISPR_screen_genes_top1pct,
            "top10_genes": ", ".join(z[m].nlargest(10).index),
            "top3_pathways": "; ".join(p3.pathway), "top3_pathways_padj": ", ".join(sci(x) for x in p3.padj),
            "CellMarker_top": top_hit("cellmarker", m), "Azimuth_top": top_hit("azimuth", m), "cell_type": CELL_TYPE[m],
            "ARCHS4_top1pct_samples": ", ".join(f"{a} {v:.1f}%" for a, v in zip(c.category, c.pct)),
            "traits_FDR_0.05": "; ".join(f"{g}: " + ", ".join(f"{d} FDR {sci(q)}" for d, q in zip(x.phenotype_desc, x.fdr))
                                         for g, x in tr.groupby("trait_group", sort=False)),
            "genes_SMultiXcan_FDR_0.05": ", ".join(f"{a} {b}" for a, b in zip(rg.gene, rg.disease)),
        })
    readme = [
        "Supplementary Table 4. Lipid-screen meta-modules (metaLVs) across ten CLAMP-full canonical ARCHS4 models and their selection.",
        "metaLVs: loading vectors of the ten models (seed 123, the reference model used in the other analyses, and seeds 456, 789 and 1001-1007), after removing LVs flagged by the technical-noise screen (Supplementary Table 3), were L2-normalised and clustered by complete-linkage agglomerative clustering on cosine distance, cut at cosine similarity 0.5. Clusters with LVs from at least 6 of the 10 seeds are metaLVs; a metaLV's loadings are the gene-wise median of its members' loadings.",
        "GSEA columns: fgsea (scoreType = 'pos', ten repetitions) of each metaLV's loadings against the HepG2 lipid-screen high-confidence sets; the smallest adjusted P across repetitions and the leading edge of that repetition.",
        "CRISPR_screen_genes_top1pct: every HepG2 lipid-screen gene (any confidence tier) among the metaLV's top 1% of loadings (185 genes), whether or not it is in a GSEA leading edge, ordered by loading rank, with the screen's lipid effect (increase/decrease; *** high-confidence sets used for GSEA, ** and * lower tiers).",
        "Selection flags (cumulative): pass_gsea, adjusted P < 0.05 for either set; pass_trait, at least one PhenoPLIER GLS trait association at FDR < 0.05 (Supplementary Table 5); pass_canonical, at least one canonical pathway at ORA FDR < 0.05 (Supplementary Table 6); pass_screen_gene_top, a leading-edge screen gene among the metaLV's top 1% of loadings (185 genes); pass_missed_by_twas, that gene has single-gene S-MultiXcan P >= 0.05 for an associated trait; pass_lipid_pathway, at least one lipid-related canonical pathway (keywords: lipid, cholesterol, sterol, fatty, lipoprotein, triglyceride, eicosanoid, phospholipid, bile, prostaglandin, leukotriene, arachidonic, sphingo, PPAR, LDL, HDL). Copy-number amplicon programs were not excluded.",
        "figure5_metaLVs sheet: summary of the four metaLVs shown in Figure 5. CRISPR_GSEA_leading_edge: leading edge of each significant screen set (same order as CRISPR_gene_set); CRISPR_screen_genes: leading-edge genes in the top 1% that pass the S-MultiXcan filter; CRISPR_screen_genes_top1pct: as in the metaLVs sheet. CellMarker_top and Azimuth_top: most significant ORA term; cell_type: label used in Figure 5. ARCHS4_top1pct_samples: sample types among the 1% of ARCHS4 samples with the highest projected metaLV activity (B), excluding 'other'. traits_FDR_0.05: traits with at least 1,000 cases grouped by disease. genes_SMultiXcan_FDR_0.05: displayed genes with S-MultiXcan FDR < 0.05 for a linked disease (circled in Figure 5a).",
    ]
    return readme, {"metaLVs": t, "selection_counts": counts, "figure5_metaLVs": pd.DataFrame(rows)}


def table_technical_noise():
    flags = pd.read_csv(ARCHS4 / "06_technical_noise/seeds/all_flags.csv")
    assert flags.model.nunique() == 10 and len(flags) == 17280
    flags["technical_joint"] = flags.technical_joint.round(4)
    flags = flags.rename(columns={"model": "seed_model", "technical_joint": "technical_joint_R2",
                                  "any_flag": "technically_noisy"})
    flags["seed_model"] = pd.Categorical(flags.seed_model, SEEDS)
    flags["LV_num"] = flags.LV.str[2:].astype(int)
    flags = flags.sort_values(["seed_model", "LV_num"]).drop(columns="LV_num")
    counts = (flags.groupby("seed_model", observed=True)
              .agg(n_LVs=("LV", "size"), gene_flag=("gene_flag", "sum"), sample_flag=("sample_flag", "sum"),
                   technically_noisy=("technically_noisy", "sum")).reset_index())
    assert counts.set_index("seed_model").at["S123", "technically_noisy"] == 116
    readme = [
        f"Supplementary Table 3. Technical-noise screen of the ten CLAMP-full canonical ARCHS4 models: {int(flags.technically_noisy.sum()):,} of {len(flags):,} LVs flagged and removed before metaLV clustering.",
        "gene_flag: gene loadings correlated with log10 cDNA length or transcript GC content (Ensembl 107), over all genes or over non-zero loadings, with |r| >= 0.2, BH FDR across the model's LVs < 0.2 and |robust z| >= 4.",
        "technical_joint_R2: cross-validated (5-fold) fraction of variance of the LV's sample scores (B) across ARCHS4 samples explained jointly by technical covariates (platform, instrument model, library selection, library source, molecule, submission year, log10 library size, log10 detected genes, mitochondrial, ribosomal and protein-coding read fractions, single-cell probability). sample_flag: technical_joint_R2 >= 0.25.",
        "technically_noisy: gene_flag or sample_flag. S123 is the reference model used in the other analyses (116 LVs flagged).",
    ]
    return readme, {"per_seed_counts": counts, "all_LVs": flags}


def table_metalv_ora():
    rows = []
    for db, label in ORA_LIBRARIES.items():
        e = pd.read_csv(METALV / f"ora/{db}/enrichment.csv.gz")
        e = e[e.LV.isin(FIG5_LVS) & (e["p.adjust"] < FDR)].copy()
        if db == "bp":
            with open(ROOT / "data/pathways/go_bp.Hs.symbols.gmt") as fh:
                go = dict(line.split("\t", 2)[:2] for line in fh)
            e["Description"] = e.ID.map(go).fillna(e.Description)
        e["library"] = label
        rows.append(e)
    o = pd.concat(rows)
    top = pd.read_csv(METALV / "pathways/canonical_bp_top.csv")
    shown = {(m, p) for m, x in top.groupby("metaLV") for p in x.pathway.head(3)}
    o["top3_pathway_in_figure5"] = [(m, d) in shown for m, d in zip(o.LV, o.Description)]
    o["genes"] = o.geneID.str.replace("/", ", ")
    o = o.rename(columns={"LV": "metaLV", "Description": "term", "ID": "term_id", "p.adjust": "padj",
                          "Count": "n_overlap"})
    o = o[["metaLV", "library", "term", "term_id", "GeneRatio", "BgRatio", "FoldEnrichment", "n_overlap", "pvalue",
           "padj", "top3_pathway_in_figure5", "genes"]]
    o["metaLV"] = pd.Categorical(o.metaLV, FIG5_LVS)
    o["library"] = pd.Categorical(o.library, list(ORA_LIBRARIES.values()))
    o = o.sort_values(["metaLV", "library", "padj"])
    readme = [
        "Supplementary Table 6. Over-representation analysis for the four Figure 5 metaLVs.",
        f"Query: each metaLV's top 1% of gene loadings ({TOP_N} of 18,423 genes); background: all model genes. clusterProfiler enricher with Benjamini-Hochberg FDR within metaLV and library; all terms with FDR < 0.05 are listed.",
        "Libraries: MSigDB C2 canonical pathways without Reactome (10-500 genes), GO Biological Process (10-500 genes), CellMarker 2.0 human markers pooled by cell name (10-500 genes), Azimuth 2023 (Enrichr; 5-3,000 genes) and ARCHS4 cell lines (Enrichr; 10-3,000 genes).",
        "top3_pathway_in_figure5: the three pathways shown in Figure 5 (canonical and GO:BP terms with FDR < 0.05 and at most 100 genes, merged by Jaccard similarity >= 0.4 and ranked by FDR).",
    ]
    return readme, {"ora": o}


def table_metalv_gls():
    g = pd.read_csv(METALV / "selection/gls_metalv.csv", dtype={"phenotype": str})
    g = g[g.metaLV.isin(FIG5_LVS)].copy()
    info = pd.read_csv(PHENOMEXCAN / "phenotypes_info.tsv", sep="\t", dtype={"short_code": str}).set_index("short_code")
    g["variable_type"] = g.phenotype.map(info.type)
    g["n_samples"] = g.phenotype.map(info.n)
    g["n_cases"] = g.phenotype.map(info.n_cases)
    g["n_controls"] = g.phenotype.map(info.n_controls)
    shown = pd.read_csv(METALV / "fig5/module_traits_all.csv", dtype={"phenotype": str})
    grp = shown.drop_duplicates("phenotype").set_index("phenotype").trait_group
    g["disease_group"] = g.phenotype.map(grp)
    g["significant"] = g.fdr < FDR
    g["at_least_1000_cases"] = g.n_cases.isna() | (g.n_cases >= 1000)
    g["in_figure5"] = [(m, p) in set(zip(shown.metaLV, shown.phenotype)) for m, p in zip(g.metaLV, g.phenotype)]
    g = g.rename(columns={"phenotype": "trait_code", "phenotype_desc": "trait", "pvalue": "gls_pvalue", "fdr": "gls_fdr"})
    g = g[["metaLV", "trait_code", "trait", "variable_type", "n_samples", "n_cases", "n_controls", "gls_pvalue",
           "gls_fdr", "significant", "at_least_1000_cases", "disease_group", "in_figure5"]]
    g["metaLV"] = pd.Categorical(g.metaLV, FIG5_LVS)
    g = g.sort_values(["metaLV", "gls_pvalue"])
    assert len(g) == 4 * 2366
    readme = [
        "Supplementary Table 5. PhenoPLIER GLS trait associations of the four Figure 5 metaLVs across 2,366 PhenomeXcan traits.",
        "GLS regression of each trait's S-MultiXcan gene associations on the metaLV's top 1% of gene loadings (binarized), accounting for gene-gene correlations and gene size and density covariates (PhenoPLIER v0.5.2, biomedical trait filter). gls_fdr: Benjamini-Hochberg FDR across all traits and the 27 GSEA-significant metaLVs.",
        "n_samples, n_cases, n_controls: UK Biobank sample sizes (PhenomeXcan); continuous traits have no cases or controls. at_least_1000_cases: traits used for Figure 5 (binary traits with at least 1,000 cases, or continuous traits). disease_group / in_figure5: disease to which the trait was merged in Figure 5.",
    ]
    return readme, {"gls": g}


def table_metalv_smultixcan():
    pv = pd.read_pickle(PHENOMEXCAN / "gene_assoc/smultixcan-mashr-pvalues.pkl")
    zs = pd.read_pickle(PHENOMEXCAN / "gene_assoc/smultixcan-mashr-zscores.pkl")
    name = pd.read_pickle(PHENOMEXCAN / "genes_metadata/genes_mapping_name_to_id.pkl")
    col = {c.split("-")[0]: c for c in pv.columns}
    zcol = {c.split("-")[0]: c for c in zs.columns}
    z = pd.read_csv(METALV / "metaZ.csv", index_col=0)[FIG5_LVS]
    shown = pd.read_csv(METALV / "fig5/module_genes.csv")
    traits = pd.read_csv(METALV / "fig5/module_traits_all.csv", dtype={"phenotype": str})
    fdr = {}

    def bh(p):
        s = p.dropna().sort_values()
        n = len(s)
        return pd.Series(np.minimum.accumulate((s * n / np.arange(1, n + 1))[::-1])[::-1].clip(upper=1).values, index=s.index)

    rows = []
    for m in FIG5_LVS:
        rank = z[m].rank(ascending=False, method="min").astype(int)
        screen = list(shown[(shown.metaLV == m) & shown.crispr].gene)
        genes = list(dict.fromkeys(screen + list(z[m].nlargest(TOP_N).index)))
        displayed = set(shown[shown.metaLV == m].gene)
        for t in traits[traits.metaLV == m].itertuples():
            c = col[t.phenotype]
            if c not in fdr:
                fdr[c] = bh(pv[c])
            for g in genes:
                gid = name.get(g)
                rows.append({"metaLV": m, "gene": g, "rank_in_metaLV": int(rank[g]), "loading": round(float(z.at[g, m]), 4),
                             "screen_gene": g in screen, "in_figure5": g in displayed, "disease_group": t.trait_group,
                             "trait_code": t.phenotype, "trait": t.phenotype_desc,
                             "smultixcan_z": zs.at[gid, zcol[t.phenotype]] if gid in zs.index else np.nan,
                             "smultixcan_pvalue": pv.at[gid, c] if gid in pv.index else np.nan,
                             "smultixcan_fdr": fdr[c].get(gid, np.nan)})
    s = pd.DataFrame(rows)
    readme = [
        "Supplementary Table 8. Single-gene S-MultiXcan associations (PhenomeXcan) of the four Figure 5 metaLVs' genes with their Figure 5 traits.",
        f"Genes: each metaLV's top 1% of gene loadings ({TOP_N} genes) and its screen genes. Traits: the metaLV's traits shown in Figure 5 (GLS FDR < 0.05, at least 1,000 cases).",
        "smultixcan_fdr: Benjamini-Hochberg FDR across all genes for that trait. in_figure5: gene drawn in Figure 5a.",
    ]
    return readme, {"smultixcan": s}


L1000_DRUG = {"metaLV91": "Clomifene", "metaLV588": "Clomifene", "metaLV175": "Cortisone acetate", "metaLV335": "Calcitriol"}


def table_metalv_l1000():
    r = pd.read_csv(METALV / "l1000/l1000_fig5_ranking.csv")
    cols = ["rank", "drugbank_id", "drug", "atc_codes", "score", "mean_z", "perm_z", "perm_p", "fdr", "n_genes"]
    sheets = {m: r[r.metaLV == m][cols] for m in ["metaLV91", "metaLV588", "metaLV175", "metaLV335"]}
    assert all(len(d) == 1170 for d in sheets.values())
    picks = pd.concat([r[(r.metaLV == m) & (r.drug == d)] for m, d in L1000_DRUG.items()])
    picks = picks.assign(n_drugs_fdr_0_05=picks.metaLV.map(r[r.fdr < FDR].groupby("metaLV").size()))
    readme = [
        "Supplementary Table 7. LINCS L1000 drug signatures scored on the four Figure 5 metaLVs.",
        f"Drugs: the 1,170 DrugBank-mapped LINCS L1000 consensus signatures used for drug-disease prediction. Genes: each metaLV's top 1% of gene loadings ({TOP_N} of 18,423) that are measured in LINCS (n_genes), mapped one-to-one from symbols to Ensembl identifiers.",
        "mean_z: mean LINCS z-score of the drug over these genes. score: mean_z standardized per drug across the 27 GSEA-significant metaLVs (the PhenoPLIER set), so drugs that change many programs are not favoured. rank: descending score (program induced).",
        "perm_p: one-sided permutation P from 10,000 random sets of the same number of LINCS genes, scored in the same way; perm_z: (score - null mean) / null SD; fdr: Benjamini-Hochberg across the 1,170 drugs within each metaLV. Random gene sets ignore co-expression among module genes, so perm_p and fdr are lenient; rank and perm_z are the more informative quantities.",
        "figure5_drugs sheet: the drug shown for each metaLV in Figure 5, the highest-ranked drug with an established mechanism on the program (clomifene, sterol-synthesis inhibitor; cortisone acetate, glucocorticoid; calcitriol, vitamin D); n_drugs_fdr_0_05: drugs at FDR < 0.05 for that metaLV.",
    ]
    return readme, {**sheets, "figure5_drugs": picks[["metaLV", *cols, "n_drugs_fdr_0_05"]]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    for n, builder in enumerate([table_traits, table_projection_targets, table_technical_noise, table_metalv_selection,
                                 table_metalv_gls, table_metalv_ora, table_metalv_l1000, table_metalv_smultixcan], start=1):
        readme, sheets = builder()
        write_workbook(OUT / f"supplementary_table_{n}.xlsx", readme, sheets)


if __name__ == "__main__":
    main()

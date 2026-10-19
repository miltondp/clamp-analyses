#!/usr/bin/env python
# Summarize screen enrichment and select interpretable metaLVs.
import argparse
import sys
import pandas as pd
from pathlib import Path
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


def gsea_best_main(args):
    ap = argparse.ArgumentParser()
    ap.add_argument("--fgsea", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)

    fg = pd.read_csv(a.fgsea)
    best = fg.groupby(["lv", "pathway"]).agg(pval=("pval", "min"), padj=("padj", "min"), NES=("NES", "max")).reset_index()
    le = fg.sort_values("padj").drop_duplicates(["lv", "pathway"])[["lv", "pathway", "leadingEdge"]]
    best = best.merge(le, on=["lv", "pathway"])
    best.to_csv(a.out, index=False)
    print(f"metaLVs: {best.lv.nunique()}; significant (best padj < 0.05): {best[best.padj < 0.05].lv.nunique()}")


def select_main(args):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cluster-summary", required=True)
    ap.add_argument("--metaz", required=True)
    ap.add_argument("--gsea-best", required=True)
    ap.add_argument("--gls", required=True)
    ap.add_argument("--canonical-ora", required=True)
    ap.add_argument("--smultixcan-zscores", required=True)
    ap.add_argument("--symbol-ensembl", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--top-n", type=int, default=185)
    ap.add_argument("--lipid-regex", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(args)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    summ = pd.read_csv(a.cluster_summary).set_index("metaLV")
    z = pd.read_csv(a.metaz, index_col=0)
    gsea = pd.read_csv(a.gsea_best)
    sig_sets = gsea[gsea.padj < a.fdr]
    gsea_sig = list(dict.fromkeys(sig_sets.sort_values("padj").lv))
    le = {m: list(dict.fromkeys(g.strip() for s in x.leadingEdge for g in str(s).split(",")))
          for m, x in sig_sets.groupby("lv")}

    gls = pd.read_csv(a.gls, sep="\t", dtype={"phenotype": str})
    gls = gls[gls.metaLV.isin(gsea_sig)].copy()
    assert gls.metaLV.nunique() == len(gsea_sig), "GLS summary is missing GSEA-significant metaLVs"
    gls["fdr"] = multipletests(gls.pvalue, method="fdr_bh")[1]
    gls.to_csv(out / "gls_metalv.csv", index=False)
    traits = gls[gls.fdr < a.fdr]

    can = pd.read_csv(a.canonical_ora)
    can = can[can["p.adjust"] < a.fdr]
    tw = pd.read_pickle(a.smultixcan_zscores)
    col = {c.split("-")[0]: c for c in tw.columns}
    sym = pd.read_csv(a.symbol_ensembl).dropna().drop_duplicates("SYMBOL").set_index("SYMBOL").ENSEMBL

    pairs = []
    for m in gsea_sig:
        rank = z[m].rank(ascending=False, method="min")
        for t in traits[traits.metaLV == m].itertuples():
            c = col.get(t.phenotype)
            for g in le[m]:
                e = sym.get(g)
                p = float(2 * norm.sf(abs(tw.at[e, c]))) if (c is not None and e in tw.index) else np.nan
                pairs.append({"metaLV": m, "gene": g, "gene_rank": int(rank[g]), "trait_code": t.phenotype,
                              "trait": t.phenotype_desc, "trait_fdr": t.fdr, "twas_p": p})
    pairs = pd.DataFrame(pairs)
    pairs["in_top"] = pairs.gene_rank <= a.top_n
    pairs["missed_by_twas"] = pairs.twas_p >= 0.05
    pairs.to_csv(out / "screen_gene_trait_pairs.csv", index=False)

    lipid = can[can.Description.str.contains(a.lipid_regex, case=False, regex=True)]
    f = pd.DataFrame({"metaLV": summ.index})
    f["n_models"] = f.metaLV.map(summ.n_models)
    f["n_lvs"] = f.metaLV.map(summ.n_lvs)
    best = gsea.groupby("lv").padj.min()
    f["gsea_padj"] = f.metaLV.map(best)
    f["crispr_leading_edge"] = f.metaLV.map(lambda m: ", ".join(le.get(m, [])))
    f["pass_gsea"] = f.metaLV.isin(gsea_sig)
    f["pass_trait"] = f.pass_gsea & f.metaLV.isin(traits.metaLV)
    f["pass_canonical"] = f.pass_trait & f.metaLV.isin(can.LV)
    f["pass_screen_gene_top"] = f.pass_canonical & f.metaLV.isin(pairs[pairs.in_top].metaLV)
    f["pass_missed_by_twas"] = f.pass_screen_gene_top & f.metaLV.isin(pairs[pairs.in_top & pairs.missed_by_twas].metaLV)
    f["pass_lipid_pathway"] = f.pass_missed_by_twas & f.metaLV.isin(lipid.LV)
    f["n_lipid_pathways"] = f.metaLV.map(lipid.groupby("LV").size()).fillna(0).astype(int)
    f.to_csv(out / "funnel.csv", index=False)

    steps = ["pass_gsea", "pass_trait", "pass_canonical", "pass_screen_gene_top", "pass_missed_by_twas", "pass_lipid_pathway"]
    print("metaLVs", len(f), "|", " -> ".join(f"{s.replace('pass_', '')} {int(f[s].sum())}" for s in steps))
    sel = f[f.pass_lipid_pathway][["metaLV", "n_models"]]
    sel.to_csv(out / "selected_metalvs.csv", index=False)
    print("selected:", ", ".join(sel.metaLV))


def main():
    commands = {
        "gsea-best": gsea_best_main,
        "select": select_main,
    }
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=commands)
    args = parser.parse_args(sys.argv[1:2])
    commands[args.command](sys.argv[2:])


if __name__ == "__main__":
    main()

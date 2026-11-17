#!/usr/bin/env python
# Summarize metaLV activity and prepare Figure 5 tables.
import argparse
import sys
import pandas as pd
import itertools
import shutil
from pathlib import Path
import numpy as np


MODULES = {
    "metaLV120": {"crispr": ["DGAT2"], "cells": "Hepatocyte",
                  "pathways": ["Complement and coagulation", "Acylglycerol homeostasis", "Bile acid metabolism"]},
    "metaLV281": {"crispr": ["DGAT2"], "cells": "Macrophage",
                  "pathways": ["Macrophage activation", "Acute inflammatory response", "Amyloid-beta clearance"]},
    "metaLV70": {"crispr": ["ACACA", "PCYT2"], "cells": "Adrenocortical, oligodendrocyte",
                 "pathways": ["Cholesterol biosynthesis", "SREBP signaling", "Omega-9 fatty acid synthesis"]},
    "metaLV270": {"crispr": ["TCF7L2"], "cells": "Primitive endoderm",
                  "pathways": ["Endoderm development", "Pancreas development", "FOXA2 (HNF3B) pathway"]},
}


TRAITS = {
    "20002_1162": ("Gallstones", "Gallstones"),
    "K80": ("Cholelithiasis (K80)", "Gallstones"),
    "K81": ("Cholecystitis (K81)", "Gallstones"),
    "20110_10": ("Mother AD", "Alzheimer's"),
    "20107_10": ("Father AD", "Alzheimer's"),
    "20111_10": ("Siblings AD", "Alzheimer's"),
    "20002_1473": ("High cholesterol", "High cholesterol"),
    "I80": ("Phlebitis (I80)", "Venous thrombosis"),
    "6152_5": ("DVT (leg)", "Venous thrombosis"),
    "2443": ("Diabetes (doctor)", "Diabetes"),
    "20002_1220": ("Diabetes (self-rep.)", "Diabetes"),
    "20110_9": ("Mother diabetes", "Diabetes"),
    "20107_9": ("Father diabetes", "Diabetes"),
}


GROUP = ({k: v[1] for k, v in TRAITS.items()}
         | {k: "Venous thrombosis" for k in ["I9_PHLETHROMBDVTLOW", "20002_1094", "I9_VTE", "I9_DVTANDPULM", "20002_1093", "6152_7"]}
         | {"K11_GALLBILPANC": "Gallstones"})


SCREEN = ["DGAT2", "ACACA", "PCYT2", "TCF7L2"]




SHORT = ["Non-cancer illness code, self-reported: ", "Diagnoses - main ICD10: ", "Illnesses of ",
         "Blood clot, DVT, bronchitis, emphysema, asthma, rhinitis, eczema, allergy diagnosed by doctor: "]


def bh(p):
    s = p.dropna().sort_values()
    n = len(s)
    return pd.Series(np.minimum.accumulate((s * n / np.arange(1, n + 1))[::-1])[::-1].clip(upper=1).values, index=s.index)


def figure_main(args):
    ap = argparse.ArgumentParser()
    ap.add_argument("--metaz", required=True)
    ap.add_argument("--cluster-summary", required=True)
    ap.add_argument("--selected", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--pathway-members", required=True)
    ap.add_argument("--pathway-top", required=True)
    ap.add_argument("--gls", required=True)
    ap.add_argument("--phenotypes-info", required=True)
    ap.add_argument("--phenomexcan", required=True)
    ap.add_argument("--composition", required=True)
    ap.add_argument("--l1000", required=True)
    ap.add_argument("--min-cases", type=int, default=1000)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--top-n", type=int, default=185)
    ap.add_argument("--reference-genes", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(args)
    reference = [g for g in a.reference_genes.split(",") if g]
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    selected = list(pd.read_csv(a.selected).metaLV)
    assert set(selected) == set(MODULES), f"selected metaLVs {selected} differ from the Figure 5 modules"
    pairs = pd.read_csv(a.pairs)
    hit = pairs[pairs.in_top]
    for m, d in MODULES.items():
        assert set(hit[hit.metaLV == m].gene) <= set(d["crispr"]), m
    summ = pd.read_csv(a.cluster_summary).set_index("metaLV")
    z = pd.read_csv(a.metaz, index_col=0)
    rank = z.rank(ascending=False, method="min").astype(int)
    for m, d in MODULES.items():
        assert all(rank.at[g, m] <= a.top_n and z.at[g, m] > 0 for g in d["crispr"]), m
    mem = pd.read_csv(a.pathway_members)
    l1000 = pd.read_csv(a.l1000)
    assert set(l1000.metaLV) <= set(selected) and l1000.fdr.lt(a.fdr).all()
    assert l1000.direction_agrees_with_null.all()
    l1000["label"] = np.where(l1000.direction.eq("Activation-like"), "↑ ", "↓ ") + l1000.drug
    l1000.to_csv(out / "module_drugs.csv", index=False)
    top = pd.read_csv(a.pathway_top)

    mods, genes, edges = [], [], []
    for m, d in MODULES.items():
        nodes = d["crispr"] + list(z[m].nlargest(10).index)
        for g in nodes:
            genes.append({"metaLV": m, "gene": g, "rank": int(rank.at[g, m]), "crispr": g in d["crispr"]})
        sets = [set(s.split("/")) for s in mem[mem.metaLV == m].geneID]
        for g1, g2 in itertools.combinations(nodes, 2):
            if any(g1 in s and g2 in s for s in sets):
                edges.append({"metaLV": m, "gene1": g1, "gene2": g2})
        mods.append({"metaLV": m, "models": f"{summ.at[m, 'n_models']}/10", "cells": d["cells"],
                     "pathways": "; ".join(d["pathways"]),
                     "pathway_terms": "; ".join(top[top.metaLV == m].pathway.head(3)),
                     "drug": "; ".join(l1000.loc[l1000.metaLV.eq(m), "label"])})
    genes = pd.DataFrame(genes)

    info = pd.read_csv(a.phenotypes_info, sep="\t", dtype={"short_code": str}).set_index("short_code")
    gls = pd.read_csv(a.gls, dtype={"phenotype": str})
    gls = gls[gls.metaLV.isin(MODULES) & (gls.fdr < a.fdr)].copy()
    gls["n_cases"] = gls.phenotype.map(info.n_cases)
    gls = gls[gls.n_cases.isna() | (gls.n_cases >= a.min_cases)].sort_values("pvalue").drop_duplicates(["metaLV", "phenotype"])
    gls["trait_group"] = gls.phenotype.map(GROUP)
    lt = gls.dropna(subset=["trait_group"])
    lt.to_csv(out / "module_traits_all.csv", index=False)
    mt = lt.groupby(["metaLV", "trait_group"]).agg(fdr=("fdr", "min"), n_traits=("phenotype", "nunique"),
                                                  traits=("phenotype_desc", lambda s: "; ".join(s))).reset_index()

    px = Path(a.phenomexcan)
    pv = pd.read_pickle(px / "gene_assoc/smultixcan-mashr-pvalues.pkl")
    rcp = pd.read_pickle(px / "gene_assoc/fastenloc-torus-rcp.pkl")
    name = pd.read_pickle(px / "genes_metadata/genes_mapping_name_to_id.pkl")
    col = {c.split("-")[0]: c for c in pv.columns}
    rcol = {c.split("-")[0]: c for c in rcp.columns}
    fq = {}

    def fdr_of(code):
        if code not in fq:
            fq[code] = bh(pv[col[code]])
        return fq[code]

    tw = []
    for g in SCREEN + reference:
        gid = name[g]
        for code, (label, group) in TRAITS.items():
            tw.append({"gene": g, "screen": g in SCREEN, "trait_code": code, "trait": label, "trait_group": group,
                       "pvalue": pv.at[gid, col[code]] if gid in pv.index else np.nan,
                       "rcp": rcp.at[gid, rcol[code]] if gid in rcp.index and code in rcol else np.nan,
                       "fdr": fdr_of(code).get(gid, np.nan)})
    tw = pd.DataFrame(tw)
    own = genes[genes.crispr].merge(lt[["metaLV", "phenotype"]], on="metaLV")
    tw["own"] = [((own.gene == g) & (own.phenotype == c)).any() for g, c in zip(tw.gene, tw.trait_code)]

    ge, ring = [], []
    for m in MODULES:
        mtr = lt[lt.metaLV == m]
        for g in genes[genes.metaLV == m].gene:
            gid = name.get(g)
            hits = {}
            for r in mtr.itertuples():
                c = col.get(r.phenotype)
                if gid not in pv.index or c is None:
                    continue
                p = pv.at[gid, c]
                if p < 0.05:
                    hits[r.trait_group] = min(hits.get(r.trait_group, 1), p)
                    ring.append({"metaLV": m, "gene": g, "disease": r.trait_group, "trait": r.phenotype_desc,
                                 "p": p, "fdr": fdr_of(r.phenotype)[gid]})
            for k, v in hits.items():
                ge.append({"metaLV": m, "gene": g, "trait_group": k, "pvalue": v})
    ring = pd.DataFrame(ring).sort_values("p").drop_duplicates(["metaLV", "gene", "disease"]).sort_values(["metaLV", "p"])
    for s in SHORT:
        ring["trait"] = ring.trait.str.replace(s, "", regex=False)
    ring["p"] = ring.p.map(lambda x: f"{x:.1e}")
    ring["fdr"] = ring.fdr.map(lambda x: f"{x:.1e}")

    pd.DataFrame(mods).to_csv(out / "modules.csv", index=False)
    genes.to_csv(out / "module_genes.csv", index=False)
    pd.DataFrame(edges).to_csv(out / "pathway_edges.csv", index=False)
    mt.to_csv(out / "module_trait_groups.csv", index=False)
    pd.DataFrame(ge).to_csv(out / "gene_trait_edges.csv", index=False)
    ring.to_csv(out / "ring_pvalues.csv", index=False)
    tw.to_csv(out / "twas_rcp_matrix.csv", index=False)
    shutil.copyfile(a.composition, out / "projection_composition.csv")
    print(mt.to_string())


def main():
    commands = {
        "figure": figure_main,
    }
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=commands)
    args = parser.parse_args(sys.argv[1:2])
    commands[args.command](sys.argv[2:])


if __name__ == "__main__":
    main()

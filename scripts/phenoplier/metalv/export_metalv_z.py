#!/usr/bin/env python
# Prepare screen-associated metaLV loadings for PhenoPLIER.
import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metaz", required=True)
    ap.add_argument("--gsea-best", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--out-pkl", required=True)
    ap.add_argument("--out-mapping", required=True)
    a = ap.parse_args()

    z = pd.read_csv(a.metaz, index_col=0)
    g = pd.read_csv(a.gsea_best)
    sig = sorted(set(g[g.padj < a.fdr].lv), key=lambda m: int(m.replace("metaLV", "")))
    z = z[sig]
    mapping = pd.DataFrame({"phenoplier_lv": [f"LV{i}" for i in range(1, len(sig) + 1)], "metaLV": sig})
    z.columns = mapping.phenoplier_lv
    z.to_pickle(a.out_pkl)
    mapping.to_csv(a.out_mapping, index=False)
    print(f"{z.shape[1]} metaLVs x {z.shape[0]} genes -> {a.out_pkl}")


if __name__ == "__main__":
    main()

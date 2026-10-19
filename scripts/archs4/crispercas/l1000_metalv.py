#!/usr/bin/env python
# Score L1000 drug signatures against metaLV genes.
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def spec(m, i, v):

    x = m.copy()
    x[i] = v
    return (v - x.mean(0)) / x.std(0, ddof=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metaz", required=True)
    ap.add_argument("--gsea-best", required=True)
    ap.add_argument("--selected", required=True)
    ap.add_argument("--lincs", required=True)
    ap.add_argument("--drugbank", required=True)
    ap.add_argument("--symbol-ensembl", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--top-n", type=int, default=185)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    z = pd.read_csv(a.metaz, index_col=0)
    g = pd.read_csv(a.gsea_best)
    lvs = sorted(set(g[g.padj < a.fdr].lv), key=lambda m: int(m.replace("metaLV", "")))
    fig5 = list(pd.read_csv(a.selected).metaLV)
    assert set(fig5) <= set(lvs)

    mp = pd.read_csv(a.symbol_ensembl).dropna()
    mp = mp[mp.SYMBOL.isin(z.index)]
    mp = mp[~mp.SYMBOL.duplicated(keep=False) & ~mp.ENSEMBL.duplicated(keep=False)].set_index("SYMBOL").ENSEMBL

    lincs = pd.read_pickle(a.lincs)
    assert not lincs.isna().any().any()
    db = pd.read_csv(a.drugbank, sep="\t").set_index("drugbank_id").reindex(lincs.columns)

    top = {}
    for lv in lvs:
        s = z[lv].nlargest(a.top_n).index
        top[lv] = [mp[x] for x in s if x in mp.index and mp[x] in lincs.index]
    m = np.vstack([lincs.loc[top[lv]].values.mean(0) for lv in lvs])
    score = (m - m.mean(0)) / m.std(0, ddof=1)
    pd.DataFrame(score, index=lvs, columns=lincs.columns).to_csv(out / "l1000_scores.csv")

    rng = np.random.default_rng(a.seed)
    xv = lincs.values
    rows = []
    for lv in fig5:
        i = lvs.index(lv)
        n = len(top[lv])
        obs = score[i]
        s1 = np.zeros_like(obs)
        s2 = np.zeros_like(obs)
        ge = np.zeros_like(obs)
        for _ in range(a.n_perm):
            b = spec(m, i, xv[rng.choice(xv.shape[0], n, replace=False)].mean(0))
            s1 += b
            s2 += b * b
            ge += b >= obs
        mu = s1 / a.n_perm
        p = (ge + 1) / (a.n_perm + 1)
        r = pd.DataFrame({"metaLV": lv, "drugbank_id": lincs.columns, "drug": db.name.values,
                          "atc_codes": db.atc_codes.values, "categories": db.categories.values,
                          "mean_z": m[i], "score": obs, "perm_z": (obs - mu) / np.sqrt(s2 / a.n_perm - mu ** 2),
                          "perm_p": p, "fdr": multipletests(p, method="fdr_bh")[1], "n_genes": n})
        r = r.sort_values("score", ascending=False, kind="stable")
        r.insert(1, "rank", np.arange(1, len(r) + 1))
        rows.append(r)
        print(f"{lv}: {n} of {a.top_n} top genes in LINCS; {(r.fdr < a.fdr).sum()} drugs at FDR < {a.fdr}")
    pd.concat(rows).to_csv(out / "l1000_fig5_ranking.csv", index=False)


if __name__ == "__main__":
    main()

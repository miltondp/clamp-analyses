#!/usr/bin/env python
# Filter technical noise and cluster modules across models.
import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster


def cache_main(args):
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", required=True)
    ap.add_argument("--z", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)
    assert len(a.models) == len(a.z)

    mats, models, lvs, genes = [], [], [], None
    for model, path in zip(a.models, a.z):
        z = pd.read_csv(path, index_col=0)
        if genes is None:
            genes = z.index
        assert set(z.index) == set(genes), model
        z = z.reindex(genes)
        mats.append(z.to_numpy(dtype=np.float32).T)
        models += [model] * z.shape[1]
        lvs += list(z.columns)
        print(model, z.shape, flush=True)
    np.savez(a.out, X=np.concatenate(mats), models=np.asarray(models), lvs=np.asarray(lvs), genes=np.asarray(genes))


COLS = ["model", "LV", "gene_flag", "technical_joint", "sample_flag", "any_flag"]


def flags(gene_path, sample_path, threshold):
    g = pd.read_csv(gene_path)[["LV", "gene_flag"]]
    sm = pd.read_csv(sample_path)
    tj = sm[sm.covariate == "technical_joint"][["LV", "value"]].rename(columns={"value": "technical_joint"})
    d = g.merge(tj, on="LV", how="inner")
    assert len(d) == len(g), (gene_path, len(d), len(g))
    d["sample_flag"] = d.technical_joint.notna() & (d.technical_joint >= threshold)
    d["any_flag"] = d.gene_flag | d.sample_flag
    return d


def noise_main(args):
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-dir", required=True)
    ap.add_argument("--seeds-dir", required=True)
    ap.add_argument("--seeds", nargs="+", type=int, required=True)
    ap.add_argument("--threshold", type=float, default=0.25)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)

    ref_dir, seeds_dir = Path(a.reference_dir), Path(a.seeds_dir)
    ref = flags(ref_dir / "lv_gene_properties_canonical.csv", ref_dir / "lv_sample_metadata_canonical.csv", a.threshold)
    cls = pd.read_csv(ref_dir / "lv_technical_classification.csv").set_index("LV").any_flag
    assert (ref.set_index("LV").any_flag == cls.loc[ref.LV].values).all(), "S123 flags differ from the classification"
    rows = [ref.assign(model="S123")[COLS]]
    for s in a.seeds:
        d = seeds_dir / f"seed{s}"
        rows.append(flags(d / "lv_gene_properties_canonical.csv", d / "lv_sample_metadata_canonical.csv",
                          a.threshold).assign(model=f"S{s}")[COLS])
    out = pd.concat(rows)
    out.to_csv(a.out, index=False)
    print(out.groupby("model", sort=False)[["gene_flag", "sample_flag", "any_flag"]].sum())


def flags_main(args):
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", nargs="+", required=True)
    ap.add_argument("--threshold", type=float, default=0.25)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)

    rows = []
    for spec in a.model:
        name, path = spec.split("=", 1)
        d = Path(path)
        rows.append(flags(d / "lv_gene_properties_canonical.csv", d / "lv_sample_metadata_canonical.csv",
                          a.threshold).assign(model=name)[COLS])
    out = pd.concat(rows)
    out.to_csv(a.out, index=False)
    print(out.groupby("model", sort=False)[["gene_flag", "sample_flag", "any_flag"]].sum())


def cluster_main(args):
    ap = argparse.ArgumentParser()
    ap.add_argument("--z-cache", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--min-cos", type=float, default=0.5)
    ap.add_argument("--min-models", type=int, default=6)
    ap.add_argument("--reference", default="S123")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(args)
    ref = a.reference
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    d = np.load(a.z_cache, allow_pickle=True)
    X, models, lvs, genes = d["X"], d["models"], d["lvs"], pd.Index(d["genes"])
    flags = pd.read_csv(a.flags)
    noisy = set(zip(flags.loc[flags.any_flag, "model"], flags.loc[flags.any_flag, "LV"]))
    keep = np.array([(m, l) not in noisy for m, l in zip(models, lvs)])
    removed = pd.Series(models[~keep]).value_counts().rename("noisy_removed")
    X, models, lvs = X[keep], models[keep], lvs[keep]
    n = len(lvs)
    Xn = X / np.linalg.norm(X, axis=1, keepdims=True)

    dist = np.empty(n * (n - 1) // 2, dtype=np.float64)
    for start in range(0, n, 256):
        stop = min(start + 256, n)
        sim = Xn[start:stop] @ Xn.T
        for k, i in enumerate(range(start, stop)):
            off = n * i - i * (i + 1) // 2
            dist[off:off + n - i - 1] = np.clip(1 - sim[k, i + 1:], 0, 2)
    cl = fcluster(linkage(dist, method="complete"), t=1 - a.min_cos, criterion="distance")
    del dist

    allc = pd.DataFrame({"model": models, "LV": lvs, "cluster": cl})
    allc.to_csv(out / "all_lv_clusters.csv", index=False)
    stats = allc.groupby("cluster").agg(n_lvs=("LV", "size"), n_models=("model", "nunique"),
                                        **{f"has_{ref}": ("model", lambda s: (s == ref).any())})
    kept = stats[stats.n_models >= a.min_models].sort_values(["n_models", "n_lvs"], ascending=False)
    kept["metaLV"] = [f"metaLV{i + 1}" for i in range(len(kept))]
    print(f"LVs clustered: {n}; clusters: {len(stats)}; kept (>= {a.min_models}/10 models): {len(kept)}", flush=True)

    meta, members, summ = {}, [], []
    for cid, row in kept.iterrows():
        idx = np.flatnonzero(cl == cid)
        meta[row.metaLV] = np.median(X[idx], axis=0)
        C = Xn[idx] @ Xn[idx].T
        off = C[~np.eye(len(idx), dtype=bool)]
        for j, i in enumerate(idx):
            members.append({"metaLV": row.metaLV, "cluster": cid, "model": models[i], "LV": lvs[i],
                            "mean_cos_to_others": np.delete(C[j], j).mean()})
        summ.append({"metaLV": row.metaLV, "cluster": cid, "n_lvs": row.n_lvs, "n_models": row.n_models,
                     "min_cos": off.min(), "mean_cos": off.mean(),
                     f"{ref}_lvs": ";".join(lvs[idx][models[idx] == ref]),
                     "members": ";".join(f"{models[i]}:{lvs[i]}" for i in idx)})

    pd.DataFrame(meta, index=genes).to_csv(out / "metaZ.csv")
    pd.DataFrame(members).to_csv(out / "cluster_members.csv", index=False)
    pd.DataFrame(summ).to_csv(out / "cluster_summary.csv", index=False)
    removed.to_csv(out / "noisy_removed_per_model.csv")
    stats.n_models.value_counts().sort_index().to_csv(out / "cluster_model_count_distribution.csv")


def main():
    commands = {
        "cache": cache_main,
        "noise": noise_main,
        "flags": flags_main,
        "cluster": cluster_main,
    }
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=commands)
    args = parser.parse_args(sys.argv[1:2])
    commands[args.command](sys.argv[2:])


if __name__ == "__main__":
    main()

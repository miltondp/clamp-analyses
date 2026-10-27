#!/usr/bin/env python
import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve
from statsmodels.stats.multitest import multipletests


def empirical_statistics(observed, upper, lower, n_perm):
    p_upper = (upper + 1) / (n_perm + 1)
    p_lower = (lower + 1) / (n_perm + 1)
    p = np.minimum(1.0, 2 * np.minimum(p_upper, p_lower))
    fdr = multipletests(p.ravel(), method="fdr_bh")[1].reshape(p.shape)
    direction = np.where(observed > 0, "Activation-like", np.where(observed < 0, "Suppression-like", "Neutral"))
    return p_upper, p_lower, p, fdr, direction


def null_moments(w, v):
    mean = w.mean(axis=1)[:, None] * v.sum(axis=0)[None, :]
    variance = ((w - w.mean(axis=1)[:, None]) ** 2).sum(axis=1)[:, None] * ((v - v.mean(axis=0)) ** 2).sum(axis=0)[None, :] / (w.shape[1] - 1)
    return mean, np.sqrt(variance)


def permutation_tails(w, v, observed, n_perm, seed):
    rng = np.random.default_rng(seed)
    upper = np.zeros_like(observed, dtype=np.int64)
    lower = np.zeros_like(observed, dtype=np.int64)
    for start in range(0, n_perm, 32):
        count = min(32, n_perm - start)
        wb = np.concatenate([w[:, rng.permutation(w.shape[1])] for _ in range(count)], axis=0)
        draws = (wb @ v).reshape(count, w.shape[0], v.shape[1])
        upper += (draws >= observed[None, :, :]).sum(axis=0)
        lower += (draws <= observed[None, :, :]).sum(axis=0)
        if start == 0 or start + count == n_perm:
            print(f"{start + count}/{n_perm} permutations", flush=True)
    return upper, lower


def signature_groups(lincs, db):
    groups = {}
    for j, h in enumerate(pd.util.hash_pandas_object(lincs.T, index=False)):
        groups.setdefault(int(h), []).append(j)
    rows = []
    for group, js in enumerate(groups.values(), 1):
        assert all(np.array_equal(lincs.iloc[:, js[0]], lincs.iloc[:, j]) for j in js)
        ids = lincs.columns[js]
        names = " / ".join(sorted(set(db.loc[ids, "name"].dropna()), key=lambda s: (len(s), s)))
        rows.extend({"drugbank_id": drug, "signature_group": group, "signature_names": names,
                     "signature_drugbank_ids": ";".join(ids)} for drug in ids)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    for name in ["metaz", "selected", "reference-model", "lincs", "drugbank", "symbol-ensembl", "shortlist", "out-dir"]:
        ap.add_argument(f"--{name}", required=True)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20261027)
    a = ap.parse_args()
    assert a.n_perm > 0
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    z = pd.read_csv(a.metaz, index_col=0)
    selected = list(pd.read_csv(a.selected).metaLV)
    assert len(selected) == 4 and len(set(selected)) == 4 and set(selected) <= set(z.columns)
    mp = pd.read_csv(a.symbol_ensembl).dropna()
    mp = mp[mp.SYMBOL.isin(z.index)]
    mp = mp[~mp.SYMBOL.duplicated(keep=False) & ~mp.ENSEMBL.duplicated(keep=False)].set_index("SYMBOL").ENSEMBL
    z = z.loc[mp.index].set_axis(mp.values)
    lincs = pd.read_pickle(a.lincs)
    assert z.index.is_unique and lincs.index.is_unique and lincs.columns.is_unique
    assert np.isfinite(z.values).all() and np.isfinite(lincs.values).all()
    db = pd.read_csv(a.drugbank, sep="\t").set_index("drugbank_id").reindex(lincs.columns)
    assert db.name.notna().all()
    l2 = float(subprocess.check_output(["Rscript", "-e", 'cat(sprintf("%.17g", readRDS(commandArgs(TRUE)[1])$L2))', a.reference_model], text=True))
    assert np.isfinite(l2) and l2 > 0
    zz = z.to_numpy(dtype=float)
    gram = zz.T @ zz + l2 * np.eye(z.shape[1])
    projected = solve(gram, zz.T @ lincs.reindex(z.index).fillna(0).to_numpy(dtype=float), assume_a="pos")
    pd.DataFrame(projected, index=z.columns, columns=lincs.columns).to_csv(out / "l1000_scores.csv")
    indices = z.columns.get_indexer(selected)
    projector = solve(gram, np.eye(z.shape[1])[:, indices], assume_a="pos").T @ zz.T
    shared = z.index.isin(lincs.index)
    w = np.ascontiguousarray(projector[:, shared])
    v = np.ascontiguousarray(lincs.reindex(z.index[shared]).to_numpy(dtype=float))
    observed = w @ v
    assert np.allclose(observed, projected[indices], atol=1e-12, rtol=1e-9)
    mu, sd = null_moments(w, v)
    assert np.isfinite(sd).all() and (sd > 0).all()
    null_z = (observed - mu) / sd
    upper, lower = permutation_tails(w, v, observed, a.n_perm, a.seed)
    p_upper, p_lower, p, fdr, direction = empirical_statistics(observed, upper, lower, a.n_perm)
    groups = signature_groups(lincs, db)
    rows = []
    for i, m in enumerate(selected):
        rows.append(pd.DataFrame({"metaLV": m, "drugbank_id": lincs.columns, "drug": db.name.values,
                                  "atc_codes": db.atc_codes.values, "categories": db.categories.values,
                                  "projection": observed[i], "null_mean": mu[i], "null_sd": sd[i],
                                  "null_z": null_z[i], "p_upper": p_upper[i], "p_lower": p_lower[i],
                                  "p_two_sided": p[i], "fdr": fdr[i], "direction": direction[i],
                                  "direction_agrees_with_null": np.sign(observed[i]) == np.sign(null_z[i]),
                                  "n_genes": int(shared.sum())}))
    results = pd.concat(rows, ignore_index=True).merge(groups, on="drugbank_id", validate="many_to_one")
    unique = results.sort_values("drugbank_id").drop_duplicates(["metaLV", "signature_group"]).copy()
    unique = unique.sort_values("null_z", key=lambda s: s.abs(), ascending=False, kind="stable")
    eligible = unique.direction_agrees_with_null & unique.direction.ne("Neutral")
    unique.loc[eligible, "direction_rank"] = unique[eligible].groupby(["metaLV", "direction"]).cumcount() + 1
    unique["direction_rank"] = unique.direction_rank.astype("Int64")
    results = results.merge(unique[["metaLV", "signature_group", "direction_rank"]], on=["metaLV", "signature_group"], validate="many_to_one")
    results.to_csv(out / "l1000_fig5_ranking.csv", index=False)
    unique.to_csv(out / "l1000_unique_signatures.csv", index=False)
    picks = pd.read_csv(a.shortlist).merge(results, on=["metaLV", "drug", "direction"], validate="one_to_one")
    assert len(picks) == len(pd.read_csv(a.shortlist)) == 10
    assert picks.fdr.lt(a.fdr).all() and picks.direction_agrees_with_null.all()
    assert not picks.duplicated(["metaLV", "signature_group"]).any()
    picks.to_csv(out / "l1000_fig5_curated.csv", index=False)
    assert results.groupby("metaLV").size().eq(lincs.shape[1]).all()
    assert np.isfinite(results[["projection", "null_mean", "null_sd", "null_z", "p_two_sided", "fdr"]]).all().all()
    metadata = {"l2": l2, "n_mgms_projected": z.shape[1], "tested_mgms": selected,
                "mapped_model_genes": len(z), "genes_measured_in_lincs": int(shared.sum()),
                "drug_entries": lincs.shape[1], "unique_signatures": groups.signature_group.nunique(),
                "n_permutations": a.n_perm, "seed": a.seed, "bh_family_size": fdr.size,
                "two_sided_p": "min(1, 2*min((upper+1)/(N+1), (lower+1)/(N+1)))",
                "null_moments": "exact mean and SD over all gene-label permutations",
                "direction": "sign of projected coefficient; ranks require agreement with null departure"}
    (out / "l1000_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(picks[["metaLV", "drug", "direction", "null_z", "fdr"]].to_string(index=False))


if __name__ == "__main__":
    main()

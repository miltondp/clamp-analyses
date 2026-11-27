#!/usr/bin/env python
"""Consensus-LV evaluation CLI (see PLAN.md). Subcommands: flags, zcache-check, reproduction,
linkage, hier-set, hier-sets, set-check, pool-hits, evaluate, compare, label-null, aggregate."""
import argparse
import filecmp
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "archs4/crispercas"))
from cluster_metalv import COLS  # noqa: E402
import cev  # noqa: E402  (same directory as this script)

MAX_SHOW = 20  # differing metaLVs listed in a failed reproduction


def rethreshold(d, r2):
    # Same rule as cluster_metalv.flags(), applied to an already-merged flags table.
    d = d.copy()
    d["sample_flag"] = d.technical_joint.notna() & (d.technical_joint >= r2)
    d["any_flag"] = d.gene_flag | d.sample_flag
    return d[COLS]


def flags_main(args):
    ap = argparse.ArgumentParser(prog="run.py flags")
    ap.add_argument("--legacy", required=True)
    ap.add_argument("--legacy-r2", type=float, required=True, help="noise_r2 the legacy flags were written at")
    ap.add_argument("--r2", type=float, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)

    legacy = pd.read_csv(a.legacy)
    assert list(legacy.columns) == COLS, legacy.columns
    assert legacy.technical_joint.notna().all(), "technical_joint has NAs"
    again = rethreshold(legacy, a.legacy_r2)
    for c in ["sample_flag", "any_flag"]:
        assert (again[c] == legacy[c]).all(), f"legacy {c} is not technical_joint >= {a.legacy_r2}"
    out = rethreshold(legacy, a.r2)
    out.to_csv(a.out, index=False)
    print(out.groupby("model", sort=False)[["gene_flag", "sample_flag", "any_flag"]].sum())
    print(f"r2={a.r2}: {int(out.any_flag.sum())} flagged, {int((~out.any_flag).sum())} kept")


def zcache_check_main(args):
    ap = argparse.ArgumentParser(prog="run.py zcache-check")
    ap.add_argument("--z-cache", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--n-genes", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)

    X, models, lvs, genes = cev.load_cache(a.z_cache)
    keys = list(zip(models, lvs))
    flags = pd.read_csv(a.flags)
    fkeys = list(zip(flags.model, flags.LV))
    problems = []
    if len(set(keys)) != len(keys):
        problems.append("duplicate (model, LV) keys in the cache")
    if set(keys) != set(fkeys):
        problems.append(f"keys differ: {len(set(keys) - set(fkeys))} only in cache, "
                        f"{len(set(fkeys) - set(keys))} only in flags")
    if X.shape != (len(keys), a.n_genes):
        problems.append(f"X shape {X.shape}, expected ({len(keys)}, {a.n_genes})")
    if len(set(genes)) != a.n_genes:
        problems.append(f"{len(set(genes))} unique genes, expected {a.n_genes}")
    if not np.isfinite(X).all():
        problems.append("X has non-finite values")
    if problems:
        sys.exit("zcache check failed:\n  " + "\n  ".join(problems))
    msg = f"{len(keys)} LVs from {len(pd.unique(models))} models x {len(genes)} genes; keys equal {a.flags}"
    print(msg)
    Path(a.out).write_text(msg + "\n")


def member_sets(df):
    return {r.metaLV: frozenset(r.members.split(";")) for r in df.itertuples()}


def reproduction_main(args):
    ap = argparse.ArgumentParser(prog="run.py reproduction")
    ap.add_argument("--set", required=True)
    ap.add_argument("--metaz", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--external", required=True)
    ap.add_argument("--members", required=True)
    ap.add_argument("--tol", type=float, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)

    ours, ext = pd.read_csv(a.metaz, index_col=0), pd.read_csv(a.external, index_col=0)
    names_equal = list(ours.columns) == list(ext.columns)
    genes_equal = list(ours.index) == list(ext.index)
    max_dz = float(np.abs(ours.to_numpy() - ext.to_numpy()).max()) if names_equal and genes_equal else np.nan
    z_ok = max_dz < a.tol  # False when max_dz is NaN
    print(f"Z: ours {ours.shape}, external {ext.shape}; names equal {names_equal}, genes equal {genes_equal}, "
          f"max|dZ| {max_dz}")
    if not names_equal:
        print(f"  metaLVs only ours: {len(set(ours.columns) - set(ext.columns))}, "
              f"only external: {len(set(ext.columns) - set(ours.columns))}")
    elif genes_equal and not z_ok:
        bad = (np.abs(ours - ext) >= a.tol).any()
        print(f"  {int(bad.sum())} metaLVs differ, first: {list(bad[bad].index[:MAX_SHOW])}")
    byte_identical = filecmp.cmp(a.metaz, a.external, shallow=False)
    print(f"metaZ.csv byte-identical to external: {byte_identical}")

    s_ours, s_leg = pd.read_csv(a.summary), pd.read_csv(a.members)
    m_ours, m_leg = member_sets(s_ours), member_sets(s_leg)
    counts_ours = s_ours.set_index("metaLV")[["n_models", "n_lvs"]]
    counts_leg = s_leg.set_index("metaLV")[["n_models", "n_lvs"]]
    common = [k for k in m_leg if k in m_ours]
    diff = [k for k in common if m_ours[k] != m_leg[k]
            or not (counts_ours.loc[k] == counts_leg.loc[k]).all()]
    members_ok = set(m_ours) == set(m_leg) and not diff
    print(f"members: ours {len(m_ours)} metaLVs, legacy {len(m_leg)}; {len(diff)} of {len(common)} common differ")
    for k in diff[:MAX_SHOW]:
        print(f"  {k}: only ours {sorted(m_ours[k] - m_leg[k])}; only legacy {sorted(m_leg[k] - m_ours[k])}")
    # Same clusters under different numbering?
    if diff or set(m_ours) != set(m_leg):
        same_clusters = set(m_ours.values()) == set(m_leg.values())
        print(f"  same member sets ignoring metaLV names: {same_clusters}")

    passed = z_ok and members_ok
    pd.DataFrame([{"set": a.set, "n_ours": ours.shape[1], "n_external": ext.shape[1],
                   "names_equal": names_equal, "genes_equal": genes_equal, "max_abs_dZ": max_dz,
                   "byte_identical": byte_identical, "n_member_diff": len(diff),
                   "members_equal": members_ok, "pass": passed}]).to_csv(a.out, index=False)
    print("PASS" if passed else "FAIL")
    if not passed:
        sys.exit(1)


def add_pool_args(ap):
    ap.add_argument("--z-cache", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--runs", nargs="+", default=None, help="keep only these runs (split-half / LORO subset)")
    ap.add_argument("--gene-perm-seed", type=int, default=None,
                    help="shuffle the gene order within every pool LV (gene-permutation null)")


def subset_pool(a):
    return cev.load_pool(a.z_cache, a.flags, runs=a.runs, gene_perm_seed=a.gene_perm_seed)


def linkage_main(args):
    ap = argparse.ArgumentParser(prog="run.py linkage")
    add_pool_args(ap)
    ap.add_argument("--method", required=True, choices=["complete", "average"])
    ap.add_argument("--dist", required=True, choices=["cos", "jac"],
                    help="cos: 1 - cosine of raw loadings; jac: 1 - Jaccard of the top-n positive genes")
    ap.add_argument("--top-n", type=int, required=True, help="top genes per LV (used by --dist jac)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)

    pool = subset_pool(a)
    t0 = time.time()
    if a.dist == "cos":
        L = cev.hier_linkage(pool.Xn, a.method)
    else:
        B =cev.top_genes(pool.X, a.top_n)
        size = B.sum(1)
        print(f"top-{a.top_n} sets: {(size < a.top_n).sum()} LVs with fewer genes, {(size == 0).sum()} empty; "
              f"size median {np.median(size):.0f}, min {size.min()}")
        L = cev.hier_linkage_jac(B, a.method)
    np.savez(a.out, L=L, models=pool.models, lvs=pool.lvs, dist=a.dist)
    print(f"{a.dist} {a.method} linkage on {len(pool.lvs)} LVs in {time.time() - t0:.1f} s -> {a.out}")


def write_set(out_dir, Z, lvs, members, meta, z_link=None):
    """z_link: an identical set's Z.npy to symlink instead of writing Z again."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if z_link is None:
        np.save(out / "Z.npy", Z)
    else:
        (out / "Z.npy").symlink_to(Path("..") / Path(z_link).parent.name / "Z.npy")
    lvs.to_csv(out / "lvs.csv", index=False)
    members.to_csv(out / "members.csv", index=False)
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def read_set(set_dir):
    sd = Path(set_dir)
    return (np.load(sd / "Z.npy"), pd.read_csv(sd / "lvs.csv"), pd.read_csv(sd / "members.csv"),
            json.loads((sd / "meta.json").read_text()))


def load_linkage(path, pool):
    d = np.load(path, allow_pickle=True)
    assert (d["models"] == pool.models).all() and (d["lvs"] == pool.lvs).all(), "linkage keys differ from the pool"
    return d


def build_hier_set(pool, d, name, subset, p, out_dir, built=None):
    """Cut the linkage d at the set's min_<dist>, keep clusters spanning enough runs, write the set.

    built: {(cut, min_models, combine): (out_dir, set)} shared across the sets of one linkage. In 5-run
    halves min_frac 0.3/0.4 (and 0.5/0.6) resolve to the same min_models, so such a set is built once
    and its Z.npy symlinked."""
    n_runs = len(pd.unique(pool.models))
    min_models = p["min_models"] if "min_models" in p else cev.resolve_min_models(p["min_frac"], n_runs)
    dist = p["dist"]
    assert str(d["dist"]) == dist, f"linkage dist {d['dist']} differs from the set's {dist}"
    key = (p[f"min_{dist}"], min_models, p["combine"])
    if built is not None and key in built:
        z_dir, (Z, lvs, members, n_clusters) = built[key]
    else:
        labels = cev.hier_partition(d["L"], p[f"min_{dist}"])
        z_dir, (Z, lvs, members, n_clusters) = None, cev.build_set(pool, labels, min_models, p["combine"])
        if built is not None:
            built[key] = (out_dir, (Z, lvs, members, n_clusters))
    meta = {"set": name, "subset": subset, **p, "n_runs": n_runs,
            "min_models_resolved": min_models, "n_pool": len(pool.lvs), "n_clusters": n_clusters,
            "n_metalvs": len(lvs), "n_members": len(members)}
    write_set(out_dir, Z, lvs, members, meta, z_link=None if z_dir is None else Path(z_dir) / "Z.npy")
    return meta


def hier_set_main(args):
    ap = argparse.ArgumentParser(prog="run.py hier-set")
    ap.add_argument("--set", required=True)
    ap.add_argument("--subset", required=True)
    add_pool_args(ap)
    ap.add_argument("--linkage", required=True)
    ap.add_argument("--params", required=True,
                    help="JSON set definition (dist, min_cos or min_jac, min_models or min_frac, combine, ...), copied to meta.json")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(args)
    pool = subset_pool(a)
    meta = build_hier_set(pool, load_linkage(a.linkage, pool), a.set, a.subset, json.loads(a.params), a.out_dir)
    print(json.dumps(meta))


def hier_sets_main(args):
    """Every set cut from one linkage, built from one pool load (split-half / LORO / null subsets)."""
    ap = argparse.ArgumentParser(prog="run.py hier-sets")
    ap.add_argument("--subset", required=True)
    add_pool_args(ap)
    ap.add_argument("--linkage", required=True)
    ap.add_argument("--params", required=True, help="JSON {set name: set definition}")
    ap.add_argument("--out-dir", required=True, help="one sub-directory per set")
    a = ap.parse_args(args)
    pool = subset_pool(a)
    d = load_linkage(a.linkage, pool)
    built = {}
    for name, p in json.loads(a.params).items():
        meta = build_hier_set(pool, d, name, a.subset, p, Path(a.out_dir) / name, built)
        print(json.dumps(meta))


def set_check_main(args):
    """Our builder vs the unchanged cluster_metalv.py output for the same set."""
    ap = argparse.ArgumentParser(prog="run.py set-check")
    ap.add_argument("--set-dir", required=True)
    ap.add_argument("--reproduce-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)

    sd, rd = Path(a.set_dir), Path(a.reproduce_dir)
    Z, lvs, members, _ = read_set(sd)
    ref_z, ref_s = pd.read_csv(rd / "metaZ.csv", index_col=0), pd.read_csv(rd / "cluster_summary.csv")
    problems = cev.compare_to_metalv(Z, lvs, members, ref_z, ref_s)
    if problems:
        sys.exit(f"set check failed for {sd}:\n  " + "\n  ".join(problems))
    msg = f"{sd} equals {rd} exactly ({Z.dtype}): {len(lvs)} metaLVs, {len(members)} members"
    print(msg)
    Path(a.out).write_text(msg + "\n")


def pool_hits_main(args):
    ap = argparse.ArgumentParser(prog="run.py pool-hits")
    ap.add_argument("--z-cache", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--top-n", type=int, required=True)
    ap.add_argument("--reproducible-cos", type=float, required=True)
    ap.add_argument("--reproducible-frac", type=float, required=True)
    ap.add_argument("--jac-null-q", type=float, required=True)
    ap.add_argument("--null-pairs", type=int, required=True)
    ap.add_argument("--null-seed", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--null-out", required=True)
    a = ap.parse_args(args)

    pool = cev.load_pool(a.z_cache, a.flags)
    runs = list(pd.unique(pool.models))
    B = cev.top_genes(pool.X, a.top_n)
    hits, arg, jac = cev.pool_best_hits(pool.Xn, B, pool.models, runs)
    recip = cev.reciprocal(arg, pool.models, runs)
    thr, quant, pairs = cev.jaccard_null(B, pool.models, runs, a.null_pairs, a.null_seed, a.jac_null_q)

    # Reproducible LV: a reciprocal best hit with cos >= c and Jaccard > the null threshold in
    # >= min_runs - 1 other runs. Counts at each stage show what each condition removes.
    min_runs = cev.resolve_min_models(a.reproducible_frac, len(runs))
    stages = cev.hit_stages(hits, recip, jac, a.reproducible_cos, thr)
    repro = cev.reproducible(stages["cos_reciprocal_jaccard"], min_runs)
    np.savez(a.out, hits=hits, arg=arg, jac=jac, repro=repro,
             runs=np.asarray(runs), models=pool.models, lvs=pool.lvs)
    real = stages["cos_reciprocal"]
    info = {"n_pool": len(pool.lvs), "n_runs": len(runs), "top_n": a.top_n,
            "cos": a.reproducible_cos, "min_frac": a.reproducible_frac, "min_runs": min_runs,
            "jac_null": {"q": a.jac_null_q, "threshold": thr, "seed": a.null_seed, "pairs": pairs, **quant},
            "observed_jac_reciprocal_cos": cev.quantiles(jac[real], [0.01, 0.1, 0.25, 0.5, 0.75, 0.9]),
            "frac_reciprocal_cos_within_null": float((jac[real] <= thr).mean()),
            "n_reproducible": {k: int(cev.reproducible(v, min_runs).sum()) for k, v in stages.items()}}
    Path(a.null_out).write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2))


def evaluate_main(args):
    ap = argparse.ArgumentParser(prog="run.py evaluate")
    ap.add_argument("--set-dir", required=True)
    ap.add_argument("--z-cache", required=True)
    ap.add_argument("--best-hits", required=True, help="pool-hits output; defines the evaluation pool")
    ap.add_argument("--pool-r2", type=float, required=True)
    ap.add_argument("--top-n", type=int, required=True)
    ap.add_argument("--cos-thresholds", type=float, nargs="+", required=True)
    ap.add_argument("--single-gene", type=float, required=True,
                    help="top1_share above which an LV (consensus or pool) counts as single-gene dominated")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(args)

    sd, out = Path(a.set_dir), Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    Z, lvs, members, meta = read_set(sd)
    cache = cev.load_pool(a.z_cache)  # all LVs: a set may hold LVs flagged at pool_r2 (built at another r2)
    assert Z.shape == (len(cache.genes), len(lvs)), (Z.shape, len(cache.genes), len(lvs))

    row_of = {k: i for i, k in enumerate(cache.keys)}
    member_rows = np.array([row_of[k] for k in cev.lv_keys(members.model, members.LV)])
    codes = pd.Categorical(members.metaLV, categories=lvs.metaLV).codes
    assert (codes >= 0).all()
    per = cev.cluster_metrics(cache, member_rows, codes, Z, a.top_n)
    assert (per.n_lvs.to_numpy() == lvs.n_lvs.to_numpy()).all() and (per.n_models.to_numpy() == lvs.n_models.to_numpy()).all()
    per.insert(0, "metaLV", lvs.metaLV)

    h = np.load(a.best_hits, allow_pickle=True)
    pool_rows = np.array([row_of[k] for k in cev.lv_keys(h["models"], h["lvs"])])
    runs = list(h["runs"])
    best = cev.best_cos(cache.Xn[pool_rows], cev.l2norm(Z.T))
    multi = cev.top1_share(cache.X[pool_rows]) <= a.single_gene   # pool LVs not single-gene dominated
    cov = cev.coverage(best, h["models"], runs, h["repro"], multi, a.cos_thresholds)

    params = {k: meta[k] for k in meta if k not in ("n_metalvs", "n_members")}
    row = {**params, "pool_r2": a.pool_r2, "n_pool_eval": len(pool_rows),
           **cev.summarize(per, cov, a.cos_thresholds, a.single_gene)}
    per.to_csv(out / "per_cluster.csv", index=False)
    cov.to_csv(out / "coverage.csv", index=False)
    pd.DataFrame([row]).to_csv(out / "metrics.csv", index=False)
    print(pd.Series(row).to_string())


def compare_main(args):
    """Split-half (Hungarian, half a vs half b), LORO recovery and the null counts of one set."""
    ap = argparse.ArgumentParser(prog="run.py compare")
    ap.add_argument("--full-dir", required=True)
    ap.add_argument("--half-a", nargs="+", required=True, help="half-a set dirs, paired in order with --half-b")
    ap.add_argument("--half-b", nargs="+", required=True)
    ap.add_argument("--loro", nargs="+", required=True, help="leave-one-run-out set dirs")
    ap.add_argument("--geneperm", nargs="+", required=True, help="gene-permutation null set dirs")
    ap.add_argument("--label-null", required=True, help="label-null table of the set's full linkage")
    ap.add_argument("--top-n", type=int, required=True)
    ap.add_argument("--match-cos", type=float, required=True)
    ap.add_argument("--match-jac", type=float, required=True, help="top-gene Jaccard cut of the Jaccard-matched scores")
    ap.add_argument("--single-gene", type=float, required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(args)
    assert len(a.half_a) == len(a.half_b)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    Z, lvs, _, meta = read_set(a.full_dir)
    t, tj = f"c{a.match_cos:g}", f"j{a.match_jac:g}"

    rows = []
    for i, (da, db) in enumerate(zip(a.half_a, a.half_b)):
        Za, _, _, ma = read_set(da)
        Zb, _, _, mb = read_set(db)
        assert ma["set"] == mb["set"] == meta["set"], (ma["set"], mb["set"], meta["set"])
        assert ma["min_models_resolved"] == mb["min_models_resolved"]
        rows.append({"split": i, "subset_a": ma["subset"], "subset_b": mb["subset"],
                     "min_models": ma["min_models_resolved"],
                     **cev.split_half_scores(Za, Zb, a.top_n, a.match_cos, a.match_jac, a.single_gene)})
    sh = pd.DataFrame(rows)
    sh.to_csv(out / "split_half.csv", index=False)

    loro = lvs[["metaLV", "n_lvs", "n_models"]].copy()
    loro["top1_share"] = cev.top1_share(Z.T)
    best, best_j, loro_sizes = [], [], []
    for d in a.loro:   # one subset set in memory at a time
        Zs, _, _, ms = read_set(d)
        c, j = cev.loro_recovery(Z, Zs, a.top_n)
        loro[ms["subset"]] = c
        best.append(c)
        best_j.append(j)
        loro_sizes.append(ms["n_metalvs"])
    loro["recovery"] = (np.column_stack(best) >= a.match_cos).mean(1)
    loro["recovery_jac"] = (np.column_stack(best_j) >= a.match_jac).mean(1)
    loro.to_csv(out / "loro.csv", index=False)

    row = {"set": meta["set"], "sh_n_splits": len(sh), "sh_n_empty": int(sh.n_empty.sum()),
           "sh_min_models": int(sh.min_models.iloc[0]), "sh_n_lvs": float(sh[["n_a", "n_b"]].to_numpy().mean())}
    mask = cev.strata(Z, a.single_gene)
    for prefix, rec_col, thr in [("", "recovery", t), ("j", "recovery_jac", tj)]:
        for st in cev.STRATA:
            row[f"sh_{prefix}share_{thr}_{st}"] = sh[f"{prefix}share_{st}"].mean()
            row[f"sh_{prefix}share_{thr}_{st}_sd"] = sh[f"{prefix}share_{st}"].std()
            row[f"sh_{prefix}pair_cos_{st}"] = sh[f"{prefix}pair_cos_{st}"].mean()
            row[f"sh_{prefix}pair_jac_{st}"] = sh[f"{prefix}pair_jac_{st}"].mean()
            rec = loro[rec_col][mask[st]]
            row[f"loro_recovery_{thr}_{st}"] = rec.mean() if len(rec) else np.nan
            row[f"loro_frac_stable_{thr}_{st}"] = (rec == 1).mean() if len(rec) else np.nan
    row["loro_n_lvs"] = float(np.mean(loro_sizes))
    row["geneperm_n_lvs"] = float(np.mean([json.loads((Path(d) / "meta.json").read_text())["n_metalvs"]
                                           for d in a.geneperm]))

    ln = pd.read_csv(a.label_null)
    cut = meta[f"min_{meta['dist']}"]
    ln = ln[np.isclose(ln.cut, cut) & (ln.m == meta["min_models_resolved"])]
    assert len(ln), f"no label-null rows for cut {cut}, m {meta['min_models_resolved']}"
    row.update({"labelnull_observed": int(ln["observed"].iloc[0]), "labelnull_mean": ln["null"].mean(),
                "labelnull_max": int(ln["null"].max())})
    pd.DataFrame([row]).to_csv(out / "robust.csv", index=False)
    print(sh.to_string(index=False))
    print(pd.Series(row).to_string())


def label_null_main(args):
    """Clusters spanning >= m runs at each cut of one linkage, observed and with run labels shuffled."""
    ap = argparse.ArgumentParser(prog="run.py label-null")
    ap.add_argument("--linkage", required=True)
    ap.add_argument("--cuts", type=float, nargs="+", required=True, help="similarity cuts (min_cos or min_jac)")
    ap.add_argument("--n-perm", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)
    d = np.load(a.linkage, allow_pickle=True)
    tabs = []
    for cut in a.cuts:
        labels = cev.hier_partition(d["L"], cut)
        tab = cev.label_perm_null(labels, d["models"], a.n_perm, a.seed)
        tab.insert(0, "cut", cut)
        tabs.append(tab)
    res = pd.concat(tabs, ignore_index=True)
    res.to_csv(a.out, index=False)
    print(res.groupby(["cut", "m"]).agg(observed=("observed", "first"), null_mean=("null", "mean"),
                                        null_max=("null", "max")).to_string())


def aggregate_main(args):
    ap = argparse.ArgumentParser(prog="run.py aggregate")
    ap.add_argument("--metrics", nargs="+", required=True)
    ap.add_argument("--params", nargs="+", required=True,
                    help="set-definition columns that lead the headline table (absent ones are skipped)")
    ap.add_argument("--headline", nargs="+", required=True, help="metric columns of the headline table")
    ap.add_argument("--robust", nargs="*", default=[], help="compare robust.csv files (grid sets only)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--headline-out", required=True)
    ap.add_argument("--robust-out", default=None)
    a = ap.parse_args(args)
    d = pd.concat([pd.read_csv(p) for p in a.metrics], ignore_index=True)
    d.to_csv(a.out, index=False)
    if a.robust:
        rob = pd.concat([pd.read_csv(p) for p in a.robust], ignore_index=True)
        d = d.merge(rob, on="set", how="left")
        d.loc[d.set.isin(rob.set), [c for c in a.params if c in d.columns] + ["n_lvs", *rob.columns[1:]]].to_csv(
            a.robust_out, index=False)
    head = d[[c for c in a.params if c in d.columns] + a.headline]
    head.to_csv(a.headline_out, index=False)
    print(f"{len(d)} sets -> {a.out}, {a.headline_out}")
    print(head.to_string(index=False))


def main():
    commands = {"flags": flags_main, "zcache-check": zcache_check_main, "reproduction": reproduction_main,
                "linkage": linkage_main, "hier-set": hier_set_main, "hier-sets": hier_sets_main,
                "set-check": set_check_main, "pool-hits": pool_hits_main, "evaluate": evaluate_main,
                "compare": compare_main, "label-null": label_null_main, "aggregate": aggregate_main}
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=commands)
    args = parser.parse_args(sys.argv[1:2])
    commands[args.command](sys.argv[2:])


if __name__ == "__main__":
    main()

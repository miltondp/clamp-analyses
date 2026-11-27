#!/usr/bin/env python
"""Consensus-LV evaluation CLI (see PLAN.md). Subcommands: flags, zcache-check, reproduction,
linkage, hier-set, set-check, pool-hits, evaluate, aggregate."""
import argparse
import filecmp
import json
import sys
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


def linkage_main(args):
    ap = argparse.ArgumentParser(prog="run.py linkage")
    ap.add_argument("--z-cache", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--method", required=True, choices=["complete", "average"])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(args)

    pool = cev.load_pool(a.z_cache, a.flags)
    L = cev.hier_linkage(pool.Xn, a.method)
    np.savez(a.out, L=L, models=pool.models, lvs=pool.lvs)
    print(f"{a.method} linkage on {len(pool.lvs)} LVs -> {a.out}")


def write_set(out_dir, Z, lvs, members, meta):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "Z.npy", Z)
    lvs.to_csv(out / "lvs.csv", index=False)
    members.to_csv(out / "members.csv", index=False)
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


def read_set(set_dir):
    sd = Path(set_dir)
    return (np.load(sd / "Z.npy"), pd.read_csv(sd / "lvs.csv"), pd.read_csv(sd / "members.csv"),
            json.loads((sd / "meta.json").read_text()))


def hier_set_main(args):
    ap = argparse.ArgumentParser(prog="run.py hier-set")
    ap.add_argument("--set", required=True)
    ap.add_argument("--subset", required=True)
    ap.add_argument("--z-cache", required=True)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--linkage", required=True)
    ap.add_argument("--params", required=True,
                    help="JSON set definition (min_cos, min_models or min_frac, combine, ...), copied to meta.json")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args(args)
    p = json.loads(a.params)

    pool = cev.load_pool(a.z_cache, a.flags)
    d = np.load(a.linkage, allow_pickle=True)
    assert (d["models"] == pool.models).all() and (d["lvs"] == pool.lvs).all(), "linkage keys differ from the pool"
    n_runs = len(pd.unique(pool.models))
    min_models = p["min_models"] if "min_models" in p else cev.resolve_min_models(p["min_frac"], n_runs)
    labels = cev.hier_partition(d["L"], p["min_cos"])
    Z, lvs, members, n_clusters = cev.build_set(pool, labels, min_models, p["combine"])
    meta = {"set": a.set, "subset": a.subset, **p, "n_runs": n_runs,
            "min_models_resolved": min_models, "n_pool": len(pool.lvs), "n_clusters": n_clusters,
            "n_metalvs": len(lvs), "n_members": len(members)}
    write_set(a.out_dir, Z, lvs, members, meta)
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
                    help="top1_share above which a consensus LV counts as single-gene dominated")
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
    cov = cev.coverage(best, h["models"], runs, h["repro"], a.cos_thresholds)

    params = {k: meta[k] for k in meta if k not in ("n_metalvs", "n_members")}
    row = {**params, "pool_r2": a.pool_r2, "n_pool_eval": len(pool_rows),
           **cev.summarize(per, cov, a.cos_thresholds, a.single_gene)}
    per.to_csv(out / "per_cluster.csv", index=False)
    cov.to_csv(out / "coverage.csv", index=False)
    pd.DataFrame([row]).to_csv(out / "metrics.csv", index=False)
    print(pd.Series(row).to_string())


def aggregate_main(args):
    ap = argparse.ArgumentParser(prog="run.py aggregate")
    ap.add_argument("--metrics", nargs="+", required=True)
    ap.add_argument("--params", nargs="+", required=True,
                    help="set-definition columns that lead the headline table (absent ones are skipped)")
    ap.add_argument("--headline", nargs="+", required=True, help="metric columns of the headline table")
    ap.add_argument("--out", required=True)
    ap.add_argument("--headline-out", required=True)
    a = ap.parse_args(args)
    d = pd.concat([pd.read_csv(p) for p in a.metrics], ignore_index=True)
    d.to_csv(a.out, index=False)
    head = d[[c for c in a.params if c in d.columns] + a.headline]
    head.to_csv(a.headline_out, index=False)
    print(f"{len(d)} sets -> {a.out}, {a.headline_out}")
    print(head.to_string(index=False))


def main():
    commands = {"flags": flags_main, "zcache-check": zcache_check_main, "reproduction": reproduction_main,
                "linkage": linkage_main, "hier-set": hier_set_main, "set-check": set_check_main,
                "pool-hits": pool_hits_main, "evaluate": evaluate_main, "aggregate": aggregate_main}
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=commands)
    args = parser.parse_args(sys.argv[1:2])
    commands[args.command](sys.argv[2:])


if __name__ == "__main__":
    main()

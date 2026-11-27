#!/usr/bin/env python
"""Consensus-LV evaluation CLI (see PLAN.md). Subcommands: flags, zcache-check, reproduction."""
import argparse
import filecmp
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "archs4/crispercas"))
from cluster_metalv import COLS  # noqa: E402

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

    d = np.load(a.z_cache, allow_pickle=True)
    X, models, lvs, genes = d["X"], d["models"], d["lvs"], d["genes"]
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


def main():
    commands = {"flags": flags_main, "zcache-check": zcache_check_main, "reproduction": reproduction_main}
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=commands)
    args = parser.parse_args(sys.argv[1:2])
    commands[args.command](sys.argv[2:])


if __name__ == "__main__":
    main()

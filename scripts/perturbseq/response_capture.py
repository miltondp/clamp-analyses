#!/usr/bin/env python
"""Perturbation-response capture of one model's LVs, against a gene-label permutation null.

For each Perturb-seq perturbation (build_regulons.py's logfc.tsv.gz: universe genes x perturbations
with >= min_de DE genes), Pearson r between its log2FC vector and every LV's loadings over the
universe; the perturbation's statistic is max |r| over LVs (|r|, so the LV sign does not matter).

* The perturbed gene's own row is left out of that perturbation's correlations (exactly, by removing
  its terms from the sums): the knockdown itself is not a downstream response.
* A zero-variance LV has r = 0.

Null: --n-perm permutations of the Z gene labels (rows) from --seed, the same stream as
ora_null.py, shared by every LV. The logFC side and the self-mask stay fixed.

A perturbation is captured when its permutation p, (1 + #{perms >= observed}) / (1 + n_perm), is
<= --alpha. For the null of the capture rate, every member of {observed, perms} is ranked against
the other n_perm, so under H0 the rate is about alpha.

Outputs in --out-dir: perturbations.csv (observed per perturbation: best LV, its signed r, max |r|,
null mean/95th percentile, perm p, captured), null.csv.gz (run statistics per permutation; perm 0
is observed) and summary.csv (ora_null.summarize's columns).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_regulons import log, self_cells  # noqa: E402
from ora_null import summarize  # noqa: E402
from signed_z import read_z  # noqa: E402

STATISTICS = ("median_max_abs_r", "mean_max_abs_r", "pct_captured")


class Response:
    """Correlations of fixed perturbation responses (genes x P) with many genes x LVs matrices."""

    def __init__(self, logfc, self_rows):
        """`self_rows`: each perturbation's own row in `logfc`, or -1 where it is not in the universe."""
        x = np.asarray(logfc, float)
        if not np.isfinite(x).all():
            raise ValueError("logFC has non-finite values")
        self.n_genes = x.shape[0]
        self.x = x - x.mean(axis=0)  # r is shift-invariant; centering keeps the sums well conditioned
        self.scale_x = (self.n_genes * (x ** 2).sum(axis=0))[:, None]
        self.rows = np.asarray(self_rows)
        self.masked = self.rows >= 0
        cols = np.flatnonzero(self.masked)
        self.xs = np.zeros(x.shape[1])
        self.xs[cols] = self.x[self.rows[cols], cols]
        self.n = (self.n_genes - self.masked)[:, None].astype(float)
        self.sx = (self.x.sum(axis=0) - self.xs)[:, None]
        self.sxx = ((self.x ** 2).sum(axis=0) - self.xs ** 2)[:, None]

    def r(self, z):
        """perturbations x LVs Pearson r, each perturbation without its own gene."""
        z = np.asarray(z, float)
        if not np.isfinite(z).all():
            raise ValueError("Z has non-finite values")
        scale_z = (z.shape[0] * (z ** 2).sum(axis=0))[None, :]
        z = z - z.mean(axis=0)
        zs = np.zeros((len(self.rows), z.shape[1]))
        zs[self.masked] = z[self.rows[self.masked]]
        sz = z.sum(axis=0)[None, :] - zs
        szz = (z ** 2).sum(axis=0)[None, :] - zs ** 2
        sxz = self.x.T @ z - self.xs[:, None] * zs
        num = self.n * sxz - self.sx * sz
        var_x = np.maximum(self.n * self.sxx - self.sx ** 2, 0)
        var_z = np.maximum(self.n * szz - sz ** 2, 0)
        # Zero-variance columns get r = 0. A constant column centres to rounding noise, not exact
        # zeros, so "zero" is relative to the uncentred sum of squares.
        flat = (var_x <= 1e-20 * self.scale_x) | (var_z <= 1e-20 * scale_z)
        den = np.sqrt(var_x * var_z)
        return np.divide(num, den, out=np.zeros_like(num), where=~flat & (den > 0))


def best(r):
    """max |r| over LVs, the LV index reaching it and its signed r, per perturbation."""
    idx = np.abs(r).argmax(axis=1)
    signed = r[np.arange(len(idx)), idx]
    return np.abs(signed), idx, signed


def null_max_abs_r(response, z, n_perm, seed):
    """(1 + n_perm) x perturbations max |r|; row 0 is observed, then gene-label permutations of Z."""
    rng = np.random.default_rng(seed)
    rows = [best(response.r(z))[0]]
    for _ in range(n_perm):
        rows.append(best(response.r(z[rng.permutation(z.shape[0])]))[0])
    return np.vstack(rows)


def leave_one_out_p(values):
    """Per member and perturbation, (1 + #{other members >= it}) / n_members.

    Row 0 (observed) is thereby ranked against the n_perm permutations, and each permutation against
    the other permutations plus the observed, which keeps the null rows exchangeable under H0.
    """
    count_ge = (values[None, :, :] >= values[:, None, :]).sum(axis=1)  # includes itself
    return count_ge / values.shape[0]


def run_statistics(values, alpha):
    p = leave_one_out_p(values)
    return pd.DataFrame({
        "perm": np.arange(values.shape[0]),
        "median_max_abs_r": np.median(values, axis=1),
        "mean_max_abs_r": values.mean(axis=1),
        "pct_captured": 100 * (p <= alpha).mean(axis=1),
    }), p[0]


def read_logfc(path):
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame


def run(args):
    z = read_z(args.z)
    logfc = read_logfc(args.logfc)
    if set(z.index) != set(logfc.index):
        raise ValueError(f"Z ({len(z)}) and logFC ({len(logfc)}) gene sets differ")
    logfc = logfc.loc[z.index]

    self_rows = np.full(logfc.shape[1], -1)
    rows, cols = self_cells(logfc)
    self_rows[cols] = rows
    response = Response(logfc.to_numpy(), self_rows)
    log(f"{z.shape[1]} LVs, {len(z)} genes, {logfc.shape[1]} perturbations "
        f"({int(response.masked.sum())} with their own gene masked)")

    values = z.to_numpy(float)
    max_abs, idx, signed = best(response.r(values))
    null = null_max_abs_r(response, values, args.n_perm, args.seed)
    stats, perm_p = run_statistics(null, args.alpha)

    perturbations = pd.DataFrame({
        "perturbation": logfc.columns,
        "self_masked": response.masked,
        "best_LV": np.asarray(z.columns)[idx],
        "best_r": signed,
        "max_abs_r": max_abs,
        "null_mean": null[1:].mean(axis=0),
        "null_p95": np.quantile(null[1:], 0.95, axis=0),
        "perm_p": perm_p,
        "captured": perm_p <= args.alpha,
    })
    if args.sizes:
        sizes = pd.read_csv(args.sizes, usecols=["perturbation", "n_de_universe"])
        perturbations.insert(1, "n_de", perturbations.perturbation.map(
            sizes.set_index("perturbation")["n_de_universe"]))

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    perturbations.to_csv(out / "perturbations.csv", index=False)
    stats.to_csv(out / "null.csv.gz", index=False)
    summary = summarize(stats, {
        "strategy": args.strategy, "model": args.model, "latent_variables": z.shape[1],
        "perturbations": logfc.shape[1], "self_masked": int(response.masked.sum()),
        "n_perm": args.n_perm, "seed": args.seed, "alpha": args.alpha,
    }, statistics=STATISTICS)
    summary.to_csv(out / "summary.csv", index=False)
    log(summary[["statistic", "observed", "null_mean", "null_p95", "perm_p"]].to_string(index=False))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--z", required=True, help="Z restricted to the regulon universe (signed_z.py's pos.csv)")
    parser.add_argument("--logfc", required=True, help="build_regulons.py's logfc.tsv.gz")
    parser.add_argument("--sizes", help="build_regulons.py's sizes.csv, to add each perturbation's DE count")
    parser.add_argument("--strategy", required=True)
    parser.add_argument("--model", required=True)
    # No defaults: the values live in workflow/config/perturbseq.yaml.
    parser.add_argument("--n-perm", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    run(parse_args(argv))


if __name__ == "__main__":
    main()

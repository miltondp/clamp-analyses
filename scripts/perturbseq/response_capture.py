#!/usr/bin/env python
"""Perturbation-response capture of one model's LVs, against a gene-label permutation null.

For each Perturb-seq perturbation (build_regulons.py's logfc.tsv.gz: universe genes x perturbations
with >= min_de DE genes), Pearson r between its log2FC vector and every LV's loadings over the
universe; the perturbation's statistic is max |r| over LVs (|r|, so the LV sign does not matter).

* The perturbed gene's own row is left out of that perturbation's correlations (exactly, by removing
  its terms from the sums): the knockdown itself is not a downstream response.
* A zero-variance LV has r = 0.
* K562 knockdowns share a strong common response (the first logFC component explains ~23% of the
  variance), which one hub LV per method tends to absorb. For each --remove-pcs n, the top-n shared
  components are projected out of the logFC first (residualize()); n = 0 is the raw response.

Null: --n-perm permutations of the Z gene labels (rows) from --seed, the same stream as
ora_null.py, shared by every LV. The logFC side and the self-mask stay fixed.

A perturbation is captured when its permutation p, (1 + #{perms >= observed}) / (1 + n_perm), is
<= --alpha. For the null of the capture rate, every member of {observed, perms} is ranked against
the other n_perm, so under H0 the rate is about alpha.

Run statistics: median and mean max |r|, % captured, and top_lv_share (% of perturbations whose best
LV is the single most common one: the hub concentration).

Outputs in --out-dir, stacked over --remove-pcs (column remove_pcs): perturbations.csv (observed per
perturbation: best LV, its signed r, max |r|, null mean/95th percentile, perm p, captured),
null.csv.gz (run statistics per permutation; perm 0 is observed) and summary.csv
(ora_null.summarize's columns).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_regulons import log, read_gene_level, self_cells  # noqa: E402
from ora_null import permuted, summarize  # noqa: E402
from signed_z import read_z  # noqa: E402

STATISTICS = ("median_max_abs_r", "mean_max_abs_r", "pct_captured", "top_lv_share")


class Response:
    """Correlations of fixed perturbation responses (genes x P) with many genes x LVs matrices."""

    def __init__(self, logfc, self_rows):
        """`self_rows`: each perturbation's own row in `logfc`, or -1 where it is not in the universe."""
        x = np.asarray(logfc, float)
        if not np.isfinite(x).all():
            raise ValueError("logFC has non-finite values")
        n_genes = x.shape[0]
        self.x = x - x.mean(axis=0)  # r is shift-invariant; centering keeps the sums well conditioned
        self.scale_x = (n_genes * (x ** 2).sum(axis=0))[:, None]
        self.rows = np.asarray(self_rows)
        self.masked = self.rows >= 0
        cols = np.flatnonzero(self.masked)
        self.xs = np.zeros(x.shape[1])
        self.xs[cols] = self.x[self.rows[cols], cols]
        self.n = (n_genes - self.masked)[:, None].astype(float)
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


def residualize(logfc, self_rows, n):
    """logFC without its top-n shared response components, and the variance share they explained.

    Each perturbation column is centred and its own-gene entry zeroed (so a strong knockdown does not
    steer the shared axes); U holds the top-n left singular vectors (gene space) of that matrix, and
    the result is the centred logFC minus its projection on U. Own-gene entries stay in the output
    but are masked in the correlations anyway.
    """
    x = np.asarray(logfc, float)
    if n == 0:
        return x, 0.0
    cols = np.flatnonzero(np.asarray(self_rows) >= 0)
    own = (np.asarray(self_rows)[cols], cols)
    masked = x.copy()
    masked[own] = np.nan
    xc = x - np.nanmean(masked, axis=0)  # the column mean leaves the own gene out too
    x0 = xc.copy()
    x0[own] = 0.0
    u, s, _ = np.linalg.svd(x0, full_matrices=False)
    u = u[:, :n]
    return xc - u @ (u.T @ x0), float((s[:n] ** 2).sum() / (s ** 2).sum())


def best(r):
    """max |r| over LVs, the LV index reaching it and its signed r, per perturbation."""
    idx = np.abs(r).argmax(axis=1)
    signed = r[np.arange(len(idx)), idx]
    return np.abs(signed), idx, signed


def null_best(response, z, n_perm, seed):
    """(1 + n_perm) x perturbations max |r| and best-LV index (row 0 is observed, then gene-label
    permutations of Z), and the observed best LV's signed r."""
    max_abs, idx, signed = best(response.r(z))
    values, lvs = [max_abs], [idx]
    for _, zp in permuted(z, n_perm, seed):
        max_abs, idx, _ = best(response.r(zp))
        values.append(max_abs)
        lvs.append(idx)
    return np.vstack(values), np.vstack(lvs), signed


def leave_one_out_p(values):
    """Per member and perturbation, (1 + #{other members >= it}) / n_members.

    Row 0 (observed) is thereby ranked against the n_perm permutations, and each permutation against
    the other permutations plus the observed, which keeps the null rows exchangeable under H0.
    """
    count_ge = (values[None, :, :] >= values[:, None, :]).sum(axis=1)  # includes itself
    return count_ge / values.shape[0]


def top_lv_share(lvs):
    """Per member, % of perturbations whose best LV is the most common best LV (hub concentration)."""
    return np.array([100 * np.bincount(row).max() / len(row) for row in lvs])


def run_statistics(values, lvs, alpha):
    p = leave_one_out_p(values)
    return pd.DataFrame({
        "perm": np.arange(values.shape[0]),
        "median_max_abs_r": np.median(values, axis=1),
        "mean_max_abs_r": values.mean(axis=1),
        "pct_captured": 100 * (p <= alpha).mean(axis=1),
        "top_lv_share": top_lv_share(lvs),
    }), p[0]


def capture(z, logfc, self_rows, n_pcs, args):
    """Per-perturbation table, per-permutation statistics and summary for one --remove-pcs value."""
    x, removed = residualize(logfc.to_numpy(), self_rows, n_pcs)
    response = Response(x, self_rows)
    log(f"remove_pcs={n_pcs}: removed {removed:.1%} of the logFC variance")

    null, lvs, signed = null_best(response, z.to_numpy(float), args.n_perm, args.seed)
    stats, perm_p = run_statistics(null, lvs, args.alpha)

    perturbations = pd.DataFrame({
        "remove_pcs": n_pcs,
        "perturbation": logfc.columns,
        "self_masked": response.masked,
        "best_LV": np.asarray(z.columns)[lvs[0]],
        "best_r": signed,
        "max_abs_r": null[0],
        "null_mean": null[1:].mean(axis=0),
        "null_p95": np.quantile(null[1:], 0.95, axis=0),
        "perm_p": perm_p,
        "captured": perm_p <= args.alpha,
    })
    stats.insert(0, "remove_pcs", n_pcs)
    summary = summarize(stats, {
        "strategy": args.strategy, "model": args.model, "remove_pcs": n_pcs,
        "removed_variance": removed, "latent_variables": z.shape[1],
        "perturbations": logfc.shape[1], "self_masked": int(response.masked.sum()),
        "n_perm": args.n_perm, "seed": args.seed, "alpha": args.alpha,
    }, columns=STATISTICS)
    log(summary[["statistic", "observed", "null_mean", "null_p95", "perm_p"]].to_string(index=False))
    return perturbations, stats, summary


def run(args):
    z = read_z(args.z)
    logfc = read_gene_level(args.logfc)
    if set(z.index) != set(logfc.index):
        raise ValueError(f"Z ({len(z)}) and logFC ({len(logfc)}) gene sets differ")
    logfc = logfc.loc[z.index]

    self_rows = np.full(logfc.shape[1], -1)
    rows, cols = self_cells(logfc)
    self_rows[cols] = rows
    log(f"{z.shape[1]} LVs, {len(z)} genes, {logfc.shape[1]} perturbations "
        f"({int((self_rows >= 0).sum())} with their own gene masked)")

    results = [capture(z, logfc, self_rows, n, args) for n in args.remove_pcs]
    perturbations, stats, summary = (pd.concat(parts, ignore_index=True) for parts in zip(*results))
    sizes = pd.read_csv(args.sizes, usecols=["perturbation", "n_de_universe"])
    perturbations.insert(2, "n_de", perturbations.perturbation.map(
        sizes.set_index("perturbation")["n_de_universe"]))

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    perturbations.to_csv(out / "perturbations.csv", index=False)
    stats.to_csv(out / "null.csv.gz", index=False)
    summary.to_csv(out / "summary.csv", index=False)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--z", required=True, help="Z restricted to the regulon universe (signed_z.py's pos.csv)")
    parser.add_argument("--logfc", required=True, help="build_regulons.py's logfc.tsv.gz")
    parser.add_argument("--sizes", required=True, help="build_regulons.py's sizes.csv, to add each perturbation's DE count")
    parser.add_argument("--strategy", required=True)
    parser.add_argument("--model", required=True)
    # No defaults: the values live in workflow/config/perturbseq.yaml.
    parser.add_argument("--remove-pcs", type=int, nargs="+", required=True,
                        help="numbers of shared logFC components to remove; each gives its own results")
    parser.add_argument("--n-perm", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    run(parse_args(argv))


if __name__ == "__main__":
    main()

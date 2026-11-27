#!/usr/bin/env python
"""Permutation null for the Perturb-seq regulon ORA of one model's LVs.

The observed ORA comes from `scripts/archs4/coverage/run_ora.R`, run on Z and on -Z (signed_z.py).
That runs clusterProfiler::enricher once per LV, which is too slow for 100 permutations x 2 tails.
This script reimplements the same test, vectorized over LVs (DOSE:::enricher_internal):

* query: the LV's top ceil(top_pct x universe) positive loadings (ties in row order);
* N: the annotated universe (genes in >= 1 eligible regulon), n: the query genes in it;
* M: regulon size, k: query genes in the regulon; p = P(X >= k) for X ~ Hypergeom(N, M, n);
* only regulons with k >= 1 are tested, and BH runs over those, within each LV;
* an LV hits a regulon when p.adjust < fdr (run_ora.R's recovery rule).

Validation: with the identity permutation this must reproduce run_ora.R's hits for both tails (and
its p-values, to rtol 1e-6), or the job fails.

Null: --n-perm permutations of the Z gene labels from --seed. Each permutation is shared by every
LV and both tails (and, with the same seed, by every model of a strategy), so each method keeps its
k, its loading distribution and the correlation between its LVs.

Outputs in --out-dir: null.csv.gz (one row per permutation; perm 0 is observed), lv_hits.csv
(observed hits per LV and tail) and summary.csv (observed, null mean/sd/95th percentile,
observed - null mean, and the permutation p-value of each statistic).
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import hypergeom

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_regulons import bh_columns  # noqa: E402

TAILS = ("pos", "neg")
STATISTICS = ("lv_pct_pos", "lv_pct_neg", "lv_pct_any", "regulon_pct_pos", "regulon_pct_any")


def log(message):
    print(message, file=sys.stderr, flush=True)


def read_gmt(path):
    """term -> genes, with run_ora.R::read_gmt's parsing (name, description, genes...)."""
    terms = {}
    for line in Path(path).read_text().splitlines():
        record = line.split("\t")
        if len(record) >= 3:
            terms.setdefault(record[0], []).extend(record[2:])
    return terms


class Ora:
    """Regulon ORA of many LVs at once, over a fixed gene universe (the Z rows)."""

    def __init__(self, genes, terms, top_pct, min_size, max_size, fdr):
        self.genes = pd.Index(genes)
        if self.genes.has_duplicates:
            raise ValueError("Z contains duplicate genes")
        self.fdr = fdr
        self.query_target = max(1, math.ceil(len(self.genes) * top_pct))

        rows, cols, names = [], [], []
        for term, members in terms.items():
            idx = np.unique(self.genes.get_indexer(members))
            idx = idx[idx >= 0]
            if min_size <= len(idx) <= max_size:
                rows.append(idx)
                cols.append(np.full(len(idx), len(names)))
                names.append(term)
        if not names:
            raise ValueError("No eligible regulons after universe and size filtering")
        self.terms = np.asarray(names)
        rows, cols = np.concatenate(rows), np.concatenate(cols)
        self.membership = sparse.csr_matrix(
            (np.ones(len(rows), np.int32), (rows, cols)), shape=(len(self.genes), len(names)))
        self.annotated = np.asarray(self.membership.sum(axis=1)).ravel() > 0
        self.N = int(self.annotated.sum())
        self.M = np.asarray(self.membership.sum(axis=0)).ravel().astype(int)
        # P(X >= k) for every (regulon size, n, k) a query can produce: one table, then lookups.
        self.sizes, self.size_idx = np.unique(self.M, return_inverse=True)
        q = self.query_target
        n, k = np.meshgrid(np.arange(q + 1), np.arange(q + 1), indexing="ij")
        self.sf = hypergeom.sf(k[None] - 1, self.N, self.sizes[:, None, None], n[None])

    def query(self, values):
        """Sparse LVs x genes indicator of each LV's top positive loadings."""
        values = np.where(np.isfinite(values) & (values > 0), values, -np.inf)
        order = np.argsort(-values, axis=0, kind="stable")[: self.query_target]
        lv = np.broadcast_to(np.arange(values.shape[1]), order.shape)
        keep = np.isfinite(np.take_along_axis(values, order, axis=0))
        return sparse.csr_matrix(
            (np.ones(int(keep.sum()), np.int32), (lv[keep], order[keep])), shape=values.shape[::-1])

    def test(self, values):
        """p and BH p.adjust (NaN where untested) for a genes x LVs loading matrix."""
        query = self.query(values)
        k = (query @ self.membership).toarray()
        n = np.asarray(query @ self.annotated.astype(np.int32)).ravel()
        tested = k >= 1
        p = np.full(k.shape, np.nan)
        lv, term = np.nonzero(tested)
        p[lv, term] = self.sf[self.size_idx[term], n[lv], k[lv, term]]
        return {
            "p": p,
            "p_adjust": bh_columns(p.T).T,
            "tested_lvs": int((np.asarray(query.sum(axis=1)).ravel() > 0).sum()),
        }

    def hits(self, values):
        with np.errstate(invalid="ignore"):
            return self.test(values)["p_adjust"] < self.fdr


def tail_values(z):
    return {"pos": z, "neg": -z}


def statistics(hits):
    n_lv, n_terms = hits["pos"].shape
    lv = {tail: hits[tail].any(axis=1) for tail in TAILS}
    reg = {tail: hits[tail].any(axis=0) for tail in TAILS}
    return {
        "lv_pct_pos": 100 * lv["pos"].sum() / n_lv,
        "lv_pct_neg": 100 * lv["neg"].sum() / n_lv,
        "lv_pct_any": 100 * (lv["pos"] | lv["neg"]).sum() / n_lv,
        "regulon_pct_pos": 100 * reg["pos"].sum() / n_terms,
        "regulon_pct_any": 100 * (reg["pos"] | reg["neg"]).sum() / n_terms,
    }


def validate(ora, z, tail, results, ora_dir):
    """Fail unless this implementation reproduces run_ora.R's output for one tail."""
    ora_dir = Path(ora_dir)
    summary = pd.read_csv(ora_dir / "summary.csv").iloc[0]
    expected = {
        "universe_size": len(ora.genes),
        "effective_annotated_universe_size": ora.N,
        "query_target": ora.query_target,
        "latent_variables": z.shape[1],
        "eligible_pathways": len(ora.terms),
        "tested_lvs": results["tested_lvs"],
    }
    for key, value in expected.items():
        if int(summary[key]) != value:
            raise AssertionError(f"{ora_dir}: {key} is {summary[key]} in run_ora.R, {value} here")

    # With no result rows, run_ora.R's fwrite(compress = "gzip") leaves a truncated header-only gzip.
    if int(summary["result_rows"]) == 0:
        r = pd.DataFrame({"LV": [], "ID": [], "pvalue": [], "p.adjust": []})
    else:
        r = pd.read_csv(ora_dir / "enrichment.csv.gz")
    lv_index = pd.Index(z.columns).get_indexer(r["LV"].astype(str))
    term_index = pd.Index(ora.terms).get_indexer(r["ID"].astype(str))
    if (lv_index < 0).any() or (term_index < 0).any():
        raise AssertionError(f"{ora_dir}: run_ora.R reports LVs or regulons this script does not know")
    np.testing.assert_allclose(results["p"][lv_index, term_index], r["pvalue"], rtol=1e-6,
                               err_msg=f"{ora_dir}: p-values differ from run_ora.R")
    np.testing.assert_allclose(results["p_adjust"][lv_index, term_index], r["p.adjust"], rtol=1e-6,
                               err_msg=f"{ora_dir}: p.adjust differs from run_ora.R")
    hit = (r["p.adjust"] < ora.fdr).to_numpy()
    r_hits = set(zip(lv_index[hit], term_index[hit]))
    with np.errstate(invalid="ignore"):
        ours = set(zip(*np.nonzero(results["p_adjust"] < ora.fdr)))
    if r_hits != ours:
        raise AssertionError(
            f"{ora_dir}: {len(ours - r_hits)} hits only here, {len(r_hits - ours)} only in run_ora.R")
    if int(summary["recovered_pathways"]) != len({t for _, t in ours}):
        raise AssertionError(f"{ora_dir}: recovered regulon count differs from run_ora.R")
    log(f"{tail}: reproduces run_ora.R ({len(ours)} LV-regulon hits, {len(r)} reported rows)")


def lv_hits(z, ora, observed):
    rows = []
    for tail in TAILS:
        padj = observed[tail]["p_adjust"]
        with np.errstate(invalid="ignore"):
            hit = padj < ora.fdr
        best = np.where(np.isnan(padj), np.inf, padj).argmin(axis=1)
        for i, lv in enumerate(z.columns):
            rows.append({
                "LV": lv,
                "tail": tail,
                "n_regulons": int(hit[i].sum()),
                "best_regulon": ora.terms[best[i]] if hit[i].any() else "",
                "best_p_adjust": float(padj[i, best[i]]) if hit[i].any() else np.nan,
            })
    return pd.DataFrame(rows)


def permutation_rows(ora, values, n_perm, seed):
    """Statistics for permutations 1..n_perm of the gene labels (rows) of a genes x LVs matrix."""
    rng = np.random.default_rng(seed)
    rows = []
    for perm in range(1, n_perm + 1):
        shuffled = values[rng.permutation(values.shape[0])]
        rows.append({"perm": perm, **statistics({t: ora.hits(v) for t, v in tail_values(shuffled).items()})})
    return rows


def summarize(null, n_perm):
    observed = null.loc[null.perm == 0].iloc[0]
    perms = null.loc[null.perm > 0]
    rows = []
    for stat in STATISTICS:
        values = perms[stat].to_numpy()
        rows.append({
            "statistic": stat,
            "observed": observed[stat],
            "null_mean": values.mean(),
            "null_sd": values.std(ddof=1) if len(values) > 1 else np.nan,
            "null_p95": np.quantile(values, 0.95),
            "observed_minus_null": observed[stat] - values.mean(),
            "perm_p": (1 + (values >= observed[stat]).sum()) / (1 + n_perm),
        })
    return pd.DataFrame(rows)


def run(args):
    z = pd.read_csv(args.z_pos, index_col=0, float_precision="round_trip")
    z.index = z.index.astype(str)
    z_neg = pd.read_csv(args.z_neg, index_col=0, float_precision="round_trip")
    if not (z_neg.index.astype(str).equals(z.index) and z_neg.columns.equals(z.columns)
            and np.array_equal(z_neg.to_numpy(), -z.to_numpy(), equal_nan=True)):
        raise ValueError("--z-neg is not -1 x --z-pos")
    ora = Ora(z.index, read_gmt(args.gmt), args.top_pct, args.min_size, args.max_size, args.fdr)
    log(f"{z.shape[1]} LVs, {len(ora.genes)} genes ({ora.N} annotated), {len(ora.terms)} regulons, "
        f"query {ora.query_target}")

    values = z.to_numpy(float)
    observed = {tail: ora.test(v) for tail, v in tail_values(values).items()}
    for tail, ora_dir in zip(TAILS, (args.ora_pos, args.ora_neg)):
        validate(ora, z, tail, observed[tail], ora_dir)
    with np.errstate(invalid="ignore"):
        rows = [{"perm": 0, **statistics({t: r["p_adjust"] < ora.fdr for t, r in observed.items()})}]

    null = pd.DataFrame(rows + permutation_rows(ora, values, args.n_perm, args.seed))

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    null.to_csv(out / "null.csv.gz", index=False)
    lv_hits(z, ora, observed).to_csv(out / "lv_hits.csv", index=False)
    summary = summarize(null, args.n_perm)
    for key, value in reversed(list({
        "strategy": args.strategy, "model": args.model, "latent_variables": z.shape[1],
        "eligible_regulons": len(ora.terms), "n_perm": args.n_perm, "seed": args.seed,
    }.items())):
        summary.insert(0, key, value)
    summary.to_csv(out / "summary.csv", index=False)
    log(summary[["statistic", "observed", "null_mean", "null_p95", "perm_p"]].to_string(index=False))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--z-pos", required=True)
    parser.add_argument("--z-neg", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--ora-pos", required=True, help="run_ora.R output directory for --z-pos")
    parser.add_argument("--ora-neg", required=True, help="run_ora.R output directory for --z-neg")
    parser.add_argument("--strategy", required=True)
    parser.add_argument("--model", required=True)
    # No defaults: the values live in workflow/config/perturbseq.yaml.
    parser.add_argument("--top-pct", type=float, required=True)
    parser.add_argument("--min-size", type=int, required=True)
    parser.add_argument("--max-size", type=int, required=True)
    parser.add_argument("--fdr", type=float, required=True)
    parser.add_argument("--n-perm", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    run(parse_args(argv))


if __name__ == "__main__":
    main()

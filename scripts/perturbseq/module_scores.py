#!/usr/bin/env python
"""aucodds and F-aucodds (Saelens et al. 2018) of one model's LVs against the Perturb-seq regulons.

Saelens, Cannoodt & Saeys (Nat Commun 2018) score modules against an incomplete regulatory gold
standard (their human data) with `modbindevalscorer` and `cal_faucodds`. This is a numpy port:

* Modules: each LV gives one module per tail (Z and -Z, as in the ORA), defined two ways:
  - `top`: the tail's top ceil(top_pct x universe) positive loadings (ora_null.Ora.query, the ORA
    query);
  - `zscore`: loadings above mean + zscore_sd x SD of the tail (numpy SD, ddof 0). On -Z that is
    z < mean - zscore_sd x SD. The paper's ICA rule thresholds |s| > zscore_sd x SD without
    centering; centering matters for nonnegative methods.
  Genes outside the annotated universe (genes in >= 1 eligible regulon) are dropped, and modules with
  fewer than min_size genes are removed.
* For each module x regulon: a 2x2 table over the annotated universe, OR = tp*tn / (fp*fn) (inf when
  fp*fn = 0 and tp > 0), and the Fisher right tail p = P(X >= tp) for X ~ Hypergeom(N, M, n).
* Holm-Sidak within each module, across all eligible regulons; ORs with q > q are set to 0.
* aucodds over regulons: each regulon's max OR over modules; aucodds over modules: each module's
  max OR over regulons. aucodds is the fraction of maxima with log10 OR >= c, integrated over
  c in linspace(0, log10_or_max, 100) and divided by log10_or_max (a max of 0 is below every c).
* F-aucodds: the harmonic mean of the two aucodds. The module half penalizes hub and uninformative
  modules, which the regulon half and the ORA's "% regulons recovered" do not see.

Null: --n-perm permutations of the Z gene labels from --seed (ora_null.permuted, the benchmark's one
permutation stream, so each model shares its permutations with the ORA and response nulls). The
summary adds fold = observed / null mean, the paper's normalization.

Outputs in --out-dir: null.csv.gz (one row per permutation; perm 0 is observed), module_hits.csv and
regulon_hits.csv (observed best matches) and summary.csv.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_regulons import log, read_gmt  # noqa: E402
from ora_null import Ora, permuted, summarize, tail_values  # noqa: E402
from signed_z import read_z  # noqa: E402

DEFINITIONS = ("top", "zscore")
SCORES = ("aucodds_regulons", "aucodds_modules", "faucodds")


def zscore_modules(values, sd):
    """LVs x genes boolean: loadings above each column's mean + sd x SD (non-finite never members)."""
    x = np.asarray(values, float)
    x = np.where(np.isfinite(x), x, np.nan)
    with np.errstate(invalid="ignore"):
        return (x > np.nanmean(x, axis=0) + sd * np.nanstd(x, axis=0)).T


def modules(values, ora, definition, zscore_sd):
    """Sparse LVs x genes membership for one tail of a genes x LVs matrix."""
    if definition == "top":
        return ora.query(values)
    if definition == "zscore":
        return sparse.csr_matrix(zscore_modules(values, zscore_sd), dtype=np.int32)
    raise ValueError(f"Unknown module definition {definition}")


def holm_sidak(p):
    """Holm-Sidak adjusted p-values within each row (statsmodels multipletests method 'hs')."""
    p = np.asarray(p, float)
    m = p.shape[1]
    order = np.argsort(p, axis=1, kind="stable")
    ranked = np.take_along_axis(p, order, axis=1)
    with np.errstate(divide="ignore"):
        adjusted = -np.expm1(np.log1p(-ranked) * np.arange(m, 0, -1))
    adjusted = np.minimum(np.maximum.accumulate(adjusted, axis=1), 1)
    out = np.empty_like(adjusted)
    np.put_along_axis(out, order, adjusted, axis=1)
    return out


def enrichment(member, ora):
    """2x2 statistics of each module (rows of `member`, annotated genes only) against each regulon.

    Unannotated genes are in no regulon, so they only need leaving out of the module size n.
    """
    member = sparse.csr_matrix(member, dtype=np.int32)
    tp = (member @ ora.membership).toarray()
    n = np.asarray(member @ ora.annotated.astype(np.int32)).ravel()[:, None]
    m = ora.M[None, :]
    fp, fn = m - tp, n - tp
    tn = ora.N - m - n + tp
    with np.errstate(divide="ignore", invalid="ignore"):
        odds = np.where(fp * fn > 0, tp * tn / (fp * fn).astype(float), np.where(tp > 0, np.inf, 0.0))
    p = np.empty(tp.shape)
    for size in np.unique(n):
        rows = np.flatnonzero(n.ravel() == size)
        p[rows] = ora.sf(int(size))[ora.size_idx[None, :], tp[rows]]
    return {"tp": tp, "n": n.ravel(), "odds": odds, "p": p}


def aucodds(max_odds, log10_or_max):
    """Area under the fraction of maxima with log10 OR >= c, c in linspace(0, log10_or_max, 100)."""
    max_odds = np.asarray(max_odds, float)
    if max_odds.size == 0 or log10_or_max <= 0:
        return 0.0
    cutoffs = np.linspace(0, log10_or_max, 100)
    with np.errstate(divide="ignore"):
        log_odds = np.log10(max_odds)  # OR 0 -> -inf, below every cutoff
    still = (log_odds[None, :] >= cutoffs[:, None]).mean(axis=1)
    return float(np.trapezoid(still, cutoffs) / log10_or_max)


def harmonic(a, b):
    return 0.0 if a == 0 or b == 0 else 2 / (1 / a + 1 / b)


def score_definition(values, ora, definition, args):
    """Filtered odds (modules x regulons) and module labels for one module definition."""
    blocks, labels = [], []
    for tail, v in tail_values(values).items():
        member = modules(v, ora, definition, args.zscore_sd)
        size = np.asarray(member @ ora.annotated.astype(np.int32)).ravel()
        keep = np.flatnonzero(size >= args.module_min_size)
        blocks.append(enrichment(member[keep], ora))
        labels.append(pd.DataFrame({"lv_index": keep, "tail": tail}))
    res = {key: np.concatenate([b[key] for b in blocks]) for key in ("tp", "n", "odds", "p")}
    res["q"] = holm_sidak(res["p"])
    res["odds"] = np.where(res["q"] <= args.q, res["odds"], 0.0)
    return res, pd.concat(labels, ignore_index=True)


def scores(odds, log10_or_max):
    regulons = aucodds(odds.max(axis=0, initial=0.0), log10_or_max)
    mods = aucodds(odds.max(axis=1, initial=0.0), log10_or_max)
    return {"aucodds_regulons": regulons, "aucodds_modules": mods, "faucodds": harmonic(regulons, mods)}


def statistics(values, ora, args):
    row, details = {}, {}
    for definition in DEFINITIONS:
        res, labels = score_definition(values, ora, definition, args)
        for key, value in scores(res["odds"], args.log10_or_max).items():
            row[f"{definition}_{key}"] = value
        details[definition] = (res, labels)
    return row, details


def best_match(res, axis):
    """Index of the best match along `axis` of the filtered odds, whether it is significant (odds > 0)
    and a getter for any modules x regulons array at those cells."""
    odds = res["odds"]
    idx = np.expand_dims(odds.argmax(axis=axis), axis)
    at = lambda a: np.take_along_axis(a, idx, axis=axis).squeeze(axis)  # noqa: E731
    return idx.squeeze(axis), at(odds) > 0, at


def module_hits(z, ora, details):
    frames = []
    for definition, (res, labels) in details.items():
        best, found, at = best_match(res, axis=1)
        frames.append(pd.DataFrame({
            "definition": definition,
            "LV": z.columns[labels["lv_index"].to_numpy(int)],
            "tail": labels["tail"].to_numpy(),
            "size": res["n"].astype(int),
            "n_regulons": (res["odds"] > 0).sum(axis=1),
            "best_regulon": np.where(found, ora.terms[best], ""),
            "best_odds": np.where(found, at(res["odds"]), np.nan),
            "best_q": np.where(found, at(res["q"]), np.nan),
            "best_tp": np.where(found, at(res["tp"]), 0).astype(int),
        }))
    return pd.concat(frames, ignore_index=True)


def regulon_hits(z, ora, details):
    frames = []
    for definition, (res, labels) in details.items():
        if len(labels) == 0:  # argmax over no modules
            continue
        best, found, at = best_match(res, axis=0)
        frames.append(pd.DataFrame({
            "definition": definition,
            "regulon": ora.terms,
            "size": ora.M,
            "n_modules": (res["odds"] > 0).sum(axis=0),
            "best_LV": np.where(found, z.columns[labels["lv_index"].to_numpy(int)][best], ""),
            "best_tail": np.where(found, labels["tail"].to_numpy()[best], ""),
            "best_odds": np.where(found, at(res["odds"]), np.nan),
            "best_tp": np.where(found, at(res["tp"]), 0).astype(int),
        }))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def summary_table(null, details, meta):
    frames = []
    for definition, (_, labels) in details.items():
        columns = [f"{definition}_{key}" for key in SCORES]
        frame = summarize(null, {**meta, "definition": definition, "modules": len(labels)}, columns)
        frame["statistic"] = frame["statistic"].str.removeprefix(f"{definition}_")
        frames.append(frame)
    out = pd.concat(frames, ignore_index=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["fold"] = out["observed"] / out["null_mean"]
    return out


def run(args):
    z = read_z(args.z)
    # Query size and regulon set as in the ORA (Ora's fdr only applies to Ora.test).
    ora = Ora(z.index, read_gmt(args.gmt), args.top_pct, args.min_size, args.max_size, fdr=None)
    log(f"{z.shape[1]} LVs, {len(ora.genes)} genes ({ora.N} annotated), {len(ora.terms)} regulons, "
        f"top query {ora.query_target}")

    values = z.to_numpy(float)
    observed, details = statistics(values, ora, args)
    for definition, (res, labels) in details.items():
        sizes = res["n"]
        log(f"{definition}: {len(labels)} modules"
            + (f", size median {np.median(sizes):.0f} (range {sizes.min():.0f}-{sizes.max():.0f})"
               if len(sizes) else ""))
    rows = [{"perm": 0, **observed}]
    for perm, shuffled in permuted(values, args.n_perm, args.seed):
        rows.append({"perm": perm, **statistics(shuffled, ora, args)[0]})
    null = pd.DataFrame(rows)

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    null.to_csv(out / "null.csv.gz", index=False)
    module_hits(z, ora, details).to_csv(out / "module_hits.csv", index=False)
    regulon_hits(z, ora, details).to_csv(out / "regulon_hits.csv", index=False)
    summary = summary_table(null, details, {
        "strategy": args.strategy, "model": args.model, "latent_variables": z.shape[1],
        "eligible_regulons": len(ora.terms), "n_perm": args.n_perm, "seed": args.seed,
    })
    summary.to_csv(out / "summary.csv", index=False)
    log(summary[["definition", "statistic", "observed", "null_mean", "fold", "perm_p"]].to_string(index=False))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--z", required=True, help="Z restricted to the regulon universe (signed_z.py's pos.csv)")
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--strategy", required=True)
    parser.add_argument("--model", required=True)
    # No defaults: the values live in workflow/config/perturbseq.yaml.
    parser.add_argument("--top-pct", type=float, required=True)
    parser.add_argument("--min-size", type=int, required=True, help="minimum regulon size")
    parser.add_argument("--max-size", type=int, required=True, help="maximum regulon size")
    parser.add_argument("--module-min-size", type=int, required=True)
    parser.add_argument("--zscore-sd", type=float, required=True)
    parser.add_argument("--q", type=float, required=True, help="Holm-Sidak q cutoff for keeping an OR")
    parser.add_argument("--log10-or-max", type=float, required=True)
    parser.add_argument("--n-perm", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    run(parse_args(argv))


if __name__ == "__main__":
    main()

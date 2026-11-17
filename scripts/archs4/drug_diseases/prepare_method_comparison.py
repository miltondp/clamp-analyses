#!/usr/bin/env python3
# Compare drug-disease prediction methods against the gold standard.
from __future__ import annotations

import argparse
import itertools
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import average_precision_score, roc_auc_score

METRICS = {"auroc": roc_auc_score, "auprc": average_precision_score}


def bh(pvalues: list[float]) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / (np.arange(len(values)) + 1)
    out = np.empty_like(ranked)
    out[order] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    return out


def parse_thresholds(raw: str) -> list[float]:
    return [-1.0 if tok == "all" else float(tok) for tok in raw.split(",")]


def read_method(path: Path, label: str, gold: pd.DataFrame, n_top: list[float]) -> pd.DataFrame:
    frames = []
    for h5 in sorted(path.glob("*.h5")):
        meta = pd.read_hdf(h5, "metadata")
        k = float(meta["n_top_genes"].iloc[0])
        if k not in n_top:
            continue
        pred = pd.read_hdf(h5, "prediction")
        pred["score"] = pred["score"].rank()
        pred = pred.merge(gold, on=["trait", "drug"], how="inner")
        tissue = str(meta["data"].iloc[0]).removeprefix("spredixcan-mashr-zscores-")
        tissue = re.sub(r"-projection-(archs4|recount2|gtex)$|-data$", "", tissue)
        pred["method"], pred["tissue"], pred["k"] = label, tissue, k
        frames.append(pred[["trait", "drug", "score", "true_class", "method", "tissue", "k"]])
    result = pd.concat(frames, ignore_index=True)
    if result["tissue"].nunique() != 49:
        raise ValueError(f"{label}: expected 49 tissues, got {result['tissue'].nunique()}")
    counts = result.groupby("tissue")["k"].nunique()
    if not (counts == len(n_top)).all():
        raise ValueError(f"{label}: expected thresholds {n_top} in every tissue")
    result["trait"], result["drug"] = result["trait"].astype(str), result["drug"].astype(str)
    return (result.groupby(["trait", "drug", "method", "tissue"], as_index=False)
            .agg(score=("score", "mean"), true_class=("true_class", "first")))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--method", action="append", required=True, help="LABEL=DIR=lv|gene")
    p.add_argument("--gold-standard", required=True, type=Path)
    p.add_argument("--n-top-lvs", default="5,10,25,50", type=parse_thresholds)
    p.add_argument("--n-top-genes", default="50,100,250,500", type=parse_thresholds)
    p.add_argument("--n-boot", default=10000, type=int)
    p.add_argument("--seed", default=42, type=int)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    gold = pd.read_pickle(args.gold_standard)
    specs = [m.rsplit("=", 2) for m in args.method]
    order = [lab for lab, _, _ in specs]
    pred = pd.concat([read_method(Path(d), lab, gold, args.n_top_lvs if kind == "lv" else args.n_top_genes)
                      for lab, d, kind in specs], ignore_index=True)
    per_tissue = pred.groupby(["method", "tissue"], as_index=False).apply(
        lambda x: pd.Series({"n_pairs": len(x), "n_positive": int(x.true_class.sum()),
                             "auroc": roc_auc_score(x.true_class, x.score),
                             "auprc": average_precision_score(x.true_class, x.score)}), include_groups=False)
    per_tissue["method"] = pd.Categorical(per_tissue["method"], order, ordered=True)
    per_tissue.sort_values(["method", "tissue"]).to_csv(args.output_dir / "per_tissue_metrics.csv", index=False)
    paired = pred.groupby(["trait", "drug", "method"], as_index=False).agg(score=("score", "max"), true_class=("true_class", "first"))
    wide = paired.pivot(index=["trait", "drug"], columns="method", values="score").reindex(columns=order)
    y = paired.groupby(["trait", "drug"])["true_class"].first().reindex(wide.index).astype(int).to_numpy()
    aggregate = pd.DataFrame({"method": order,
                              **{metric: [func(y, wide[m]) for m in order] for metric, func in METRICS.items()}})
    aggregate.to_csv(args.output_dir / "max_aggregate_reference.csv", index=False)
    paired.rename(columns={"method": "model"}).to_csv(args.output_dir / "prediction_scores_max.csv", index=False)
    rng = np.random.default_rng(args.seed)
    scores = {m: wide[m].to_numpy() for m in order}
    traits = wide.index.get_level_values("trait").to_numpy()
    trait_ids = np.unique(traits)
    trait_rows = {t: np.flatnonzero(traits == t) for t in trait_ids}
    boot = {m: np.empty(args.n_boot) for m in order}
    cluster_boot = {m: np.full(args.n_boot, np.nan) for m in order}
    for i in range(args.n_boot):
        idx = rng.integers(0, len(y), len(y))
        cidx = np.concatenate([trait_rows[t] for t in rng.choice(trait_ids, len(trait_ids))])
        valid = 0 < y[cidx].sum() < len(cidx)
        for m in order:
            boot[m][i] = roc_auc_score(y[idx], scores[m][idx])
            if valid:
                cluster_boot[m][i] = roc_auc_score(y[cidx], scores[m][cidx])
    dominance = pd.DataFrame([{"row_m": a, "col_m": b, "value": float(np.nanmean(cluster_boot[a] > cluster_boot[b])),
                               "value_pair": float(np.mean(boot[a] > boot[b]))}
                              for a, b in itertools.permutations(order, 2)])
    dominance.to_csv(args.output_dir / "ordering_stability.csv", index=False)
    pd.DataFrame([{"method": m, "auroc_ci_low": float(np.nanpercentile(cluster_boot[m], 2.5)),
                   "auroc_ci_high": float(np.nanpercentile(cluster_boot[m], 97.5))} for m in order]
                 ).to_csv(args.output_dir / "max_aggregate_ci.csv", index=False)
    by_method = {m: per_tissue[per_tissue.method == m].set_index("tissue")["auroc"] for m in order}
    tests = []
    for group1, group2 in itertools.combinations(order, 2):
        diff = by_method[group1] - by_method[group2].reindex(by_method[group1].index)
        tests.append({"group1": group1, "group2": group2,
                      "median_difference": float(diff.median()),
                      "higher": group1 if diff.median() > 0 else group2,
                      "p_value": wilcoxon(diff, alternative="two-sided").pvalue})
    tests = pd.DataFrame(tests); tests["q_value_bh"] = bh(tests.p_value.tolist())
    tests.to_csv(args.output_dir / "paired_tissue_tests.csv", index=False)
    (args.output_dir / "figure_statistics.txt").write_text(
        f"top_lvs={','.join(f'{k:g}' for k in args.n_top_lvs)} (mean); "
        f"top_genes={','.join(f'{k:g}' for k in args.n_top_genes)} (mean); bootstrap=disease-clustered; pairs={len(y)}; positives={int(y.sum())}; "
        f"auprc_baseline={y.mean():.4f}; diseases={len(trait_ids)}; tissues=49\n" +
        "\n".join(f"{r.method}: max_AUROC={r.auroc:.4f}, max_AUPRC={r.auprc:.4f}" for r in aggregate.itertuples()) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# Combine drug-disease predictions across tissues and methods.
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

N_TISSUES = 49

def parse_thresholds(raw: str) -> list[float]:
    return [-1.0 if tok == "all" else float(tok) for tok in raw.split(",")]


def _get_tissue(data_value: str) -> str:
    prefix = "spredixcan-mashr-zscores-"
    assert data_value.startswith(prefix), data_value
    tissue = data_value[len(prefix):]
    for suffix in (
        "-projection-archs4",
        "-projection-gtex",
        "-projection-recount2",
        "-projection",
        "-data",
    ):
        if tissue.endswith(suffix):
            return tissue[: -len(suffix)]
    raise ValueError(f"Cannot extract tissue from metadata data value: {data_value}")


def load_predictions(
    prediction_dirs: dict[str, Path],
    gold_standard: pd.DataFrame,
    method_thresholds: dict[str, list[float]],
) -> pd.DataFrame:
    expected_methods = tuple(method_thresholds)
    predictions = []
    for method_name, d in prediction_dirs.items():
        for f in sorted(d.glob("*.h5")):
            metadata = pd.read_hdf(f, key="metadata")
            if float(metadata["n_top_genes"].values[0]) not in method_thresholds[method_name]:
                continue
            prediction_data = pd.read_hdf(f, key="prediction")
            prediction_data["score"] = prediction_data["score"].rank()
            prediction_data = pd.merge(prediction_data, gold_standard, on=["trait", "drug"], how="inner")
            prediction_data["trait"] = prediction_data["trait"].astype("category")
            prediction_data["drug"] = prediction_data["drug"].astype("category")
            prediction_data = prediction_data.assign(method=method_name)
            prediction_data["method"] = pd.Categorical(prediction_data["method"], categories=expected_methods, ordered=True)
            prediction_data = prediction_data.assign(n_top_genes=metadata["n_top_genes"].values[0])
            data_value = metadata["data"].values[0]
            prediction_data = prediction_data.assign(data=data_value)
            prediction_data["data"] = prediction_data["data"].astype("category")
            prediction_data = prediction_data.assign(tissue=_get_tissue(data_value))
            predictions.append(prediction_data)

    predictions = pd.concat(predictions, ignore_index=True)
    assert not predictions.isna().any().any()

    method_counts = predictions["method"].value_counts().reindex(expected_methods)
    n_predictions = predictions[["drug", "trait"]].drop_duplicates().shape[0]
    print(f"Unique drug-disease pairs: {n_predictions}")

    for method_name, thresholds in method_thresholds.items():
        expected = N_TISSUES * len(thresholds) * n_predictions
        actual = int(method_counts.loc[method_name])
        assert actual == expected, f"{method_name}: expected {expected}, got {actual}"

    threshold_counts = predictions.groupby(["method", "n_top_genes"], observed=True).size().rename("n")
    for method_name, thresholds in method_thresholds.items():
        counts = threshold_counts.loc[method_name]
        assert sorted(float(x) for x in counts.index) == sorted(thresholds), method_name
        assert np.all(counts.values == N_TISSUES * n_predictions), method_name

    return predictions


def _reduce_mean(x: pd.DataFrame) -> pd.Series:
    return pd.Series({"score": x["score"].mean(), "true_class": x["true_class"].unique()[0]})


def _reduce_max(x: pd.DataFrame) -> pd.Series:
    return pd.Series({"score": x["score"].max(), "true_class": x["true_class"].unique()[0]})


def reduce_mean_then_max(predictions: pd.DataFrame) -> pd.DataFrame:
    return (
        predictions
        .groupby(["trait", "drug", "method", "tissue"], observed=True)
        .apply(_reduce_mean, include_groups=False)
        .dropna()
        .groupby(["trait", "drug", "method"], observed=True)
        .apply(_reduce_max, include_groups=False)
        .dropna()
        .sort_index()
        .reset_index()
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", action="append", required=True, help="LABEL=DIR=lv|gene")
    parser.add_argument("--gold-standard", required=True, type=Path)
    parser.add_argument("--n-top-lvs", required=True, help="comma-separated ints, or 'all' entries")
    parser.add_argument("--n-top-genes", required=True, help="comma-separated ints, or 'all' entries")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    specs = [m.rsplit("=", 2) for m in args.method]
    n_top = {"lv": parse_thresholds(args.n_top_lvs), "gene": parse_thresholds(args.n_top_genes)}
    prediction_dirs = {label: Path(d) for label, d, _ in specs}
    method_thresholds = {label: n_top[kind] for label, _, kind in specs}
    for name, d in prediction_dirs.items():
        assert d.exists(), f"{name} predictions missing: {d}"

    gold_standard = pd.read_pickle(args.gold_standard)
    print(f"Gold standard shape: {gold_standard.shape}")

    predictions = load_predictions(prediction_dirs, gold_standard, method_thresholds)
    print(f"Loaded predictions shape: {predictions.shape}")

    predictions.to_pickle(args.output_dir / "predictions_results.pkl")

    predictions_avg = reduce_mean_then_max(predictions)
    assert predictions_avg.shape[0] == len(method_thresholds) * predictions[["drug", "trait"]].drop_duplicates().shape[0]
    predictions_avg.to_pickle(args.output_dir / "predictions_results_aggregated.pkl")

    auroc_final = predictions_avg.groupby("method", observed=True).apply(
        lambda x: roc_auc_score(x["true_class"], x["score"]), include_groups=False
    ).rename("AUROC")
    ap_final = predictions_avg.groupby("method", observed=True).apply(
        lambda x: average_precision_score(x["true_class"], x["score"]), include_groups=False
    ).rename("AvgPrecision")
    summary = pd.concat([auroc_final, ap_final], axis=1)
    print(summary)
    summary.to_csv(args.output_dir / "performance_summary.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Build Perturb-seq regulons (one target-gene set per CRISPRi knockdown) over an LV gene universe.

Inputs are the Ota et al. 2025 gene-level K562GW matrices (targets in Ensembl rows x perturbations
named by symbol): limma log2FC and raw P versus non-targeting controls. P is BH-adjusted within each
perturbation over every measured target, then:

* targets are mapped to symbols (unmapped and ambiguous symbols are dropped);
* the universe is the Z genes (identical across every --z file) that the screen measured;
* a regulon is a perturbation's q < --q targets, minus the perturbed gene itself, within the
  universe, kept when its size is in [--min-size, --max-size];
* near-duplicate regulons (Jaccard > --jaccard; e.g. ribosomal and other essential knockdowns) are
  collapsed greedily onto the largest one.

Outputs in --out-dir: regulons.gmt (representatives), clusters.csv, sizes.csv, logfc.tsv.gz
(universe x perturbations with >= --min-de DE genes in the universe) and report.json.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def log(message):
    print(message, file=sys.stderr, flush=True)


def bh_columns(p):
    """Benjamini-Hochberg q-values for every column of a genes x perturbations p matrix.

    Each column's denominator is its count of finite p; NaN p stays NaN.
    """
    p = np.asarray(p, float)
    n = np.isfinite(p).sum(axis=0)
    order = np.argsort(p, axis=0)  # NaN sorts last
    ps = np.take_along_axis(p, order, axis=0)
    ranks = np.arange(1, p.shape[0] + 1)[:, None]
    # fmin ignores the trailing NaN, so they do not poison the running minimum.
    qs = np.fmin.accumulate((ps * n / ranks)[::-1], axis=0)[::-1]
    q = np.empty_like(p)
    np.put_along_axis(q, order, np.minimum(qs, 1), axis=0)
    return q


def read_gene_level(path, dtype=float):
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.astype(dtype)


def read_id_map(path):
    ids = pd.read_csv(path, sep="\t", header=None, names=["ensembl", "symbol"], dtype=str)
    return ids.drop_duplicates("ensembl").set_index("ensembl")["symbol"]


def target_symbols(ensembl, id_map):
    """Symbol per Ensembl row, NaN where unmapped or where the symbol is shared by several rows."""
    symbols = id_map.reindex(ensembl)
    ambiguous = symbols.duplicated(keep=False) & symbols.notna()
    return pd.Series(np.where(ambiguous, np.nan, symbols.values), index=ensembl, dtype=object), {
        "unmapped": int(symbols.isna().sum()),
        "ambiguous_rows": int(ambiguous.sum()),
    }


def read_universe(z_paths):
    """The gene set every Z file shares; fails if the files disagree."""
    gene_sets = []
    for path in z_paths:
        genes = pd.read_csv(path, usecols=[0]).iloc[:, 0].astype(str)
        gene_sets.append(set(genes))
    first = gene_sets[0]
    for path, genes in zip(z_paths[1:], gene_sets[1:]):
        if genes != first:
            raise ValueError(f"{path} has a different gene set from {z_paths[0]}")
    return first


def drop_self(de):
    """Clear each perturbation's own target row (the knockdown itself is not a regulon member)."""
    values = de.to_numpy(copy=True)
    rows = de.index.get_indexer(de.columns)
    cols = np.flatnonzero(rows >= 0)
    values[rows[cols], cols] = False
    return pd.DataFrame(values, index=de.index, columns=de.columns)


def collapse(sets, threshold):
    """Greedily cluster regulons whose Jaccard with a representative exceeds `threshold`.

    `sets` is a genes x regulons boolean frame. Regulons are visited by size (largest first, then
    name); each joins the representative it overlaps most if that Jaccard is > threshold, else it
    starts a new cluster. Returns a frame: perturbation, representative, size, jaccard_to_rep.
    """
    names = np.asarray(sets.columns)
    m = sets.values.astype(np.float32)
    size = m.sum(axis=0)
    inter = m.T @ m
    jaccard = inter / (size[:, None] + size[None, :] - inter)
    order = sorted(range(len(names)), key=lambda i: (-size[i], names[i]))
    reps, rep_of, jac_of = [], {}, {}
    for i in order:
        if reps:
            scores = jaccard[i, reps]
            best = int(np.argmax(scores))
            if scores[best] > threshold:
                rep_of[i], jac_of[i] = reps[best], float(scores[best])
                continue
        reps.append(i)
        rep_of[i], jac_of[i] = i, 1.0
    return pd.DataFrame({
        "perturbation": names[order],
        "representative": [names[rep_of[i]] for i in order],
        "size": size[order].astype(int),
        "jaccard_to_rep": [round(jac_of[i], 4) for i in order],
    })


def write_gmt(path, sets, clusters):
    cluster_size = clusters.groupby("representative").size()
    with open(path, "w") as handle:
        for rep in clusters.loc[clusters.perturbation == clusters.representative, "perturbation"]:
            genes = sets.index[sets[rep].values]
            handle.write("\t".join([rep, f"cluster_size={cluster_size[rep]}", *genes]) + "\n")


def quantiles(values):
    values = np.asarray(values)
    if not len(values):
        return {}
    return {f"q{int(q * 100)}": float(np.quantile(values, q)) for q in (0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)}


def build(args):
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    log(f"Reading {args.pvalues}")
    pvalues = read_gene_level(args.pvalues)
    log(f"Reading {args.logfc}")
    logfc = read_gene_level(args.logfc, np.float32)
    if not (pvalues.index.equals(logfc.index) and pvalues.columns.equals(logfc.columns)):
        raise ValueError("logFC and P matrices are not aligned")
    log(f"Screen: {pvalues.shape[0]} targets x {pvalues.shape[1]} perturbations")
    n_nan_p = int(np.isnan(pvalues.values).sum())

    q = pd.DataFrame(bh_columns(pvalues.values), index=pvalues.index, columns=pvalues.columns)

    symbols, id_stats = target_symbols(pvalues.index, read_id_map(args.id_map))
    keep = symbols.notna().values
    q = q.loc[keep]
    logfc = logfc.loc[keep]
    q.index = logfc.index = pd.Index(symbols[keep].values, name="gene")
    measured = set(q.index)

    z_genes = read_universe(args.z)
    universe = sorted(measured & z_genes)
    log(f"Universe: {len(universe)} genes ({len(universe) / len(measured):.1%} of mapped screen targets)")

    de_all = drop_self(q < args.q)
    own = q.index.get_indexer(q.columns)
    own_q = np.where(own >= 0, q.values[np.maximum(own, 0), np.arange(q.shape[1])], np.nan)
    own_logfc = np.where(own >= 0, logfc.values[np.maximum(own, 0), np.arange(q.shape[1])], np.nan)
    de = de_all.loc[universe]

    sizes = pd.DataFrame({
        "perturbation": q.columns,
        "n_de_screen": de_all.sum().values,
        "n_de_universe": de.sum().values,
        "self_measured": own >= 0,
        "self_de": own_q < args.q,
        "self_logfc": own_logfc,
    })
    in_range = sizes.n_de_universe.between(args.min_size, args.max_size)
    sizes["regulon"] = in_range.values
    sizes.to_csv(out / "sizes.csv", index=False)

    regulons = de.loc[:, in_range.values]
    clusters = collapse(regulons, args.jaccard)
    clusters.to_csv(out / "clusters.csv", index=False)
    write_gmt(out / "regulons.gmt", regulons, clusters)
    n_reps = int((clusters.perturbation == clusters.representative).sum())
    log(f"Regulons: {int(in_range.sum())} in [{args.min_size}, {args.max_size}] -> {n_reps} representatives")

    response = sizes.perturbation[sizes.n_de_universe >= args.min_de]
    logfc.loc[universe, response.values].to_csv(out / "logfc.tsv.gz", sep="\t", float_format="%.5g")
    log(f"Response matrix: {len(universe)} genes x {len(response)} perturbations")

    cluster_size = clusters.groupby("representative").size().sort_values(ascending=False)
    report = {
        "params": {k: v for k, v in vars(args).items() if k not in {"out_dir"}},
        "screen": {
            "targets": int(pvalues.shape[0]),
            "perturbations": int(pvalues.shape[1]),
            "nan_p": n_nan_p,
            "targets_unmapped": id_stats["unmapped"],
            "targets_ambiguous_symbol": id_stats["ambiguous_rows"],
            "targets_mapped": len(measured),
            "perturbations_self_measured": int(sizes.self_measured.sum()),
            "perturbations_self_de": int(sizes.self_de.sum()),
        },
        "universe": {
            "z_genes": len(z_genes),
            "size": len(universe),
            "fraction_of_mapped_targets": round(len(universe) / len(measured), 4),
            "fraction_of_z_genes": round(len(universe) / len(z_genes), 4),
        },
        "regulons": {
            "with_any_de": int((sizes.n_de_universe > 0).sum()),
            "below_min_size": int((sizes.n_de_universe.between(1, args.min_size - 1)).sum()),
            "above_max_size": int((sizes.n_de_universe > args.max_size).sum()),
            "in_size_range": int(in_range.sum()),
            "representatives": n_reps,
            "size_quantiles": quantiles(sizes.n_de_universe[in_range]),
            "multi_member_clusters": int((cluster_size > 1).sum()),
            "largest_clusters": {str(k): int(v) for k, v in cluster_size[cluster_size > 1].head(15).items()},
        },
        "response": {"min_de": args.min_de, "perturbations": int(len(response))},
    }
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--logfc", required=True)
    parser.add_argument("--pvalues", required=True)
    parser.add_argument("--id-map", required=True, help="Ensembl<TAB>symbol, no header")
    parser.add_argument("--z", required=True, nargs="+", help="Z.csv files defining the gene universe")
    parser.add_argument("--q", type=float, default=0.05)
    parser.add_argument("--min-size", type=int, default=10)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--jaccard", type=float, default=0.5)
    parser.add_argument("--min-de", type=int, default=20)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    build(parse_args(argv))


if __name__ == "__main__":
    main()

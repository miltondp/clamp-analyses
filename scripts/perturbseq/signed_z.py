#!/usr/bin/env python
"""Restrict an LV loading matrix to the Perturb-seq gene universe and split it into two signed tails.

`scripts/archs4/coverage/run_ora.R` uses every Z row as its universe and queries only positive
loadings, so it needs Z restricted to the screen's genes, once as Z (`pos`) and once as -Z (`neg`).
For nonnegative methods the `neg` tail has no positive loadings, and run_ora.R tests no LV in it.

Inputs: --z (genes x LVs CSV, gene in the first column), --universe (a file whose first column lists
the universe genes, such as build_regulons.py's logfc.tsv.gz) and --b (LVs x samples CSV, for the
sample count). Outputs: --out-pos, --out-neg and --out-manifest (JSON with `samples`, read by
run_ora.R's --model-manifest).
"""
import argparse
import json
from pathlib import Path

import pandas as pd


def read_z(path):
    z = pd.read_csv(path, index_col=0, float_precision="round_trip")
    z.index = z.index.astype(str)
    if z.index.has_duplicates:
        raise ValueError(f"{path} has duplicate genes")
    return z


def read_universe(path):
    sep = "\t" if ".tsv" in Path(path).name else ","
    return pd.read_csv(path, sep=sep, usecols=[0]).iloc[:, 0].astype(str)


def restrict(z, universe):
    """Z rows in the universe, in Z's own row order; fails if a universe gene is missing from Z."""
    universe = set(universe)
    missing = universe.difference(z.index)
    if missing:
        raise ValueError(f"{len(missing)} universe genes are missing from Z, e.g. {sorted(missing)[:5]}")
    return z.loc[z.index.isin(universe)]


def count_samples(b_path):
    """Number of samples in a B.csv (LVs x samples, LV names in the first column)."""
    return len(pd.read_csv(b_path, index_col=0, nrows=0).columns)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--z", required=True)
    parser.add_argument("--universe", required=True)
    parser.add_argument("--b", required=True)
    parser.add_argument("--out-pos", required=True)
    parser.add_argument("--out-neg", required=True)
    parser.add_argument("--out-manifest", required=True)
    args = parser.parse_args(argv)

    z = read_z(args.z)
    signed = restrict(z, read_universe(args.universe))
    for path in (args.out_pos, args.out_neg, args.out_manifest):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    signed.rename_axis("gene").to_csv(args.out_pos)
    (-signed).rename_axis("gene").to_csv(args.out_neg)
    manifest = {
        "z": args.z,
        "samples": count_samples(args.b),
        "z_genes": int(z.shape[0]),
        "universe_genes": int(signed.shape[0]),
        "latent_variables": int(signed.shape[1]),
    }
    Path(args.out_manifest).write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

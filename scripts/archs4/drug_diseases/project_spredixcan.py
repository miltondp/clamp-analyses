#!/usr/bin/env python3
# Project tissue gene-trait associations into canonical models.
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from drug_disease_utils import prepare_clamp_projector, prepare_z_projector, project_to_clamp


def load_clean_tissue(path: Path) -> pd.DataFrame:
    data = pd.read_pickle(path)
    n_dups = data.index.duplicated().sum()
    if n_dups > 0:
        print(f"  dropping {n_dups} duplicate gene rows")
        data = data[~data.index.duplicated(keep="first")]
    assert data.index.is_unique
    assert data.columns.is_unique
    data = data.dropna(how="any")
    assert not data.isna().any().any()
    return data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", type=Path)
    parser.add_argument("--model-name")
    parser.add_argument("--z", type=Path)
    parser.add_argument("--l2", type=float)
    parser.add_argument("--spredixcan-dir", required=True, type=Path)
    parser.add_argument("--compendium", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()


    raw_dir = args.output_dir / "raw"
    proj_dir = args.output_dir / "proj"
    raw_dir.mkdir(parents=True, exist_ok=True)
    proj_dir.mkdir(parents=True, exist_ok=True)

    if args.z is not None:
        assert args.z.exists() and args.l2 is not None, args.z
        CLAMP, clamp_sub, mapped_ensembl, lv_names = prepare_z_projector(args.z, args.l2)
    else:
        clamp_model_file = args.model_dir / f"{args.model_name}.rds"
        assert clamp_model_file.exists(), clamp_model_file
        CLAMP, clamp_sub, mapped_ensembl, lv_names = prepare_clamp_projector(clamp_model_file)
        
    input_file_list = sorted(args.spredixcan_dir.glob("spredixcan-mashr-zscores-*.pkl"))
    n_files = len(input_file_list)
    assert n_files == 49, f"Expected 49 tissue files, found {n_files}"

    for input_file in input_file_list:
        print(input_file.name)
        data = load_clean_tissue(input_file)
        print(f"  shape (no NaN, no dups): {data.shape}")

        output_raw = raw_dir / f"{input_file.stem}-data.pkl"
        data.to_pickle(output_raw)
        print(f"  saved raw to: {output_raw}")

        print("  projecting through CLAMP...")
        projection = project_to_clamp(data, CLAMP, clamp_sub, mapped_ensembl, lv_names)
        assert not projection.isna().any().any()
        print(f"    projection shape: {projection.shape}")

        output_proj = proj_dir / f"{input_file.stem}-projection-{args.compendium}.pkl"
        projection.to_pickle(output_proj)
        print(f"    saved projection to: {output_proj}")


if __name__ == "__main__":
    main()

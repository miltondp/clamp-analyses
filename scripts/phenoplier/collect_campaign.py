#!/usr/bin/env python3
"""Expose a delivered campaign unit without fitting models or rebuilding stores."""
import argparse
import csv
import json
import os
from pathlib import Path


def collect(source, dest, row):
    source = source.resolve(strict=True)
    dest = dest.resolve()
    stem = "summary/gls-summary-phenomexcan_rapid_gwas-phenomexcan-curated"
    files = {
        "curated.tsv.gz": f"{stem}.tsv.gz",
        "curated.pkl.gz": f"{stem}.pkl.gz",
        "study-phenomexcan-curated.h5": "study-phenomexcan-curated.h5",
        "curated_excluded.tsv": "summary/trait-filter-excluded-phenomexcan-curated.tsv",
    }
    # Validate all sources before publishing any links. Never replace source data.
    for name, relative in files.items():
        path = source / relative
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty delivered file: {path}")
        if (dest / name).resolve() == path.resolve():
            if not (dest / name).is_symlink():
                raise ValueError(f"Destination overlaps source: {dest / name}")
    dest.mkdir(parents=True, exist_ok=True)
    for name, relative in files.items():
        target = dest / name
        tmp = dest / f".{name}.{os.getpid()}.tmp"
        try:
            tmp.symlink_to(source / relative)
            tmp.replace(target)
        finally:
            tmp.unlink(missing_ok=True)
    provenance = dict(model=row, source=str(source), mode="recorded",
                      files={name: str(source / path) for name, path in files.items()})
    tmp = dest / f".source.{os.getpid()}.tmp"
    try:
        tmp.write_text(json.dumps(provenance, indent=2) + "\n")
        tmp.replace(dest / "source.json")
    finally:
        tmp.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--model", required=True)
    args = parser.parse_args()
    with args.manifest.open() as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["model"] == args.model]
    if len(rows) != 1:
        parser.error(f"Expected one manifest entry for {args.model}, found {len(rows)}")
    collect(args.source, args.dest, rows[0])


if __name__ == "__main__":
    main()

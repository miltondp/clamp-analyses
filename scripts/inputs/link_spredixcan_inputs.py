#!/usr/bin/env python3
"""Stage verified tissue pickles, downloading/converting public TSVs if needed."""
import argparse
from pathlib import Path
import os
import json

from spredixcan_inputs import MANIFEST, recover, verify_pickle

INVENTORY = Path(__file__).resolve().parents[2] / "references/phenoplier/spredixcan_inventory.json"


def validate(path, inventory=None):
    expected = inventory if inventory is not None else json.loads(INVENTORY.read_text())['files']
    actual = {file.name: file.stat().st_size for file in path.glob('spredixcan-mashr-zscores-*.pkl')}
    if inventory is not None:
        if actual != expected:
            raise ValueError(f'{path}: tissue filenames/sizes do not match supplied inventory')
        return
    manifest = json.loads(MANIFEST.read_text())['tissues']
    if set(actual) != set(expected):
        raise ValueError(f'{path}: tissue filenames do not match the recorded 49-file inventory')
    for name in actual:
        verify_pickle(path / name, manifest[name])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    source_root = args.source_root or (Path(os.environ['PHENOPLIER_DATA_ROOT']) if os.environ.get('PHENOPLIER_DATA_ROOT') else None)
    candidates = [] if source_root is None else [
        source_root / 'data/phenomexcan/gene_assoc/spredixcan/pkl',
        source_root / 'phenomexcan/gene_assoc/spredixcan/pkl',
    ]
    for source in candidates:
        if source.is_dir():
            validate(source)
            args.destination.mkdir(parents=True, exist_ok=True)
            for name, asset in manifest['tissues'].items():
                destination = args.destination / name
                if destination.exists():
                    verify_pickle(destination, asset)
                else:
                    destination.symlink_to((source / name).resolve())
            return
    recover(args.destination, manifest)
    validate(args.destination)


if __name__ == '__main__':
    main()

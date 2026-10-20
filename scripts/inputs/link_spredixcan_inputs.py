#!/usr/bin/env python3
"""Link the 49 processed S-PrediXcan tissue pickles from a reference workspace."""
import argparse
from pathlib import Path
import os
import json

INVENTORY = Path(__file__).resolve().parents[2] / "references/phenoplier/spredixcan_inventory.json"


def validate(path, inventory=None):
    expected = inventory or json.loads(INVENTORY.read_text())['files']
    actual = {file.name: file.stat().st_size for file in path.glob('spredixcan-mashr-zscores-*.pkl')}
    if actual != expected:
        raise ValueError(f'{path}: processed tissue filenames/sizes do not match the recorded 49-file inventory')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path)
    parser.add_argument('--destination', type=Path, required=True)
    args = parser.parse_args()
    if args.destination.exists():
        validate(args.destination)
        return
    source_root = args.source_root or (Path(os.environ['PHENOPLIER_DATA_ROOT']) if os.environ.get('PHENOPLIER_DATA_ROOT') else None)
    if source_root is None:
        raise SystemExit('Set PHENOPLIER_DATA_ROOT to an existing workspace/extracted bundle. The old raw/HDF5 archive does not directly supply these processed pickles.')
    candidates = [source_root / 'data/phenomexcan/gene_assoc/spredixcan/pkl', source_root / 'phenomexcan/gene_assoc/spredixcan/pkl']
    for source in candidates:
        if source.is_dir():
            validate(source)
            args.destination.mkdir(parents=True, exist_ok=False)
            for file in sorted(source.glob('spredixcan-mashr-zscores-*.pkl')):
                (args.destination / file.name).symlink_to(file.resolve())
            return
    raise SystemExit('Processed S-PrediXcan pickles not found in the supplied workspace/bundle.')


if __name__ == '__main__':
    main()

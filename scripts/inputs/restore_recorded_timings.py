#!/usr/bin/env python3
"""Restore author-generated timing tables without fitting any models."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def restore(dataset, output):
    manifest = json.loads((ROOT / 'references/runtime/recorded_timings.json').read_text())
    output.mkdir(parents=True, exist_ok=True)
    for entry in manifest['datasets'][dataset]:
        data = (ROOT / entry['source']).read_bytes()
        if hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError(f"Recorded timing evidence changed: {entry['source']}")
        destination = output / entry['filename']
        if destination.exists():
            if destination.read_bytes() != data:
                raise ValueError(f'Refusing to replace different timing results: {destination}')
            continue
        with tempfile.TemporaryDirectory(prefix='.recorded-', dir=output) as temporary:
            path = Path(temporary) / entry['filename']
            path.write_bytes(data)
            os.link(path, destination)
        print(f'Restored {destination}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', choices=['pseudobulk', 'gtex'], required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    restore(args.dataset, args.out_dir)


if __name__ == '__main__':
    main()

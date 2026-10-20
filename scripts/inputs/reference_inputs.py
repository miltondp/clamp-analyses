#!/usr/bin/env python3
"""Fetch/verify pinned small reference inputs; never overwrite an existing file."""
import argparse
import gzip
import hashlib
import json
import os
import shutil
from pathlib import Path
import subprocess
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / 'workflow/config/reference_inputs.json'


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def verify(path, asset):
    actual = sha256(path)
    if actual != asset['sha256']:
        raise ValueError(f"Checksum mismatch for {path}: expected {asset['sha256']}, found {actual}. Existing files are never replaced; inspect this input before proceeding.")


def fetch(key, asset, destination_root, check=False, upstream=False, offline=False, destination_override=None, source_root=None):
    destination = destination_override.resolve() if destination_override else destination_root / asset['path']
    if destination.exists():
        verify(destination, asset)
        print(f'OK {key}: {destination}')
        return
    if check:
        raise FileNotFoundError(f'Missing {destination}')
    local_source = asset.get('local_source')
    candidates = [ROOT / local_source] if local_source else []
    if source_root:
        candidates.append(source_root / asset['path'])
        candidates.extend(source_root / candidate for candidate in asset.get('source_candidates', []))
    for candidate in candidates:
        if candidate.is_file():
            verify(candidate, asset)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.symlink_to(candidate.resolve())
            print(f'LINKED {key}: {destination}')
            return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.reference-', dir=destination.parent) as scratch:
        temporary = Path(scratch) / 'input'
        recovered = False
        if asset['bundled'] and not upstream:
            with temporary.open('wb') as handle:
                result = subprocess.run(['git', 'show', f"HEAD:{asset['path']}"], cwd=ROOT, stdout=handle, stderr=subprocess.PIPE)
            recovered = result.returncode == 0
        if not recovered:
            if offline or not asset.get('url'):
                raise RuntimeError(f"{key}: no verified live download is available. Restore its bundled Git snapshot or stage a matching PhenoPLIER workspace with --source-root.")
            request = urllib.request.Request(asset['url'], headers={'User-Agent': 'CLAMP-reference-inputs/1.0'})
            with urllib.request.urlopen(request, timeout=60) as response, temporary.open('wb') as output:
                if asset.get('compression') == 'gzip':
                    with gzip.GzipFile(fileobj=response) as decompressed:
                        shutil.copyfileobj(decompressed, output, length=1024 * 1024)
                else:
                    shutil.copyfileobj(response, output, length=1024 * 1024)
        verify(temporary, asset)
        # Same-filesystem link publishes atomically and refuses to replace a
        # destination created by another process while this download was running.
        os.link(temporary, destination)
        print(f'FETCHED {key}: {destination}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=MANIFEST)
    parser.add_argument('--source-root', type=Path, help='Existing PhenoPLIER workspace/bundle root; matching files are linked read-only by convention')
    parser.add_argument('--only', nargs='+', help='Manifest keys; default: all small reference inputs')
    parser.add_argument('--output-root', type=Path, default=ROOT)
    parser.add_argument('--destination', type=Path, help='Override output path for exactly one --only key')
    parser.add_argument('--check', action='store_true', help='Check presence and SHA-256 without downloads')
    parser.add_argument('--upstream', action='store_true', help='Test upstream downloads instead of recovering bundled Git snapshots')
    parser.add_argument('--offline', action='store_true', help='Recover bundled snapshots from Git only')
    args = parser.parse_args()
    assets = json.loads(args.manifest.read_text())['assets']
    keys = args.only or list(assets)
    if args.destination and len(keys) != 1:
        parser.error('--destination requires exactly one --only key')
    for key in keys:
        if key not in assets:
            parser.error(f'Unknown reference key: {key}')
        fetch(key, assets[key], args.output_root.resolve(), args.check, args.upstream, args.offline, args.destination, args.source_root)


if __name__ == '__main__':
    main()

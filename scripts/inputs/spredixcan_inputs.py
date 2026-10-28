#!/usr/bin/env python3
"""Recover the exact S-PrediXcan matrices from the pinned public TSV archive."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tarfile
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / 'references/phenoplier/spredixcan_download.json'


def checksum(path, algorithm='sha256'):
    digest = hashlib.new(algorithm)
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(frame):
    """v1: SHA256 of pandas axis hashes followed by row-major float64 bytes."""
    import numpy as np
    import pandas as pd
    digest = hashlib.sha256()
    digest.update(pd.util.hash_pandas_object(frame.index).to_numpy(dtype='<u8').tobytes())
    digest.update(pd.util.hash_pandas_object(frame.columns).to_numpy(dtype='<u8').tobytes())
    digest.update(np.ascontiguousarray(frame.to_numpy(), dtype='<f8').tobytes())
    return digest.hexdigest()


def verify_pickle(path, asset):
    if checksum(path) == asset['pickle_sha256']:
        return
    # Pickle serialization can differ across pandas releases. Accept only the
    # identical axes and floating-point values, never a numerical tolerance.
    import pandas as pd
    frame = pd.read_pickle(path)
    if (list(frame.shape) != asset['shape'] or frame.index.name != asset['index_name']
            or frame.columns.name != asset['columns_name']
            or fingerprint(frame) != asset['semantic_sha256']):
        raise ValueError(f'S-PrediXcan matrix mismatch: {path}')


def convert(tsv, destination, asset):
    import pandas as pd
    if destination.exists():
        verify_pickle(destination, asset)
        return
    if checksum(tsv) != asset['tsv_sha256']:
        raise ValueError(f'S-PrediXcan TSV checksum mismatch: {tsv}')
    frame = pd.read_csv(tsv, sep='\t', index_col=0)
    if (list(frame.shape) != asset['shape'] or frame.index.name != asset['index_name']
            or frame.columns.name != asset['columns_name']
            or fingerprint(frame) != asset['semantic_sha256']):
        raise ValueError(f'Converted S-PrediXcan matrix mismatch: {tsv}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.spredixcan-', dir=destination.parent) as scratch:
        temporary = Path(scratch) / 'matrix.pkl'
        frame.to_pickle(temporary)
        # No replacement if another process publishes the same path first.
        os.link(temporary, destination)
    print(f'CONVERTED {destination.name}', flush=True)


def obtain_archive(path, asset):
    if path.exists():
        if path.stat().st_size != asset['bytes'] or checksum(path, 'md5') != asset['md5']:
            raise ValueError(f'Input archive checksum mismatch: {path}')
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.download-', dir=path.parent) as scratch:
        temporary = Path(scratch) / 'archive.tar'
        request = urllib.request.Request(asset['url'], headers={'User-Agent': 'CLAMP-inputs/1.0'})
        with urllib.request.urlopen(request, timeout=90) as response, temporary.open('wb') as output:
            shutil.copyfileobj(response, output, length=1024 * 1024)
        if temporary.stat().st_size != asset['bytes'] or checksum(temporary, 'md5') != asset['md5']:
            raise ValueError('Downloaded Input archive checksum mismatch')
        os.link(temporary, path)


def recover(destination, manifest, archive=None, tsv_root=None):
    assets = manifest['tissues']
    destination.mkdir(parents=True, exist_ok=True)
    # Keep completed conversions outside Snakemake's directory output, which
    # Snakemake may clean after a failed job. Hard links avoid duplicate storage.
    cache = destination.parent / 'downloads/recovered'
    pending = {}
    for filename, asset in assets.items():
        path = destination / filename
        cached = cache / filename
        if path.exists():
            verify_pickle(path, asset)
        elif cached.exists():
            verify_pickle(cached, asset)
            os.link(cached, path)
        else:
            pending[filename] = asset
    if not pending:
        return
    cache.mkdir(parents=True, exist_ok=True)
    if tsv_root:
        for filename, asset in pending.items():
            cached = cache / filename
            convert(tsv_root / Path(asset['archive_member']).name, cached, asset)
            os.link(cached, destination / filename)
        return
    archive = archive or destination.parent / 'downloads/spredixcan-mashr-zscores.tar'
    obtain_archive(archive, manifest['archive'])
    with tarfile.open(archive, 'r:') as handle:
        for filename, asset in pending.items():
            member = handle.getmember(asset['archive_member'])
            if not member.isfile() or member.size != asset['tsv_bytes']:
                raise ValueError(f'Invalid archive member: {member.name}')
            # Never extract archive paths: only copy the expected regular member
            # into a private temporary file, then publish a verified matrix.
            with tempfile.TemporaryDirectory(prefix='.spredixcan-tsv-', dir=destination.parent) as scratch:
                tsv = Path(scratch) / 'matrix.tsv.gz'
                with handle.extractfile(member) as source, tsv.open('wb') as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
                cached = cache / filename
                convert(tsv, cached, asset)
                os.link(cached, destination / filename)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, default=MANIFEST)
    parser.add_argument('--archive', type=Path, help='Reuse a checksum-verified local archive')
    parser.add_argument('--tsv-root', type=Path, help='Reuse the original downloaded tissue TSV directory')
    args = parser.parse_args()
    recover(args.destination, json.loads(args.manifest.read_text()), args.archive, args.tsv_root)


if __name__ == '__main__':
    main()

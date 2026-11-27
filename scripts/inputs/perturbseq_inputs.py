#!/usr/bin/env python3
"""Fetch the pinned Ota et al. 2025 release and stage its K562 genome-wide Perturb-seq members.

Only the manifest's members are extracted from the 7.7 GB Zenodo ZIP. With --source-root, an
already-extracted copy of the release is linked instead, after every file passes its checksum.
Existing destinations are verified and never replaced.
"""
import argparse
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tempfile
import zipfile

import reference_inputs
from reference_inputs import ROOT, link_verified
from spredixcan_inputs import obtain_archive

MANIFEST = ROOT / 'workflow/config/perturbseq_inputs.json'


def verify(path, member):
    # Size first: a mismatch fails fast without hashing a multi-GB file.
    if path.stat().st_size != member['bytes']:
        raise ValueError(f"Checksum mismatch for {path}: expected {member['bytes']} bytes. Existing files are never replaced; inspect this input before proceeding.")
    reference_inputs.verify(path, member)


def safe_member(name):
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or not path.parts:
        raise ValueError(f'Unsafe archive member: {name}')
    return name


def stage(manifest, root=ROOT, archive=None, source_root=None, check=False):
    pending = {}
    for key, member in manifest['members'].items():
        destination = root / member['path']
        if destination.exists():
            verify(destination, member)
            print(f'OK {key}: {destination}')
            continue
        # --check never links: it only reports what is missing.
        if not check and source_root and link_verified(key, member, destination, [source_root / safe_member(member['member'])]):
            continue
        pending[key] = (member, destination)
    if not pending:
        return
    if check:
        raise FileNotFoundError('Missing: ' + ', '.join(str(destination) for _, destination in pending.values()))
    archive = archive or root / manifest['archive']['path']
    obtain_archive(archive, manifest['archive'])
    with zipfile.ZipFile(archive) as handle:
        for key, (member, destination) in pending.items():
            try:
                info = handle.getinfo(safe_member(member['member']))
            except KeyError:
                raise ValueError(f"Archive {archive} has no member {member['member']}") from None
            if info.file_size != member['bytes']:
                raise ValueError(f"Unexpected size for archive member {member['member']}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(prefix='.perturbseq-', dir=destination.parent) as scratch:
                temporary = Path(scratch) / 'input'
                with handle.open(info) as source, temporary.open('wb') as target:
                    shutil.copyfileobj(source, target, length=1024 * 1024)
                verify(temporary, member)
                os.link(temporary, destination)
            print(f'EXTRACTED {key}: {destination}')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--manifest', type=Path, default=MANIFEST)
    parser.add_argument('--output-root', type=Path, default=ROOT)
    parser.add_argument('--archive', type=Path, help='Reuse a local copy of the pinned ZIP (md5-verified)')
    parser.add_argument('--source-root', type=Path, help='Already-extracted release root; verified files are symlinked')
    parser.add_argument('--check', action='store_true', help='Verify presence and checksums without downloading')
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    stage(manifest, args.output_root.resolve(), args.archive, args.source_root, args.check)


if __name__ == '__main__':
    main()

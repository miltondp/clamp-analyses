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

from spredixcan_inputs import checksum, obtain_archive

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / 'workflow/config/perturbseq_inputs.json'


def verify(path, member):
    if path.stat().st_size != member['bytes'] or checksum(path) != member['sha256']:
        raise ValueError(f'Checksum mismatch for {path}. Existing files are never replaced; inspect this input before proceeding.')


def safe_member(name):
    path = PurePosixPath(name)
    if path.is_absolute() or '..' in path.parts or not path.parts:
        raise ValueError(f'Unsafe archive member: {name}')
    return name


def stage(manifest, root=ROOT, archive=None, source_root=None, check=False):
    members = manifest['members']
    pending = {}
    for key, member in members.items():
        destination = root / member['path']
        if destination.exists():
            verify(destination, member)
            print(f'OK {key}: {destination}')
        else:
            pending[key] = member
    if not pending:
        return
    if check:
        raise FileNotFoundError('Missing: ' + ', '.join(str(root / m['path']) for m in pending.values()))
    if source_root:
        for key, member in list(pending.items()):
            candidate = source_root / safe_member(member['member'])
            if candidate.is_file():
                verify(candidate, member)
                destination = root / member['path']
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.symlink_to(candidate.resolve())
                print(f'LINKED {key}: {destination} -> {candidate.resolve()}')
                del pending[key]
    if not pending:
        return
    archive = archive or root / manifest['archive']['path']
    obtain_archive(archive, manifest['archive'])
    with zipfile.ZipFile(archive) as handle:
        for key, member in pending.items():
            try:
                info = handle.getinfo(safe_member(member['member']))
            except KeyError:
                raise ValueError(f"Archive {archive} has no member {member['member']}") from None
            if info.file_size != member['bytes']:
                raise ValueError(f"Unexpected size for archive member {member['member']}")
            destination = root / member['path']
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

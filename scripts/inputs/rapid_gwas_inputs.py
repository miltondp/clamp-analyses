#!/usr/bin/env python3
"""Download/stage the pinned 4,049 rapid-GWAS inputs in a PhenoPLIER workspace."""
import argparse
import hashlib
from pathlib import Path
import os
import shutil
import tempfile
import zipfile

from spredixcan_inputs import obtain_archive

ARCHIVE = {
    'url': 'https://zenodo.org/records/13274496/files/rapid_gwas.zip?download=1',
    'bytes': 9894462383,
    'md5': '97b5938db1c3508fdd5aecea2f555a80',
    'files': 4049,
}


def digest(handle):
    result = hashlib.sha256()
    for chunk in iter(lambda: handle.read(1024 * 1024), b''):
        result.update(chunk)
    return result.hexdigest()


def stage(workspace, archive=None, asset=ARCHIVE):
    if not (workspace / 'config.toml').is_file():
        raise ValueError('Initialize this PhenoPLIER workspace first: phenoplier workspace init')
    archive = archive or workspace / '.clamp-input-downloads/rapid_gwas.zip'
    obtain_archive(archive, asset)
    destination = workspace / 'data/smultixcan/rapid_gwas'
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as handle:
        members = [member for member in handle.infolist() if not member.is_dir()]
        names = [member.filename for member in members]
        if len(names) != asset['files'] or len(set(names)) != len(names):
            raise ValueError('Unexpected rapid-GWAS archive inventory')
        for member in members:
            path = Path(member.filename)
            if (len(path.parts) != 2 or path.parts[0] != 'rapid_gwas'
                    or not path.name.startswith('smultixcan_') or not path.name.endswith('.tsv.gz')):
                raise ValueError(f'Unexpected archive path: {member.filename}')
            output = destination / path.name
            if output.exists():
                with output.open('rb') as existing, handle.open(member) as published:
                    if output.stat().st_size != member.file_size or digest(existing) != digest(published):
                        raise ValueError(f'Existing input differs from published rapid-GWAS file: {output}')
                continue
            with tempfile.TemporaryDirectory(prefix='.rapid-gwas-', dir=destination) as scratch:
                temporary = Path(scratch) / 'input.tsv.gz'
                with handle.open(member) as source, temporary.open('wb') as target:
                    shutil.copyfileobj(source, target, length=1024 * 1024)
                os.link(temporary, output)
    print(f'Verified {len(members)} rapid-GWAS inputs in {destination}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--archive', type=Path, help='Reuse the pinned local rapid_gwas.zip')
    args = parser.parse_args()
    stage(args.workspace.resolve(), args.archive)


if __name__ == '__main__':
    main()

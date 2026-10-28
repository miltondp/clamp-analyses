#!/usr/bin/env python3
"""Extract the original recount2 project metadata without replacing a destination."""
import argparse
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile


def extract(archive, destination):
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix='.recount2-') as tmp:
        root = Path(tmp) / 'metadata'
        root.mkdir()
        count = 0
        with tarfile.open(archive, 'r|gz') as source:
            for member in source:
                path = PurePosixPath(member.name)
                if path.is_absolute() or '..' in path.parts or member.issym() or member.islnk():
                    raise ValueError(f'Unsafe archive member: {member.name}')
                if member.isdir():
                    continue
                if not member.isfile():
                    raise ValueError(f'Unsupported archive member: {member.name}')
                if path.parts[:2] != ('recount2', 'data'):
                    raise ValueError(f'Unexpected archive layout: {member.name}')
                target = root.joinpath(*path.parts[2:])
                target.parent.mkdir(parents=True, exist_ok=True)
                with source.extractfile(member) as incoming, target.open('xb') as outgoing:
                    shutil.copyfileobj(incoming, outgoing)
                count += target.name == 'rse_gene.Rdata'
        if not count:
            raise ValueError('Archive contains no per-project rse_gene.Rdata files')
        if destination.exists():
            raise FileExistsError(destination)
        root.rename(destination)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive')
    parser.add_argument('destination')
    args = parser.parse_args()
    extract(args.archive, args.destination)

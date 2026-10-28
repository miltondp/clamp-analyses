#!/usr/bin/env python3
"""Acquire pinned raw and original lab-prepared pseudobulk inputs."""
import argparse
import io
import json
import os
from pathlib import Path
import tempfile
import urllib.request
import zipfile

from reference_inputs import ROOT, fetch, verify

MANIFEST = ROOT / 'workflow/config/pseudobulk_inputs.json'


class RemoteZipFile(io.RawIOBase):
    """Seekable HTTP Range reader: selected members, never a whole-archive fallback."""
    def __init__(self, url, size):
        self.url = url
        self.size = size
        self.position = 0
        self.cache_start = 0
        self.cache = b''

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=0):
        self.position = offset if whence == 0 else self.position + offset if whence == 1 else self.size + offset
        if self.position < 0:
            raise ValueError('Negative ZIP seek')
        return self.position

    def read(self, size=-1):
        if size < 0:
            size = self.size - self.position
        size = min(size, self.size - self.position)
        if size <= 0:
            return b''
        offset = self.position - self.cache_start
        if 0 <= offset and offset + size <= len(self.cache):
            result = self.cache[offset:offset + size]
        else:
            # Read ahead for sequential compressed members, with bounded buffers.
            count = min(max(size, 4 * 1024 * 1024), self.size - self.position)
            start, stop = self.position, self.position + count - 1
            request = urllib.request.Request(self.url, headers={
                'Range': f'bytes={start}-{stop}', 'User-Agent': 'CLAMP-pseudobulk-inputs/1.0'})
            with urllib.request.urlopen(request, timeout=90) as response:
                if response.status != 206 or response.headers.get('Content-Range') != f'bytes {start}-{stop}/{self.size}':
                    raise ValueError('Archive server did not honor the requested HTTP byte range')
                self.cache = response.read(count + 1)
            if len(self.cache) != count:
                raise ValueError('Incomplete archive range response')
            self.cache_start = start
            result = self.cache[:size]
        self.position += len(result)
        return result


def fetch_member(key, asset, destination_root, check=False):
    destination = destination_root / asset['path']
    if destination.exists():
        verify(destination, asset)
        print(f'OK {key}: {destination}')
        return
    if check:
        raise FileNotFoundError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    remote = RemoteZipFile(asset['url'], asset['archive_bytes'])
    with tempfile.TemporaryDirectory(prefix='.pseudobulk-', dir=destination.parent) as scratch:
        temporary = Path(scratch) / 'input.csv'
        with zipfile.ZipFile(remote) as archive:
            member = archive.getinfo(asset['archive_member'])
            if member.is_dir() or member.file_size != asset['bytes']:
                raise ValueError(f'Unexpected prepared input member: {member.filename}')
            with archive.open(member) as source, temporary.open('wb') as output:
                while chunk := source.read(1024 * 1024):
                    output.write(chunk)
        verify(temporary, asset)
        os.link(temporary, destination)
    print(f'FETCHED {key}: {destination}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=MANIFEST)
    parser.add_argument('--only', nargs='+')
    parser.add_argument('--output-root', type=Path, default=ROOT)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    assets = json.loads(args.manifest.read_text())['assets']
    for key in args.only or assets:
        if key not in assets:
            parser.error(f'Unknown input: {key}')
        asset = assets[key]
        if asset.get('archive_member'):
            fetch_member(key, asset, args.output_root.resolve(), args.check)
        else:
            fetch(key, asset, args.output_root.resolve(), check=args.check)


if __name__ == '__main__':
    main()

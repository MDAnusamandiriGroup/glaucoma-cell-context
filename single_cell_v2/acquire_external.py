#!/usr/bin/env python3
"""Retrieve only the RNA members of the public GEO tar by validated byte ranges."""
from pathlib import Path
import urllib.request
import tarfile
import hashlib
import json
import concurrent.futures
import time

ROOT = Path(__file__).resolve().parent / 'data/raw'
TAR_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196235/suppl/GSE196235_RAW.tar'


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def ranged(start, end):
    req = urllib.request.Request(TAR_URL + '?offset=' + str(start), headers={'Range': f'bytes={start}-{end}'})
    response = urllib.request.urlopen(req, timeout=60)
    expected = f'bytes {start}-{end}/10720030720'
    if response.status != 206 or response.headers.get('Content-Range') != expected:
        response.close()
        raise ValueError('Unverified tar byte range')
    return response


def main():
    offset = 0; selected = []
    for line in (ROOT / 'wang_filelist.txt').read_text().splitlines():
        fields = line.split('\t')
        if fields[0] != 'File':
            continue
        name, size = fields[1], int(fields[3])
        if name.startswith('GSM5866081_'):
            selected.append((name, size, offset))
        offset += 512 + ((size + 511) // 512) * 512
    if len(selected) != 3:
        raise ValueError('Expected three RNA matrix members')
    manifests = []
    for name, size, offset in selected:
        with ranged(offset, offset + 511) as response:
            header = response.read()
        member = tarfile.TarInfo.frombuf(header, encoding='utf-8', errors='strict')
        if member.name != name or member.size != size:
            raise ValueError('Tar member disagrees with public file list')
        target = ROOT / name
        if not target.exists():
            partial = target.with_name(target.name + '.partial')
            if partial.exists():
                raise ValueError('Interrupted download requires explicit recovery: ' + name)
            with ranged(offset + 512, offset + 511 + size) as response, partial.open('xb') as handle:
                copied = 0
                while True:
                    block = response.read(4 * 1024 * 1024)
                    if not block:
                        break
                    copied += len(block); handle.write(block)
                    if copied % (32 * 1024 * 1024) == 0:
                        print(name + ': ' + str(copied // (1024 * 1024)) + ' MiB', flush=True)
            if copied != size:
                raise ValueError('Truncated download')
            partial.rename(target)
        if target.stat().st_size != size:
            raise ValueError('Downloaded size mismatch')
        manifests.append({'name': name, 'size': size, 'tar_header_offset': offset,
            'sha256': digest(target), 'source': TAR_URL,
            'validation': 'HTTP 206 Content-Range, tar checksum/header name and size, gzip and matrix checked downstream'})
        print('DOWNLOADED ' + name, flush=True)
    (ROOT / 'wang_RNA_acquisition.json').write_text(json.dumps({'files': manifests}, indent=2) + '\n')


if __name__ == '__main__':
    main()

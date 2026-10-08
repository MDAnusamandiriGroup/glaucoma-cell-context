#!/usr/bin/env python3
"""Fetch frozen, checksummed public source files; retain incompatible files."""
from pathlib import Path
import json
import hashlib
import os
import urllib.request
import zipfile
import anndata as ad


def sha256(path):
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def main():
    root=Path(__file__).resolve().parent
    raw=root/'data/raw';raw.mkdir(parents=True,exist_ok=True)
    registry=json.loads((root/'data/source_manifest.json').read_text())
    for source in registry['files']:
        if source['method']!='direct_https':continue
        target=raw/source['name']
        if target.exists():
            if sha256(target)!=source['sha256']:raise ValueError('Preserving changed source: '+str(target))
            print('VERIFIED '+source['name'],flush=True);continue
        temporary=target.with_name(target.name+'.download.partial')
        if temporary.exists():raise ValueError('Preserving interrupted download: '+str(temporary))
        with urllib.request.urlopen(source['url'],timeout=60) as response,temporary.open('xb') as handle:
            while True:
                block=response.read(4*1024*1024)
                if not block:break
                handle.write(block)
            handle.flush();os.fsync(handle.fileno())
        if sha256(temporary)!=source['sha256']:raise ValueError('Remote source changed: '+source['name'])
        temporary.rename(target);print('DOWNLOADED '+source['name'],flush=True)
    snapshot=raw/'ensembl_index_variants_GRCh38.json'
    if not snapshot.exists():
        snapshot.write_bytes((root/'inputs/ensembl_index_variants_GRCh38.json').read_bytes())
    # GEO publishes the three RNA members within a 10 GiB tar. Validated byte
    # ranges retrieve only 268 MiB of required members, never the full archive.
    import acquire_external
    acquire_external.main()
    archive=raw/'yuan_sciadv.adv9162_tables_s1_to_s15.zip'
    with zipfile.ZipFile(archive) as handle:
        target=raw/'yuan_table_S6.tsv'
        if not target.exists():target.write_bytes(handle.read('Tables_revision1/table S6.tsv'))
    data=ad.read_h5ad(raw/'wang_cellxgene.h5ad',backed='r')
    metadata=raw/'wang_metadata.csv'
    if not metadata.exists():data.obs.to_csv(metadata,index_label='cell_id')
    data.file.close()
    for source in registry['files']:
        target=raw/source['name']
        if not target.exists() or sha256(target)!=source['sha256']:
            raise ValueError('Frozen acquisition contract changed: '+source['name'])
    print('ALL V2 INPUTS VERIFIED',flush=True)


if __name__=='__main__':main()

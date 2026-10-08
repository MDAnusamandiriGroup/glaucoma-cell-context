from contextlib import contextmanager
from pathlib import Path
import csv
import gzip
import io
import os
import re
import tarfile
import time
import urllib.request

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .common import file_hash, publish, save_adata, save_csv, save_json


def acquire(root, config, offline=False):
    raw = Path(root) / 'data/raw'
    raw.mkdir(parents=True, exist_ok=True)
    sources = {}
    for source in config['sources']:
        path = raw / source['file']
        if path.exists():
            if file_hash(path, 'md5') != source['md5']:
                raise ValueError('Input checksum mismatch; preserving input: ' + str(path))
            print('REUSE input: ' + path.name, flush=True)
        else:
            if offline:
                raise FileNotFoundError('Offline input is missing: ' + str(path))
            def download(temp):
                request = urllib.request.Request(source['url'], headers={'User-Agent': 'MilaResearchPortfolio/1.0'})
                with urllib.request.urlopen(request, timeout=60) as response, temp.open('wb') as handle:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        handle.write(block)
                if file_hash(temp, 'md5') != source['md5']:
                    raise ValueError('Downloaded input failed the published checksum: ' + source['file'])
            print('DOWNLOAD input: ' + path.name, flush=True)
            publish(path, download)
        sources[path.name] = {'sha256': file_hash(path), 'md5': source['md5'], 'bytes': path.stat().st_size, 'url': source['url']}
    return sources


@contextmanager
def matrix_stream(path):
    """The source .csv.gz is actually a gzip-compressed TAR with one CSV."""
    with gzip.open(path, 'rb') as probe:
        header = probe.read(512)
    if len(header) >= 262 and header[257:262] == b'ustar':
        # TextIOWrapper requires a seekable() capability; ExFileObject from
        # tarfile's pipe mode lacks it. r:gz still reads the CSV incrementally
        # and never extracts the 881 MB member onto the filesystem.
        with tarfile.open(path, mode='r:gz') as archive:
            selected = False
            for member in archive:
                if member.isfile() and Path(member.name).name == 'lukowski_embo2019_raw_count_matrix.csv':
                    if selected:
                        raise ValueError('Multiple count-matrix entries in source archive')
                    selected = True
                    binary = archive.extractfile(member)
                    with io.TextIOWrapper(binary, encoding='utf-8') as handle:
                        yield handle
                    return
            if not selected:
                raise ValueError('Expected CSV member is absent from source archive')
    else:
        with gzip.open(path, 'rt', encoding='utf-8') as handle:
            yield handle


def read_sparse_counts(path):
    genes = []
    arrays = []
    indices = []
    indptr = [0]
    with matrix_stream(path) as handle:
        cells = next(csv.reader([handle.readline().rstrip('\r\n')]))
        if cells and cells[0] == '':
            cells = cells[1:]
        if len(cells) != len(set(cells)):
            raise ValueError('Duplicate cell barcodes in source count matrix')
        n_cells = len(cells)
        for row_number, line in enumerate(handle, start=1):
            gene, separator, values = line.rstrip('\r\n').partition(',')
            if separator != ',' or not gene:
                raise ValueError('Invalid count row: ' + str(row_number))
            counts = np.fromstring(values, sep=',', dtype=np.int32)
            if counts.size != n_cells or np.any(counts < 0):
                raise ValueError('Invalid count vector at row ' + str(row_number))
            nonzero = np.flatnonzero(counts).astype(np.int32)
            genes.append(gene.strip('"'))
            arrays.append(counts[nonzero])
            indices.append(nonzero)
            indptr.append(indptr[-1] + len(nonzero))
            if row_number % 5000 == 0:
                print('INGEST genes: ' + str(row_number), flush=True)
    if len(genes) != len(set(genes)):
        raise ValueError('Duplicate gene symbols: explicit resolution is required')
    gene_by_cell = sparse.csr_matrix((np.concatenate(arrays), np.concatenate(indices), np.array(indptr)), shape=(len(genes), n_cells))
    return gene_by_cell.T.tocsr(), genes, cells


def broad_label(label):
    label = re.sub(r' C(?:CA)?\d+$', '', str(label))
    mapping = {'Rod PR': 'Rod', 'Cone PR': 'Cone', 'MG': 'Muller glia', 'Astrocytes': 'Astrocyte', 'Others': 'Unassigned'}
    return mapping.get(label, label)


def attach_metadata(cells, meta, labels):
    if meta['cell.bc'].duplicated().any() or labels['cell.bc'].duplicated().any():
        raise ValueError('Duplicate cell IDs in source metadata')
    meta = meta.set_index('cell.bc')
    labels = labels.set_index('cell.bc')
    if set(cells) != set(meta.index) or set(cells) != set(labels.index):
        raise ValueError('Count-matrix and metadata cells do not match exactly')
    meta = meta.loc[cells]
    labels = labels.loc[cells]
    obs = pd.DataFrame(index=pd.Index(cells, name='cell_barcode'))
    obs['donor'] = meta['sample.group'].astype(str).values
    obs['library'] = meta['batch'].astype(str).values
    obs['published_cell_type'] = labels['cell.id.orig'].map(broad_label).values
    obs['published_cca_type'] = labels['cell.id.cca'].map(broad_label).values
    obs['published_original_cluster'] = labels['cluster.id.orig'].astype(str).values
    obs['published_cca_cluster'] = labels['cluster.id.cca'].astype(str).values
    obs['donor_age'] = meta['Age'].astype(int).values
    obs['donor_sex'] = meta['Gender'].astype(str).values
    obs['collection_hours'] = meta['Collection'].astype(float).values
    obs['source_n_counts'] = meta['nCount_RNA'].astype(int).values
    obs['source_n_genes'] = meta['nFeature_RNA'].astype(int).values
    obs['source_mito_percent'] = meta['percent.mito'].astype(float).values * 100.0
    if obs['donor'].nunique() != 3 or obs['library'].nunique() != 5:
        raise ValueError('Unexpected donor/library structure')
    return obs


def run(root, directory, config):
    raw = Path(root) / 'data/raw'
    counts, genes, cells = read_sparse_counts(raw / 'lukowski_embo2019_raw_count_matrix.csv.gz')
    meta = pd.read_csv(raw / 'lukowski_embo2019_CCA_metadata.csv')
    labels = pd.read_csv(raw / 'lukowski_embo2019_cellbc_cellid.csv')
    obs = attach_metadata(cells, meta, labels)
    data = ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name='gene_symbol')))
    totals = np.asarray(counts.sum(axis=1)).ravel()
    detected = counts.getnnz(axis=1)
    count_disagreement = totals != obs['source_n_counts'].to_numpy()
    gene_disagreement = detected != obs['source_n_genes'].to_numpy()
    # Validate alignment using independently published count/QC metadata.
    if count_disagreement.any() or gene_disagreement.any():
        raise ValueError('Counts or detected genes disagree with the corresponding published per-cell QC metadata')
    if counts.data.dtype.kind not in 'iu' or np.any(counts.data < 0):
        raise ValueError('Non-integer or negative UMI counts')
    save_adata(directory / 'counts.h5ad', data)
    save_csv(directory / 'donor_cell_type_counts.csv', obs.groupby(['donor', 'published_cell_type'], observed=True).size().rename('n_cells').reset_index())
    save_json(directory / 'summary.json', {
        'cells': int(data.n_obs), 'genes': int(data.n_vars), 'nonzero_counts': int(counts.nnz),
        'donors': sorted(obs['donor'].unique().tolist()), 'libraries': sorted(obs['library'].unique().tolist()),
        'counts_match_published_QC': True, 'detected_genes_match_published_QC': True,
        'source_container': 'gzip-compressed TAR containing CSV, despite .csv.gz filename',
        'input_status': 'Published filtered UMI count matrix; not FASTQ or an unfiltered droplet matrix',
        'annotation_status': 'Published original and CCA labels are retained explicitly'
    })

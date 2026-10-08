from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from .common import save_adata, save_csv, save_json


def qc_mask(obs, rule):
    return (obs['n_genes_by_counts'].to_numpy() >= rule['min_genes']) & (obs['pct_counts_mt'].to_numpy() < rule['max_mito_percent'])


def run(ingest_dir, directory, config):
    data = ad.read_h5ad(Path(ingest_dir) / 'counts.h5ad')
    data.var['mt'] = data.var_names.str.startswith('MT-')
    if int(data.var['mt'].sum()) < 10:
        raise ValueError('Mitochondrial gene identification is unexpectedly incomplete')
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    data.obs['qc_primary'] = qc_mask(data.obs, config['qc']['primary'])
    data.obs['qc_strict'] = qc_mask(data.obs, config['qc']['strict'])
    audit = data.obs.reset_index()
    save_csv(directory / 'cell_qc_audit.csv', audit)
    donor_counts = data.obs.groupby('donor', observed=True)[['qc_primary', 'qc_strict']].sum().reset_index()
    donor_counts['input_cells'] = donor_counts['donor'].map(data.obs['donor'].value_counts())
    save_csv(directory / 'qc_by_donor.csv', donor_counts)
    if int(data.obs['qc_primary'].sum()) < 1000 or data.obs.loc[data.obs['qc_primary'], 'donor'].nunique() != 3:
        raise ValueError('Primary QC leaves insufficient cells or donors')
    normalized = data[data.obs['qc_primary'].to_numpy()].copy()
    normalized.layers['counts'] = normalized.X.copy()
    normalized.X = normalized.X.astype(np.float32)
    sc.pp.normalize_total(normalized, target_sum=config['normalization_target'])
    sc.pp.log1p(normalized)
    normalized.var['detected_cells'] = normalized.layers['counts'].getnnz(axis=0)
    normalized.var['model_eligible'] = (normalized.var['detected_cells'] >= config['hvg_min_cells']) & ~normalized.var['mt']
    if not np.isfinite(normalized.X.data).all():
        raise ValueError('Non-finite normalized expression')
    save_adata(directory / 'normalized.h5ad', normalized)
    save_json(directory / 'summary.json', {
        'input_cells': int(data.n_obs), 'primary_cells': int(normalized.n_obs),
        'strict_cells': int(data.obs['qc_strict'].sum()), 'mitochondrial_genes': int(data.var['mt'].sum()),
        'model_eligible_genes': int(normalized.var['model_eligible'].sum()),
        'normalization': 'Total-count scaling to 10,000 followed by log1p; raw UMI layer retained',
        'doublet_status': 'No newly fitted doublet classifier; published filtering and Unassigned labels retained',
        'rules': config['qc']
    })

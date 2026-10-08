from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from .common import save_adata, save_csv, save_json


def pca(preprocess_dir, directory, config):
    normalized = ad.read_h5ad(Path(preprocess_dir) / 'normalized.h5ad')
    eligible = normalized[:, normalized.var['model_eligible'].to_numpy()].copy()
    eligible.layers.clear()
    sc.pp.highly_variable_genes(eligible, n_top_genes=config['hvg_n_genes'], batch_key=config['hvg_batch_key'], flavor='seurat')
    save_csv(directory / 'highly_variable_gene_audit.csv', eligible.var.reset_index())
    selected = eligible.var['highly_variable'].to_numpy()
    model = ad.AnnData(X=eligible.X[:, selected].copy(), obs=eligible.obs.copy(), var=eligible.var.loc[selected].copy())
    del eligible, normalized
    sc.pp.pca(model, n_comps=config['n_pcs'], zero_center=True, svd_solver='arpack', random_state=config['seed'])
    if not np.isfinite(model.obsm['X_pca']).all():
        raise ValueError('Invalid PCA embedding')
    save_adata(directory / 'pca.h5ad', model)
    save_json(directory / 'summary.json', {
        'cells': int(model.n_obs), 'highly_variable_genes': int(model.n_vars), 'PCs': int(model.obsm['X_pca'].shape[1]),
        'variance_ratio': model.uns['pca']['variance_ratio'].astype(float).tolist(),
        'batch_handling': 'Library-aware highly-variable-gene selection; no batch-corrected expression or integrated embedding'
    })


def graph(pca_dir, directory, config):
    model = ad.read_h5ad(Path(pca_dir) / 'pca.h5ad')
    sc.pp.neighbors(model, n_neighbors=config['n_neighbors'], n_pcs=config['n_pcs'], random_state=config['seed'])
    for resolution in config['leiden_resolutions']:
        key = 'leiden_' + str(resolution)
        sc.tl.leiden(model, resolution=resolution, key_added=key, flavor='igraph', directed=False, n_iterations=2, random_state=config['seed'])
    save_adata(directory / 'graph.h5ad', model)
    save_json(directory / 'summary.json', {
        'clusters': {str(resolution): int(model.obs['leiden_' + str(resolution)].nunique()) for resolution in config['leiden_resolutions']},
        'seed': config['seed'], 'neighbors': config['n_neighbors'], 'leiden_implementation': 'igraph, undirected, two iterations'
    })


def embedding(graph_dir, directory, config):
    model = ad.read_h5ad(Path(graph_dir) / 'graph.h5ad')
    sc.tl.umap(model, min_dist=config['umap_min_dist'], maxiter=config['umap_iterations'], random_state=config['seed'])
    if not np.isfinite(model.obsm['X_umap']).all() or model.obsm['X_umap'].shape != (model.n_obs, 2):
        raise ValueError('Invalid UMAP embedding')
    frame = model.obs.reset_index()
    frame['UMAP1'] = model.obsm['X_umap'][:, 0]
    frame['UMAP2'] = model.obsm['X_umap'][:, 1]
    save_csv(directory / 'cell_embeddings.csv', frame)
    save_adata(directory / 'embedding.h5ad', model)
    save_json(directory / 'summary.json', {
        'cells': int(model.n_obs), 'dimensions': 2,
        'note': 'UMAP is visualization; it is not evidence of gene causality, ageing effects or therapeutic efficacy'
    })

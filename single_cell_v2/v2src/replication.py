"""Cross-study RNA cell-context concordance, using published labels and donors."""
import gzip
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits
from .cache import save_json, winners


TYPE_MAP = {'Rod': 'Rod', 'Cone': 'Cone', 'OFF-cone bipolar': 'Bipolar',
    'ON-cone bipolar': 'Bipolar', 'Rod bipolar': 'Bipolar', 'GABA-amacrine': 'Amacrine',
    'Gly-amacrine': 'Amacrine', 'AII-amacrine': 'Amacrine', 'Horizontal': 'Horizontal',
    'Retinal ganglion cell': 'RGC', 'Muller glia': 'Muller glia',
    'Astrocyte': 'Astrocyte', 'Microglia': 'Microglia'}


def run(root, stages, directory, config):
    raw_root = root / 'data/raw'
    metadata = pd.read_csv(raw_root / 'wang_metadata.csv', index_col='cell_id')
    barcodes = pd.read_csv(raw_root / 'GSM5866081_barcodes.tsv.gz', header=None)[0].astype(str)
    features = pd.read_csv(raw_root / 'GSM5866081_features.tsv.gz', header=None, sep='\t')
    if metadata.index.duplicated().any() or barcodes.duplicated().any() or features[0].duplicated().any():
        raise ValueError('Duplicate external identifiers')
    if len(metadata) != 51645 or set(metadata.index) != set(barcodes):
        raise ValueError('Published labels do not exactly cover GEO barcodes')
    metadata = metadata.loc[barcodes].copy()
    if set(metadata['cell_type__custom']) != set(TYPE_MAP):
        raise ValueError('Unexpected external annotation labels')
    if set(metadata['donor_id']) != {'LGS1', 'LGS2', 'LGS3', 'LVG1'} or not metadata['disease'].eq('normal').all():
        raise ValueError('Expected four published healthy donors')
    metadata['harmonized_type'] = metadata['cell_type__custom'].map(TYPE_MAP)
    donors = sorted(metadata['donor_id'].unique())
    donor_details = metadata.groupby('donor_id', observed=True).agg(
        n_cells=('harmonized_type', 'size'), age=('age', 'first'), sex=('sex', 'first'),
        n_eyes=('biosample_id', 'nunique')).reset_index()
    donor_details.to_csv(directory / 'external_donor_metadata.csv', index=False)
    counts_by_type = pd.crosstab(metadata['harmonized_type'], metadata['donor_id']).reindex(config['all_types'], fill_value=0)
    counts_by_type.to_csv(directory / 'external_donor_type_counts.csv', index_label='cell_type')
    for key, threshold in [('primary_types', config['primary_min_cells']), ('extended_types', config['extended_min_cells'])]:
        if (counts_by_type.loc[config[key]] < threshold).any().any():
            raise ValueError('External fixed comparison set lacks eligible donor coverage: ' + key)
    candidates = pd.read_csv(stages['expanded'] / 'candidate_mapping.csv')
    # This GEO distribution repeats gene symbols in both feature columns.
    # They must not be interpreted as Ensembl IDs.
    feature_ids = features[0].astype(str).to_numpy()
    feature_symbols = features[1].astype(str).to_numpy()
    id_to_index = {g: i for i, g in enumerate(feature_ids)}
    positions = np.asarray([id_to_index.get(g, -1) for g in candidates['gene']])
    mapped = positions >= 0
    audit = candidates[['gene_id', 'gene', 'exact_symbol_in_matrix']].copy()
    audit['external_exact_symbol'] = mapped
    audit['external_symbol'] = [feature_symbols[i] if i >= 0 else '' for i in positions]
    audit['external_symbol_agrees'] = audit['external_symbol'].eq(audit['gene']) & mapped
    audit['mapping_policy'] = 'Exact gene symbol in the original GEO features; unique symbols; no inferred aliases or genomic proximity. Curated Ensembl IDs are used separately for count alignment checks.'
    audit.to_csv(directory / 'external_mapping_audit.csv', index=False)
    print('EXTERNAL ingest GEO matrix: 89,102,830 integer nonzeros', flush=True)
    with threadpool_limits(limits=4):
        with gzip.open(raw_root / 'GSM5866081_matrix.mtx.gz', 'rb') as handle:
            matrix = mmread(handle).T.tocsr()
    if matrix.shape != (51645, 36601) or matrix.nnz != 89102830:
        raise ValueError('External matrix dimension/nonzero mismatch')
    if (matrix.data < 0).any() or not np.equal(matrix.data, np.floor(matrix.data)).all():
        raise ValueError('Expected nonnegative integer external UMIs')
    # Confirm independent count and metadata distributions align cell by cell.
    curated = ad.read_h5ad(raw_root / 'wang_cellxgene.h5ad', backed='r')
    validation_genes = ['NPC2', 'PSMC3', 'DGKG', 'SLC2A12', 'COL8A2']
    checks = []
    for gene in validation_genes:
        row = candidates.loc[candidates['gene'].eq(gene)].iloc[0]
        gene_id = row['gene_id']
        if gene_id not in curated.raw.var_names or gene not in id_to_index:
            raise ValueError('Required count validation gene missing')
        curated_values = curated.raw.X[:, curated.raw.var_names.get_loc(gene_id)].toarray().ravel()
        curated_series = pd.Series(curated_values, index=curated.obs_names).reindex(barcodes)
        geo_values = matrix[:, id_to_index[gene]].toarray().ravel()
        exact = np.array_equal(curated_series.to_numpy(), geo_values)
        checks.append({'gene': gene, 'n_cells_checked': len(barcodes), 'cellwise_count_exact_match': exact})
        if not exact:
            raise ValueError('GEO/CellxGene count mismatch: ' + gene)
    curated.file.close()
    pd.DataFrame(checks).to_csv(directory / 'external_count_alignment_validation.csv', index=False)
    rows, group_sums, group_info = [], [], []
    for donor in donors:
        for cell_type in config['all_types']:
            selected = metadata['donor_id'].eq(donor).to_numpy() & metadata['harmonized_type'].eq(cell_type).to_numpy()
            sums = np.asarray(matrix[selected].sum(axis=0)).ravel().astype(np.int64)
            total = int(sums.sum()); n = int(selected.sum())
            group_sums.append(sums); group_info.append((donor, cell_type, n, total))
            for j, gene in enumerate(candidates['gene']):
                count = int(sums[positions[j]]) if mapped[j] else np.nan
                cpm = 1e6 * count / total if mapped[j] and total else np.nan
                rows.append({'gene': gene, 'donor': donor, 'cell_type': cell_type, 'n_cells': n,
                    'sum_gene_UMIs': count, 'sum_all_UMIs': total, 'pseudobulk_CPM': cpm,
                    'log1p_pseudobulk_CPM': np.log1p(cpm), 'external_gene_mapped': bool(mapped[j])})
    frame = pd.DataFrame(rows)
    frame.to_csv(directory / 'external_donor_candidate_expression.csv', index=False)
    info = pd.DataFrame(group_info, columns=['donor', 'cell_type', 'n_cells', 'total'])
    np.savez_compressed(directory / 'external_all_gene_pseudobulks.npz', sums=np.asarray(group_sums),
        genes=feature_symbols, donor=info['donor'].to_numpy(dtype=str),
        cell_type=info['cell_type'].to_numpy(dtype=str), n_cells=info['n_cells'].to_numpy(), total=info['total'].to_numpy())
    del matrix
    original = pd.read_csv(stages['expanded'] / 'expanded_robustness.csv').set_index('gene')
    original_scores = pd.read_csv(stages['expanded'] / 'balanced_candidate_expression.csv')
    results = []
    for j, gene in enumerate(candidates['gene']):
        record = {'gene': gene, 'gene_id': candidates.iloc[j]['gene_id'], 'external_exact_symbol_mapped': bool(mapped[j]),
            'V1_exact_symbol_mapped': bool(candidates.iloc[j]['exact_symbol_in_matrix'])}
        for key, label in [('primary_types', 'primary'), ('extended_types', 'extended')]:
            types = config[key]
            subset = frame.loc[frame['gene'].eq(gene) & frame['cell_type'].isin(types)]
            values = subset.groupby('cell_type')['log1p_pseudobulk_CPM'].mean().reindex(types)
            winner = int(winners(values.to_numpy()[:, None])[0]) if mapped[j] else -3
            top = types[winner] if winner >= 0 else ''
            prior = original.loc[gene, label + '_top_type']
            prior = prior if isinstance(prior, str) else ''
            evaluable = winner >= 0 and original.loc[gene, label + '_status'] == 'unique'
            same = evaluable and top == prior
            lodo_agree = lodo_valid = 0
            for omitted in donors:
                reduced = subset.loc[~subset['donor'].eq(omitted)].groupby('cell_type')['log1p_pseudobulk_CPM'].mean().reindex(types)
                reduced_top = int(winners(reduced.to_numpy()[:, None])[0]) if mapped[j] else -3
                lodo_valid += int(winner >= 0 and reduced_top >= 0)
                lodo_agree += int(winner >= 0 and reduced_top == winner)
            prior_scores = original_scores.loc[original_scores['gene'].eq(gene) & original_scores['comparison'].eq(label)].set_index('cell_type')['equal_donor_mean_log1p_CPM'].reindex(types)
            rank_corr = float(spearmanr(values, prior_scores).statistic) if evaluable and np.std(values) > 0 and np.std(prior_scores) > 0 else np.nan
            record.update({label + '_external_top': top, label + '_V1_top': prior,
                label + '_external_status': 'unique' if winner >= 0 else {-1: 'no_UMIs', -2: 'tie', -3: 'missing_symbol'}[winner],
                label + '_comparison_evaluable': evaluable, label + '_context_concordant': same,
                label + '_rank_spearman': rank_corr, label + '_external_LODO_agree': lodo_agree,
                label + '_external_LODO_evaluable': lodo_valid})
        results.append(record)
    summary = pd.DataFrame(results)
    summary.to_csv(directory / 'cross_study_context_concordance.csv', index=False)
    save_json(directory / 'summary.json', {'external_cells': 51645, 'external_donors': donors, 'external_n_donors': 4,
        'external_eyes': int(donor_details['n_eyes'].sum()), 'features': 36601,
        'external_exact_symbol_mapped_candidates': int(mapped.sum()),
        'primary_evaluable': int(summary['primary_comparison_evaluable'].sum()),
        'primary_context_concordant': int(summary['primary_context_concordant'].sum()),
        'extended_evaluable': int(summary['extended_comparison_evaluable'].sum()),
        'extended_context_concordant': int(summary['extended_context_concordant'].sum()),
        'barcode_alignment': 'All 51,645 GEO barcodes exactly match CellxGene labels; five genes agree cell by cell across both distributions',
        'source': 'Wang et al. 2022; DOI 10.1016/j.xgen.2022.100164; GSE196235; CellxGene collection 348da6dc-5bf6-435d-adc5-37747b9ae38a',
        'scope': 'Cross-study descriptive RNA context; V1 whole-cell RNA versus Wang nuclear multiome RNA; not disease replication, causal validation, or pooled-cell inference',
        'independence': 'Distinct studies and published donor metadata/ages; no shared global person identifiers are available for formal cross-study identity linkage'})

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from .cache import save_json, winners


def run(root, directory, config):
    data = ad.read_h5ad(root / config['source_v1_normalized'])
    candidates = pd.read_csv(root / config['source_candidates'])
    if len(candidates) != 228 or candidates['gene'].duplicated().any():
        raise ValueError('Candidate set must have 228 unique symbols')
    if not data.var_names.is_unique or not data.obs_names.is_unique:
        raise ValueError('Duplicate matrix identifiers')
    raw = data.layers['counts']
    if (raw.data < 0).any() or not np.equal(raw.data, np.floor(raw.data)).all():
        raise ValueError('Expected nonnegative integer UMIs')
    donors, types = config['donors'], config['all_types']
    obs = data.obs
    gene_indices = data.var_names.get_indexer(candidates['gene'])
    mapped = gene_indices >= 0
    mapping = candidates.copy()
    mapping['exact_symbol_in_matrix'] = mapped
    mapping['mapping_status'] = np.where(mapped, 'exact_symbol', 'missing_symbol')
    mapping.to_csv(directory / 'candidate_mapping.csv', index=False)
    groups, all_sums, rows = [], [], []
    for qc in ['primary', 'strict']:
        qc_mask = np.ones(data.n_obs, bool) if qc == 'primary' else obs['qc_strict'].to_numpy(bool)
        for donor in donors:
            for cell_type in types:
                selected = qc_mask & obs['donor'].eq(donor).to_numpy() & obs['published_cell_type'].eq(cell_type).to_numpy()
                n = int(selected.sum())
                sums = np.asarray(raw[selected].sum(axis=0)).ravel().astype(np.int64)
                total = int(sums.sum())
                groups.append((qc, donor, cell_type, n, total))
                all_sums.append(sums)
                for j, gene in enumerate(candidates['gene']):
                    count = int(sums[gene_indices[j]]) if mapped[j] else np.nan
                    cpm = 1e6 * count / total if total and mapped[j] else np.nan
                    rows.append({'qc': qc, 'donor': donor, 'cell_type': cell_type, 'gene': gene,
                        'n_cells': n, 'sum_all_UMIs': total, 'sum_gene_UMIs': count,
                        'log1p_pseudobulk_CPM': np.log1p(cpm), 'pseudobulk_CPM': cpm,
                        'gene_in_matrix': bool(mapped[j])})
    group_frame = pd.DataFrame(groups, columns=['qc', 'donor', 'cell_type', 'n_cells', 'sum_all_UMIs'])
    group_frame.to_csv(directory / 'donor_type_counts.csv', index=False)
    np.savez_compressed(directory / 'all_gene_pseudobulks.npz', sums=np.asarray(all_sums),
        genes=data.var_names.to_numpy(dtype=str), qc=group_frame['qc'].to_numpy(dtype=str),
        donor=group_frame['donor'].to_numpy(dtype=str), cell_type=group_frame['cell_type'].to_numpy(dtype=str),
        n_cells=group_frame['n_cells'].to_numpy(), total=group_frame['sum_all_UMIs'].to_numpy())
    frame = pd.DataFrame(rows)
    frame.to_csv(directory / 'donor_candidate_expression.csv', index=False)
    for key, minimum in [('primary_types', config['primary_min_cells']), ('extended_types', config['extended_min_cells'])]:
        subset = group_frame.loc[group_frame['qc'].eq('primary') & group_frame['cell_type'].isin(config[key])]
        if len(subset) != len(donors) * len(config[key]) or (subset['n_cells'] < minimum).any():
            raise ValueError('Frozen comparison-set eligibility failed: ' + key)
    strict_groups = group_frame.loc[group_frame['qc'].eq('strict') & group_frame['cell_type'].isin(config['primary_types'])]
    if (strict_groups['n_cells'] < config['primary_min_cells']).any():
        raise ValueError('Strict QC lacks frozen primary groups')
    summaries, balanced, folds = [], [], []
    for j, gene in enumerate(candidates['gene']):
        record = {'gene': gene, 'gene_in_matrix': bool(mapped[j])}
        for qc, key, label in [('primary', 'primary_types', 'primary'), ('strict', 'primary_types', 'strict'), ('primary', 'extended_types', 'extended')]:
            group = frame.loc[frame['gene'].eq(gene) & frame['qc'].eq(qc) & frame['cell_type'].isin(config[key])]
            scores = group.groupby('cell_type')['log1p_pseudobulk_CPM'].mean().reindex(config[key])
            top = int(winners(scores.to_numpy()[:, None])[0]) if mapped[j] else -3
            record[label + '_top_type'] = config[key][top] if top >= 0 else ''
            record[label + '_status'] = 'unique' if top >= 0 else {-1: 'no_UMIs', -2: 'tie', -3: 'missing_symbol'}[top]
            for cell_type, score in scores.items():
                balanced.append({'gene': gene, 'comparison': label, 'cell_type': cell_type,
                    'equal_donor_mean_log1p_CPM': score})
            if label == 'primary' and mapped[j]:
                agree = valid = 0
                for omitted in donors:
                    reduced = group.loc[~group['donor'].eq(omitted)].groupby('cell_type')['log1p_pseudobulk_CPM'].mean().reindex(config[key])
                    fold_top = int(winners(reduced.to_numpy()[:, None])[0])
                    evaluable = fold_top >= 0 and top >= 0
                    same = evaluable and fold_top == top
                    valid += int(evaluable); agree += int(same)
                    folds.append({'gene': gene, 'omitted_donor': omitted,
                        'top_type': config[key][fold_top] if fold_top >= 0 else '',
                        'evaluable': evaluable, 'agrees_with_primary': same,
                        'comparison_types': ';'.join(config[key])})
                record['LODO_agree'] = agree
                record['LODO_evaluable'] = valid
                record['LODO_all3'] = valid == 3 and agree == 3
        record['QC_same_top'] = record['primary_status'] == record['strict_status'] == 'unique' and record['primary_top_type'] == record['strict_top_type']
        assigned = frame.loc[frame['gene'].eq(gene) & frame['qc'].eq('primary')]
        record['assigned_cell_UMIs'] = assigned['sum_gene_UMIs'].sum() if mapped[j] else np.nan
        record['assigned_donors_detected'] = int((assigned.groupby('donor')['sum_gene_UMIs'].sum() > 0).sum()) if mapped[j] else 0
        record['sparse_UMI_flag_lt50'] = bool(mapped[j] and record['assigned_cell_UMIs'] < 50)
        summaries.append(record)
    summary = pd.DataFrame(summaries)
    summary.to_csv(directory / 'expanded_robustness.csv', index=False)
    pd.DataFrame(balanced).to_csv(directory / 'balanced_candidate_expression.csv', index=False)
    pd.DataFrame(folds).to_csv(directory / 'leave_one_donor_out.csv', index=False)
    save_json(directory / 'summary.json', {'cells': data.n_obs, 'genes': data.n_vars,
        'donors': donors, 'candidates': len(candidates), 'exact_mapped': int(mapped.sum()),
        'missing_symbols': int((~mapped).sum()), 'LODO_stable_all3': int(summary['LODO_all3'].fillna(False).sum()),
        'QC_same_top': int(summary['QC_same_top'].sum()), 'primary_types': config['primary_types'],
        'extended_types': config['extended_types'], 'scope': 'Descriptive healthy-donor cell context, not disease or causal effects'})

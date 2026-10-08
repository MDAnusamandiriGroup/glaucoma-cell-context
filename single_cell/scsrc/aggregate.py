"""Descriptive donor-level expression; no cell-level significance tests."""
from pathlib import Path
import json
import anndata as ad
import numpy as np
import pandas as pd

from .common import save_csv, save_json


def group_expression(counts, normalized, obs, genes, targets, types, minimum):
    """Sum raw UMIs within each donor/type; retain small groups as descriptive."""
    index = {gene: i for i, gene in enumerate(genes)}
    available = [gene for gene in targets if gene in index]
    raw = counts[:, [index[g] for g in available]].toarray()
    lognorm = normalized[:, [index[g] for g in available]].toarray()
    column = {gene: i for i, gene in enumerate(available)}
    totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    donors = sorted(obs['donor'].astype(str).unique())
    rows = []
    donor_values = obs['donor'].astype(str).to_numpy()
    type_values = obs['published_cell_type'].astype(str).to_numpy()
    for qc in ['primary', 'strict']:
        qc_mask = np.ones(len(obs), dtype=bool) if qc == 'primary' else obs['qc_strict'].to_numpy(dtype=bool)
        for donor in donors:
            for cell_type in types:
                selected = qc_mask & (donor_values == donor) & (type_values == cell_type)
                n = int(selected.sum())
                denominator = float(totals[selected].sum())
                for gene in targets:
                    present = gene in column
                    gene_sum = float(raw[selected, column[gene]].sum()) if n and present else np.nan
                    cpm = 1e6 * gene_sum / denominator if denominator > 0 and present else np.nan
                    rows.append({
                        'qc': qc, 'donor': donor, 'cell_type': cell_type, 'gene': gene,
                        'n_cells': n, 'sum_all_UMIs': denominator,
                        'sum_gene_UMIs': gene_sum, 'pseudobulk_CPM': cpm,
                        'log1p_pseudobulk_CPM': float(np.log1p(cpm)) if np.isfinite(cpm) else np.nan,
                        'mean_log1p_CP10k': float(lognorm[selected, column[gene]].mean()) if n and present else np.nan,
                        'pct_cells_detected': float(100 * np.mean(raw[selected, column[gene]] > 0)) if n and present else np.nan,
                        'gene_in_matrix': present, 'eligible': n >= minimum and present,
                        'status': 'gene_not_in_matrix' if not present else ('no_cells' if n == 0 else ('below_minimum' if n < minimum else 'eligible'))
                    })
    return pd.DataFrame(rows)


def balanced_expression(frame, minimum):
    rows = []
    for (qc, cell_type, gene), group in frame.groupby(['qc', 'cell_type', 'gene'], sort=False):
        eligible = group.loc[(group['n_cells'] >= minimum) & group['gene_in_matrix']]
        record = {'qc': qc, 'cell_type': cell_type, 'gene': gene, 'minimum_cells': minimum,
                  'eligible_donors': int(len(eligible)), 'total_cells': int(group['n_cells'].sum())}
        for source, destination in [
            ('log1p_pseudobulk_CPM', 'donor_mean_log1p_pseudobulk_CPM'),
            ('mean_log1p_CP10k', 'donor_mean_log1p_CP10k'),
            ('pct_cells_detected', 'donor_mean_pct_detected')
        ]:
            record[destination] = float(eligible[source].mean()) if len(eligible) else np.nan
        rows.append(record)
    return pd.DataFrame(rows)


def shared_types(frame, qc, types, minimum, donors):
    counts = frame.loc[frame['qc'].eq(qc)].drop_duplicates(['donor', 'cell_type'])
    counts = counts.pivot(index='cell_type', columns='donor', values='n_cells').reindex(index=types, columns=donors).fillna(0)
    return [cell_type for cell_type in types if bool((counts.loc[cell_type] >= minimum).all())]


def top_type(frame, gene, qc, types, donors):
    group = frame.loc[frame['gene'].eq(gene) & frame['qc'].eq(qc) & frame['cell_type'].isin(types) & frame['donor'].isin(donors)]
    if group.empty or not bool(group['gene_in_matrix'].all()):
        return '', 'gene_not_in_matrix', np.nan
    scores = group.groupby('cell_type')['log1p_pseudobulk_CPM'].mean().reindex(types)
    if scores.isna().any() or len(scores) == 0:
        return '', 'insufficient_groups', np.nan
    maximum = float(scores.max())
    if maximum == 0:
        return '', 'no_detected_UMIs', 0.0
    winners = scores.index[np.isclose(scores.to_numpy(), maximum, rtol=0, atol=1e-10)].tolist()
    return ';'.join(winners), 'unique' if len(winners) == 1 else 'tie', maximum


def robustness(frame, targets, types, minimum, sensitivity_minimum):
    donors = sorted(frame['donor'].unique())
    shared = {qc: shared_types(frame, qc, types, minimum, donors) for qc in ['primary', 'strict']}
    common_qc = [t for t in types if t in shared['primary'] and t in shared['strict']]
    sensitivity = shared_types(frame, 'primary', types, sensitivity_minimum, donors)
    summaries, folds, donor_tops = [], [], []
    for gene in targets:
        baseline, status, score = top_type(frame, gene, 'primary', shared['primary'], donors)
        agreements = 0
        valid_folds = 0
        for omitted in donors:
            included = [donor for donor in donors if donor != omitted]
            winner, fold_status, fold_score = top_type(frame, gene, 'primary', shared['primary'], included)
            valid = status == 'unique' and fold_status == 'unique'
            same = bool(valid and winner == baseline)
            agreements += int(same)
            valid_folds += int(valid)
            folds.append({'gene': gene, 'omitted_donor': omitted, 'included_donors': ';'.join(included),
                          'baseline_top_type': baseline, 'fold_top_type': winner, 'fold_status': fold_status,
                          'same_as_baseline': same, 'top_score': fold_score,
                          'comparison_types': ';'.join(shared['primary']),
                          'operation': 'Reaggregation of saved donor pseudobulks; no reclustering or model refit'})
            winner_d, status_d, score_d = top_type(frame, gene, 'primary', shared['primary'], [omitted])
            donor_tops.append({'gene': gene, 'donor': omitted, 'top_type': winner_d, 'status': status_d, 'top_score': score_d})
        primary_qc, primary_status, _ = top_type(frame, gene, 'primary', common_qc, donors)
        strict_qc, strict_status, _ = top_type(frame, gene, 'strict', common_qc, donors)
        sens_top, sens_status, _ = top_type(frame, gene, 'primary', sensitivity, donors)
        qc_evaluable = primary_status == 'unique' and strict_status == 'unique'
        summaries.append({'gene': gene, 'primary_top_type': baseline, 'primary_status': status,
                          'primary_top_score': score, 'common_types_n': len(shared['primary']),
                          'LODO_agree_folds': agreements, 'LODO_evaluable_folds': valid_folds,
                          'LODO_all_three_agree': agreements == len(donors) and valid_folds == len(donors),
                          'primary_QC_comparable_top': primary_qc, 'strict_QC_comparable_top': strict_qc,
                          'QC_evaluable': qc_evaluable,
                          'QC_top_agrees': bool(qc_evaluable and primary_qc == strict_qc),
                          'sensitivity_min10_top': sens_top, 'sensitivity_min10_status': sens_status,
                          'common_QC_types': ';'.join(common_qc), 'min10_types': ';'.join(sensitivity)})
    return pd.DataFrame(summaries), pd.DataFrame(folds), pd.DataFrame(donor_tops), shared, common_qc, sensitivity


def expression_panels(data, labels, panels):
    """Pooled descriptive marker means; never used for donor-level inference."""
    genes = list(dict.fromkeys(gene for values in panels.values() for gene in values))
    available = [g for g in genes if g in data.var_names]
    norm = data.X[:, data.var_names.get_indexer(available)].toarray()
    raw = data.layers['counts'][:, data.var_names.get_indexer(available)].toarray()
    column = {g: i for i, g in enumerate(available)}
    rows = []
    for label in sorted(set(labels), key=str):
        mask = np.asarray(labels) == label
        for expected, panel in panels.items():
            for gene in panel:
                rows.append({'group': str(label), 'expected_cell_type': expected, 'gene': gene,
                             'n_cells': int(mask.sum()), 'gene_in_matrix': gene in column,
                             'mean_log1p_CP10k': float(norm[mask, column[gene]].mean()) if gene in column else np.nan,
                             'pct_detected': float(100 * np.mean(raw[mask, column[gene]] > 0)) if gene in column else np.nan})
    return pd.DataFrame(rows)


def run(root, preprocess_dir, graph_dir, directory, config):
    data = ad.read_h5ad(Path(preprocess_dir) / 'normalized.h5ad')
    graph = ad.read_h5ad(Path(graph_dir) / 'graph.h5ad')
    if not data.obs_names.equals(graph.obs_names):
        raise ValueError('Cluster and expression cell order does not match')
    candidates = pd.read_csv(Path(root) / 'inputs/published_candidates.csv')
    if len(candidates) != 228 or candidates['gene'].duplicated().any():
        raise ValueError('Published candidate set changed; review before analysing')
    mapping = candidates[['gene_id', 'gene', 'retina_eQTL', 'strict_retina_eQTL', 'highest_mean_expression_cell']].copy()
    mapping['exact_symbol_in_count_matrix'] = mapping['gene'].isin(data.var_names)
    mapping['targeted_in_this_extension'] = mapping['gene'].isin(config['target_genes'])
    mapping['mapping_policy'] = 'Exact source symbol; absent symbols retained as missing; no inferred aliases'
    save_csv(directory / 'candidate_mapping_audit.csv', mapping)
    frame = group_expression(data.layers['counts'], data.X, data.obs, data.var_names.tolist(), config['target_genes'], config['cell_types'], config['min_cells_per_donor_type'])
    save_csv(directory / 'donor_gene_expression.csv', frame)
    balanced = balanced_expression(frame, config['min_cells_per_donor_type'])
    save_csv(directory / 'balanced_gene_expression.csv', balanced)
    save_csv(directory / 'balanced_gene_expression_min10.csv', balanced_expression(frame, config['min_cells_sensitivity']))
    counts = frame.drop_duplicates(['qc', 'donor', 'cell_type'])[['qc', 'donor', 'cell_type', 'n_cells']].copy()
    counts['minimum_cells'] = config['min_cells_per_donor_type']
    counts['eligible_group'] = counts['n_cells'] >= config['min_cells_per_donor_type']
    save_csv(directory / 'donor_cell_type_counts.csv', counts)
    summary, folds, donor_tops, shared, common_qc, sensitivity = robustness(frame, config['target_genes'], config['cell_types'], config['min_cells_per_donor_type'], config['min_cells_sensitivity'])
    summary['prior_table_top_type'] = summary['gene'].map(mapping.set_index('gene')['highest_mean_expression_cell'])
    summary['prior_table_type_in_primary_comparison'] = summary['prior_table_top_type'].isin(shared['primary'])
    # The old pooled class means and new equal-donor log-CPM have different
    # estimands; this comparison is descriptive, not independent replication.
    save_csv(directory / 'robustness_summary.csv', summary)
    save_csv(directory / 'leave_one_donor_out.csv', folds)
    save_csv(directory / 'individual_donor_top_type.csv', donor_tops)
    markers = expression_panels(data, data.obs['published_cell_type'].astype(str).to_numpy(), config['markers'])
    save_csv(directory / 'marker_expression_published_types.csv', markers)
    cluster_audits = []
    for resolution in config['leiden_resolutions']:
        key = 'leiden_' + str(resolution)
        clusters = graph.obs[key].astype(str)
        contingency = pd.crosstab(clusters, data.obs['published_cell_type'], dropna=False)
        donor_counts = pd.crosstab(clusters, data.obs['donor'], dropna=False)
        save_csv(directory / ('cluster_cell_type_contingency_' + str(resolution) + '.csv'), contingency.reset_index())
        save_csv(directory / ('cluster_donor_contingency_' + str(resolution) + '.csv'), donor_counts.reset_index())
        marker_frame = expression_panels(data, clusters.to_numpy(), config['markers'])
        save_csv(directory / ('marker_expression_clusters_' + str(resolution) + '.csv'), marker_frame)
        wide = marker_frame.pivot(index='group', columns='gene', values='mean_log1p_CP10k')
        # Per-gene z scores describe a panel's relative enrichment across
        # clusters. Provisional labels are an audit aid, not expert annotation.
        standard = (wide - wide.mean(axis=0)) / wide.std(axis=0, ddof=0).replace(0, np.nan)
        for cluster, row in contingency.iterrows():
            scores = {cell_type: float(standard.loc[str(cluster), [g for g in panel if g in standard]].mean()) for cell_type, panel in config['markers'].items()}
            finite = {label: score for label, score in scores.items() if np.isfinite(score)}
            provisional = max(finite, key=finite.get) if finite else 'Unassigned'
            n = int(row.sum())
            dominant = str(row.idxmax())
            cluster_audits.append({'resolution': resolution, 'cluster': str(cluster), 'n_cells': n,
                                   'published_dominant_type': dominant, 'published_dominant_fraction': float(row.max() / n),
                                   'provisional_marker_panel': provisional,
                                   'provisional_marker_score': finite.get(provisional, np.nan),
                                   'marker_panel_matches_published_dominant': provisional == dominant,
                                   'largest_donor_fraction': float(donor_counts.loc[cluster].max() / n),
                                   'interpretation': 'Descriptive cluster/marker audit; published labels remain primary'})
    save_csv(directory / 'cluster_annotation_audit.csv', pd.DataFrame(cluster_audits))
    matched = int(mapping['exact_symbol_in_count_matrix'].sum())
    save_json(directory / 'summary.json', {
        'target_genes': config['target_genes'],
        'target_genes_in_count_matrix': [g for g in config['target_genes'] if g in data.var_names],
        'published_candidate_genes_exactly_mapped': matched, 'published_candidate_genes_total': len(mapping),
        'primary_shared_cell_types': shared['primary'], 'strict_shared_cell_types': shared['strict'],
        'QC_comparable_cell_types': common_qc, 'min10_sensitivity_shared_cell_types': sensitivity,
        'minimum_cells': config['min_cells_per_donor_type'],
        'primary_estimand': 'Equal-donor mean of log1p pseudobulk CPM across cell types with >=20 cells in each of all three donors',
        'LODO_stable_genes': summary.loc[summary['LODO_all_three_agree'], 'gene'].tolist(),
        'QC_stable_genes': summary.loc[summary['QC_top_agrees'], 'gene'].tolist(),
        'LODO_operation': 'Reaggregate existing donor-by-cell-type pseudobulks; no reclustering or model refit',
        'statistical_scope': 'Descriptive; no significance tests, no cell-level pseudo-replication and no disease/age effect estimates',
        'annotation_scope': 'Published original labels for donor summaries; independent Leiden clusters with provisional marker-panel audit'
    })

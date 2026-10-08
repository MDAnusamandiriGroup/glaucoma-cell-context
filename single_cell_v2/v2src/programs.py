import numpy as np
import pandas as pd
from .cache import save_json, unit_rng, winners


def run(root, expanded_dir, directory, config):
    bulk = np.load(expanded_dir / 'all_gene_pseudobulks.npz')
    genes = bulk['genes'].astype(str)
    matrix = np.log1p(1e6 * bulk['sums'] / bulk['total'][:, None])
    primary = (bulk['qc'] == 'primary') & np.isin(bulk['cell_type'], config['primary_types'])
    reference_mean = matrix[primary].mean(axis=0)
    gene_to_index = {g: i for i, g in enumerate(genes)}
    sets = {}
    for line in (root / 'data/raw/hallmark2025.gmt').read_text().splitlines():
        fields = line.split('\t')
        sets[fields[0]] = list(dict.fromkeys(fields[2:]))
    panels = {name: [g for g in sets[name] if not g.startswith('MT-')] for name in config['programs']}
    excluded = set(g for panel in panels.values() for g in panel)
    expressed_nuclear = np.flatnonzero((reference_mean > 0) & ~np.char.startswith(genes, 'MT-'))
    order = expressed_nuclear[np.argsort(reference_mean[expressed_nuclear], kind='stable')]
    bins = np.full(len(genes), -1, dtype=int)
    bins[order] = np.minimum(config['reference_bins'] - 1,
        np.arange(len(order)) * config['reference_bins'] // len(order))
    control_pool = {b: [i for i in order if bins[i] == b and genes[i] not in excluded]
                    for b in range(config['reference_bins'])}
    score_rows, coverage, references, membership, summaries = [], [], [], [], []
    candidates = pd.read_csv(expanded_dir / 'candidate_mapping.csv')
    for name, panel in panels.items():
        available = [g for g in panel if g in gene_to_index and reference_mean[gene_to_index[g]] > 0]
        fraction = len(available) / len(panel)
        eligible = len(available) >= config['program_min_genes'] and fraction >= config['program_min_fraction']
        coverage.append({'program': name, 'source_release': config['hallmark_release'],
            'nuclear_panel_size': len(panel), 'expressed_exact_mapped': len(available),
            'coverage_fraction': fraction, 'eligible': eligible,
            'excluded_mt_genes': ';'.join(g for g in sets[name] if g.startswith('MT-'))})
        for gene in candidates['gene']:
            membership.append({'gene': gene, 'program': name, 'member_exact_symbol': gene in panel,
                'in_matrix_and_expressed': gene in available})
        if not eligible:
            continue
        reference_means = []
        for gene in available:
            pool = control_pool[bins[gene_to_index[gene]]]
            n = min(config['references_per_gene'], len(pool))
            if n < 5:
                raise ValueError('Insufficient matched-reference genes')
            rng = unit_rng(config['seed'], 'program_controls', name, gene)
            selected = rng.choice(pool, n, replace=False)
            reference_means.append(matrix[:, selected].mean(axis=1))
            for idx in selected:
                references.append({'program': name, 'panel_gene': gene, 'reference_gene': genes[idx],
                    'expression_bin': int(bins[idx]), 'weight_within_panel': 1 / (len(available) * n)})
        score = matrix[:, [gene_to_index[g] for g in available]].mean(axis=1) - np.stack(reference_means).mean(axis=0)
        group_frame = pd.DataFrame({'qc': bulk['qc'], 'donor': bulk['donor'], 'cell_type': bulk['cell_type'],
            'n_cells': bulk['n_cells'], 'score': score})
        for record in group_frame.to_dict('records'):
            score_rows.append({'program': name, **record})
        for qc, key, label in [('primary', 'primary_types', 'primary'), ('strict', 'primary_types', 'strict'), ('primary', 'extended_types', 'extended')]:
            frame = group_frame.loc[group_frame['qc'].eq(qc) & group_frame['cell_type'].isin(config[key])]
            balanced = frame.groupby('cell_type')['score'].mean().reindex(config[key])
            # Program scores can be negative; top is a relative maximum, not activation.
            top = str(balanced.idxmax())
            lodo = []
            for donor in config['donors']:
                reduced = frame.loc[~frame['donor'].eq(donor)].groupby('cell_type')['score'].mean().reindex(config[key])
                lodo.append(str(reduced.idxmax()))
            summaries.append({'program': name, 'comparison': label, 'relative_top_type': top,
                'LODO_same_top_folds': sum(t == top for t in lodo), 'LODO_top_types': ';'.join(lodo),
                'comparison_types': ';'.join(config[key]),
                **{'score_' + ct: float(balanced[ct]) for ct in config[key]}})
    pd.DataFrame(coverage).to_csv(directory / 'program_coverage.csv', index=False)
    pd.DataFrame(references).to_csv(directory / 'matched_reference_genes.csv', index=False)
    pd.DataFrame(score_rows).to_csv(directory / 'donor_program_scores.csv', index=False)
    pd.DataFrame(summaries).to_csv(directory / 'program_context_summary.csv', index=False)
    pd.DataFrame(membership).to_csv(directory / 'candidate_program_membership.csv', index=False)
    mito = pd.read_excel(root / 'data/raw/Human.MitoCarta3.0.xls', sheet_name='A Human MitoCarta3.0')
    if len(mito) != 1136 or mito['Symbol'].duplicated().any():
        raise ValueError('Unexpected MitoCarta3.0 inventory')
    mito_map = mito.set_index('Symbol')
    annotation = candidates[['gene', 'exact_symbol_in_matrix']].copy()
    annotation['MitoCarta3_member'] = annotation['gene'].isin(mito_map.index)
    annotation['MitoCarta3_subcompartment'] = annotation['gene'].map(mito_map['MitoCarta3.0_SubMitoLocalization'])
    annotation['MitoCarta3_pathways'] = annotation['gene'].map(mito_map['MitoCarta3.0_MitoPathways'])
    annotation.to_csv(directory / 'candidate_mitocarta_annotation.csv', index=False)
    save_json(directory / 'summary.json', {'programs_tested': len(panels), 'eligible_programs': sum(c['eligible'] for c in coverage),
        'gene_sets': 'MSigDB Hallmark 2025.1.Hs; exact symbols; mitochondrial-genome genes excluded',
        'reference_strategy': '20 mean-expression rank bins, up to 20 controls per gene, all six panels excluded from controls',
        'candidate_MitoCarta3_members': int(annotation['MitoCarta3_member'].sum()),
        'score_scope': 'Descriptive transcriptional program contrasts; not metabolic flux, mitochondrial function, disease effects or pathway enrichment p-values',
        'inference': 'No population significance tests from three healthy donors'})

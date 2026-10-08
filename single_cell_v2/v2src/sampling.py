import anndata as ad
import numpy as np
import pandas as pd
from .cache import save_json, unit_rng, winners


def thin_group(counts, total, fraction, rng):
    counts = np.asarray(counts, dtype=np.int64)
    other = int(total) - int(counts.sum())
    if other < 0:
        raise ValueError('Selected gene UMIs exceed library UMIs')
    selected = rng.binomial(counts, fraction)
    denominator = int(selected.sum()) + int(rng.binomial(other, fraction))
    return selected, denominator


def run(root, expanded_dir, directory, config):
    data = ad.read_h5ad(root / config['source_v1_normalized'])
    mapping = pd.read_csv(expanded_dir / 'candidate_mapping.csv')
    genes = mapping.loc[mapping['exact_symbol_in_matrix'], 'gene'].tolist()
    if len(genes) != len(set(genes)):
        raise ValueError('Duplicate selected genes')
    raw = data.layers['counts'][:, data.var_names.get_indexer(genes)].toarray().astype(np.int64)
    totals = np.asarray(data.layers['counts'].sum(axis=1)).ravel().astype(np.int64)
    donors = config['donors']; obs = data.obs
    groups = {}
    for donor in donors:
        for cell_type in config['extended_types']:
            idx = np.flatnonzero(obs['donor'].eq(donor).to_numpy() & obs['published_cell_type'].eq(cell_type).to_numpy())
            groups[(donor, cell_type)] = idx
    scenarios = [('equal_cells_20', 'cell_n', 20, 'primary_types'), ('equal_cells_10_extended', 'cell_n', 10, 'extended_types')]
    scenarios += [('cell_fraction_' + str(p), 'cell_fraction', p, 'primary_types') for p in config['cell_fraction_scenarios']]
    scenarios += [('UMI_fraction_' + str(p), 'umi_fraction', p, 'primary_types') for p in config['umi_fraction_scenarios']]
    base = pd.read_csv(expanded_dir / 'expanded_robustness.csv').set_index('gene')
    summaries, group_design, arrays = [], [], {}
    repetitions = config['replicates']
    for name, operation, parameter, key in scenarios:
        types = config[key]
        baseline_column = 'extended_top_type' if key == 'extended_types' else 'primary_top_type'
        baseline = base.loc[genes, baseline_column].fillna('').to_numpy()
        scores = np.zeros((repetitions, len(types), len(genes)), dtype=np.float32)
        for d, donor in enumerate(donors):
            for t, cell_type in enumerate(types):
                idx = groups[(donor, cell_type)]
                n = int(parameter) if operation == 'cell_n' else max(1, int(np.floor(len(idx) * parameter))) if operation == 'cell_fraction' else len(idx)
                if n > len(idx):
                    raise ValueError('Cannot sample more cells than available')
                group_design.append({'scenario': name, 'donor': donor, 'cell_type': cell_type,
                    'available_cells': len(idx), 'selected_cells': n, 'operation': operation,
                    'parameter': parameter, 'comparison_types': ';'.join(types)})
                all_sums = raw[idx].sum(axis=0)
                all_total = int(totals[idx].sum())
                for rep in range(repetitions):
                    rng = unit_rng(config['seed'], name, rep, donor, cell_type)
                    if operation == 'umi_fraction':
                        counts, denominator = thin_group(all_sums, all_total, parameter, rng)
                    else:
                        chosen = rng.choice(idx, size=n, replace=False)
                        counts = raw[chosen].sum(axis=0); denominator = int(totals[chosen].sum())
                    if denominator <= 0:
                        raise ValueError('Empty sampled UMI library')
                    scores[rep, t] += (np.log1p(1e6 * counts / denominator) / len(donors)).astype(np.float32)
        top_draws = np.stack([winners(s) for s in scores])
        arrays[name + '_top_indices'] = top_draws
        arrays[name + '_types'] = np.asarray(types)
        arrays[name + '_scores'] = scores
        for g, gene in enumerate(genes):
            counts = np.bincount(top_draws[:, g][top_draws[:, g] >= 0], minlength=len(types))
            valid_draws = int(counts.sum())
            mode = int(np.argmax(counts)) if valid_draws else -1
            baseline_id = types.index(baseline[g]) if baseline[g] in types else -1
            summaries.append({'gene': gene, 'scenario': name, 'n_draws': repetitions,
                'comparison_types': ';'.join(types), 'baseline_top_type': baseline[g],
                'modal_top_type': types[mode] if mode >= 0 else '', 'valid_draws': valid_draws,
                'baseline_agreement_fraction': float(counts[baseline_id] / repetitions) if baseline_id >= 0 else np.nan,
                'modal_top_fraction': float(counts[mode] / repetitions) if mode >= 0 else 0.0,
                'zero_UMI_fraction': float(np.mean(top_draws[:, g] == -1)),
                'tie_fraction': float(np.mean(top_draws[:, g] == -2)),
                **{'top_fraction_' + ct: float(counts[t] / repetitions) for t, ct in enumerate(types)}})
        print('SAMPLING ' + name + ': ' + str(repetitions) + ' draws complete', flush=True)
    pd.DataFrame(summaries).to_csv(directory / 'sampling_stability.csv', index=False)
    pd.DataFrame(group_design).to_csv(directory / 'sampling_design.csv', index=False)
    arrays['genes'] = np.asarray(genes)
    np.savez_compressed(directory / 'sampling_draws.npz', **arrays)
    save_json(directory / 'summary.json', {'replicates_per_scenario': repetitions, 'scenarios': [s[0] for s in scenarios],
        'mapped_genes': len(genes), 'seed': config['seed'], 'number_scenario_draws': repetitions * len(scenarios),
        'scope': config['sampling_scope'], 'unit': 'Cell draws without replacement within donor and published type; depth draws thin selected genes and other UMIs jointly'})

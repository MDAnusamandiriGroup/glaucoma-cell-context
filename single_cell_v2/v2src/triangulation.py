"""Use cached aggregates to combine evidence without refitting upstream models."""
import json
import os
import numpy as np
import pandas as pd
from .cache import save_json
from .regulatory import tokens


def run(root, stages, directory, config):
    original = np.load(stages['expanded'] / 'all_gene_pseudobulks.npz')
    # This is our own manifest-verified stage, never an untrusted uploaded NPZ.
    # Pandas 3 produced an object string array for symbols. Convert that saved
    # representation to portable Unicode without reparsing/refitting the atlas.
    saved_external = np.load(stages['replication'] / 'external_all_gene_pseudobulks.npz', allow_pickle=True)
    external = {}
    for key in saved_external.files:
        value = saved_external[key]
        if value.dtype.kind == 'O':
            if not all(isinstance(x, (str, np.str_)) for x in value.ravel()):
                raise ValueError('Expected strings in the local cached object array')
            value = value.astype('U')
        external[key] = value
    saved_external.close()
    portable = directory / 'external_pseudobulks_portable.npz'
    temporary = directory / 'external_pseudobulks_portable.npz.partial'
    with temporary.open('xb') as handle:
        np.savez_compressed(handle, **external)
        handle.flush(); os.fsync(handle.fileno())
    if temporary.stat().st_size == 0:
        raise ValueError('Empty portable pseudobulk archive')
    temporary.rename(portable)
    refs = pd.read_csv(stages['programs'] / 'matched_reference_genes.csv')
    source_genes = original['genes'].astype(str); target_genes = external['genes'].astype(str)
    source_index = {g: i for i, g in enumerate(source_genes)}
    target_index = {g: i for i, g in enumerate(target_genes)}
    if len(target_index) != len(target_genes):
        raise ValueError('External symbols are not unique')
    source_matrix = np.log1p(1e6 * original['sums'] / original['total'][:, None])
    target_matrix = np.log1p(1e6 * external['sums'] / external['total'][:, None])
    coverage_rows, score_rows, program_summaries = [], [], []
    for program in config['programs']:
        program_refs = refs.loc[refs['program'].eq(program)]
        shared_pairs = []
        for gene, group in program_refs.groupby('panel_gene', sort=True):
            controls = [g for g in group['reference_gene'] if g in target_index and g in source_index]
            if gene in source_index and gene in target_index and len(controls) >= 5:
                shared_pairs.append((gene, controls))
        fraction = len(shared_pairs) / program_refs['panel_gene'].nunique()
        coverage_rows.append({'program': program, 'V1_panel_genes': program_refs['panel_gene'].nunique(),
            'common_panel_genes': len(shared_pairs), 'common_gene_fraction': fraction,
            'method': 'Identical intersected panel/control pairs and weights in both studies; no new control selection'})
        if len(shared_pairs) < config['program_min_genes'] or fraction < config['program_min_fraction']:
            continue
        for study, bulk, matrix, index in [('Lukowski2019', original, source_matrix, source_index), ('Wang2022', external, target_matrix, target_index)]:
            gene_contrasts = [matrix[:, index[g]] - matrix[:, [index[c] for c in controls]].mean(axis=1)
                              for g, controls in shared_pairs]
            score = np.stack(gene_contrasts).mean(axis=0)
            qc = bulk['qc'].astype(str) if study == 'Lukowski2019' else np.repeat('primary', len(score))
            frame = pd.DataFrame({'study': study, 'program': program, 'qc': qc, 'donor': bulk['donor'],
                'cell_type': bulk['cell_type'], 'n_cells': bulk['n_cells'], 'score': score})
            score_rows.extend(frame.to_dict('records'))
            donors = sorted(frame['donor'].unique())
            for key, label in [('primary_types', 'primary'), ('extended_types', 'extended')]:
                selected = frame.loc[frame['qc'].eq('primary') & frame['cell_type'].isin(config[key])]
                balanced = selected.groupby('cell_type')['score'].mean().reindex(config[key])
                top = str(balanced.idxmax())
                lodo_tops = [str(selected.loc[~selected['donor'].eq(donor)].groupby('cell_type')['score'].mean().reindex(config[key]).idxmax()) for donor in donors]
                program_summaries.append({'study': study, 'program': program, 'comparison': label,
                    'relative_top_type': top, 'LODO_agree': sum(t == top for t in lodo_tops),
                    'LODO_folds': len(donors), 'comparison_types': ';'.join(config[key]),
                    **{'score_' + ct: float(balanced[ct]) for ct in config[key]}})
    pd.DataFrame(coverage_rows).to_csv(directory / 'shared_program_mapping.csv', index=False)
    pd.DataFrame(score_rows).to_csv(directory / 'harmonized_donor_program_scores.csv', index=False)
    program_frame = pd.DataFrame(program_summaries)
    program_frame.to_csv(directory / 'cross_study_program_context.csv', index=False)
    concordance = []
    for (program, comparison), group in program_frame.groupby(['program', 'comparison']):
        a = group.set_index('study')
        concordance.append({'program': program, 'comparison': comparison,
            'V1_relative_top': a.loc['Lukowski2019', 'relative_top_type'],
            'external_relative_top': a.loc['Wang2022', 'relative_top_type'],
            'same_relative_top': a.loc['Lukowski2019', 'relative_top_type'] == a.loc['Wang2022', 'relative_top_type']})
    pd.DataFrame(concordance).to_csv(directory / 'program_concordance.csv', index=False)
    evidence = pd.read_csv(stages['expanded'] / 'expanded_robustness.csv')
    rep = pd.read_csv(stages['replication'] / 'cross_study_context_concordance.csv')
    evidence = evidence.merge(rep.drop(columns='gene_id'), on='gene', validate='one_to_one')
    sampling = pd.read_csv(stages['sampling'] / 'sampling_stability.csv')
    primary_sampling = sampling.loc[~sampling['scenario'].eq('equal_cells_10_extended')].groupby('gene')['baseline_agreement_fraction'].min()
    extended_sampling = sampling.loc[sampling['scenario'].eq('equal_cells_10_extended')].set_index('gene')['baseline_agreement_fraction']
    evidence['primary_min_sampling_agreement'] = evidence['gene'].map(primary_sampling)
    evidence['extended_equal10_sampling_agreement'] = evidence['gene'].map(extended_sampling)
    evidence['descriptive_primary_context_flag'] = (evidence['LODO_all3'].fillna(False)
        & evidence['QC_same_top'] & evidence['primary_min_sampling_agreement'].ge(0.9)
        & evidence['primary_context_concordant'] & evidence['primary_external_LODO_agree'].eq(4)
        & evidence['primary_external_LODO_evaluable'].eq(4))
    # Keep DNA/RNA/chromatin evidence separate; do not invent a combined causal score.
    regulatory = pd.read_csv(stages['regulatory'] / 'variant_candidate_regulatory_evidence.csv')
    for source, label in [('any_published_gene_link', 'n_indices_with_published_gene_link'),
        ('published_HiChIP_target_match', 'n_indices_with_published_HiChIP_link'),
        ('published_predicted_target_match', 'n_indices_with_published_predicted_target'),
        ('published_retina_eQTL_target_match', 'n_indices_with_published_retina_eQTL'),
        ('extended_RNA_type_in_index_ATAC', 'n_indices_ATAC_matches_extended_RNA')]:
        counts = regulatory.loc[regulatory[source]].groupby('gene')['rsid'].nunique()
        evidence[label] = evidence['gene'].map(counts).fillna(0).astype(int)
    memberships = pd.read_csv(stages['programs'] / 'candidate_program_membership.csv')
    membership = memberships.loc[memberships['member_exact_symbol']].groupby('gene')['program'].agg(lambda x: ';'.join(sorted(set(x))))
    evidence['predefined_program_memberships'] = evidence['gene'].map(membership).fillna('')
    mito = pd.read_csv(stages['programs'] / 'candidate_mitocarta_annotation.csv')
    evidence = evidence.merge(mito.drop(columns='exact_symbol_in_matrix'), on='gene', validate='one_to_one')
    evidence.to_csv(directory / 'candidate_evidence_matrix.csv', index=False)
    evidence.loc[evidence['descriptive_primary_context_flag']].to_csv(directory / 'descriptive_context_shortlist.csv', index=False)
    cases = ['COL8A2', 'DGKG', 'NPC2', 'PSMC3', 'RP11-466F5.8', 'SLC2A12']
    evidence.set_index('gene').reindex(cases).reset_index().to_csv(directory / 'six_gene_case_studies.csv', index=False)
    # Compute index-variant/loop/gene-body overlap from the published HiChIP table.
    loops = pd.read_excel(root / 'data/raw/wang_mmc8.xlsx', header=1)
    coordinates = pd.read_csv(stages['regulatory'] / 'index_variant_coordinate_audit.csv')
    coord = coordinates.loc[coordinates['status'].eq('unique_GRCh38_position')].set_index('rsid')
    loop_rows = []
    for pair in regulatory.to_dict('records'):
        if pair['rsid'] not in coord.index:
            continue
        chrom = coord.loc[pair['rsid'], 'chromosome']; pos = int(coord.loc[pair['rsid'], 'GRCh38_position_1based']) - 1
        for left, right in [(1, 2), (2, 1)]:
            selected = loops.loc[loops[f'anchor_{left}_chr'].eq(chrom)
                & loops[f'anchor_{left}_start'].le(pos) & loops[f'anchor_{left}_end'].gt(pos)]
            for loop_id, loop in selected.iterrows():
                if pair['gene'] in tokens(loop[f'anchor_{right}_intersecting_genes']):
                    loop_rows.append({'rsid': pair['rsid'], 'gene': pair['gene'], 'published_loop_row': int(loop_id),
                        'variant_anchor': left, 'gene_anchor': right,
                        'scope': 'Variant overlaps one loop anchor; the other intersects the gene body. Does not establish promoter-specific regulation or causality.'})
    pd.DataFrame(loop_rows, columns=['rsid', 'gene', 'published_loop_row', 'variant_anchor', 'gene_anchor', 'scope']).to_csv(directory / 'computed_index_HiChIP_gene_anchor_overlaps.csv', index=False)
    # The recent atlas does not specify the interval endpoint convention in its TSV.
    # Validate that reported intersections survive plausible one-base conventions.
    recent = pd.read_csv(stages['regulatory'] / 'Yuan2026_index_cCRE_overlaps.csv')
    robust = []
    for row in recent.to_dict('records'):
        chrom, start, end = row['peak'].rsplit('-', 2)
        pos = int(coord.loc[row['rsid'], 'GRCh38_position_1based']) - 1
        row['interior_overlap_robust_to_1bp_conventions'] = int(start) <= pos < int(end)
        robust.append(row)
    recent['interior_overlap_robust_to_1bp_conventions'] = [r['interior_overlap_robust_to_1bp_conventions'] for r in robust]
    recent.to_csv(directory / 'Yuan2026_interval_convention_audit.csv', index=False)
    primary_programs = pd.DataFrame(concordance).loc[lambda d:d['comparison'].eq('primary')]
    save_json(directory / 'summary.json', {'descriptive_primary_context_flag_genes': int(evidence['descriptive_primary_context_flag'].sum()),
        'primary_program_top_concordance': int(primary_programs['same_relative_top'].sum()),
        'primary_programs_evaluable': len(primary_programs), 'computed_HiChIP_gene_anchor_pairs': len(loop_rows),
        'Yuan2026_boundary_robust_intersections': int(recent['interior_overlap_robust_to_1bp_conventions'].sum()),
        'context_flag_definition': 'All 3 V1 LODO folds agree; strict QC agrees; >=90% baseline agreement in all seven primary technical scenarios; external RNA context agrees in all four external LODO folds',
        'context_flag_scope': 'Descriptive screening rule, not a p-value, causal posterior or clinical decision threshold',
        'shared_program_scope': 'Same intersected gene/control pairs in both studies; relative transcriptional contrasts, not metabolic flux or disease-specific activation'})

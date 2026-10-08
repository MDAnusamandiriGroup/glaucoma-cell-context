"""Check coverage, use nominal df, and align literal regional variants."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

from coloc_core import derived_standard_error


def allele_key(chromosome, position, effect, baseline):
    effect, baseline = str(effect).upper(), str(baseline).upper()
    if not effect or not baseline or effect == baseline:
        return None
    if any(base not in 'ACGT' for allele in [effect, baseline] for base in allele):
        return None
    return str(chromosome).removeprefix('chr') + ':' + str(int(position)) + ':' + ':'.join(sorted([effect, baseline]))


def preprocess(root, plan, output):
    root, output = Path(root), Path(output)
    inputs = root / 'data' / 'processed'
    qtl = pd.read_csv(inputs / 'eyegex_target_nominal_raw.tsv', sep='\t')
    permutations = pd.read_csv(inputs / 'eyegex_target_permutations.tsv', sep='\t')
    source_audit, coverage_audit, exclusions = [], [], []
    for gene, gene_id in plan['genes'].items():
        source = qtl.loc[qtl.V1 == gene_id].copy()
        permutation = permutations.loc[permutations.V1 == gene_id]
        if len(permutation) != 1 or len(source) == 0:
            raise RuntimeError('Missing/duplicate source gene: ' + gene)
        declared = source.V6.unique()
        if len(declared) != 1 or len(source) != int(declared[0]) or len(source) != int(permutation.V6.iloc[0]):
            raise RuntimeError('Incomplete nominal gene statistics: ' + gene)
        if not source.V9.eq('chr' + plan['index_chromosome']).all():
            raise RuntimeError('Unexpected QTL chromosome')
        nominal_df = float(permutation.V12.iloc[0])
        if not np.isfinite(nominal_df) or nominal_df <= 0:
            raise RuntimeError('Missing nominal-test df')
        source_audit.append({'gene': gene, 'gene_id': gene_id, 'nominal_source_rows': len(source),
            'reported_n_var_in_cis': int(declared[0]), 'nominal_df_dof1': nominal_df,
            'optimized_df_dof2_NOT_USED': float(permutation.V13.iloc[0]),
            'source_nominal_p_min': float(source.V12.min()),
            'source_nominal_p_max': float(source.V12.max()),
            'source_nominal_p_gt_005': int((source.V12 > 0.05).sum()),
            'source_nominal_p_gt_05': int((source.V12 > 0.5).sum()),
            'gene_start_source_GRCh38': int(source.V3.iloc[0]),
            'source_full_row_count_validated': True})
        source['gene'] = gene
        source['variant_key'] = [allele_key(c, p, a, b) for c, p, a, b in zip(source.V9, source.V10, source.Effect_allele, source.Baseline_allele)]
        source['eqtl_beta'] = pd.to_numeric(source.V13)
        source['eqtl_p'] = pd.to_numeric(source.V12)
        source['nominal_df'] = nominal_df
        valid = source.variant_key.notna() & np.isfinite(source.eqtl_beta) & source.eqtl_beta.ne(0) & source.eqtl_p.between(0, 1, inclusive='neither')
        for _, row in source.loc[~valid].iterrows():
            exclusions.append({'gene': gene, 'study': 'EyeGEx', 'source_variant_id': row.V8, 'position': row.V10, 'reason': 'invalid_alleles_or_nonreconstructible_beta_p'})
        source = source.loc[valid].copy()
        if source.variant_key.duplicated().any():
            raise RuntimeError('Duplicate QTL coordinate/allele keys: ' + gene)
        source['eqtl_se_DERIVED'] = derived_standard_error(source.eqtl_beta, source.eqtl_p, source.nominal_df)
        reconstructed_p = 2 * t.sf(np.abs(source.eqtl_beta / source.eqtl_se_DERIVED), nominal_df)
        source['derived_se_p_roundtrip_relative_error'] = np.abs(reconstructed_p - source.eqtl_p) / source.eqtl_p
        if source.derived_se_p_roundtrip_relative_error.max() > 1e-7:
            raise RuntimeError('SE derivation p-value roundtrip failed: ' + gene)
        for study in [plan['primary_gwas'], plan['sensitivity_gwas']]:
            gwas = pd.read_csv(inputs / ('gwas_' + study + '_NPC2_region.tsv'), sep='\t')
            if not gwas.chromosome.eq(int(plan['index_chromosome'])).all():
                raise RuntimeError('Unexpected GWAS chromosome')
            gwas['variant_key'] = [allele_key(c, p, a, b) for c, p, a, b in zip(gwas.chromosome, gwas.base_pair_location, gwas.effect_allele, gwas.other_allele)]
            gvalid = gwas.variant_key.notna() & np.isfinite(gwas.beta) & np.isfinite(gwas.standard_error) & gwas.standard_error.gt(0) & gwas.p_value.between(0, 1, inclusive='both')
            gwas = gwas.loc[gvalid].copy()
            if gwas.variant_key.duplicated().any():
                raise RuntimeError('Duplicate GWAS coordinate/allele keys')
            low = plan['index_position'] - plan['primary_half_window_bp']
            high = plan['index_position'] + plan['primary_half_window_bp']
            region = source.loc[source.V10.between(low, high)].copy()
            position_set = set(gwas.base_pair_location)
            key_set = set(gwas.variant_key)
            for _, row in region.loc[~region.variant_key.isin(key_set)].iterrows():
                reason = 'literal_allele_mismatch' if row.V10 in position_set else 'position_absent_in_GWAS'
                exclusions.append({'gene': gene, 'study': study, 'source_variant_id': row.V8, 'position': row.V10, 'reason': reason})
            merged = region.merge(gwas, on='variant_key', how='inner', validate='one_to_one')
            if len(merged) < 1000:
                raise RuntimeError('Insufficient dense shared coverage: ' + gene + '/' + study)
            if not merged.V10.eq(merged.base_pair_location).all():
                raise RuntimeError('Position mismatch after allele-key match')
            same = merged.Effect_allele.eq(merged.effect_allele) & merged.Baseline_allele.eq(merged.other_allele)
            swap = merged.Effect_allele.eq(merged.other_allele) & merged.Baseline_allele.eq(merged.effect_allele)
            if not (same | swap).all():
                raise RuntimeError('Nonliteral strand match reached analysis')
            merged['gwas_beta_aligned'] = np.where(swap, -merged.beta, merged.beta)
            merged['gwas_beta_published'] = merged.beta
            merged['gwas_se_published'] = merged.standard_error
            merged['gwas_p_published'] = merged.p_value
            merged['eqtl_effect_allele'] = merged.Effect_allele
            merged['eqtl_other_allele'] = merged.Baseline_allele
            merged['position_GRCh38'] = merged.V10.astype(int)
            merged['eqtl_source_variant_id'] = merged.V8
            merged['gwas_rsid'] = merged.rsid
            merged['palindromic'] = [set([a, b]) in [{'A', 'T'}, {'C', 'G'}] for a, b in zip(merged.Effect_allele, merged.Baseline_allele)]
            merged['effect_orientation_ambiguous'] = merged.palindromic
            merged['gwas_effect_allele_swapped'] = swap
            columns = ['gene','variant_key','position_GRCh38','eqtl_source_variant_id','gwas_rsid','eqtl_effect_allele','eqtl_other_allele','eqtl_beta','eqtl_p','nominal_df','eqtl_se_DERIVED','gwas_beta_published','gwas_beta_aligned','gwas_se_published','gwas_p_published','palindromic','effect_orientation_ambiguous','gwas_effect_allele_swapped','derived_se_p_roundtrip_relative_error']
            result = merged[columns].sort_values(['position_GRCh38', 'variant_key'])
            if not result.gwas_rsid.eq(plan['index_variant']).any():
                raise RuntimeError('Prespecified index variant absent')
            result.to_csv(output / (study + '_' + gene + '.tsv'), sep='\t', index=False)
            coverage_audit.append({'gene': gene, 'study': study, 'qtl_eligible_region_rows': len(region),
                'gwas_eligible_region_rows': len(gwas), 'shared_literal_variants': len(result),
                'qtl_matching_fraction': len(result) / len(region), 'palindromic_rows': int(result.palindromic.sum()),
                'allele_swaps': int(result.gwas_effect_allele_swapped.sum()),
                'shared_min_position': int(result.position_GRCh38.min()), 'shared_max_position': int(result.position_GRCh38.max()),
                'max_p_roundtrip_error': float(result.derived_se_p_roundtrip_relative_error.max()),
                'shared_eqtl_min_p': float(result.eqtl_p.min()), 'shared_gwas_min_p': float(result.gwas_p_published.min())})
    pd.DataFrame(source_audit).to_csv(output / 'eqtl_source_audit.csv', index=False)
    pd.DataFrame(coverage_audit).to_csv(output / 'shared_coverage_audit.csv', index=False)
    pd.DataFrame(exclusions, columns=['gene','study','source_variant_id','position','reason']).to_csv(output / 'excluded_variants.csv', index=False)
    summary = {'genes': len(source_audit), 'gene_study_pairs': len(coverage_audit), 'all_nominal_gene_row_counts_match': True,
               'eqtl_SE': 'derived from beta, nominal p and published dof1; not a supplied native SE',
               'p_value_significance_filter': False, 'exclusion_rows': len(exclusions)}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary

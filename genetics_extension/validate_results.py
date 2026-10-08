"""Validate saved inputs and outputs without recomputing the fitted analysis."""
import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DEFAULT_INDEX = 'results/run_index_b5a8e45b71c2f2b4.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(index_path=DEFAULT_INDEX):
    index_path = Path(index_path)
    if not index_path.is_absolute():
        index_path = ROOT / index_path
    index = json.loads(index_path.read_text())
    registry = json.loads((ROOT / 'data/source_registry.json').read_text())
    for record in registry['derived_inputs']:
        assert sha(ROOT / record['file']) == record['sha256'], record['file']
    for stage in ['preprocess', 'fit']:
        folder = ROOT / index[stage]
        assert sha(folder / 'manifest.json') == index[stage + '_manifest_sha256']
        manifest = json.loads((folder / 'manifest.json').read_text())
        assert manifest['status'] == 'complete'
        for relative, digest in manifest['outputs'].items():
            assert sha(folder / relative) == digest, relative
    prep, fit = ROOT / index['preprocess'], ROOT / index['fit']
    primary = pd.read_csv(fit / 'primary_settings_posteriors.csv')
    sensitivity = pd.read_csv(fit / 'all_sensitivity_posteriors.csv')
    columns = ['PP_H' + str(i) for i in range(5)]
    probabilities = sensitivity[columns].to_numpy()
    assert len(primary) == 6 and len(sensitivity) == 216
    assert np.isfinite(probabilities).all()
    assert ((probabilities >= 0) & (probabilities <= 1)).all()
    np.testing.assert_allclose(probabilities.sum(axis=1), 1, atol=2e-12, rtol=0)
    assert sensitivity.groupby(['gene', 'study']).size().eq(36).all()
    numerical_errors = []
    for row in primary.itertuples():
        matched = pd.read_csv(prep / (row.study + '_' + row.gene + '.tsv'), sep='\t')
        variants = pd.read_csv(fit / (row.study + '_' + row.gene + '_snp_posteriors.csv'))
        assert len(matched) == len(variants) == row.shared_variants
        assert matched.variant_key.is_unique and variants.variant_key.is_unique
        assert matched.gwas_rsid.eq('rs754458').any()
        assert np.isfinite(matched.eqtl_se_DERIVED).all() and matched.eqtl_se_DERIVED.gt(0).all()
        assert matched.nominal_df.eq(404).all()
        assert matched.derived_se_p_roundtrip_relative_error.max() < 1e-7
        expected_aligned = np.where(matched.gwas_effect_allele_swapped,
                                    -matched.gwas_beta_published, matched.gwas_beta_published)
        np.testing.assert_allclose(matched.gwas_beta_aligned, expected_aligned, rtol=1e-14)
        mass = variants.SNP_PP_H4_CONDITIONAL.to_numpy()
        np.testing.assert_allclose(mass.sum(), 1, atol=2e-12, rtol=0)
        assert np.all(np.diff(mass) <= 1e-15)
        chosen = variants.loc[variants.in_95pct_set_CONDITIONAL_ON_H4, 'SNP_PP_H4_CONDITIONAL']
        assert len(chosen) == row.conditional_95pct_set_variants
        assert chosen.sum() >= 0.95 and chosen.iloc[:-1].sum() < 0.95
        # Independent arithmetic reference: ordinary long-double Bayes factors
        # and the distinct-pair aggregate, rather than the production log sums.
        # The five small-fixture tests separately exercise explicit enumeration
        # and extreme cases where subtraction would lose precision.
        a = np.exp(variants.GWAS_log_ABF.to_numpy().astype(np.longdouble))
        b = np.exp(variants.eQTL_log_ABF.to_numpy().astype(np.longdouble))
        shared = np.dot(a, b)
        different = a.sum() * b.sum() - shared
        assert different > 0
        weights = np.array([1, row.p1*a.sum(), row.p2*b.sum(),
                            row.p1*row.p2*different, row.p12*shared], dtype=np.longdouble)
        reference = weights / weights.sum()
        actual = np.array([getattr(row, col) for col in columns])
        np.testing.assert_allclose(reference.astype(float), actual, atol=2e-12, rtol=1e-10)
        numerical_errors.append(float(np.max(np.abs(reference.astype(float) - actual))))
    audit = pd.read_csv(prep / 'eqtl_source_audit.csv')
    assert audit.nominal_source_rows.eq(audit.reported_n_var_in_cis).all()
    assert audit.source_nominal_p_gt_005.gt(0).all()
    assert audit.nominal_source_rows.sum() == 17507
    for source in ROOT.rglob('*.py'):
        ast.parse(source.read_text(), filename=str(source))
    return {'status': 'passed', 'gene_study_pairs': 6, 'sensitivity_settings': 216,
            'source_nominal_gene_rows': 17507, 'all_input_and_stage_output_hashes_match': True,
            'posterior_normalization': 'passed', 'conditional_H4_sets': 'passed',
            'effect_allele_orientation': 'passed', 'native_vs_derived_SE_labels': 'passed',
            'independent_primary_posterior_reference_max_absolute_error': max(numerical_errors),
            'R_coloc_package_comparison_performed': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', default=DEFAULT_INDEX)
    parser.add_argument('--write-json')
    args = parser.parse_args()
    result = validate(args.index)
    if args.write_json:
        target = Path(args.write_json)
        target.parent.mkdir(parents=True, exist_ok=True)
        data = json.dumps(result, indent=2) + '\n'
        if target.exists() and target.read_text() != data:
            raise RuntimeError('Preserve existing incompatible validation: ' + str(target))
        if not target.exists():
            target.write_text(data)
    print(json.dumps(result))

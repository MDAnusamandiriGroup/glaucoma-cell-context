"""Frozen ABF posterior and region/prior/allele-class sensitivities."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from coloc_core import log_abf, combine_log_abf


def fit(preprocessed, plan, output):
    preprocessed, output = Path(preprocessed), Path(output)
    rows, primary_snp_tables = [], []
    for study in [plan['primary_gwas'], plan['sensitivity_gwas']]:
        for gene in plan['genes']:
            source = pd.read_csv(preprocessed / (study + '_' + gene + '.tsv'), sep='\t')
            for half_window in [plan['primary_half_window_bp'], plan['sensitivity_half_window_bp']]:
                for allele_policy in ['exact_all', 'exclude_palindromic']:
                    subset = source.loc[source.position_GRCh38.between(plan['index_position'] - half_window, plan['index_position'] + half_window)].copy()
                    if allele_policy == 'exclude_palindromic':
                        subset = subset.loc[~subset.palindromic].copy()
                    if len(subset) < 1000:
                        raise RuntimeError('Insufficient variants for frozen sensitivity')
                    l1 = log_abf(subset.gwas_beta_aligned, subset.gwas_se_published, plan['gwas_prior_sd_log_odds'])
                    for eqtl_prior_sd in plan['eqtl_prior_sd_grid']:
                        l2 = log_abf(subset.eqtl_beta, subset.eqtl_se_DERIVED, eqtl_prior_sd)
                        for p12 in plan['p12_grid']:
                            posterior, conditional_snp_h4 = combine_log_abf(l1, l2, plan['primary_p1'], plan['primary_p2'], p12)
                            primary_setting = half_window == plan['primary_half_window_bp'] and allele_policy == 'exact_all' and eqtl_prior_sd == plan['primary_eqtl_prior_sd'] and p12 == plan['primary_p12']
                            winner = int(np.argmax(posterior))
                            conditional_ratio = posterior[4] / (posterior[3] + posterior[4])
                            rows.append({'gene':gene,'study':study,'primary_GWAS':study==plan['primary_gwas'],'primary_settings':primary_setting,
                                'half_window_bp':half_window,'allele_policy':allele_policy,'shared_variants':len(subset),'p1':plan['primary_p1'],
                                'p2':plan['primary_p2'],'p12':p12,'eqtl_prior_sd':eqtl_prior_sd,'gwas_prior_sd':plan['gwas_prior_sd_log_odds'],
                                **{'PP_H'+str(i):float(posterior[i]) for i in range(5)},'PP_H4_given_H3_or_H4':float(conditional_ratio),
                                'highest_posterior_hypothesis':'H'+str(winner),'H4_ge_08':bool(posterior[4]>=plan['decision_reference_H4'])})
                            if primary_setting:
                                variants=subset.copy()
                                variants['GWAS_log_ABF']=l1
                                variants['eQTL_log_ABF']=l2
                                variants['SNP_PP_H4_CONDITIONAL']=conditional_snp_h4
                                variants=variants.sort_values('SNP_PP_H4_CONDITIONAL',ascending=False)
                                variants['conditional_H4_cumulative_probability']=variants.SNP_PP_H4_CONDITIONAL.cumsum()
                                before=variants.conditional_H4_cumulative_probability-variants.SNP_PP_H4_CONDITIONAL
                                variants['in_95pct_set_CONDITIONAL_ON_H4']=before.lt(0.95)
                                variants.to_csv(output/(study+'_'+gene+'_snp_posteriors.csv'),index=False)
                                primary_snp_tables.append({'gene':gene,'study':study,'top_H4_conditional_variant':variants.variant_key.iloc[0],
                                    'top_H4_conditional_rsid':variants.gwas_rsid.iloc[0],
                                    'top_H4_conditional_probability':float(variants.SNP_PP_H4_CONDITIONAL.iloc[0]),
                                    'conditional_95pct_set_variants':int(variants.in_95pct_set_CONDITIONAL_ON_H4.sum())})
    results=pd.DataFrame(rows)
    assert len(results)==216
    primary=results.loc[results.primary_settings].copy()
    primary=primary.merge(pd.DataFrame(primary_snp_tables),on=['gene','study'],validate='one_to_one')
    for gene in plan['genes']:
        for study in [plan['primary_gwas'],plan['sensitivity_gwas']]:
            selection=(primary.gene==gene)&(primary.study==study)
            sensitivity=results.loc[(results.gene==gene)&(results.study==study)]
            primary.loc[selection,'H4_sensitivity_min']=float(sensitivity.PP_H4.min())
            primary.loc[selection,'H4_sensitivity_max']=float(sensitivity.PP_H4.max())
            primary.loc[selection,'H4_ge_08_all_sensitivities']=bool(sensitivity.H4_ge_08.all())
    results.to_csv(output/'all_sensitivity_posteriors.csv',index=False)
    primary.to_csv(output/'primary_settings_posteriors.csv',index=False)
    summary={'gene_study_pairs':6,'posterior_sensitivity_settings':len(results),
        'method':'Python unsigned Wakefield ABF / coloc uniform-prior equations',
        'SuSiE_coloc_run':False,'MR_run':False,
        'strong_primary_European_H4_genes':primary.loc[primary.primary_GWAS&primary.H4_ge_08,'gene'].tolist(),
        'robust_primary_European_H4_genes':primary.loc[primary.primary_GWAS&primary.H4_ge_08_all_sensitivities,'gene'].tolist()}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

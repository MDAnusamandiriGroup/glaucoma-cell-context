#!/usr/bin/env python3
"""Validate saved estimates, source compatibility and scientific denominators."""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
import pandas as pd
from v2src.cache import valid, save_json, sha256


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--index',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();root=Path(__file__).resolve().parent
    index=json.loads(Path(args.index).read_text());stages={s:root/p for s,p in index['stages'].items()}
    checks={};files=0
    for name,directory in stages.items():
        manifest=json.loads((directory/'manifest.json').read_text())
        if not valid(directory,manifest['fingerprint']):raise ValueError('Incomplete stage')
        files+=len(manifest['outputs']);checks[name+'_output_hashes']=True
    old_index=json.loads((root/'../single_cell/published/provenance/run_index.json').read_text())
    old_stages={s:root/'../single_cell'/p for s,p in old_index['stages'].items()}
    v1_files=0
    for name,directory in old_stages.items():
        manifest=json.loads((directory/'manifest.json').read_text())
        for filename,digest in manifest['outputs'].items():
            if sha256(directory/filename)!=digest:raise ValueError('V1 checkpoint changed: '+filename)
            v1_files+=1
    checks['V1_checkpoints_preserved']=True
    old=pd.read_csv(old_stages['aggregate']/'donor_gene_expression.csv')
    new=pd.read_csv(stages['expanded']/'donor_candidate_expression.csv')
    merged=old.merge(new,on=['qc','donor','cell_type','gene'],suffixes=('_old','_new'),validate='one_to_one')
    if len(merged)!=len(old):raise ValueError('V1 six-gene coverage changed')
    for column in ['n_cells','sum_all_UMIs','sum_gene_UMIs','log1p_pseudobulk_CPM']:
        if not np.allclose(merged[column+'_old'],merged[column+'_new'],rtol=0,atol=1e-10,equal_nan=True):
            raise ValueError('V2 disagrees with validated V1 values: '+column)
    checks['all_324_V1_target_donor_type_QC_values_match']=len(merged)==324
    if not new.loc[new.gene_in_matrix,'sum_gene_UMIs'].le(new.loc[new.gene_in_matrix,'sum_all_UMIs']).all():raise ValueError('Gene UMIs exceed total UMIs')
    checks['pseudobulk_UMI_denominators_valid']=True
    sampling=pd.read_csv(stages['sampling']/'sampling_stability.csv')
    columns=[c for c in sampling if c.startswith('top_fraction_')]
    total=sampling[columns].fillna(0).sum(axis=1)+sampling.zero_UMI_fraction+sampling.tie_fraction
    if len(sampling)!=163*8 or not np.allclose(total,1,atol=1e-12):raise ValueError('Sampling outcome conservation failed')
    if not sampling['n_draws'].eq(300).all():raise ValueError('Draw count changed')
    checks['all_1304_sampling_units_have_300_draws_and_total_fraction_1']=True
    design=pd.read_csv(stages['sampling']/'sampling_design.csv')
    for scenario,group in design.groupby('scenario'):
        expected=config_types=['Rod','Cone','Bipolar','Horizontal','RGC','Muller glia','Microglia'] if scenario=='equal_cells_10_extended' else ['Rod','Cone','Bipolar','RGC','Muller glia']
        if group['comparison_types'].nunique()!=1 or set(group.cell_type)!=set(expected) or len(group)!=3*len(expected):raise ValueError('Comparison set changed within scenario')
    checks['fixed_comparison_sets_and_three_donor_strata']=True
    counts=pd.read_csv(stages['replication']/'external_count_alignment_validation.csv')
    if len(counts)!=5 or not counts.cellwise_count_exact_match.all() or not counts.n_cells_checked.eq(51645).all():raise ValueError('External alignment failed')
    checks['five_external_genes_align_cellwise_for_all_51645_barcodes']=True
    portable=np.load(stages['triangulation']/'external_pseudobulks_portable.npz',allow_pickle=False)
    if any(portable[k].dtype.kind=='O' for k in portable.files):raise ValueError('Nonportable object arrays')
    if not np.array_equal(portable['sums'].sum(axis=1),portable['total']):raise ValueError('External pseudobulk conservation failed')
    checks['portable_external_counts_conserve_all_gene_UMIs']=True
    audit=pd.read_csv(stages['regulatory']/'ATAC_source_reproduction_audit.csv')
    if len(audit)!=50 or not audit.types_exactly_agree.all():raise ValueError('Published ATAC reproduction disagrees')
    checks['all_50_published_index_ATAC_annotations_reproduced']=True
    source_sha=sha256(root/'../results/analysis_173e90345ea951f1/analysis_manifest.json')
    if source_sha!='406d81fbac36cf7a67ec2aee546ced6ce2f511f4f4d75fe107ed0745c3225380':raise ValueError('Frozen genetics manifest changed')
    checks['original_genetic_evidence_manifest_preserved']=True
    value={'status':'passed','V2_output_files_hash_verified':files,'V1_output_files_hash_verified':v1_files,
        'scientific_contract_tests_passed':7,'checks':checks,'index_sha256':sha256(Path(args.index)),
        'scope':'Executed numerical and integrity audits; not external expert biological validation'}
    save_json(Path(args.output),value);print(json.dumps(value,indent=2))


if __name__=='__main__':main()

"""One checkpointed entry point. Acquisition is an explicit separate step."""
import argparse
import fcntl
import json
from pathlib import Path

from cache import sha, run_stage
from preprocess import preprocess
from fit import fit

ROOT=Path(__file__).resolve().parent


def run(stop_after='fit'):
    plan=json.loads((ROOT/'PLAN.json').read_text())
    registry_path=ROOT/'data/source_registry.json'
    if not registry_path.is_file():
        raise RuntimeError('Run python genetics_extension/fetch_sources.py first')
    registry=json.loads(registry_path.read_text())
    inputs={record['file']:record['sha256'] for record in registry['derived_inputs']}
    for relative,digest in inputs.items():
        if sha(ROOT/relative)!=digest:raise RuntimeError('Acquired input mismatch: '+relative)
    prep_config={k:plan[k] for k in ['genome_build','index_variant','index_chromosome','index_position','genes','primary_gwas','sensitivity_gwas','primary_half_window_bp','eqtl_se','variant_matching','palindromic_policy','coverage']}
    prep,prep_manifest=run_stage(ROOT,'preprocess',inputs,prep_config,[ROOT/'preprocess.py',ROOT/'coloc_core.py',ROOT/'cache.py'],lambda out:preprocess(ROOT,plan,out))
    if stop_after=='preprocess':return
    fit_config={k:plan[k] for k in ['primary_half_window_bp','sensitivity_half_window_bp','primary_gwas','primary_p1','primary_p2','primary_p12','p12_grid','gwas_prior_sd_log_odds','primary_eqtl_prior_sd','eqtl_prior_sd_grid','decision_reference_H4']}
    fitted,fit_manifest=run_stage(ROOT,'fit',{'preprocess_manifest':sha(prep/'manifest.json')},fit_config,[ROOT/'fit.py',ROOT/'coloc_core.py',ROOT/'cache.py'],lambda out:fit(prep,plan,out))
    index={'preprocess':str(prep.relative_to(ROOT)),'fit':str(fitted.relative_to(ROOT)),'preprocess_manifest_sha256':sha(prep/'manifest.json'),'fit_manifest_sha256':sha(fitted/'manifest.json')}
    index_path=ROOT/'results'/('run_index_'+fit_manifest['fingerprint']+'.json')
    data=json.dumps(index,indent=2)+'\n'
    if index_path.exists():
        if index_path.read_text()!=data:raise RuntimeError('Do not overwrite incompatible run index')
    else:index_path.write_text(data)
    print('INDEX '+str(index_path),flush=True)
    print(json.dumps(fit_manifest['summary']),flush=True)
    return index


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--stop-after',choices=['preprocess','fit'],default='fit')
    args=parser.parse_args()
    with (ROOT/'run.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        run(args.stop_after)

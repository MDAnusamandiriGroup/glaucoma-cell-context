#!/usr/bin/env python3
"""Execute V2 downstream modules without changing validated V1 checkpoints."""
from pathlib import Path
import argparse
import json
import os
import sys
import importlib
import fcntl

for variable in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ.setdefault(variable, '4')

from v2src.cache import checkpoint, save_json, sha256, software


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stop-after', choices=['expanded', 'sampling', 'programs', 'replication', 'regulatory', 'triangulation', 'report'], default='report')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    config = json.loads((root / 'config.json').read_text())
    results = root / 'results'; results.mkdir(exist_ok=True)
    lock = root / '.v2.lock'
    lock_handle = lock.open('a+')
    fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    lock_handle.seek(0); lock_handle.truncate(); lock_handle.write(str(os.getpid())); lock_handle.flush()
    stages = {}
    source = {'counts': sha256(root / config['source_v1_normalized']),
              'candidates': sha256(root / config['source_candidates']),
              'v1_index': sha256(root / config['source_v1_index'])}
    environment = software()
    try:
        for name in ['expanded', 'sampling', 'programs', 'replication', 'regulatory', 'triangulation', 'report']:
            module_path = root / 'v2src' / (name + '.py')
            module = importlib.import_module('v2src.' + name)
            provenance = {'config': config, 'source': source, 'software': environment,
                'module_sha256': sha256(module_path), 'cache_code_sha256': sha256(root / 'v2src/cache.py'),
                'dependencies': {s: sha256(p / 'manifest.json') for s, p in stages.items()}}
            if name == 'programs':
                provenance['hallmark_sha256'] = sha256(root / 'data/raw/hallmark2025.gmt')
                provenance['mitocarta_sha256'] = sha256(root / 'data/raw/Human.MitoCarta3.0.xls')
            if name == 'replication':
                provenance['external_inputs'] = {n: sha256(root / 'data/raw' / n) for n in [
                    'wang_metadata.csv', 'wang_cellxgene.h5ad', 'GSM5866081_barcodes.tsv.gz',
                    'GSM5866081_features.tsv.gz', 'GSM5866081_matrix.mtx.gz']}
            if name == 'regulatory':
                provenance['regulatory_inputs'] = {n: sha256(root / 'data/raw' / n) for n in [
                    'wang_mmc3.xlsx', 'wang_mmc4.xlsx', 'wang_mmc7.xlsx', 'wang_mmc9.xlsx',
                    'ensembl_index_variants_GRCh38.json', 'yuan_table_S6.tsv']}
                provenance['regulatory_plan_sha256'] = sha256(root / 'REGULATORY_PLAN.json')
                provenance['genetic_pairs_sha256'] = sha256(root / '../results/analysis_173e90345ea951f1/tables/colocalization_records.csv')
            if name == 'triangulation':
                provenance['HiChIP_loops_sha256'] = sha256(root / 'data/raw/wang_mmc8.xlsx')
            if name == 'report':
                provenance['gene_annotation_sha256'] = sha256(root / 'data/raw/gencode.v32.annotation.gtf.gz')
                provenance['rendering_versions'] = {n: __import__('importlib').metadata.version(n) for n in ['matplotlib', 'reportlab', 'Pillow']}
            if name == 'expanded':
                compute = lambda directory: module.run(root, directory, config)
            elif name in ['sampling', 'programs']:
                compute = lambda directory: module.run(root, stages['expanded'], directory, config)
            else:
                compute = lambda directory: module.run(root, stages, directory, config)
            stages[name] = checkpoint(results, name, provenance, compute)
            if name == args.stop_after:
                break
        index = {'config_sha256': sha256(root / 'config.json'), 'v1_inputs_sha256': source,
            'software': environment, 'stages': {s: str(p.relative_to(root)) for s, p in stages.items()},
            'scope': 'Executed Python extension; V1 checkpoints are read-only inputs'}
        digest = __import__('hashlib').sha256(json.dumps(index, sort_keys=True).encode()).hexdigest()[:16]
        target = results / ('run_index_' + digest + '.json')
        if not target.exists():
            save_json(target, index)
        print('INDEX ' + str(target), flush=True)
    finally:
        lock_handle.seek(0); lock_handle.truncate(); lock_handle.write('completed\n'); lock_handle.flush()
        fcntl.flock(lock_handle, fcntl.LOCK_UN)
        lock_handle.close()


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Run/resume the cell-level retina extension without changing the earlier pilot."""
from pathlib import Path
import argparse
import hashlib
import inspect
import json
import os

for variable in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ.setdefault(variable, '4')
os.environ.setdefault('MPLBACKEND', 'Agg')

from scsrc import ingest, preprocess, model
from scsrc.common import checkpoint, file_hash, fingerprint, save_json, software


def function_hash(function):
    return hashlib.sha256(inspect.getsource(function).encode()).hexdigest()


def dependency(directory, filename):
    manifest = json.loads((Path(directory) / 'manifest.json').read_text())
    return {'stage_fingerprint': manifest['fingerprint'], 'output_sha256': file_hash(Path(directory) / filename)}


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true', help='Use cached source files only; no network calls.')
    parser.add_argument('--stop-after', choices=['ingest', 'preprocess', 'pca', 'graph', 'embedding', 'aggregate', 'report'], default='report')
    parser.add_argument('--rebuild-report', action='store_true', help='Render a new report attempt from saved valid analysis; preserve earlier figures.')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    config = json.loads((root / 'config.json').read_text())
    versions = software()
    sources = ingest.acquire(root, config, args.offline)
    results = root / 'results'
    results.mkdir(exist_ok=True)
    lock = results / '.pipeline.lock'
    try:
        lock.mkdir()
    except FileExistsError:
        raise RuntimeError('Another run may be active. Inspect single_cell/results/.pipeline.lock/owner.json and confirm its process has stopped before removing the lock.')
    stages = {}
    try:
        save_json(lock / 'owner.json', {'pid': os.getpid(), 'entry_point': str(Path(__file__).resolve())})
        general = {'schema': config['schema_version'], 'software': versions}
        stages['ingest'] = checkpoint(results, '01_ingest', {
            **general, 'inputs': sources, 'code_sha256': file_hash(root / 'scsrc/ingest.py')
        }, lambda directory: ingest.run(root, directory, config))
        if args.stop_after == 'ingest':
            return
        stages['preprocess'] = checkpoint(results, '02_preprocess', {
            **general, 'input': dependency(stages['ingest'], 'counts.h5ad'),
            'qc': config['qc'], 'normalization_target': config['normalization_target'], 'hvg_min_cells': config['hvg_min_cells'],
            'code_sha256': file_hash(root / 'scsrc/preprocess.py')
        }, lambda directory: preprocess.run(stages['ingest'], directory, config))
        if args.stop_after == 'preprocess':
            return
        stages['pca'] = checkpoint(results, '03_pca', {
            **general, 'input': dependency(stages['preprocess'], 'normalized.h5ad'),
            'parameters': {key: config[key] for key in ['hvg_n_genes', 'hvg_batch_key', 'n_pcs', 'seed']},
            'function_sha256': function_hash(model.pca)
        }, lambda directory: model.pca(stages['preprocess'], directory, config))
        if args.stop_after == 'pca':
            return
        stages['graph'] = checkpoint(results, '04_graph', {
            **general, 'input': dependency(stages['pca'], 'pca.h5ad'),
            'parameters': {key: config[key] for key in ['n_neighbors', 'n_pcs', 'leiden_resolutions', 'seed']},
            'function_sha256': function_hash(model.graph)
        }, lambda directory: model.graph(stages['pca'], directory, config))
        if args.stop_after == 'graph':
            return
        stages['embedding'] = checkpoint(results, '05_embedding', {
            **general, 'input': dependency(stages['graph'], 'graph.h5ad'),
            'parameters': {key: config[key] for key in ['umap_min_dist', 'umap_iterations', 'seed']},
            'function_sha256': function_hash(model.embedding)
        }, lambda directory: model.embedding(stages['graph'], directory, config))
        if args.stop_after == 'embedding':
            return
        from scsrc import aggregate
        stages['aggregate'] = checkpoint(results, '06_aggregate', {
            **general, 'expression': dependency(stages['preprocess'], 'normalized.h5ad'),
            'clusters': dependency(stages['graph'], 'graph.h5ad'),
            'candidates_sha256': file_hash(root / 'inputs/published_candidates.csv'),
            'parameters': {key: config[key] for key in ['target_genes', 'cell_types', 'markers', 'min_cells_per_donor_type', 'min_cells_sensitivity']},
            'code_sha256': file_hash(root / 'scsrc/aggregate.py')
        }, lambda directory: aggregate.run(root, stages['preprocess'], stages['graph'], directory, config))
        if args.stop_after == 'aggregate':
            return
        from scsrc import reporting
        stages['report'] = checkpoint(results, '07_report', {
            'schema': config['schema_version'],
            'inputs': {name: dependency(path, 'manifest.json') for name, path in stages.items()},
            'reporting_sha256': file_hash(root / 'scsrc/reporting.py'),
            'matplotlib': versions['matplotlib']
        }, lambda directory: reporting.run(stages, directory, config), rebuild=args.rebuild_report)
        index = {'status': 'complete', 'stages': {name: str(path.relative_to(root)) for name, path in stages.items()}, 'software': versions}
        index_path = results / ('run_index_' + fingerprint(index) + '.json')
        if not index_path.exists():
            save_json(index_path, index)
        print('COMPLETE report: ' + str(stages['report'] / 'REPORT.md'), flush=True)
    finally:
        (lock / 'owner.json').unlink(missing_ok=True)
        lock.rmdir()


if __name__ == '__main__':
    run()

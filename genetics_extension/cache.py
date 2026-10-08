"""Immutable, validated stage checkpoints. No expensive work on import."""
import hashlib
import json
import os
from pathlib import Path
import platform

import numpy
import pandas
import scipy


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def versions():
    return {'python': platform.python_version(), 'numpy': numpy.__version__,
            'pandas': pandas.__version__, 'scipy': scipy.__version__}


def run_stage(root, name, inputs, configuration, code_paths, producer):
    contract = {'schema_version': 1, 'stage': name, 'inputs': inputs,
                'configuration': configuration, 'software': versions(),
                'code': {Path(p).name: sha(p) for p in code_paths}}
    fingerprint = hashlib.sha256(json.dumps(contract, sort_keys=True).encode()).hexdigest()[:16]
    directory = Path(root) / 'results' / (name + '_' + fingerprint)
    manifest_path = directory / 'manifest.json'
    if directory.exists():
        if not manifest_path.exists():
            raise RuntimeError('Preserve incomplete checkpoint for inspection: ' + str(directory))
        manifest = json.loads(manifest_path.read_text())
        if manifest.get('status') != 'complete' or manifest.get('contract') != contract:
            raise RuntimeError('Checkpoint contract mismatch: ' + str(directory))
        for relative, digest in manifest['outputs'].items():
            file = directory / relative
            if not file.is_file() or sha(file) != digest:
                raise RuntimeError('Checkpoint checksum mismatch: ' + str(file))
        print('REUSE ' + directory.name, flush=True)
        return directory, manifest
    directory.mkdir(parents=True, exist_ok=False)
    try:
        summary = producer(directory)
        files = sorted(p for p in directory.rglob('*') if p.is_file())
        if not files or any('.partial' in p.name for p in files):
            raise RuntimeError('Incomplete stage outputs')
        outputs = {str(p.relative_to(directory)): sha(p) for p in files}
        manifest = {'status': 'complete', 'fingerprint': fingerprint,
                    'contract': contract, 'outputs': outputs, 'summary': summary}
        temporary = manifest_path.with_name('manifest.json.partial')
        with temporary.open('x') as stream:
            stream.write(json.dumps(manifest, indent=2) + '\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.rename(temporary, manifest_path)
        print('COMPLETE ' + directory.name, flush=True)
        return directory, manifest
    except Exception as error:
        failure_path = directory / 'FAILED.json'
        if not failure_path.exists():
            failure_path.write_text(json.dumps({'status': 'failed', 'type': type(error).__name__, 'message': str(error), 'contract': contract}, indent=2) + '\n')
        raise

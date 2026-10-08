from pathlib import Path
import hashlib
import importlib.metadata
import json
import os
import platform
import tempfile
import time


def file_hash(path, algorithm='sha256'):
    digest = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        while True:
            block = handle.read(4 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def fingerprint(value):
    payload = json.dumps(value, sort_keys=True, separators=(',', ':')).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def publish(path, writer):
    """Atomic publication which cannot replace an existing file."""
    path = Path(path)
    if path.exists():
        raise FileExistsError('Preserving existing output: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=path.stem + '_', suffix='.partial' + path.suffix, dir=path.parent)
    os.close(fd)
    temp = Path(name)
    try:
        writer(temp)
        if not temp.is_file() or temp.stat().st_size == 0:
            raise ValueError('Writer produced an empty output: ' + str(path))
        os.link(temp, path)
        temp.unlink()
    except Exception:
        # Preserve incomplete bytes for diagnosis; never publish them as valid.
        raise


def save_json(path, value):
    publish(path, lambda temp: temp.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8'))


def save_csv(path, frame):
    publish(path, lambda temp: frame.to_csv(temp, index=False))


def save_adata(path, data):
    publish(path, lambda temp: data.write_h5ad(temp, compression='gzip'))


def validate_checkpoint(directory, expected):
    directory = Path(directory)
    manifest_path = directory / 'manifest.json'
    if not manifest_path.exists():
        return False
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('status') != 'complete' or manifest.get('fingerprint') != expected:
        raise ValueError('Checkpoint provenance mismatch: ' + str(directory))
    if not manifest.get('outputs'):
        raise ValueError('Checkpoint has no outputs: ' + str(directory))
    for name, digest in manifest['outputs'].items():
        path = directory / name
        if not path.is_file() or file_hash(path) != digest:
            raise ValueError('Corrupt/changed checkpoint; preserving it for diagnosis: ' + str(path))
    return True


def checkpoint(results, name, provenance, compute, rebuild=False):
    fp = fingerprint(provenance)
    base = Path(results) / (name + '_' + fp)
    if not rebuild:
        for attempt in range(100):
            directory = base if attempt == 0 else Path(str(base) + '_attempt' + str(attempt + 1))
            if directory.exists() and validate_checkpoint(directory, fp):
                print('REUSE ' + name + ': ' + directory.name, flush=True)
                return directory
    for attempt in range(100):
        directory = base if attempt == 0 else Path(str(base) + '_attempt' + str(attempt + 1))
        if not directory.exists():
            break
    else:
        raise RuntimeError('Too many incomplete stage attempts: ' + str(base))
    directory.mkdir(parents=True)
    print('RUN ' + name + ': ' + directory.name, flush=True)
    started = time.time()
    try:
        compute(directory)
        outputs = {}
        for path in sorted(directory.rglob('*')):
            if path.is_file() and not path.name.endswith('.partial'):
                outputs[str(path.relative_to(directory))] = file_hash(path)
        if not outputs:
            raise RuntimeError('Stage did not produce outputs: ' + name)
        save_json(directory / 'manifest.json', {
            'status': 'complete', 'fingerprint': fp, 'stage': name,
            'provenance': provenance, 'outputs': outputs,
            'elapsed_seconds': round(time.time() - started, 3)
        })
        validate_checkpoint(directory, fp)
    except Exception as exc:
        save_json(directory / 'failed.json', {'status': 'failed', 'error': repr(exc), 'stage': name})
        raise
    print('DONE ' + name + ': ' + str(round(time.time() - started, 1)) + ' seconds', flush=True)
    return directory


def software():
    result = {'python': platform.python_version()}
    for name in ['scanpy', 'anndata', 'numpy', 'scipy', 'pandas', 'scikit-learn', 'numba', 'umap-learn', 'igraph', 'leidenalg', 'h5py', 'matplotlib']:
        result[name] = importlib.metadata.version(name)
    return result

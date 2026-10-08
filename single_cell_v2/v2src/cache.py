from pathlib import Path
import hashlib
import json
import os
import platform
import importlib.metadata
import time


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def save_json(path, value):
    path = Path(path)
    if path.exists():
        raise FileExistsError(str(path))
    temporary = path.with_name(path.name + '.partial')
    with temporary.open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')
    os.rename(temporary, path)


def valid(directory, fingerprint):
    path = Path(directory) / 'manifest.json'
    if not path.exists():
        return False
    manifest = json.loads(path.read_text())
    if manifest['status'] != 'complete' or manifest['fingerprint'] != fingerprint:
        raise ValueError('Checkpoint provenance mismatch: ' + str(directory))
    for filename, digest in manifest['outputs'].items():
        if sha256(Path(directory) / filename) != digest:
            raise ValueError('Checkpoint content changed: ' + filename)
    return True


def checkpoint(root, name, provenance, compute):
    fingerprint = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()[:16]
    base = Path(root) / (name + '_' + fingerprint)
    for attempt in range(100):
        directory = base if attempt == 0 else Path(str(base) + '_attempt' + str(attempt + 1))
        if directory.exists():
            if valid(directory, fingerprint):
                print('REUSE ' + directory.name, flush=True)
                return directory
            continue
        directory.mkdir(parents=True)
        break
    else:
        raise RuntimeError('Too many incomplete attempts')
    print('RUN ' + directory.name, flush=True)
    start = time.monotonic()
    try:
        compute(directory)
        outputs = {str(p.relative_to(directory)): sha256(p) for p in sorted(directory.rglob('*'))
                   if p.is_file() and '.partial' not in p.name}
        if not outputs:
            raise ValueError('Empty stage')
        save_json(directory / 'manifest.json', {'status': 'complete', 'stage': name,
            'fingerprint': fingerprint, 'provenance': provenance, 'outputs': outputs,
            'elapsed_seconds': round(time.monotonic() - start, 3)})
        valid(directory, fingerprint)
    except Exception as exc:
        save_json(directory / 'failed.json', {'status': 'failed', 'error': repr(exc)})
        raise
    print('DONE ' + directory.name, flush=True)
    return directory


def software():
    return {'python': platform.python_version(), **{n: importlib.metadata.version(n)
        for n in ['anndata', 'numpy', 'pandas', 'scipy']}}


def unit_rng(seed, *labels):
    import numpy as np
    payload = json.dumps([seed, *labels]).encode()
    derived = int.from_bytes(hashlib.sha256(payload).digest()[:8], 'little')
    return np.random.default_rng(derived)


def winners(scores):
    """Unique positive maxima only. -1 denotes no UMIs, -2 a tied maximum."""
    import numpy as np
    scores = np.asarray(scores)
    maximum = scores.max(axis=0)
    ties = np.isclose(scores, maximum[None, :], rtol=0, atol=1e-10).sum(axis=0)
    result = np.argmax(scores, axis=0).astype('int16')
    result[ties > 1] = -2
    result[maximum <= 0] = -1
    return result

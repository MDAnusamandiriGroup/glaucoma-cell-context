"""Acquire checksummed EyeGEx data and indexed regional GWAS extracts."""
from __future__ import annotations

import gzip
import argparse
import hashlib
import json
import os
from pathlib import Path
import urllib.request

import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'data' / 'raw'
PROCESSED = ROOT / 'data' / 'processed'
GWAS_BASE = 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/'
GWAS_HEADER = 'chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\tvariant_id\thm_coordinate_conversion\thm_code\trsid'
EYE = [
    ('eyegex_chr14_nominal.tsv.gz', 'new_nominal.eQTL.chr14.with_thresholds.txt.gz', 'e23894c240705f7fd48b175de7eef9b6', 113356961),
    ('eyegex_permutations.tsv', 'new_perms_full.txt', 'a0cdb73cacf67dcbdd8ac61bcf65e614', 3022816),
    ('eyegex_README.txt', 'README_Retina_eQTL_files_021125.txt', '848bec73de0876d7e1df6a26f08a9c6e', 3758),
]
GWAS = {
    'GCST90011766': {'full_md5': '9d87009f4a75ba475b012ae84b0d3335', 'index_md5': '286a6c863595dd98d05633d03a815c2e', 'n_cases': 16677, 'n_controls': 199580, 'ancestry': 'European'},
    'GCST90011770': {'full_md5': 'd56e70bbaff73ada41e71fb95ea3a3b1', 'index_md5': '793ee33d62975d044a29b687cb6fb250', 'n_cases': 34179, 'n_controls': 349321, 'ancestry': 'Combined European, Asian and African'},
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_new(path, data):
    path = Path(path)
    if path.exists():
        if path.read_bytes() != data:
            raise RuntimeError('Preserve existing incompatible artifact: ' + str(path))
        return
    temporary = path.with_name(path.name + '.partial')
    if temporary.exists():
        raise RuntimeError('Preserve and inspect incomplete artifact: ' + str(temporary))
    with temporary.open('xb') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.rename(temporary, path)


def verified_download(path, url, expected_md5=None, expected_bytes=None):
    path = Path(path)
    if path.exists():
        data = path.read_bytes()
        if expected_md5 and hashlib.md5(data).hexdigest() != expected_md5:
            raise RuntimeError('Existing source checksum mismatch: ' + str(path))
        if expected_bytes is not None and len(data) != expected_bytes:
            raise RuntimeError('Existing source size mismatch: ' + str(path))
        print('REUSE source ' + path.name, flush=True)
        return
    temporary = path.with_name(path.name + '.partial')
    if temporary.exists():
        raise RuntimeError('Preserve and inspect incomplete download: ' + str(temporary))
    md5 = hashlib.md5()
    total = 0
    with urllib.request.urlopen(url, timeout=60) as response, temporary.open('xb') as stream:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            stream.write(chunk)
            md5.update(chunk)
            total += len(chunk)
        stream.flush()
        os.fsync(stream.fileno())
    if expected_md5 and md5.hexdigest() != expected_md5:
        raise RuntimeError('Downloaded source checksum mismatch: ' + str(path))
    if expected_bytes is not None and total != expected_bytes:
        raise RuntimeError('Downloaded source size mismatch: ' + str(path))
    os.rename(temporary, path)
    print('DOWNLOADED source ' + path.name, flush=True)


def acquire():
    RAW.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    plan = json.loads((ROOT / 'PLAN.json').read_text())
    registry_path = ROOT / 'data' / 'source_registry.json'
    if registry_path.exists():
        registry = json.loads(registry_path.read_text())
        for section in ['local_sources', 'derived_inputs']:
            for record in registry[section]:
                path = ROOT / record['file']
                if not path.exists() or sha(path) != record['sha256']:
                    raise RuntimeError('Cached acquisition mismatch: ' + record['file'] + '. For the filtered-data GitHub package, use --verify-filtered; do not delete validated inputs or checkpoints.')
        print('REUSE acquisition registry', flush=True)
        return registry
    registry = {'schema_version': 1, 'retrieval_date_UTC': '2026-10-06',
        'genome_build': 'GRCh38', 'local_sources': [], 'derived_inputs': [], 'gwas': {}}
    for local_name, remote_name, md5, size in EYE:
        url = 'https://zenodo.org/records/18331341/files/' + remote_name + '?download=1'
        path = RAW / local_name
        verified_download(path, url, md5, size)
        registry['local_sources'].append({'file': str(path.relative_to(ROOT)), 'url': url, 'publisher_md5': md5, 'sha256': sha(path), 'bytes': path.stat().st_size})
    nominal = PROCESSED / 'eyegex_target_nominal_raw.tsv'
    permutations = PROCESSED / 'eyegex_target_permutations.tsv'
    ids = set(plan['genes'].values())
    if not nominal.exists():
        selected = []
        with gzip.open(RAW / EYE[0][0], 'rt') as stream:
            header = stream.readline().strip().split('\t')
            for line in stream:
                if line.split('\t', 1)[0] in ids:
                    selected.append(line.rstrip().split('\t'))
        frame = pd.DataFrame(selected, columns=header)
        for column in ['V3', 'V6', 'V10', 'V12', 'V13']:
            frame[column] = pd.to_numeric(frame[column])
        atomic_new(nominal, frame.to_csv(sep='\t', index=False).encode())
    if not permutations.exists():
        frame = pd.read_csv(RAW / EYE[1][0], sep='\t')
        atomic_new(permutations, frame.loc[frame.V1.isin(ids)].to_csv(sep='\t', index=False).encode())
    for study, metadata in GWAS.items():
        base = GWAS_BASE + study + '/harmonised/'
        index = RAW / (study + '.h.tsv.gz.tbi')
        meta = RAW / (study + '.h.tsv.gz-meta.yaml')
        checksums = RAW / (study + '_md5sum.txt')
        verified_download(meta, base + meta.name)
        verified_download(checksums, base + 'md5sum.txt')
        checksum_text = checksums.read_text()
        if metadata['full_md5'] not in checksum_text or metadata['index_md5'] not in checksum_text:
            raise RuntimeError('GWAS release no longer matches frozen checksums: ' + study)
        if 'GRCh38' not in meta.read_text():
            raise RuntimeError('Wrong genomic build: ' + study)
        verified_download(index, base + index.name, metadata['index_md5'])
        for path, url in [(meta, base + meta.name), (checksums, base + 'md5sum.txt'), (index, base + index.name)]:
            registry['local_sources'].append({'file': str(path.relative_to(ROOT)), 'url': url, 'sha256': sha(path), 'bytes': path.stat().st_size})
        url = base + study + '.h.tsv.gz'
        region_path = PROCESSED / ('gwas_' + study + '_NPC2_region.tsv')
        if not region_path.exists():
            with pysam.TabixFile(url, index=str(index)) as reader:
                records = list(reader.fetch(plan['index_chromosome'], plan['index_position'] - plan['primary_half_window_bp'] - 1, plan['index_position'] + plan['primary_half_window_bp']))
            if len(records) < 1000:
                raise RuntimeError('Unexpectedly sparse GWAS regional extract: ' + study)
            atomic_new(region_path, (GWAS_HEADER + '\n' + '\n'.join(records) + '\n').encode())
        print('REGION ' + study + ' rows=' + str(len(pd.read_csv(region_path, sep='\t'))), flush=True)
        registry['gwas'][study] = {**metadata, 'url': url, 'metadata_url': base + meta.name,
            'verification': 'Published full-file checksum retained; local index MD5 verified; full GWAS file not downloaded. HTTPS Tabix regional extraction has local SHA256 below',
            'region_start_1based': plan['index_position'] - plan['primary_half_window_bp'], 'region_end_1based': plan['index_position'] + plan['primary_half_window_bp']}
    for path in [nominal, permutations] + [PROCESSED / ('gwas_' + study + '_NPC2_region.tsv') for study in GWAS]:
        registry['derived_inputs'].append({'file': str(path.relative_to(ROOT)), 'sha256': sha(path), 'bytes': path.stat().st_size})
    atomic_new(registry_path, (json.dumps(registry, indent=2) + '\n').encode())
    return registry


def verify_filtered():
    """Check the public package without downloading omitted large raw sources."""
    registry = json.loads((ROOT / 'data/source_registry.json').read_text())
    for record in registry['derived_inputs']:
        path = ROOT / record['file']
        if not path.is_file() or sha(path) != record['sha256']:
            raise RuntimeError('Filtered input checksum mismatch: ' + record['file'])
    print('VERIFIED packaged filtered inputs; omitted full raw files not re-downloaded', flush=True)
    return registry


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify-filtered', action='store_true',
                        help='Validate the packaged filtered data without requiring omitted full raw files')
    args = parser.parse_args()
    if args.verify_filtered:
        verify_filtered()
    else:
        acquire()

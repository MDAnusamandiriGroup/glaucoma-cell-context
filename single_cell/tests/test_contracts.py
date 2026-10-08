from pathlib import Path
import gzip
import io
import json
import tempfile
import tarfile
import unittest
import numpy as np
import pandas as pd
from scipy import sparse

from scsrc.common import checkpoint, fingerprint, save_json
from scsrc.ingest import read_sparse_counts, attach_metadata
from scsrc.aggregate import group_expression, balanced_expression, robustness


class CheckpointContracts(unittest.TestCase):
    def test_reuse_never_recomputes_or_rewrites(self):
        with tempfile.TemporaryDirectory() as root:
            provenance = {'input': 'fixture', 'parameter': 7}
            calls = []
            def compute(directory):
                calls.append(True)
                save_json(directory / 'answer.json', {'value': 42})
            first = checkpoint(root, 'fixture', provenance, compute)
            before = {p.name: p.stat().st_mtime_ns for p in first.iterdir()}
            second = checkpoint(root, 'fixture', provenance, compute)
            self.assertEqual(first, second)
            self.assertEqual(len(calls), 1)
            self.assertEqual(before, {p.name: p.stat().st_mtime_ns for p in second.iterdir()})

    def test_corrupted_checkpoint_fails_closed(self):
        with tempfile.TemporaryDirectory() as root:
            provenance = {'input': 'fixture'}
            folder = checkpoint(root, 'fixture', provenance, lambda d: save_json(d / 'answer.json', {'v': 1}))
            (folder / 'answer.json').write_text('changed')
            with self.assertRaisesRegex(ValueError, 'Corrupt/changed'):
                checkpoint(root, 'fixture', provenance, lambda d: self.fail('must not compute'))
            self.assertEqual((folder / 'answer.json').read_text(), 'changed')

    def test_partial_stage_preserved_and_new_attempt_used(self):
        with tempfile.TemporaryDirectory() as root:
            provenance = {'input': 'fixture'}
            partial = Path(root) / ('fixture_' + fingerprint(provenance))
            partial.mkdir()
            (partial / 'diagnostic.txt').write_text('preserve me')
            complete = checkpoint(root, 'fixture', provenance, lambda d: save_json(d / 'answer.json', {'v': 1}))
            self.assertNotEqual(partial, complete)
            self.assertEqual((partial / 'diagnostic.txt').read_text(), 'preserve me')

    def test_publish_cannot_overwrite(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / 'output.json'
            save_json(path, {'v': 1})
            before = path.read_bytes()
            with self.assertRaises(FileExistsError):
                save_json(path, {'v': 2})
            self.assertEqual(path.read_bytes(), before)


class ScientificContracts(unittest.TestCase):
    def test_mislabeled_gzip_tar_count_matrix(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / 'matrix.csv.gz'
            payload = b'cellA,cellB\nGENE1,0,2\nGENE2,3,0\n'
            with tarfile.open(path, 'w:gz') as archive:
                member = tarfile.TarInfo('lukowski_embo2019_raw_count_matrix.csv')
                member.size = len(payload)
                archive.addfile(member, io.BytesIO(payload))
            counts, genes, cells = read_sparse_counts(path)
            self.assertEqual(genes, ['GENE1', 'GENE2'])
            self.assertEqual(cells, ['cellA', 'cellB'])
            np.testing.assert_array_equal(counts.toarray(), [[0, 3], [2, 0]])

    def test_duplicate_gene_symbols_rejected(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / 'matrix.csv.gz'
            with gzip.open(path, 'wt') as handle:
                handle.write('cellA,cellB\nGENE,0,2\nGENE,3,0\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate gene'):
                read_sparse_counts(path)

    def metadata_fixture(self):
        cells = ['a', 'b', 'c', 'd', 'e']
        meta = pd.DataFrame({'cell.bc': cells, 'sample.group': ['R1','R2','R2','R3','R3'],
                             'batch': ['L1','L2','L3','L4','L5'], 'Age': [80,42,42,53,53],
                             'Gender': ['F','M','M','F','F'], 'Collection': [11.5,6.2,6.2,14.5,14.5],
                             'nCount_RNA': [1,2,3,4,5], 'nFeature_RNA': [1]*5, 'percent.mito': [0.01]*5})
        labels = pd.DataFrame({'cell.bc': cells, 'cell.id.orig': ['Rod PR C0']*5,
                               'cell.id.cca': ['Rod PR CCA0']*5, 'cluster.id.orig': [0]*5, 'cluster.id.cca': [0]*5})
        return cells, meta, labels

    def test_metadata_alignment_uses_exact_barcodes(self):
        cells, meta, labels = self.metadata_fixture()
        obs = attach_metadata(list(reversed(cells)), meta, labels)
        self.assertEqual(obs.index.tolist(), list(reversed(cells)))
        self.assertEqual(obs['source_n_counts'].tolist(), [5,4,3,2,1])
        self.assertEqual(obs['published_cell_type'].unique().tolist(), ['Rod'])
        with self.assertRaisesRegex(ValueError, 'match exactly'):
            attach_metadata(['wrong'] + cells[1:], meta, labels)
        with self.assertRaisesRegex(ValueError, 'Duplicate cell'):
            attach_metadata(cells, pd.concat([meta,meta.iloc[:1]]), labels)

    def test_missing_groups_and_genes_are_not_zero_expression(self):
        counts = sparse.csr_matrix([[1,9],[0,10],[2,8],[0,10],[0,10],[0,10]], dtype=np.int32)
        norm = counts.astype(float).copy()
        norm.data = np.log1p(norm.data * 1000)
        obs = pd.DataFrame({'donor': ['R1','R1','R2','R2','R3','R3'], 'published_cell_type': ['A']*6,
                            'qc_strict': [True,False,True,True,True,True]})
        result = group_expression(counts, norm, obs, ['G','background'], ['G','absent'], ['A','B'], 2)
        missing = result.loc[result['cell_type'].eq('B') & result['gene'].eq('G')]
        self.assertTrue(missing['pseudobulk_CPM'].isna().all())
        self.assertTrue(missing['n_cells'].eq(0).all())
        present = result.loc[result['qc'].eq('primary') & result['gene'].eq('G') & result['cell_type'].eq('A')]
        self.assertEqual(present.set_index('donor').loc['R1','pseudobulk_CPM'], 50000.0)
        self.assertEqual(present.set_index('donor').loc['R3','pseudobulk_CPM'], 0.0)
        self.assertTrue(result.loc[result['gene'].eq('absent'), 'pseudobulk_CPM'].isna().all())
        strict_r1 = result.loc[result['qc'].eq('strict') & result['donor'].eq('R1') & result['cell_type'].eq('A') & result['gene'].eq('G')].iloc[0]
        self.assertFalse(strict_r1['eligible'])
        balanced = balanced_expression(result, 2)
        self.assertEqual(balanced.loc[balanced['qc'].eq('strict') & balanced['cell_type'].eq('A') & balanced['gene'].eq('G'), 'eligible_donors'].iloc[0], 2)

    def test_equal_donor_weighting_and_fixed_LODO_comparison(self):
        rows = []
        for qc in ['primary','strict']:
            for donor, scores in zip(['R1','R2','R3'], [[8,2],[1,4],[8,2]]):
                for cell_type, score, n in [('A',scores[0],20),('B',scores[1],2000),('Rare',99,19)]:
                    rows.append({'qc':qc,'gene':'G','donor':donor,'cell_type':cell_type,'n_cells':n,
                                 'gene_in_matrix':True,'log1p_pseudobulk_CPM':score,
                                 'mean_log1p_CP10k':score,'pct_cells_detected':10})
        frame = pd.DataFrame(rows)
        balanced = balanced_expression(frame,20)
        score = balanced.loc[balanced['qc'].eq('primary') & balanced['cell_type'].eq('A'),'donor_mean_log1p_pseudobulk_CPM'].iloc[0]
        self.assertAlmostEqual(score,17/3)
        summary, folds, _, shared, _, sensitivity = robustness(frame,['G'],['A','B','Rare'],20,10)
        self.assertEqual(shared['primary'],['A','B'])
        self.assertEqual(summary.iloc[0]['primary_top_type'],'A')
        self.assertTrue(summary.iloc[0]['LODO_all_three_agree'])
        self.assertTrue(folds['comparison_types'].eq('A;B').all())
        self.assertEqual(sensitivity,['A','B','Rare'])
        self.assertEqual(summary.iloc[0]['sensitivity_min10_top'],'Rare')


if __name__ == '__main__':
    unittest.main()

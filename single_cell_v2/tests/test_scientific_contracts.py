from pathlib import Path
import tempfile
import unittest
import json
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from v2src.cache import checkpoint, valid, unit_rng, winners
from v2src.sampling import thin_group
from v2src.regulatory import tokens, canonical_types


class ScientificContracts(unittest.TestCase):
    def test_umi_thinning_conserves_selected_and_background_counts(self):
        for seed in range(40):
            drawn, total = thin_group(np.array([1, 3, 20]), 100, .25, np.random.default_rng(seed))
            self.assertLessEqual(int(drawn.sum()), total)
            self.assertTrue(np.all(drawn <= [1, 3, 20]))
        drawn, total = thin_group(np.array([1, 3, 20]), 100, 1, np.random.default_rng(1))
        np.testing.assert_array_equal(drawn, [1, 3, 20]); self.assertEqual(total, 100)

    def test_invalid_count_partition_rejected(self):
        with self.assertRaises(ValueError):
            thin_group(np.array([3, 4]), 6, .5, np.random.default_rng(1))

    def test_ties_and_no_signal_are_not_unique_winners(self):
        result = winners(np.array([[1., 0., 3.], [1., 0., 2.]]))
        np.testing.assert_array_equal(result, [-2, -1, 0])

    def test_unit_seed_reproduces_completed_draw(self):
        a=unit_rng(42,'scenario',3,'donor','type').integers(0,100000,10)
        b=unit_rng(42,'scenario',3,'donor','type').integers(0,100000,10)
        c=unit_rng(42,'scenario',4,'donor','type').integers(0,100000,10)
        np.testing.assert_array_equal(a,b);self.assertFalse(np.array_equal(a,c))

    def test_gene_links_use_complete_names(self):
        self.assertIn('PSMC3',tokens('PSMC3, NPC2'))
        self.assertNotIn('PSMC3',tokens('PSMC3IP, NPC2'))
        self.assertEqual(tokens(float('nan')),set())

    def test_subtypes_collapse_without_inventing_labels(self):
        self.assertEqual(canonical_types('OFF.cone.bipolar, ON.cone.bipolar, Microglia'),{'Bipolar','Microglia'})
        with self.assertRaises(ValueError):canonical_types('Unknown subtype')

    def test_checkpoint_reuse_and_changed_output_rejection(self):
        with tempfile.TemporaryDirectory() as temp:
            calls=[]
            def compute(directory):
                calls.append(1);(directory/'estimate.txt').write_text('validated estimate\n')
            first=checkpoint(Path(temp),'fixture',{'source':'fixed'},compute)
            second=checkpoint(Path(temp),'fixture',{'source':'fixed'},compute)
            self.assertEqual(first,second);self.assertEqual(len(calls),1)
            fingerprint=json.loads((first/'manifest.json').read_text())['fingerprint']
            (first/'estimate.txt').write_text('changed estimate\n')
            with self.assertRaises(ValueError):valid(first,fingerprint)


if __name__=='__main__':unittest.main()

import datetime
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from reproduce_all import complete_manifest, unused_attempt, validate_manifest
from src.workflow import bh_adjust, evidence_matrix, load_expression, matched_reference, write_new


class ScientificContracts(unittest.TestCase):
    def test_method_agreement_requires_same_qtl_and_tissue(self):
        base = {"locus": "L1", "rsid": "rs1", "gene_id": "ENSG1", "gene": "GENE1", "qtl_type": "eQTL", "qc_pass": True, "score": 0.8, "baseline_pass": True, "strict_pass": True}
        records = pd.DataFrame([
            dict(base, method="eCAVIAR", tissue="Retina"),
            dict(base, method="enloc", tissue="Brain"),
            dict(base, method="eCAVIAR", tissue="Brain"),
            dict(base, method="enloc", tissue="Retina", qtl_type="sQTL"),
        ])
        result = evidence_matrix(records)
        self.assertEqual(int(result["baseline_same_context_agreement"].sum()), 1)
        self.assertFalse(result.loc[result["tissue"].eq("Retina"), "baseline_same_context_agreement"].any())
        self.assertTrue(result.loc[result["tissue"].eq("Retina") & result["qtl_type"].eq("eQTL"), "enloc_score"].isna().all())

    def test_excel_date_labels_are_audited_not_guessed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "atlas.xlsx"
            raw = pd.DataFrame({"gene": ["GENE1", datetime.datetime(2001, 3, 1), datetime.datetime(2001, 3, 1)], "C0": [1, 2, 3], "C1": [4, 5, 6]})
            raw.to_excel(path, sheet_name="Atlas", index=False)
            expression, audit = load_expression(path, "Atlas", {"RGC": ["C0"], "Rod": ["C1"]})
            self.assertEqual(expression.index.tolist(), ["GENE1"])
            self.assertEqual(audit["excel_row"].tolist(), [3, 4])

    def test_real_duplicate_text_gene_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "atlas.xlsx"
            pd.DataFrame({"gene": ["GENE1", "GENE1"], "C0": [1, 2]}).to_excel(path, sheet_name="Atlas", index=False)
            with self.assertRaises(ValueError):
                load_expression(path, "Atlas", {"RGC": ["C0"]})

    def test_bh_adjustment(self):
        np.testing.assert_allclose(bh_adjust(np.array([0.01, 0.04, 0.03])), [0.03, 0.04, 0.04])

    def test_matched_reference_is_deterministic_and_nonzero(self):
        atlas = pd.DataFrame({"RGC": np.arange(1, 41), "Rod": np.arange(40, 0, -1)}, index=[f"G{i}" for i in range(40)])
        config = {"expression_match_bins": 2, "n_permutations": 50, "seed": 99}
        first, audit = matched_reference(atlas, ["G2", "G30"], config)
        second, _ = matched_reference(atlas, ["G2", "G30"], config)
        pd.testing.assert_frame_equal(first, second)
        self.assertEqual(int(audit["candidate_genes"].sum()), 2)
        self.assertEqual(int(audit["available_control_genes"].sum()), 38)
        self.assertTrue(first["one_sided_empirical_p"].ge(1 / 51).all())


class CheckpointContracts(unittest.TestCase):
    def test_valid_cache_reuses_and_modified_output_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "analysis_example"
            root.mkdir()
            write_new(root / "result.csv", "gene,value\nGENE1,2\n")
            complete_manifest(root, "analysis_manifest.json", {"fingerprint": "abc"})
            selected, reused = unused_attempt(root, "analysis_manifest.json", "abc")
            self.assertTrue(reused)
            self.assertEqual(selected, root)
            with (root / "result.csv").open("a") as handle:
                handle.write("GENE2,3\n")
            with self.assertRaises(ValueError):
                validate_manifest(root, "analysis_manifest.json", "abc")

    def test_interrupted_stage_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "analysis_example"
            root.mkdir()
            write_new(root / "partial.csv", "preserve me\n")
            selected, reused = unused_attempt(root, "analysis_manifest.json", "abc")
            self.assertFalse(reused)
            self.assertEqual(selected.name, "analysis_example_attempt2")
            self.assertEqual((root / "partial.csv").read_text(), "preserve me\n")

    def test_write_never_overwrites(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            write_new(path, json.dumps({"value": 1}))
            with self.assertRaises(FileExistsError):
                write_new(path, json.dumps({"value": 2}))
            self.assertEqual(json.loads(path.read_text())["value"], 1)

    def test_explicit_figure_recovery_preserves_corrupt_attempt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "figures_example"
            root.mkdir()
            write_new(root / "figure.png", b"original placeholder")
            complete_manifest(root, "figure_manifest.json", {"fingerprint": "abc"})
            (root / "figure.png").write_bytes(b"")
            with self.assertRaises(ValueError):
                unused_attempt(root, "figure_manifest.json", "abc")
            selected, reused = unused_attempt(root, "figure_manifest.json", "abc", rebuild=True)
            self.assertFalse(reused)
            self.assertEqual(selected.name, "figures_example_attempt2")
            self.assertEqual((root / "figure.png").read_bytes(), b"")


if __name__ == "__main__":
    unittest.main()

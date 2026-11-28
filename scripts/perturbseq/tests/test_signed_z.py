"""Offline tests for signed_z.py."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import signed_z as sz


class SignedZTests(unittest.TestCase):
    def test_restrict_keeps_z_order_and_drops_other_genes(self):
        z = pd.DataFrame({"LV1": [1.0, 2.0, 3.0, 4.0]}, index=["D", "B", "C", "A"])
        out = sz.restrict(z, ["A", "D", "B"])
        self.assertEqual(out.index.tolist(), ["D", "B", "A"])

    def test_missing_universe_gene_fails(self):
        z = pd.DataFrame({"LV1": [1.0]}, index=["A"])
        with self.assertRaises(ValueError):
            sz.restrict(z, ["A", "B"])

    def test_check_genes_requires_the_exact_gene_set(self):
        z = pd.DataFrame({"LV1": [1.0, 2.0]}, index=["B", "A"])
        sz.check_genes(z, ["A", "B"])
        for genes in (["A"], ["A", "B", "C"]):
            with self.assertRaises(ValueError):
                sz.check_genes(z, genes)

    def test_main_writes_both_tails_and_sample_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            z = pd.DataFrame({"LV1": [0.1, -0.2, 1 / 3], "LV2": [0.0, 5.0, -1e-17]},
                             index=pd.Index(["A", "B", "C"], name="gene"))
            z.to_csv(tmp / "Z.csv")
            pd.DataFrame({"gene": ["C", "A"], "KD1": [0.1, 0.2]}).to_csv(
                tmp / "logfc.tsv.gz", sep="\t", index=False)
            pd.DataFrame({"mean": [1.0, 2.0, 3.0]}, index=["C", "B", "A"]).to_csv(tmp / "row_stats.csv")
            pd.DataFrame([[1, 2, 3]], index=["LV1"], columns=["s1", "s2", "s3"]).to_csv(tmp / "B.csv")
            sz.main(["--z", str(tmp / "Z.csv"), "--genes", str(tmp / "row_stats.csv"), "--universe", str(tmp / "logfc.tsv.gz"),
                     "--b", str(tmp / "B.csv"), "--out-pos", str(tmp / "out/pos.csv"),
                     "--out-neg", str(tmp / "out/neg.csv"), "--out-manifest", str(tmp / "out/manifest.json")])
            pos = pd.read_csv(tmp / "out/pos.csv", index_col=0, float_precision="round_trip")
            neg = pd.read_csv(tmp / "out/neg.csv", index_col=0, float_precision="round_trip")
            self.assertEqual(pos.index.tolist(), ["A", "C"])
            np.testing.assert_array_equal(pos.to_numpy(), z.loc[["A", "C"]].to_numpy())
            np.testing.assert_array_equal(neg.to_numpy(), -pos.to_numpy())
            manifest = json.loads((tmp / "out/manifest.json").read_text())
            self.assertEqual(manifest["samples"], 3)
            self.assertEqual(manifest["universe_genes"], 2)


if __name__ == "__main__":
    unittest.main()

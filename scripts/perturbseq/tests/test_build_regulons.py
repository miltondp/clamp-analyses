"""Offline tests for build_regulons.py on tiny synthetic Perturb-seq matrices."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import build_regulons as br


def read_gmt(path):
    """Same parsing as scripts/archs4/coverage/run_ora.R::read_gmt."""
    terms = {}
    for line in Path(path).read_text().splitlines():
        record = line.split("\t")
        if len(record) >= 3:
            terms[record[0]] = record[2:]
    return terms


class BhTests(unittest.TestCase):
    def test_matches_hand_computed_values(self):
        p = np.array([0.01, 0.04, 0.03, 0.20])
        # sorted 0.01, 0.03, 0.04, 0.20 -> raw 0.04, 0.06, 0.0533, 0.20 -> monotone 0.04, 0.0533, 0.0533, 0.20
        expected = [0.04, 0.16 / 3, 0.16 / 3, 0.20]
        np.testing.assert_allclose(br.bh_columns(p[:, None])[:, 0], expected)

    def test_nan_is_excluded_from_the_denominator(self):
        p = np.array([[0.01, 0.01], [0.02, np.nan], [0.03, 0.02]])
        q = br.bh_columns(p)
        np.testing.assert_allclose(q[:, 0], [0.03, 0.03, 0.03])
        np.testing.assert_allclose(q[[0, 2], 1], [0.02, 0.02])
        self.assertTrue(np.isnan(q[1, 1]))

    def test_columns_are_independent(self):
        rng = np.random.default_rng(0)
        p = rng.uniform(size=(50, 3))
        q = br.bh_columns(p)
        for j in range(3):
            np.testing.assert_allclose(q[:, j], br.bh_columns(p[:, [j]])[:, 0])


class HelperTests(unittest.TestCase):
    def test_drop_self_clears_only_the_own_row(self):
        de = pd.DataFrame(True, index=["A", "B", "C"], columns=["A", "C", "Z"])
        out = br.drop_self(de)
        self.assertFalse(out.loc["A", "A"])
        self.assertFalse(out.loc["C", "C"])
        self.assertEqual(int(out.values.sum()), 9 - 2)
        self.assertTrue(de.values.all(), "input must not be modified")

    def test_ambiguous_and_unmapped_symbols_are_dropped(self):
        ids = pd.Series({"E1": "A", "E2": "B", "E3": "B"})
        symbols, stats = br.target_symbols(pd.Index(["E1", "E2", "E3", "E4"]), ids)
        self.assertEqual(symbols.tolist()[0], "A")
        self.assertTrue(symbols.iloc[1:].isna().all())
        self.assertEqual(stats, {"unmapped": 1, "ambiguous_rows": 2})

    def test_collapse_assigns_near_duplicates_to_the_largest(self):
        genes = [f"g{i}" for i in range(20)]
        sets = pd.DataFrame(False, index=genes, columns=["big", "dup", "other"])
        sets.loc[genes[:10], "big"] = True       # 10 genes
        sets.loc[genes[:8], "dup"] = True        # Jaccard with big = 0.8
        sets.loc[genes[12:18], "other"] = True   # disjoint
        clusters = br.collapse(sets, 0.5).set_index("perturbation")
        self.assertEqual(clusters.loc["dup", "representative"], "big")
        self.assertAlmostEqual(clusters.loc["dup", "jaccard_to_rep"], 0.8)
        self.assertEqual(clusters.loc["other", "representative"], "other")
        self.assertEqual(clusters.loc["big", "representative"], "big")

    def test_collapse_threshold_is_strict(self):
        sets = pd.DataFrame({"a": [True, True, False], "b": [False, True, True]}, index=["x", "y", "z"])
        # Jaccard = 1/3
        reps = br.collapse(sets, 1 / 3).representative.tolist()
        self.assertEqual(sorted(reps), ["a", "b"])


class BuildTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        # 30 measured targets T0..T29 (Ensembl E0..E29); E29 is unmapped.
        self.n = 30
        symbols = [f"T{i}" for i in range(self.n)]
        ensembl = [f"E{i}" for i in range(self.n)]
        id_map = pd.DataFrame({"e": ensembl[:-1], "s": symbols[:-1]})
        id_map.to_csv(self.root / "ids.tsv", sep="\t", header=False, index=False)

        perts = ["T0", "BIG", "DUP", "SMALL", "NONE"]
        p = np.full((self.n, len(perts)), 0.9)
        p[:, 0][:12] = 1e-8   # T0: 12 hits incl. itself (T0) -> 11 others
        p[:, 1][:20] = 1e-8   # BIG: T0..T19
        p[:, 2][:18] = 1e-8   # DUP: T0..T17, Jaccard with BIG = 18/20 within the universe
        p[:, 3][:3] = 1e-8    # SMALL: 3 hits
        p[29, 1] = 1e-8       # unmapped row is never counted
        lfc = np.where(p < 0.05, -1.0, 0.0)
        self.write("p.tsv", p, ensembl, perts)
        self.write("lfc.tsv", lfc, ensembl, perts)

        # Universe: Z genes T0..T27 plus genes the screen did not measure.
        z_genes = symbols[:28] + ["ZONLY1", "ZONLY2"]
        for name in ("z1.csv", "z2.csv"):
            pd.DataFrame({"LV1": 0.0, "LV2": 1.0}, index=z_genes).to_csv(self.root / name)

    def write(self, name, values, rows, cols):
        frame = pd.DataFrame(values, index=pd.Index(rows, name="GENE"), columns=cols)
        frame.to_csv(self.root / name, sep="\t")

    def run_build(self, *extra):
        out = self.root / "out"
        br.main([
            "--logfc", str(self.root / "lfc.tsv"), "--pvalues", str(self.root / "p.tsv"),
            "--id-map", str(self.root / "ids.tsv"),
            "--z", str(self.root / "z1.csv"), str(self.root / "z2.csv"),
            "--min-size", "10", "--max-size", "500", "--jaccard", "0.5", "--min-de", "15",
            "--out-dir", str(out), *extra,
        ])
        return out

    def test_end_to_end(self):
        out = self.run_build()
        report = json.loads((out / "report.json").read_text())
        self.assertEqual(report["universe"]["size"], 28)
        self.assertEqual(report["screen"]["targets_unmapped"], 1)
        self.assertEqual(report["screen"]["perturbations_self_de"], 1)

        sizes = pd.read_csv(out / "sizes.csv").set_index("perturbation")
        self.assertEqual(sizes.loc["T0", "n_de_universe"], 11)      # own row dropped
        self.assertEqual(sizes.loc["BIG", "n_de_universe"], 20)     # unmapped row ignored
        self.assertEqual(sizes.loc["T0", "self_logfc"], -1.0)
        self.assertFalse(bool(sizes.loc["SMALL", "regulon"]))

        gmt = read_gmt(out / "regulons.gmt")
        self.assertEqual(set(gmt), {"BIG"})                         # T0 and DUP collapse onto BIG
        self.assertEqual(len(gmt["BIG"]), 20)
        clusters = pd.read_csv(out / "clusters.csv").set_index("perturbation")
        self.assertEqual(clusters.loc["DUP", "representative"], "BIG")
        self.assertEqual(clusters.loc["T0", "representative"], "BIG")
        self.assertEqual(report["regulons"]["largest_clusters"], {"BIG": 3})

        logfc = pd.read_csv(out / "logfc.tsv.gz", sep="\t", index_col=0)
        self.assertEqual(logfc.shape, (28, 2))                      # BIG and DUP have >= 15
        self.assertEqual(sorted(logfc.columns), ["BIG", "DUP"])

    def test_size_filter_applies_after_universe(self):
        out = self.run_build("--max-size", "19")
        sizes = pd.read_csv(out / "sizes.csv").set_index("perturbation")
        self.assertFalse(bool(sizes.loc["BIG", "regulon"]))
        self.assertEqual(set(read_gmt(out / "regulons.gmt")), {"DUP"})

    def test_mismatched_z_gene_sets_fail(self):
        pd.DataFrame({"LV1": [0.0]}, index=["T0"]).to_csv(self.root / "z2.csv")
        with self.assertRaises(ValueError):
            self.run_build()


if __name__ == "__main__":
    unittest.main()

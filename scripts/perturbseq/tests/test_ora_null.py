"""Offline tests for ora_null.py; the last class also checks it against run_ora.R when R is available."""
from math import comb
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import ora_null as on

REPO = Path(__file__).resolve().parents[3]
RUN_ORA = REPO / "scripts" / "archs4" / "coverage" / "run_ora.R"

GENES = list("ABCDEFGHIJ")
TERMS = {"T1": list("ABCD"), "T2": list("EFG"), "T3": list("HI"), "TINY": ["A"]}


def toy_ora(top_pct=0.3):
    return on.Ora(GENES, TERMS, top_pct=top_pct, min_size=2, max_size=100, fdr=0.05)


def loadings(*tops):
    """genes x LVs matrix whose LV j has positive loadings exactly on the genes in tops[j]."""
    z = np.full((len(GENES), len(tops)), -1.0)
    for j, top in enumerate(tops):
        for rank, gene in enumerate(top):
            z[GENES.index(gene), j] = 10.0 - rank
    return z


class OraTests(unittest.TestCase):
    def test_universe_size_and_query(self):
        ora = toy_ora()
        self.assertEqual(ora.terms.tolist(), ["T1", "T2", "T3"])  # TINY is below min_size
        self.assertEqual(ora.N, 9)  # J is in no eligible regulon
        self.assertEqual(ora.query_target, 3)

    def test_hand_computed_p_values(self):
        ora = toy_ora()
        result = ora.test(loadings("ABJ", "AEH"))
        # LV1: query {A, B, J}, n = 2 annotated; only T1 has a hit (k = 2 of M = 4).
        np.testing.assert_allclose(result["p"][0, 0], comb(4, 2) / comb(9, 2))
        self.assertTrue(np.isnan(result["p"][0, 1:]).all(), "k = 0 regulons are not tested")
        np.testing.assert_allclose(result["p_adjust"][0, 0], comb(4, 2) / comb(9, 2))
        # LV2: query {A, E, H}, n = 3; each regulon has k = 1, P(X >= 1) = 1 - C(N - M, 3) / C(N, 3).
        p = [1 - comb(9 - m, 3) / comb(9, 3) for m in (4, 3, 2)]
        np.testing.assert_allclose(result["p"][1], p)
        np.testing.assert_allclose(result["p_adjust"][1], [max(p)] * 3)  # BH over the 3 tested
        self.assertEqual(result["tested_lvs"], 2)

    def test_short_positive_tail_uses_every_positive_gene(self):
        ora = toy_ora(top_pct=0.5)  # query target 5
        result = ora.test(loadings("AB"))
        np.testing.assert_allclose(result["p"][0, 0], comb(4, 2) / comb(9, 2))

    def test_nonnegative_z_has_an_empty_negative_tail(self):
        ora = toy_ora()
        z = np.abs(loadings("ABC", "EFG"))
        self.assertEqual(ora.test(-z)["tested_lvs"], 0)
        self.assertFalse(ora.test(-z)["hits"].any())

    def test_permutations_are_seeded(self):
        rng = np.random.default_rng(0)
        genes = [f"G{i}" for i in range(200)]
        terms = {f"R{t}": list(rng.choice(genes, 20, replace=False)) for t in range(10)}
        ora = on.Ora(genes, terms, top_pct=0.05, min_size=5, max_size=100, fdr=0.05)
        z = rng.normal(size=(200, 4))
        first = on.permutation_rows(ora, z, 5, seed=1)
        self.assertEqual(first, on.permutation_rows(ora, z, 5, seed=1))
        self.assertEqual([row["perm"] for row in first], [1, 2, 3, 4, 5])

    def test_summary_statistics(self):
        null = pd.DataFrame({"perm": [0, 1, 2, 3], **{s: [50.0, 10.0, 20.0, 60.0] for s in on.STATISTICS}})
        summary = on.summarize(null).set_index("statistic").loc["lv_pct_any"]
        self.assertEqual(summary.observed, 50.0)
        self.assertAlmostEqual(summary.null_mean, 30.0)
        self.assertAlmostEqual(summary.observed_minus_null, 20.0)
        self.assertAlmostEqual(summary.perm_p, (1 + 1) / (1 + 3))


def r_has_clusterprofiler():
    if not shutil.which("Rscript"):
        return False
    probe = subprocess.run(["Rscript", "-e", "library(clusterProfiler)"], capture_output=True)
    return probe.returncode == 0


@unittest.skipUnless(r_has_clusterprofiler(), "Rscript with clusterProfiler is not available")
class MatchesRunOraTests(unittest.TestCase):
    """End to end: signed tails through run_ora.R, then ora_null.py must reproduce its hits."""

    def test_reproduces_run_ora(self):
        rng = np.random.default_rng(3)
        genes = [f"G{i:03d}" for i in range(400)]
        terms = {f"R{t}": sorted(rng.choice(genes, int(rng.integers(10, 60)), replace=False)) for t in range(30)}
        z = rng.normal(size=(400, 6))
        # Plant signal: LV0's positive tail and LV1's negative tail hit regulons; LV5 is nonnegative.
        z[[genes.index(g) for g in terms["R0"]], 0] += 6
        z[[genes.index(g) for g in terms["R1"]], 1] -= 6
        z[:, 5] = np.abs(z[:, 5])
        frame = pd.DataFrame(z, index=pd.Index(genes, name="gene"), columns=[f"LV{i + 1}" for i in range(6)])
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            gmt = tmp / "regulons.gmt"
            gmt.write_text("".join("\t".join([t, "x", *g]) + "\n" for t, g in terms.items()))
            paths = {"pos": tmp / "pos.csv", "neg": tmp / "neg.csv"}
            frame.to_csv(paths["pos"])
            (-frame).to_csv(paths["neg"])
            for tail, path in paths.items():
                subprocess.run([
                    "Rscript", str(RUN_ORA), "--z", str(path), "--out-dir", str(tmp / f"ora_{tail}"),
                    "--dataset", "toy", "--model", "toy", "--n-samples", "10", "--database", "toy",
                    "--database-type", "gmt", "--database-path", str(gmt), "--top-pct", "0.05",
                    "--min-size", "10", "--max-size", "1000", "--pvalue-cutoff", "0.05",
                    "--qvalue-cutoff", "0.2", "--fdr", "0.05",
                ], check=True, capture_output=True)
            args = ["--z", str(paths["pos"]), "--gmt", str(gmt),
                    "--ora-pos", str(tmp / "ora_pos"), "--ora-neg", str(tmp / "ora_neg"),
                    "--strategy", "toy", "--model", "toy", "--top-pct", "0.05", "--min-size", "10",
                    "--max-size", "1000", "--fdr", "0.05", "--n-perm", "5", "--seed", "1"]
            on.main(args + ["--out-dir", str(tmp / "a")])
            on.main(args + ["--out-dir", str(tmp / "b")])

            hits = pd.read_csv(tmp / "a" / "lv_hits.csv").set_index(["LV", "tail"])
            self.assertEqual(hits.loc[("LV1", "pos"), "best_regulon"], "R0")
            self.assertEqual(hits.loc[("LV2", "neg"), "best_regulon"], "R1")
            self.assertEqual(hits.loc[("LV6", "neg"), "n_regulons"], 0)
            pd.testing.assert_frame_equal(pd.read_csv(tmp / "a" / "null.csv.gz"),
                                          pd.read_csv(tmp / "b" / "null.csv.gz"))
            summary = pd.read_csv(tmp / "a" / "summary.csv")
            self.assertEqual(len(summary), len(on.STATISTICS))
            self.assertTrue((summary.model == "toy").all())


if __name__ == "__main__":
    unittest.main()

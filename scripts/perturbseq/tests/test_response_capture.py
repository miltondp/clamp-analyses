"""Offline tests for response_capture.py, on small synthetic matrices."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import response_capture as rc


def corr(a, b):
    return np.corrcoef(a, b)[0, 1]


class CorrelationTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.x = rng.normal(size=(30, 4)) + 5  # offset: r must not depend on it
        self.z = rng.normal(size=(30, 3)) - 2
        self.rows = np.array([3, -1, 0, 29])  # perturbation 1's own gene is not in the universe

    def test_matches_corrcoef_without_own_gene(self):
        r = rc.Response(self.x, self.rows).r(self.z)
        for p, row in enumerate(self.rows):
            keep = np.arange(30) != row
            for lv in range(3):
                self.assertAlmostEqual(r[p, lv], corr(self.x[keep, p], self.z[keep, lv]), places=12)

    def test_own_gene_is_excluded(self):
        x = self.x.copy()
        x[3, 0] = -50.0  # a strong knockdown of perturbation 0's own gene
        z = self.z.copy()
        z[3, 0] = -50.0
        masked = rc.Response(x, self.rows).r(z)[0, 0]
        unmasked = rc.Response(x, np.full(4, -1)).r(z)[0, 0]
        self.assertGreater(unmasked, 0.9)
        self.assertAlmostEqual(masked, corr(np.delete(x[:, 0], 3), np.delete(z[:, 0], 3)), places=12)
        self.assertLess(abs(masked), 0.9)

    def test_zero_variance_lv_has_zero_r(self):
        z = np.column_stack([self.z[:, 0], np.full(30, 0.1), np.zeros(30)])
        r = rc.Response(self.x, self.rows).r(z)
        self.assertTrue(np.isfinite(r).all())
        np.testing.assert_array_equal(r[:, 1:], 0.0)

    def test_non_finite_input_fails(self):
        z = self.z.copy()
        z[0, 0] = np.nan
        with self.assertRaises(ValueError):
            rc.Response(self.x, self.rows).r(z)

    def test_best_uses_absolute_r(self):
        max_abs, idx, signed = rc.best(np.array([[0.2, -0.9, 0.5], [0.3, 0.1, -0.2]]))
        np.testing.assert_allclose(max_abs, [0.9, 0.3])
        np.testing.assert_array_equal(idx, [1, 0])
        np.testing.assert_allclose(signed, [-0.9, 0.3])


class NullTests(unittest.TestCase):
    def test_planted_lv_is_captured_with_either_sign(self):
        rng = np.random.default_rng(1)
        x = rng.normal(size=(200, 5))
        z = np.column_stack([rng.normal(size=(200, 4)), -x[:, 2]])
        response = rc.Response(x, np.full(5, -1))
        null = rc.null_max_abs_r(response, z, n_perm=20, seed=7)
        self.assertAlmostEqual(null[0, 2], 1.0, places=12)
        _, p = rc.run_statistics(null, alpha=0.05)
        self.assertAlmostEqual(p[2], 1 / 21)

    def test_leave_one_out_p(self):
        values = np.array([[3.0], [1.0], [2.0], [2.0]])
        np.testing.assert_allclose(rc.leave_one_out_p(values)[:, 0], [1 / 4, 4 / 4, 3 / 4, 3 / 4])

    def test_noise_capture_rate_is_about_alpha(self):
        rng = np.random.default_rng(2)
        x = rng.normal(size=(100, 400))
        z = rng.normal(size=(100, 6))
        null = rc.null_max_abs_r(rc.Response(x, np.full(400, -1)), z, n_perm=99, seed=3)
        stats, _ = rc.run_statistics(null, alpha=0.05)
        self.assertAlmostEqual(stats.pct_captured.mean(), 5.0, delta=0.5)
        self.assertLess(stats.pct_captured.iloc[0], 15.0)

    def test_seed_is_reproducible(self):
        rng = np.random.default_rng(4)
        x, z = rng.normal(size=(50, 8)), rng.normal(size=(50, 3))
        response = rc.Response(x, np.full(8, -1))
        np.testing.assert_array_equal(rc.null_max_abs_r(response, z, 5, 11), rc.null_max_abs_r(response, z, 5, 11))
        self.assertFalse(np.array_equal(rc.null_max_abs_r(response, z, 5, 11), rc.null_max_abs_r(response, z, 5, 12)))


class CliTests(unittest.TestCase):
    def test_end_to_end(self):
        rng = np.random.default_rng(5)
        genes = [f"G{i}" for i in range(60)]
        perts = ["G0", "G1", "OTHER"]  # G0 and G1 are measured; OTHER is not in the universe
        logfc = pd.DataFrame(rng.normal(size=(60, 3)), index=pd.Index(genes, name="gene"), columns=perts)
        z = pd.DataFrame(rng.normal(size=(60, 2)), index=pd.Index(genes[::-1], name="gene"), columns=["LV1", "LV2"])
        z["LV2"] = logfc.loc[z.index, "OTHER"].to_numpy()
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            logfc.to_csv(tmp / "logfc.tsv.gz", sep="\t")
            z.to_csv(tmp / "pos.csv")
            pd.DataFrame({"perturbation": perts, "n_de_universe": [25, 30, 40]}).to_csv(tmp / "sizes.csv", index=False)
            rc.main(["--z", str(tmp / "pos.csv"), "--logfc", str(tmp / "logfc.tsv.gz"), "--sizes", str(tmp / "sizes.csv"),
                     "--strategy", "main", "--model", "TOY", "--n-perm", "9", "--seed", "1", "--alpha", "0.1",
                     "--out-dir", str(tmp / "out")])
            per = pd.read_csv(tmp / "out" / "perturbations.csv")
            null = pd.read_csv(tmp / "out" / "null.csv.gz")
            summary = pd.read_csv(tmp / "out" / "summary.csv")
        self.assertEqual(per.self_masked.tolist(), [True, True, False])
        self.assertEqual(per.n_de.tolist(), [25, 30, 40])
        other = per.set_index("perturbation").loc["OTHER"]
        self.assertEqual(other.best_LV, "LV2")
        self.assertAlmostEqual(other.best_r, 1.0, places=10)
        self.assertTrue(other.captured)
        self.assertEqual(null.perm.tolist(), list(range(10)))
        self.assertEqual(summary.statistic.tolist(), list(rc.STATISTICS))
        self.assertEqual(summary.self_masked.unique().tolist(), [2])


if __name__ == "__main__":
    unittest.main()

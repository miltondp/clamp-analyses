import importlib.util
import itertools
from pathlib import Path
import unittest

import numpy as np

spec = importlib.util.spec_from_file_location("l1000_metalv", Path(__file__).resolve().parents[1] / "scripts/archs4/crispercas/l1000_metalv.py")
l1000 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(l1000)


class ProjectionStatisticsTest(unittest.TestCase):
    def test_exact_null_moments_with_fixed_missing_positions(self):
        full_w = np.array([[1., 100., 2., 4.], [-1., -100., -2., -4.]])
        measured = np.array([True, False, True, True])
        v = np.array([[1., -2.], [3., 0.], [6., 2.]])
        w = full_w[:, measured]
        draws = np.array([w[:, p] @ v for p in itertools.permutations(range(3))])
        mu, sd = l1000.null_moments(w, v)
        np.testing.assert_allclose(mu, draws.mean(axis=0))
        np.testing.assert_allclose(sd, draws.std(axis=0))
        full_v = np.zeros((4, 2))
        full_v[measured] = v
        np.testing.assert_allclose(full_w @ full_v, w @ v)

    def test_ties_tails_plus_one_joint_bh_and_direction(self):
        observed = np.array([[2., -2.], [0., 1.]])
        upper = np.array([[0, 98], [50, 99]])
        lower = np.array([[99, 1], [2, 99]])
        pu, pl, p, q, direction = l1000.empirical_statistics(observed, upper, lower, 99)
        np.testing.assert_allclose(pu, [[.01, .99], [.51, 1.]])
        np.testing.assert_allclose(pl, [[1., .02], [.03, 1.]])
        np.testing.assert_allclose(p, [[.02, .04], [.06, 1.]])
        np.testing.assert_allclose(q, [[.08, .08], [.08, 1.]])
        np.testing.assert_array_equal(direction, [["Activation-like", "Suppression-like"], ["Neutral", "Activation-like"]])
        w = np.ones((1, 3))
        v = np.array([[1.], [2.], [3.]])
        upper, lower = l1000.permutation_tails(w, v, w @ v, 5, 7)
        np.testing.assert_array_equal(upper, [[5]])
        np.testing.assert_array_equal(lower, [[5]])
        self.assertEqual(l1000.empirical_statistics(w @ v, upper, lower, 5)[2].item(), 1.)

    def test_seed_reproducibility_and_both_tails(self):
        w = np.array([[1., 3., -2.], [-1., -3., 2.]])
        v = np.array([[2., 4.], [5., -1.], [1., 3.]])
        first = l1000.permutation_tails(w, v, w @ v, 65, 17)
        second = l1000.permutation_tails(w, v, w @ v, 65, 17)
        for a, b in zip(first, second):
            np.testing.assert_array_equal(a, b)
        np.testing.assert_array_equal(first[0][0], first[1][1])
        np.testing.assert_array_equal(first[1][0], first[0][1])


if __name__ == "__main__":
    unittest.main()

import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_samples

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cev = load("cev", HERE.parent / "cev.py")
# Loaded, never edited (see PLAN.md Key facts): the reference the builder must match.
metalv = load("cluster_metalv", HERE.parents[1] / "archs4/crispercas/cluster_metalv.py")


def toy_runs(n_runs=5, n_modules=8, n_genes=300, n_noise=6, seed=0):
    """Non-negative sparse modules, each present in a run with prob 0.8, plus noise LVs."""
    rng = np.random.default_rng(seed)
    base = rng.gamma(1.0, 1.0, (n_modules, n_genes)) * (rng.random((n_modules, n_genes)) < 0.15)
    X, models, lvs = [], [], []
    for r in range(n_runs):
        rows = [np.clip(base[m] * rng.uniform(0.5, 2) + rng.normal(0, 0.15, n_genes) * (base[m] > 0), 0, None)
                for m in range(n_modules) if rng.random() < 0.8]
        rows += [rng.gamma(1.0, 1.0, n_genes) * (rng.random(n_genes) < 0.15) for _ in range(n_noise)]
        X += rows
        models += [f"R{r + 1}"] * len(rows)
        lvs += [f"LV{i + 1}" for i in range(len(rows))]
    return np.asarray(X, dtype=np.float32), np.asarray(models), np.asarray(lvs), np.asarray([f"g{i}" for i in range(n_genes)])


class BuilderReproducesClusterMetalvTest(unittest.TestCase):
    def test_hier_median_equals_cluster_metalv(self):
        X, models, lvs, genes = toy_runs()
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            np.savez(tmp / "Z_all.npz", X=X, models=models, lvs=lvs, genes=genes)
            flagged = np.zeros(len(lvs), bool)
            flagged[[2, 11, 20]] = True
            pd.DataFrame({"model": models, "LV": lvs, "gene_flag": flagged, "technical_joint": 0.0,
                          "sample_flag": False, "any_flag": flagged}).to_csv(tmp / "flags.csv", index=False)
            for min_cos, min_models in [(0.5, 3), (0.7, 2), (0.3, 4)]:
                out = tmp / f"ref_{min_cos}_{min_models}"
                metalv.cluster_main(["--z-cache", str(tmp / "Z_all.npz"), "--flags", str(tmp / "flags.csv"),
                                     "--min-cos", str(min_cos), "--min-models", str(min_models),
                                     "--reference", "R1", "--out-dir", str(out)])
                ref_z = pd.read_csv(out / "metaZ.csv", index_col=0)
                ref_s = pd.read_csv(out / "cluster_summary.csv")

                pool = cev.load_pool(tmp / "Z_all.npz", tmp / "flags.csv")
                labels = cev.hier_partition(cev.hier_linkage(pool.Xn, "complete"), min_cos)
                Z, set_lvs, members, _ = cev.build_set(pool, labels, min_models)

                self.assertGreater(Z.shape[1], 0)
                self.assertEqual(list(set_lvs.metaLV), list(ref_z.columns))
                self.assertEqual(list(pool.genes), list(ref_z.index))
                np.testing.assert_array_equal(Z, ref_z.to_numpy().astype(np.float32))
                ours = members.assign(key=members.model + ":" + members.LV).groupby("metaLV", sort=False).key.apply(";".join)
                self.assertEqual(list(ours.loc[ref_s.metaLV]), list(ref_s.members))
                self.assertEqual(list(set_lvs.n_lvs), list(ref_s.n_lvs))
                self.assertEqual(list(set_lvs.n_models), list(ref_s.n_models))

    def test_condensed_distance_matches_pdist(self):
        X = cev.l2norm(toy_runs(n_runs=2)[0])
        np.testing.assert_allclose(cev.condensed_cosine_distance(X), np.clip(pdist(X, "cosine"), 0, 2), atol=1e-6)


class MinModelsTest(unittest.TestCase):
    def test_resolution(self):
        expected = {10: [3, 4, 5, 6, 7], 9: [3, 4, 5, 6, 7], 5: [2, 2, 3, 3, 4], 20: [6, 8, 10, 12, 14]}
        for n, want in expected.items():
            self.assertEqual([cev.resolve_min_models(f, n) for f in [0.3, 0.4, 0.5, 0.6, 0.7]], want, n)


class MetricsTest(unittest.TestCase):
    def test_top_genes_positive_only_and_stable(self):
        M = np.array([[0, 3, 1, 3, 0, 2], [0, 0, 5, 0, 0, 0]], dtype=np.float32)
        B = cev.top_genes(M, 2)
        self.assertEqual(B[0].nonzero()[0].tolist(), [1, 3])
        self.assertEqual(B[1].nonzero()[0].tolist(), [2])

    def test_iq_silhouette_brute_force_and_sklearn(self):
        rng = np.random.default_rng(1)
        Xm = cev.l2norm(np.abs(rng.normal(size=(9, 6))).astype(np.float32))
        labels = np.array([0, 0, 0, 1, 1, 2, 2, 2, 3])  # cluster 3 is a singleton
        cev.BLOCK, old = 4, cev.BLOCK  # exercise several blocks
        try:
            iq, sil = cev.iq_silhouette(Xm, labels)
        finally:
            cev.BLOCK = old
        S = Xm @ Xm.T
        for c in range(3):
            i, o = labels == c, labels != c
            intra = S[np.ix_(i, i)][~np.eye(i.sum(), dtype=bool)].mean()
            self.assertAlmostEqual(iq[c], intra - S[np.ix_(i, o)].mean(), places=5)
        self.assertTrue(np.isnan(iq[3]))
        np.testing.assert_allclose(sil, silhouette_samples(Xm, labels, metric="cosine"), atol=1e-5)

    def test_redundancy(self):
        Zn = cev.l2norm(np.array([[1, 0, 0], [1, 0, 0], [0, 1, 1]], dtype=np.float32))
        B = cev.top_genes(Zn, 2)
        cos, jac = cev.redundancy(Zn, B)
        np.testing.assert_allclose(cos, [1, 1, 0], atol=1e-6)
        np.testing.assert_allclose(jac, [1, 1, 0])

    def test_pool_best_hits_reproducible_coverage(self):
        Xn = cev.l2norm(np.array([[1, 0], [1, 0.1], [0, 1], [1, 0.05], [0.7, 0.7]], dtype=np.float32))
        models = np.array(["A", "A", "B", "C", "C"])
        hits = cev.pool_best_hits(Xn, models, ["A", "B", "C"])
        self.assertTrue(np.isnan(hits[0, 0]) and np.isnan(hits[2, 1]) and np.isnan(hits[4, 2]))
        self.assertAlmostEqual(hits[0, 2], 1, places=2)       # LV0 ~ LV3
        self.assertAlmostEqual(hits[2, 2], np.sqrt(0.5), places=5)
        repro = cev.reproducible(hits, 0.9, 2)
        self.assertEqual(repro.tolist(), [True, True, False, True, False])
        best = cev.best_cos(Xn, cev.l2norm(np.array([[1, 0]], dtype=np.float32)))
        cov = cev.coverage(best, models, ["A", "B", "C"], repro, [0.9])
        share = cov.set_index("scope").share
        self.assertEqual(share[["A", "B", "C", "all", "reproducible"]].tolist(), [1.0, 0.0, 0.5, 0.6, 1.0])

    def test_cluster_metrics_purity_and_cohesion(self):
        X, models, lvs, genes = toy_runs(n_runs=3, n_modules=3, n_genes=60, n_noise=0)
        cache = cev.LVTable(X, cev.l2norm(X), models, lvs, pd.Index(genes))
        rows = np.array([0, 1, 3, 2])          # two members of run 1 in cluster 0
        codes = np.array([0, 0, 0, 1])
        Z = np.stack([np.median(X[rows[:3]], 0), X[2]], axis=1)
        per, sil = cev.cluster_metrics(cache, rows, codes, Z, top_n=10)
        n_runs_c0 = len(set(models[rows[:3]]))
        self.assertAlmostEqual(per.purity[0], n_runs_c0 / 3)
        self.assertEqual(per.n_lvs.tolist(), [3, 1])
        self.assertTrue(np.isnan(per.mean_cos[1]))
        self.assertAlmostEqual(per.cos_to_consensus[1], 1, places=5)
        C = cache.Xn[rows[:3]] @ cache.Xn[rows[:3]].T
        self.assertAlmostEqual(per.min_cos[0], C[np.triu_indices(3, 1)].min(), places=6)
        self.assertEqual(sil[3], 0)


if __name__ == "__main__":
    unittest.main()

import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

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
                self.assertEqual(list(pool.genes), list(ref_z.index))
                self.assertEqual(cev.compare_to_metalv(Z, set_lvs, members, ref_z, ref_s), [])

    def test_condensed_distance_matches_pdist(self):
        X = cev.l2norm(toy_runs(n_runs=2)[0])
        np.testing.assert_allclose(cev.condensed_cosine_distance(X), np.clip(pdist(X, "cosine"), 0, 2), atol=1e-6)


class JaccardLinkageTest(unittest.TestCase):
    def test_condensed_jaccard_matches_pdist(self):
        rng = np.random.default_rng(4)
        B = rng.random((300, 80)) < 0.1   # > 256 rows: several blocks
        B[[17, 290]] = False              # empty rows: distance 1 (scipy gives 0 when both are empty)
        d = cev.condensed_jaccard_distance(B)
        ref = pdist(B, "jaccard")
        n = len(B)
        i, j = np.triu_indices(n, 1)
        empty = ~B.any(1)
        both = empty[i] & empty[j]
        np.testing.assert_allclose(d[~both], ref[~both], rtol=0, atol=1e-12)
        self.assertTrue((d[empty[i] | empty[j]] == 1).all())

    def test_planted_modules_recovered(self):
        # 6 modules x 4 runs; each run's copy keeps ~80% of the module's 40 genes (within Jaccard ~0.65)
        # and adds 8 random genes, so between-module Jaccard stays near 0.
        rng = np.random.default_rng(5)
        n_genes, mods = 2000, np.arange(6)
        genes = [rng.choice(n_genes, 40, replace=False) for _ in mods]
        rows, truth = [], []
        for _ in range(4):
            for m in mods:
                b = np.zeros(n_genes, dtype=bool)
                b[rng.choice(genes[m], 32, replace=False)] = True
                b[rng.choice(n_genes, 8, replace=False)] = True
                rows.append(b)
                truth.append(m)
        B, truth = np.asarray(rows), np.asarray(truth)
        for method in ["complete", "average"]:
            labels = cev.hier_partition(cev.hier_linkage_jac(B, method), 0.3)
            self.assertEqual(len(np.unique(labels)), len(mods), method)
            self.assertTrue(all(len(np.unique(labels[truth == m])) == 1 for m in mods), method)


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

    def test_redundancy(self):
        Zn = cev.l2norm(np.array([[1, 0, 0], [1, 0, 0], [0, 1, 1]], dtype=np.float32))
        B = cev.top_genes(Zn, 2)
        cos, jac = cev.redundancy(Zn, B)
        np.testing.assert_allclose(cos, [1, 1, 0], atol=1e-6)
        np.testing.assert_allclose(jac, [1, 1, 0])

    def test_pool_best_hits_reproducible_coverage(self):
        Xn = cev.l2norm(np.array([[1, 0], [1, 0.1], [0, 1], [1, 0.05], [0.7, 0.7]], dtype=np.float32))
        models = np.array(["A", "A", "B", "C", "C"])
        B = cev.top_genes(Xn, 1)
        hits, arg, jac = cev.pool_best_hits(Xn, B, models, ["A", "B", "C"])
        self.assertTrue(np.isnan(hits[0, 0]) and np.isnan(hits[2, 1]) and np.isnan(hits[4, 2]))
        self.assertEqual(arg[0].tolist(), [-1, 2, 3])
        self.assertAlmostEqual(hits[0, 2], 1, places=2)       # LV0 ~ LV3
        self.assertAlmostEqual(hits[2, 2], np.sqrt(0.5), places=5)
        self.assertEqual(jac[0, 2], 1)                        # both top gene 0
        self.assertEqual(jac[0, 1], 0)
        repro = cev.reproducible(hits >= 0.9, 2)
        self.assertEqual(repro.tolist(), [True, True, False, True, False])
        best = cev.best_cos(Xn, cev.l2norm(np.array([[1, 0]], dtype=np.float32)))
        # reproducible_multi: reproducible LVs 0, 1, 3; only LV1 counted as multi-gene, and it is covered
        multi = np.array([False, True, True, False, True])
        cov = cev.coverage(best, models, ["A", "B", "C"], repro, multi, [0.9]).set_index("scope")
        self.assertEqual(cov.share[["A", "B", "C", "all", "reproducible"]].tolist(), [1.0, 0.0, 0.5, 0.6, 1.0])
        self.assertEqual((cov.n["reproducible_multi"], cov.covered["reproducible_multi"]), (1, 1))

    def test_cluster_metrics_purity_and_cohesion(self):
        X, models, lvs, genes = toy_runs(n_runs=3, n_modules=3, n_genes=60, n_noise=0)
        cache = cev.LVTable(X, models, lvs, pd.Index(genes))
        rows = np.array([0, 1, 3, 2])          # two members of run 1 in cluster 0
        codes = np.array([0, 0, 0, 1])
        Z = np.stack([np.median(X[rows[:3]], 0), X[2]], axis=1)
        per = cev.cluster_metrics(cache, rows, codes, Z, top_n=10)
        n_runs_c0 = len(set(models[rows[:3]]))
        self.assertAlmostEqual(per.purity[0], n_runs_c0 / 3)
        self.assertEqual(per.n_lvs.tolist(), [3, 1])
        self.assertTrue(np.isnan(per.mean_cos[1]))
        self.assertAlmostEqual(per.cos_to_consensus[1], 1, places=5)
        self.assertAlmostEqual(per.overlap_to_consensus[1], 1)
        C = cache.Xn[rows[:3]] @ cache.Xn[rows[:3]].T
        self.assertAlmostEqual(per.min_cos[0], C[np.triu_indices(3, 1)].min(), places=6)

    def test_overlap_and_top1_share(self):
        Bm = np.array([[1, 1, 0, 0], [1, 0, 1, 1], [0, 0, 0, 0]], dtype=bool)
        bz = np.array([1, 1, 1, 0], dtype=bool)
        got = cev.overlap(Bm, bz)
        np.testing.assert_allclose(got[:2], [1.0, 2 / 3])
        self.assertTrue(np.isnan(got[2]))               # member with no positive genes
        M = np.array([[3, 4, 0], [1, 1, 1]], dtype=np.float32)
        np.testing.assert_allclose(cev.top1_share(M), [16 / 25, 1 / 3])


class ReproducibleTest(unittest.TestCase):
    def test_broad_lv_is_non_reciprocal_hit(self):
        # Run A has two specific LVs (genes 0 and 1); run B has a broad LV (genes 0+1, leaning to
        # gene 1) and a gene-0 LV. First every best hit is reciprocal; then B's second LV moves
        # towards gene 1, so A0's best hit becomes the broad LV, whose own best hit is A1.
        X = np.array([[1, 0, 0, 0],      # A0: gene 0
                      [0, 1, 0, 0],      # A1: gene 1
                      [1, 1.2, 0, 0],    # B0: broad, closer to A1
                      [0.9, 0, 0.1, 0]], dtype=np.float32)  # B1: gene 0
        models = np.array(["A", "A", "B", "B"])
        runs = ["A", "B"]
        Xn = cev.l2norm(X)
        hits, arg, jac = cev.pool_best_hits(Xn, cev.top_genes(X, 2), models, runs)
        self.assertEqual(arg[:, 1].tolist(), [3, 2, -1, -1])   # A0 -> B1, A1 -> B0
        self.assertEqual(arg[:, 0].tolist(), [-1, -1, 1, 0])   # B0 -> A1, B1 -> A0
        recip = cev.reciprocal(arg, models, runs)
        self.assertEqual(recip[:, 1].tolist(), [True, True, False, False])
        # Make B1 broad too, so that A0's best hit is B0, which reciprocates A1 instead.
        X[3] = [0.5, 1.4, 0, 0]
        hits, arg, jac = cev.pool_best_hits(cev.l2norm(X), cev.top_genes(X, 2), models, runs)
        recip = cev.reciprocal(arg, models, runs)
        self.assertEqual(arg[0, 1], 2)
        self.assertEqual(arg[2, 0], 1)
        self.assertFalse(recip[0, 1])
        # Jaccard guard: a threshold at the pair's Jaccard removes the hit.
        stages = cev.hit_stages(hits, recip, jac, 0.5, jac[0, 1])
        self.assertEqual([cev.reproducible(m, 2)[0] for m in stages.values()], [True, False, False])
        self.assertFalse(stages["cos"][0, 0])                  # own run (NaN) never qualifies

    def test_jaccard_null(self):
        rng = np.random.default_rng(0)
        n_genes, size, n_lvs = 1000, 20, 30
        modules = [rng.choice(n_genes, size, replace=False) for _ in range(n_lvs)]
        B = np.zeros((2 * n_lvs, n_genes), dtype=bool)
        for r in range(2):
            for i, m in enumerate(modules):
                keep = m[rng.random(size) < 0.8]        # shared module, ~80% of genes kept
                B[r * n_lvs + i, keep] = True
        models = np.repeat(["A", "B"], n_lvs)
        thr, quant, pairs = cev.jaccard_null(B, models, ["A", "B"], 2, seed=5, q=0.999)
        observed = np.diag(cev.jaccard(B[:n_lvs], B[n_lvs:]))
        self.assertLess(thr, 0.15)
        self.assertGreater(observed.min(), 2 * thr)
        self.assertEqual(sorted(map(tuple, pairs)), [("A", "B"), ("B", "A")])
        self.assertEqual(cev.jaccard_null(B, models, ["A", "B"], 2, seed=5, q=0.999)[0], thr)
        self.assertLessEqual(quant["p99.9"], quant["max"])

    def test_quantiles_labels(self):
        q = cev.quantiles(np.arange(101), [0.01, 0.5, 0.999])
        self.assertEqual(list(q), ["p01", "p50", "p99.9"])
        self.assertAlmostEqual(q["p50"], 50)


class RobustnessTest(unittest.TestCase):
    def test_load_pool_runs_and_gene_perm(self):
        X, models, lvs, genes = toy_runs()
        with tempfile.TemporaryDirectory() as tmp:
            npz = Path(tmp) / "Z_all.npz"
            np.savez(npz, X=X, models=models, lvs=lvs, genes=genes)
            full = cev.load_pool(npz)
            sub = cev.load_pool(npz, runs=["R4", "R2"])
            keep = np.isin(models, ["R2", "R4"])
            np.testing.assert_array_equal(sub.X, X[keep])            # cache order, not the runs order
            np.testing.assert_array_equal(sub.lvs, lvs[keep])
            perm = cev.load_pool(npz, gene_perm_seed=1)
            np.testing.assert_array_equal(np.sort(perm.X, 1), np.sort(full.X, 1))   # same values per LV
            self.assertFalse(np.array_equal(perm.X, full.X))
            np.testing.assert_array_equal(perm.X, cev.load_pool(npz, gene_perm_seed=1).X)
            with self.assertRaises(AssertionError):
                cev.load_pool(npz, runs=["R9"])

    def test_hungarian_and_split_half(self):
        rng = np.random.default_rng(0)
        Za = (rng.gamma(1.0, 1.0, (200, 12)) * (rng.random((200, 12)) < 0.3)).astype(np.float32)
        perm = rng.permutation(12)
        for by in ["cos", "jac"]:
            ia, ib, cos, jac = cev.hungarian_match(Za, Za[:, perm], top_n=20, by=by)
            np.testing.assert_array_equal(perm[ib], ia)
            np.testing.assert_allclose(cos, 1, atol=1e-5)
            np.testing.assert_allclose(jac, 1)
        # half b: the same 12 LVs permuted plus 4 unrelated ones -> 24 of 28 LVs matched
        extra = np.zeros((200, 4), dtype=np.float32)
        for j in range(4):
            extra[j * 50:(j + 1) * 50, j] = 1   # disjoint gene blocks: low cosine to everything in Za
        Zb = np.hstack([Za[:, perm], extra])
        row = cev.split_half_scores(Za, Zb, top_n=20, match_cos=0.9, match_jac=0.25, single_gene=0.5)
        self.assertEqual((row["n_a"], row["n_b"], row["n_pairs"]), (12, 16, 12))
        self.assertAlmostEqual(row["share_all"], 24 / 28)
        self.assertAlmostEqual(row["jshare_all"], 24 / 28)
        self.assertAlmostEqual(row["jpair_cos_all"], 1, places=5)
        self.assertAlmostEqual(row["pair_jac_all"], 1)
        self.assertEqual(row["n_single"] + row["n_multi"], 28)

    def test_empty_consensus_lv_never_matches(self):
        Za = np.eye(5, 3, dtype=np.float32)
        Zb = Za.copy()
        Zb[:, 2] = 0                      # all-zero consensus LV
        row = cev.split_half_scores(Za, Zb, top_n=2, match_cos=0.9, match_jac=0.25, single_gene=0.5)
        self.assertEqual(row["n_empty"], 1)
        self.assertAlmostEqual(row["share_all"], 4 / 6)
        self.assertAlmostEqual(row["jshare_all"], 4 / 6)
        best, best_j = cev.loro_recovery(Zb, [Za], top_n=2)
        np.testing.assert_allclose(best[:, 0], [1, 1, 0], atol=1e-6)
        np.testing.assert_allclose(best_j[:, 0], [1, 1, 0])

    def test_loro_recovery(self):
        Z = np.eye(4, 3, dtype=np.float32)
        best, best_j = cev.loro_recovery(Z, [Z[:, :2], np.zeros((4, 0), np.float32)], top_n=1)
        np.testing.assert_allclose(best[:, 0], [1, 1, 0], atol=1e-6)
        np.testing.assert_allclose(best_j[:, 0], [1, 1, 0])
        self.assertTrue(np.isnan(best[:, 1]).all() and np.isnan(best_j[:, 1]).all())

    def test_label_perm_null(self):
        # 6 clusters of one LV per run (span 3 runs) and 3 clusters drawn from one run (span 1)
        labels = np.r_[np.repeat(np.arange(6), 3), np.repeat([6, 7, 8], 2)]
        models = np.r_[np.tile(["A", "B", "C"], 6), ["A"] * 6]
        tab = cev.label_perm_null(labels, models, n_perm=5, seed=0)
        obs = tab.groupby("m")["observed"].first()
        self.assertEqual(obs.tolist(), [9, 6, 6])
        self.assertTrue((tab[tab.m == 1].null == 9).all())          # every cluster spans >= 1 run
        self.assertEqual(len(tab), 5 * 3)


if __name__ == "__main__":
    unittest.main()

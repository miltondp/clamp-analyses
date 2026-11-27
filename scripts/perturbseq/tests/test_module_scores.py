"""Offline tests for module_scores.py, including a check against the Saelens et al. 2018 scorer."""
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import module_scores as ms
import ora_null as on

CUTOFFS = np.linspace(0, 3, 100)


# The paper's scorer (lib/module-evaluation/code/lib/modulescomparison.py: test_enrichment,
# modbindevalscorer, cal_aucodds, cal_faucodds), kept as close to the original as numpy 2 allows.
def paper_cal_aucodds(odds, cutoffs=CUTOFFS):
    maxodds = odds.max(1)
    if len(maxodds) == 0:
        return 0
    maxodds[maxodds == 0] = 1e-10
    stillenriched = [(np.log10(maxodds) >= cutoff).sum() / len(maxodds) for cutoff in cutoffs]
    return np.trapezoid(stillenriched, cutoffs) / (cutoffs[-1] - cutoffs[0])


def paper_cal_faucodds(odds, cutoffs=CUTOFFS):
    a1, a2 = paper_cal_aucodds(odds, cutoffs), paper_cal_aucodds(odds.T, cutoffs)
    return 0 if a1 + a2 == 0 else 2 / (1 / a1 + 1 / a2)


def paper_scores(modmem, binmem, q=0.05):
    """modmem: genes x modules, binmem: genes x regulators (DataFrames of 0/1 over the same genes)."""
    modsizes, binsizes = modmem.sum(), binmem.sum()
    tps = modmem.T.dot(binmem.astype(int))
    fps = binsizes - tps
    fns = (modsizes - tps.T).T
    tns = binmem.shape[0] - tps - fps - fns
    with np.errstate(divide="ignore", invalid="ignore"):
        odds = (tps * tns) / (fps * fns)
    pvals = np.vectorize(lambda a, b, c, d: fisher_exact([[a, b], [c, d]], alternative="greater")[1])(
        tps.values, fps.values, fns.values, tns.values)
    with np.errstate(divide="ignore"):  # statsmodels takes log1p(-1) for p = 1
        qvals = np.array([multipletests(row)[1] for row in pvals])  # default method: Holm-Sidak
    filtered = odds.values.copy()
    filtered[qvals > q] = 0
    filtered = np.nan_to_num(filtered, nan=0.0, posinf=np.inf)  # 0/0 cells never pass the q filter
    regulators = paper_cal_aucodds(filtered.T.copy())  # modbindevalscorer: max over modules
    return {"aucodds_regulons": regulators, "aucodds_modules": paper_cal_aucodds(filtered.copy()),
            "faucodds": paper_cal_faucodds(filtered.copy()), "odds": filtered}


def toy(seed=0, n_genes=60, max_size=100):
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n_genes)]
    terms = {f"R{j}": list(rng.choice(genes, size=int(rng.integers(6, 15)), replace=False)) for j in range(6)}
    covered = set().union(*terms.values())
    terms["REST"] = [g for g in genes if g not in covered] + list(rng.choice(genes, 6, replace=False))
    ora = on.Ora(genes, terms, top_pct=0.1, min_size=5, max_size=max_size, fdr=0.05)
    assert ora.N == n_genes, "every toy gene is annotated, so both universes agree"
    member = np.zeros((5, n_genes), bool)
    gi = {g: i for i, g in enumerate(genes)}
    member[0, [gi[g] for g in terms["R0"]]] = True  # exactly R0: fn = 0, OR = inf
    member[1, [gi[g] for g in terms["R1"][:-2]] + [gi[g] for g in terms["R2"][:2]]] = True
    member[2, rng.choice(n_genes, 8, replace=False)] = True  # noise
    member[3, [gi[g] for g in terms["R3"]] + list(rng.choice(n_genes, 10, replace=False))] = True
    member[4, rng.choice(n_genes, 5, replace=False)] = True
    return genes, terms, ora, member


def args(**kw):
    base = dict(zscore_sd=3.0, module_min_size=5, q=0.05, log10_or_max=3.0)
    return SimpleNamespace(**{**base, **kw})


class RightTailTests(unittest.TestCase):
    def test_matches_scipy(self):
        sizes = np.array([5, 10, 82, 550, 1000])
        for n in (1, 5, 75, 300):
            ours = on.right_tail(7483, sizes, n)
            ref = hypergeom.sf(np.arange(n + 1) - 1, 7483, sizes[:, None], n)
            ok = ref > 1e-280
            np.testing.assert_allclose(ours[ok], ref[ok], rtol=1e-9)
            self.assertTrue((ours[~ok] < 1e-270).all())
            np.testing.assert_allclose(ours[:, 0], 1.0)


class HolmSidakTests(unittest.TestCase):
    def test_matches_statsmodels(self):
        rng = np.random.default_rng(1)
        p = np.concatenate([rng.uniform(size=(20, 30)) ** 4, np.ones((20, 5)), np.full((20, 1), 1e-300)], axis=1)
        p[3] = 0.01  # ties
        with np.errstate(divide="ignore"):
            expected = np.array([multipletests(row, method="hs")[1] for row in p])
        np.testing.assert_allclose(ms.holm_sidak(p), expected, rtol=1e-12)


class AgainstPaperTests(unittest.TestCase):
    def test_scores_match_paper_scorer(self):
        for seed in range(5):
            genes, _, ora, member = toy(seed)
            res = ms.enrichment(member, ora)
            q = ms.holm_sidak(res["p"])
            odds = np.where(q <= 0.05, res["odds"], 0.0)
            ours = ms.scores(odds, 3.0)
            modmem = pd.DataFrame(member.T.astype(int), index=genes)
            binmem = pd.DataFrame(ora.membership.toarray(), index=genes, columns=ora.terms)
            ref = paper_scores(modmem, binmem)
            np.testing.assert_array_equal(odds, ref["odds"])
            for key in ms.SCORES:
                self.assertAlmostEqual(ours[key], ref[key], places=12, msg=f"seed {seed} {key}")
            self.assertTrue(np.isinf(odds[0, list(ora.terms).index("R0")]), "module == regulon gives OR inf")
            self.assertGreater(ours["faucodds"], 0)

    def test_p_values_are_fisher_right_tail(self):
        genes, _, ora, member = toy(2)
        res = ms.enrichment(member, ora)
        for i, j in [(0, 0), (1, 1), (2, 3), (3, 3), (4, 5)]:
            tp = res["tp"][i, j]
            fp, fn = ora.M[j] - tp, res["n"][i] - tp
            tn = ora.N - tp - fp - fn
            expected = fisher_exact([[tp, fp], [fn, tn]], alternative="greater")[1]
            self.assertAlmostEqual(res["p"][i, j], expected, places=10)


class AucoddsTests(unittest.TestCase):
    def test_edge_cases(self):
        self.assertEqual(ms.aucodds([], 3), 0.0)
        self.assertEqual(ms.aucodds([0, 0], 3), 0.0)
        self.assertAlmostEqual(ms.aucodds([np.inf, np.inf], 3), 1.0)
        self.assertAlmostEqual(ms.aucodds([np.inf, 0], 3), 0.5)
        # log10 OR = 1.5 survives cutoffs up to 1.5: about half the area.
        self.assertAlmostEqual(ms.aucodds([10 ** 1.5], 3), 1.5 / 3, delta=0.02)
        self.assertEqual(ms.harmonic(0.0, 0.7), 0.0)
        self.assertEqual(ms.scores(np.zeros((0, 4)), 3), dict.fromkeys(ms.SCORES, 0.0))

    def test_orientation(self):
        # Two modules, three regulons: one strong pair only. Regulon half: 1 of 3 regulons enriched;
        # module half: 1 of 2 modules enriched.
        odds = np.array([[np.inf, 0, 0], [0, 0, 0]])
        s = ms.scores(odds, 3)
        self.assertAlmostEqual(s["aucodds_regulons"], 1 / 3)
        self.assertAlmostEqual(s["aucodds_modules"], 1 / 2)


class DefinitionTests(unittest.TestCase):
    def test_zscore_modules(self):
        rng = np.random.default_rng(3)
        z = rng.normal(size=(500, 2))
        z[:7, 0] = 10
        z[7:10, 1] = -10
        pos, neg = ms.zscore_modules(z, 3), ms.zscore_modules(-z, 3)
        self.assertEqual(pos.shape, (2, 500))
        self.assertTrue(pos[0, :7].all() and neg[1, 7:10].all())
        x = z[:, 0]
        np.testing.assert_array_equal(pos[0], x > x.mean() + 3 * x.std())
        x = z[:, 1]
        np.testing.assert_array_equal(neg[1], x < x.mean() - 3 * x.std())

    def test_nonnegative_z_has_no_negative_top_modules(self):
        genes, _, ora, _ = toy(0)
        z = np.abs(np.random.default_rng(4).normal(size=(len(genes), 3)))
        res, labels = ms.score_definition(z, ora, "top", args())
        self.assertEqual(labels["tail"].tolist(), ["pos"] * 3)
        self.assertEqual(res["n"].tolist(), [ora.query_target] * 3)

    def test_small_modules_are_dropped(self):
        genes, _, ora, _ = toy(0)
        z = np.zeros((len(genes), 2))
        z[:4, 0] = 10  # 4 genes beyond 3 SD: below module_min_size 5
        z[:5, 1] = 10  # 5 genes: kept
        res, labels = ms.score_definition(z, ora, "zscore", args())
        self.assertEqual(labels["lv_index"].tolist(), [1])
        self.assertEqual(res["n"].tolist(), [5])


class RunTests(unittest.TestCase):
    def test_run_writes_outputs_and_observed_is_perm_zero(self):
        genes, terms, ora, member = toy(0, n_genes=200, max_size=1000)
        rng = np.random.default_rng(5)
        z = rng.normal(size=(len(genes), 4))
        z[member[0], 0] += 20  # LV1 carries R0
        z[member[3], 2] -= 20  # LV3's negative tail carries R3
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            pd.DataFrame(z, index=genes, columns=[f"LV{i + 1}" for i in range(4)]).to_csv(tmp / "z.csv")
            with open(tmp / "r.gmt", "w") as handle:
                for t, g in terms.items():
                    handle.write("\t".join([t, "na", *g]) + "\n")
            ms.main(["--z", str(tmp / "z.csv"), "--gmt", str(tmp / "r.gmt"), "--strategy", "s",
                     "--model", "m", "--top-pct", "0.1", "--min-size", "5", "--max-size", "1000",
                     "--module-min-size", "5", "--zscore-sd", "3", "--q", "0.05",
                     "--log10-or-max", "3", "--n-perm", "4", "--seed", "7", "--out-dir", str(tmp / "out")])
            null = pd.read_csv(tmp / "out" / "null.csv.gz")
            summary = pd.read_csv(tmp / "out" / "summary.csv")
            hits = pd.read_csv(tmp / "out" / "module_hits.csv", keep_default_na=False)
            regulons = pd.read_csv(tmp / "out" / "regulon_hits.csv", keep_default_na=False)

        self.assertEqual(null["perm"].tolist(), list(range(5)))
        observed, _ = ms.statistics(z, ora, args())
        for key, value in observed.items():
            self.assertAlmostEqual(null.loc[0, key], value)
        # Permutations come from the shared stream: perm 1 is ora_null.permuted's first matrix.
        _, first = next(on.permuted(z, 4, 7))
        for key, value in ms.statistics(first, ora, args())[0].items():
            self.assertAlmostEqual(null.loc[1, key], value)
        self.assertEqual(set(summary["definition"]), set(ms.DEFINITIONS))
        self.assertEqual(set(summary["statistic"]), set(ms.SCORES))
        np.testing.assert_allclose(summary["fold"], summary["observed"] / summary["null_mean"])
        best = hits.set_index(["definition", "LV", "tail"])
        self.assertEqual(best.loc[("zscore", "LV1", "pos"), "best_regulon"], "R0")
        self.assertEqual(best.loc[("zscore", "LV3", "neg"), "best_regulon"], "R3")
        self.assertEqual(regulons.set_index(["definition", "regulon"]).loc[("zscore", "R0"), "best_LV"], "LV1")


if __name__ == "__main__":
    unittest.main()

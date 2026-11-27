"""Consensus-LV library: pool loading, hierarchical sets and computational metrics (see PLAN.md).

numpy/scipy/pandas only (no cupy/cuml), so it imports in every env. The builder copies
cluster_metalv.py (which must not be edited) and is checked against it in tests/test_cev.py
and, on real data, by the cev_set_check rule.
"""
import math
from dataclasses import dataclass
from functools import cached_property

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, linkage

BLOCK = 1024  # rows per similarity block


# ------------------------------------------------------------
# Loading and distance
# ------------------------------------------------------------

def l2norm(M):
    return M / np.linalg.norm(M, axis=1, keepdims=True)


def lv_keys(models, lvs):
    """'model:LV' keys, the member format of cluster_metalv.py's cluster_summary.csv."""
    return [f"{m}:{l}" for m, l in zip(models, lvs)]


@dataclass
class LVTable:
    X: np.ndarray       # n x genes, raw loadings (float32)
    models: np.ndarray
    lvs: np.ndarray
    genes: pd.Index

    @cached_property
    def Xn(self):
        """L2-normalized rows, built on first use (hier-set never needs them)."""
        return l2norm(self.X)

    @property
    def keys(self):
        return lv_keys(self.models, self.lvs)


def load_cache(npz):
    d = np.load(npz, allow_pickle=True)
    return d["X"], d["models"], d["lvs"], pd.Index(d["genes"])


def noisy_keys(flags_csv):
    flags = pd.read_csv(flags_csv)
    return set(zip(flags.loc[flags.any_flag, "model"], flags.loc[flags.any_flag, "LV"]))


def load_pool(npz, flags_csv=None):
    """Cache rows minus flagged LVs, in cache order (as cluster_metalv.py:101-109)."""
    X, models, lvs, genes = load_cache(npz)
    if flags_csv is not None:
        noisy = noisy_keys(flags_csv)
        keep = np.array([(m, l) not in noisy for m, l in zip(models, lvs)])
        X, models, lvs = X[keep], models[keep], lvs[keep]
    return LVTable(X, models, lvs, genes)


def condensed_cosine_distance(Xn):
    """Exact copy of cluster_metalv.py:111-117, so linkages (and tie-breaks) are identical."""
    n = len(Xn)
    dist = np.empty(n * (n - 1) // 2, dtype=np.float64)
    for start in range(0, n, 256):
        stop = min(start + 256, n)
        sim = Xn[start:stop] @ Xn.T
        for k, i in enumerate(range(start, stop)):
            off = n * i - i * (i + 1) // 2
            dist[off:off + n - i - 1] = np.clip(1 - sim[k, i + 1:], 0, 2)
    return dist


# ------------------------------------------------------------
# Hierarchical clustering and set building
# ------------------------------------------------------------

def condensed_jaccard_distance(B):
    """1 - Jaccard between boolean rows of B (e.g. top_genes), in the condensed_cosine_distance layout.

    Distance is 1 when both rows are empty (scipy's pdist gives 0 there). Sparse intersections,
    since top-gene rows hold ~185 of ~18k genes."""
    n = len(B)
    S = sparse.csr_matrix(B, dtype=np.float32)
    ST = S.T.tocsr()   # CSR @ CSR: a CSC right operand is converted back on every block
    size = np.asarray(B.sum(1), dtype=np.float64)
    dist = np.empty(n * (n - 1) // 2, dtype=np.float64)
    for start in range(0, n, 256):
        stop = min(start + 256, n)
        inter = (S[start:stop] @ ST).toarray()
        union = size[start:stop, None] + size[None, :] - inter
        with np.errstate(invalid="ignore", divide="ignore"):
            d = np.where(union > 0, 1 - inter / union, 1.0)
        for k, i in enumerate(range(start, stop)):
            off = n * i - i * (i + 1) // 2
            dist[off:off + n - i - 1] = d[k, i + 1:]
    return dist


def hier_linkage(Xn, method):
    return linkage(condensed_cosine_distance(Xn), method=method)


def hier_linkage_jac(B, method):
    return linkage(condensed_jaccard_distance(B), method=method)


def hier_partition(L, min_sim):
    """Cut at distance 1 - min_sim (min_cos for cosine linkages, min_jac for Jaccard ones)."""
    return fcluster(L, t=1 - min_sim, criterion="distance")


def resolve_min_models(min_frac, n_runs):
    """ceil(min_frac * n_runs), rounded first so 0.3 * 10 gives 3, not 4."""
    return math.ceil(round(min_frac * n_runs, 9))


def replicated(labels, models, min_models):
    """Clusters spanning >= min_models runs, ordered and named as cluster_metalv.py:121-126."""
    allc = pd.DataFrame({"model": models, "cluster": labels})
    stats = allc.groupby("cluster").agg(n_lvs=("model", "size"), n_models=("model", "nunique"))
    kept = stats[stats.n_models >= min_models].sort_values(["n_models", "n_lvs"], ascending=False).copy()
    kept["metaLV"] = [f"metaLV{i + 1}" for i in range(len(kept))]
    return kept, len(stats)


COMBINE = {"median": lambda X: np.median(X, axis=0)}


def build_set(pool, labels, min_models, combine="median"):
    """Returns Z (genes x k, float32), lvs (one row per metaLV) and members (one row per member LV)."""
    kept, n_clusters = replicated(labels, pool.models, min_models)
    f = COMBINE[combine]
    Z = np.empty((pool.X.shape[1], len(kept)), dtype=pool.X.dtype)
    members = []
    for j, (cid, row) in enumerate(kept.iterrows()):
        idx = np.flatnonzero(labels == cid)
        Z[:, j] = f(pool.X[idx])
        members += [{"metaLV": row.metaLV, "cluster": cid, "model": pool.models[i], "LV": pool.lvs[i]} for i in idx]
    lvs = kept.reset_index()[["metaLV", "cluster", "n_lvs", "n_models"]]
    return Z, lvs, pd.DataFrame(members, columns=["metaLV", "cluster", "model", "LV"]), n_clusters


def compare_to_metalv(Z, lvs, members, ref_z, ref_s):
    """Problems (empty if none) of a built set vs cluster_metalv.py's metaZ.csv / cluster_summary.csv."""
    problems = []
    if list(lvs.metaLV) != list(ref_z.columns):
        problems.append("metaLV names/order differ")
    elif Z.shape != ref_z.shape:
        problems.append(f"Z shape {Z.shape} vs {ref_z.shape}")
    else:
        # metaZ.csv holds float32 values in shortest repr, so the round trip is exact; an absolute
        # float64 tolerance would fail on large loadings (~1e-8 relative).
        ref32 = ref_z.to_numpy().astype(Z.dtype)
        if not (Z == ref32).all():
            problems.append(f"{int((Z != ref32).sum())} values differ in {Z.dtype}, "
                            f"max|dZ| {float(np.abs(Z.astype(np.float64) - ref_z.to_numpy()).max())}")
    ours = members.assign(key=lv_keys(members.model, members.LV)).groupby("metaLV", sort=False).key.apply(";".join)
    ref = ref_s.set_index("metaLV")
    if set(ours.index) != set(ref.index) or (ours.loc[ref.index] != ref.members).any():
        problems.append("member lists differ")
    if not (lvs.set_index("metaLV")[["n_lvs", "n_models"]].loc[ref.index] == ref[["n_lvs", "n_models"]]).all().all():
        problems.append("n_lvs/n_models differ")
    return problems


# ------------------------------------------------------------
# Metrics
# ------------------------------------------------------------

def top_genes(M, n):
    """Boolean rows x genes: the n largest positive values per row (fewer if fewer are positive),
    as run_ora.R selects its query."""
    B = np.zeros(M.shape, dtype=bool)
    for i, v in enumerate(M):
        pos = np.flatnonzero(np.isfinite(v) & (v > 0))
        B[i, pos[np.argsort(-v[pos], kind="stable")[:n]]] = True
    return B


def jaccard(A, B):
    """Pairwise Jaccard between boolean row sets A (a x g) and B (b x g)."""
    inter = A.astype(np.float32) @ B.T.astype(np.float32)
    union = A.sum(1)[:, None] + B.sum(1)[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(union > 0, inter / union, np.nan)


def jaccard_rows(A, B):
    """Jaccard of paired boolean rows A[i], B[i] (same NaN-on-empty rule as jaccard)."""
    inter = (A & B).sum(1)
    union = A.sum(1) + B.sum(1) - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(union > 0, inter / union, np.nan)


def quantiles(v, ps):
    """{'p50': ..., 'p99.9': ...} for probabilities ps (labels zero-padded: 0.01 -> 'p01')."""
    return {f"p{100 * p:02g}": float(x) for p, x in zip(ps, np.quantile(v, ps))}


def offdiag(C):
    return C[~np.eye(len(C), dtype=bool)]


def redundancy(Zn, B):
    """Per consensus LV: max cosine and max top-gene Jaccard to any other consensus LV."""
    if len(Zn) < 2:
        return np.full(len(Zn), np.nan), np.full(len(Zn), np.nan)
    C = Zn @ Zn.T
    J = jaccard(B, B)
    np.fill_diagonal(C, -np.inf)
    np.fill_diagonal(J, -np.inf)
    return C.max(1), J.max(1)


def best_cos(A, Bn):
    """Max cosine of each row of A (normalized) against the rows of Bn, blockwise."""
    out = np.empty(len(A), dtype=np.float32)
    for s in range(0, len(A), BLOCK):
        out[s:s + BLOCK] = (A[s:s + BLOCK] @ Bn.T).max(1)
    return out


def pool_best_hits(Xn, B, models, runs):
    """Best hit of each pool LV in each run. Returns n_pool x n_runs arrays: hits (cosine), arg
    (pool row of the hit, -1 for the LV's own run) and jac (top-gene Jaccard of the pair; B is
    top_genes of the pool). Own run: hits and jac NaN."""
    cols = {r: np.flatnonzero(models == r) for r in runs}
    hits = np.full((len(Xn), len(runs)), np.nan, dtype=np.float32)
    arg = np.full((len(Xn), len(runs)), -1, dtype=np.int64)
    jac = np.full((len(Xn), len(runs)), np.nan, dtype=np.float32)
    for s in range(0, len(Xn), BLOCK):
        e = min(s + BLOCK, len(Xn))
        S = Xn[s:e] @ Xn.T
        for j, r in enumerate(runs):
            sub = S[:, cols[r]]
            best = sub.argmax(1)
            hits[s:e, j] = sub[np.arange(e - s), best]
            arg[s:e, j] = a = cols[r][best]
            jac[s:e, j] = jaccard_rows(B[s:e], B[a])
    own = models[:, None] == np.asarray(runs)[None, :]
    hits[own], jac[own], arg[own] = np.nan, np.nan, -1
    return hits, arg, jac


def reciprocal(arg, models, runs):
    """n_pool x n_runs: LV i's best hit in run r has i as its own best hit in i's run."""
    own = pd.Index(runs).get_indexer(models)
    assert (own >= 0).all(), "a pool LV's run is not in runs"
    valid = arg >= 0
    i, j = np.nonzero(valid)
    out = np.zeros(arg.shape, dtype=bool)
    out[i, j] = arg[arg[i, j], own[i]] == i
    return out


def jaccard_null(B, models, runs, n_pairs, seed, q):
    """Best-hit top-gene Jaccard of each run-a LV against run b with gene order shuffled within
    every row (same set sizes), over n_pairs ordered run pairs drawn with seed. Returns the q
    quantile, summary quantiles and the pairs."""
    rng = np.random.default_rng(seed)
    pairs = [(a, b) for a in runs for b in runs if a != b]
    pairs = [pairs[i] for i in sorted(rng.choice(len(pairs), n_pairs, replace=False))]
    vals = []
    for a, b in pairs:
        Bb = rng.permuted(B[models == b], axis=1)
        J = jaccard(B[models == a], Bb)
        vals.append(np.where(np.isnan(J), -np.inf, J).max(1))
    v = np.concatenate(vals)
    v = v[np.isfinite(v)]
    quant = {**quantiles(v, [0.5, 0.99, 0.999]), "max": float(v.max()), "n": int(len(v))}
    return float(np.quantile(v, q)), quant, [[str(a), str(b)] for a, b in pairs]


def hit_stages(hits, recip, jac, cos, jac_thr):
    """Qualifying-hit masks (n_pool x n_runs), each adding one condition: cosine >= cos, then
    reciprocal, then top-gene Jaccard > jac_thr. NaN entries (own run) never qualify."""
    stages = {"cos": hits >= cos}
    stages["cos_reciprocal"] = stages["cos"] & recip
    stages["cos_reciprocal_jaccard"] = stages["cos_reciprocal"] & (jac > jac_thr)
    return stages


def reproducible(ok, min_runs):
    """LV whose qualifying hits (ok, n_pool x n_runs), plus its own run, span >= min_runs runs."""
    return ok.sum(1) + 1 >= min_runs


POOL_SCOPES = ["all", "reproducible", "reproducible_multi"]


def coverage(best, models, runs, repro, multi, thresholds):
    """Long table: share of each run's pool LVs, and of the POOL_SCOPES pool LVs (all, reproducible,
    reproducible and multi-gene), with best cosine to the set >= t. multi: pool mask of LVs that are
    not single-gene dominated."""
    rows = []
    pools = dict(zip(POOL_SCOPES, [np.ones(len(best), bool), repro, repro & multi]))
    scopes = [(r, models == r) for r in runs] + list(pools.items())
    for t in thresholds:
        hit = best >= t
        for scope, mask in scopes:
            n = int(mask.sum())
            rows.append({"scope": scope, "t": t, "n": n, "covered": int(hit[mask].sum()),
                         "share": hit[mask].mean() if n else np.nan})
    return pd.DataFrame(rows)


def top1_share(M):
    """Per row: largest squared value / squared norm (single-gene dominance)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.abs(M).max(1).astype(np.float64) ** 2 / np.einsum("ij,ij->i", M, M, dtype=np.float64)


def overlap(Bm, bz):
    """|top(member) ∩ top(consensus)| / |top(member)| for each row of Bm (NaN if empty)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return (Bm & bz).sum(1) / Bm.sum(1)


def cluster_metrics(cache, member_rows, codes, Z, top_n):
    """Per consensus LV metrics. member_rows index cache rows; codes are 0..k-1 metaLV codes."""
    Zn = l2norm(Z.T)
    Xm = cache.Xn[member_rows]
    Bm = top_genes(cache.X[member_rows], top_n)
    Bz = top_genes(Z.T, top_n)
    models = cache.models[member_rows]
    max_cos, max_jac = redundancy(Zn, Bz)
    share = top1_share(Z.T)
    rows = []
    for c in range(len(Zn)):
        idx = np.flatnonzero(codes == c)
        C = Xm[idx] @ Xm[idx].T
        J = jaccard(Bm[idx], Bm[idx])
        multi = len(idx) > 1
        rows.append({"n_lvs": len(idx), "n_models": len(set(models[idx])),
                     "mean_cos": offdiag(C).mean() if multi else np.nan,
                     "min_cos": offdiag(C).min() if multi else np.nan,
                     "cos_to_consensus": (Xm[idx] @ Zn[c]).mean(),
                     "jac_members": np.nanmean(offdiag(J)) if multi else np.nan,
                     "overlap_to_consensus": np.nanmean(overlap(Bm[idx], Bz[c])),
                     "top1_share": share[c],
                     "max_cos_other": max_cos[c], "max_jac_other": max_jac[c]})
    d = pd.DataFrame(rows)
    d.insert(2, "purity", d.n_models / d.n_lvs)
    return d


def summarize(per, cov, thresholds, single_gene):
    """One metrics row; single_gene is the top1_share above which a consensus LV is single-gene dominated."""
    row = {"n_lvs": len(per), "n_members": int(per.n_lvs.sum()),
           "median_n_models": per.n_models.median(), "mean_purity": per.purity.mean()}
    for c in ["mean_cos", "min_cos", "cos_to_consensus", "jac_members", "overlap_to_consensus",
              "max_cos_other", "max_jac_other"]:
        row[f"median_{c}"] = per[c].median()
    row["frac_single_gene"] = (per.top1_share > single_gene).mean()
    runs = cov[~cov.scope.isin(POOL_SCOPES)]
    by = cov.set_index(["scope", "t"])
    # recall_multi: recall of the reproducible LVs that are not single-gene dominated
    pooled = {"pool_coverage": "all", "recall": "reproducible", "recall_multi": "reproducible_multi"}
    for t in thresholds:
        row[f"frac_redundant_c{t}"] = (per.max_cos_other >= t).mean()
        row[f"coverage_c{t}"] = runs[runs.t == t].share.mean()
        for name, scope in pooled.items():
            row[f"{name}_c{t}"] = by.share[scope, t]
    for scope in ["reproducible", "reproducible_multi"]:
        row[f"n_{scope}"] = int(by.n[scope, thresholds[0]])
    return row

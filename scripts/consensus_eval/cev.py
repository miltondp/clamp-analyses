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

def hier_linkage(Xn, method):
    return linkage(condensed_cosine_distance(Xn), method=method)


def hier_partition(L, min_cos):
    return fcluster(L, t=1 - min_cos, criterion="distance")


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


def offdiag(C):
    return C[~np.eye(len(C), dtype=bool)]


def indicator(labels, k):
    H = np.zeros((len(labels), k), dtype=np.float32)
    H[np.arange(len(labels)), labels] = 1
    return H


def iq_silhouette(Xm, labels):
    """Icasso Iq per cluster and cosine silhouette per member, over the set's members only.

    labels are 0..k-1. Iq = mean intra-cluster similarity (pairs, no self) minus mean
    similarity of the cluster's members to all other members. Silhouette uses distance
    1 - cos, b = the nearest other cluster; members of singleton clusters get 0.
    """
    m, k = len(labels), labels.max() + 1
    H = indicator(labels, k)
    size = H.sum(0)
    intra = np.zeros(k)
    extra = np.zeros(k)
    sil = np.zeros(m)
    for s in range(0, m, BLOCK):
        e = min(s + BLOCK, m)
        S = (Xm[s:e] @ Xm.T).astype(np.float64)
        sums = S @ H
        rows = np.arange(e - s)
        lab = labels[s:e]
        self_sim = S[rows, rows + s]
        own = sums[rows, lab] - self_sim
        np.add.at(intra, lab, own)
        np.add.at(extra, lab, S.sum(1) - sums[rows, lab])
        n_own = size[lab] - 1
        with np.errstate(invalid="ignore", divide="ignore"):
            a = 1 - own / n_own
            dist_other = 1 - sums / size[None, :]
        dist_other[rows, lab] = np.inf
        b = dist_other.min(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            si = (b - a) / np.maximum(a, b)
        si[n_own == 0] = 0
        sil[s:e] = si
    with np.errstate(invalid="ignore", divide="ignore"):
        iq = intra / (size * (size - 1)) - extra / (size * (m - size))
    return iq, sil


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


def pool_best_hits(Xn, models, runs):
    """n_pool x n_runs: best cosine of each LV to any LV of each run; NaN for its own run."""
    cols = {r: np.flatnonzero(models == r) for r in runs}
    out = np.empty((len(Xn), len(runs)), dtype=np.float32)
    for s in range(0, len(Xn), BLOCK):
        S = Xn[s:s + BLOCK] @ Xn.T
        for j, r in enumerate(runs):
            out[s:s + BLOCK, j] = S[:, cols[r]].max(1) if len(cols[r]) else np.nan
    own = models[:, None] == np.asarray(runs)[None, :]
    out[own] = np.nan
    return out


def reproducible(best_hits, cos, min_runs):
    """LV whose hits at >= cos, plus its own run, span >= min_runs runs."""
    return (np.nan_to_num(best_hits, nan=-1) >= cos).sum(1) + 1 >= min_runs


def coverage(best, models, runs, repro, thresholds):
    """Long table: share of each run's pool LVs (and of all / reproducible pool LVs) with
    best cosine to the set >= t."""
    rows = []
    scopes = [(r, models == r) for r in runs] + [("all", np.ones(len(best), bool)), ("reproducible", repro)]
    for t in thresholds:
        hit = best >= t
        for scope, mask in scopes:
            n = int(mask.sum())
            rows.append({"scope": scope, "t": t, "n": n, "covered": int(hit[mask].sum()),
                         "share": hit[mask].mean() if n else np.nan})
    return pd.DataFrame(rows)


def cluster_metrics(cache, member_rows, codes, Z, top_n):
    """Per consensus LV metrics. member_rows index cache rows; codes are 0..k-1 metaLV codes."""
    Zn = l2norm(Z.T)
    Xm = cache.Xn[member_rows]
    Bm = top_genes(cache.X[member_rows], top_n)
    Bz = top_genes(Z.T, top_n)
    models = cache.models[member_rows]
    iq, sil = iq_silhouette(Xm, codes)
    max_cos, max_jac = redundancy(Zn, Bz)
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
                     "jac_to_consensus": np.nanmean(jaccard(Bm[idx], Bz[c:c + 1])),
                     "iq": iq[c], "silhouette": sil[idx].mean(),
                     "max_cos_other": max_cos[c], "max_jac_other": max_jac[c]})
    d = pd.DataFrame(rows)
    d.insert(2, "purity", d.n_models / d.n_lvs)
    return d, sil


def summarize(per, sil, cov, thresholds):
    row = {"n_lvs": len(per), "n_members": int(per.n_lvs.sum()),
           "median_n_models": per.n_models.median(), "mean_purity": per.purity.mean()}
    for c in ["mean_cos", "min_cos", "cos_to_consensus", "jac_members", "jac_to_consensus", "iq",
              "max_cos_other", "max_jac_other"]:
        row[f"median_{c}"] = per[c].median()
    row["mean_silhouette"] = sil.mean() if len(sil) else np.nan
    runs = cov[~cov.scope.isin(["all", "reproducible"])]
    for t in thresholds:
        row[f"frac_redundant_c{t}"] = (per.max_cos_other >= t).mean()
        row[f"coverage_c{t}"] = runs[runs.t == t].share.mean()
        row[f"pool_coverage_c{t}"] = cov[(cov.t == t) & (cov.scope == "all")].share.iloc[0]
        row[f"recall_c{t}"] = cov[(cov.t == t) & (cov.scope == "reproducible")].share.iloc[0]
    row["n_reproducible"] = int(cov[cov.scope == "reproducible"].n.iloc[0])
    return row

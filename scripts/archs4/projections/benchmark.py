#!/usr/bin/env python3
# Score projection recovery and GTEx clustering benchmarks.
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def provenance(output, inputs, parameters):
    import importlib.metadata
    versions = {'python': platform.python_version()}
    for name in ['numpy', 'pandas', 'pyreadr', 'scikit-learn', 'cuml']:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    Path(str(output) + '.json').write_text(json.dumps({
        'inputs': {str(p): digest(p) for p in inputs},
        'parameters': parameters, 'versions': versions}, indent=2))


def matrix(path):
    if str(path).endswith('.rds'):
        import pyreadr
        return next(iter(pyreadr.read_r(str(path)).values()))
    return pd.read_csv(path, index_col=0)


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    q = np.empty_like(p)
    q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    return np.minimum(q, 1)


def greedy(c):
    work = c.copy()
    out = np.full(c.shape[1], np.nan)
    lv = np.full(c.shape[1], -1, dtype=int)

    for _ in range(min(c.shape)):
        index = int(np.argmax(work.ravel(order='F')))
        row, col = index % c.shape[0], index // c.shape[0]
        if not np.isfinite(work[row, col]):
            break
        out[col], lv[col] = work[row, col], row
        work[row, :] = -np.inf
        work[:, col] = -np.inf
    return out, lv


def score(a):
    fac = matrix(a.matrix).T
    truth = pd.read_csv(a.truth, index_col=0).drop(columns=a.exclude)
    fac.index = fac.index.astype(str)
    truth.index = truth.index.astype(str)
    shared = fac.index.intersection(truth.index, sort=False)
    if len(shared) != len(fac) or len(shared) != len(truth):
        raise ValueError('Observed and permuted calculations must use the identical donor cohort')
    if fac.index.duplicated().any() or truth.index.duplicated().any():
        raise ValueError('Duplicate donor IDs')
    X = fac.loc[shared].to_numpy(dtype=np.float64)
    Y = truth.loc[shared].to_numpy(dtype=np.float64)
    if not np.isfinite(X).all() or not np.isfinite(Y).all():
        raise ValueError('Nonfinite values require explicit pairwise missingness handling')
    X -= X.mean(axis=0)
    norms = np.linalg.norm(X, axis=0)
    valid = norms > 0
    X = X[:, valid] / norms[valid]
    names = fac.columns[valid]
    Y -= Y.mean(axis=0)
    yn = np.linalg.norm(Y, axis=0)
    if np.any(yn == 0):
        raise ValueError('Constant cell-type proportions')
    Y /= yn
    cc = X.T @ Y
    if a.matching == 'max':
        assigned = cc.argmax(axis=0)
        observed = cc.max(axis=0)
    else:
        observed, assigned = greedy(cc)
    rng = np.random.default_rng(a.seed)
    null = np.full((a.permutations, Y.shape[1]), np.nan)
    for start in range(0, a.permutations, 100):
        stop = min(start+100, a.permutations)
        nb = stop-start
        yp = np.empty((len(Y), nb, Y.shape[1]))
        for j in range(Y.shape[1]):
            order = np.vstack([rng.permutation(len(Y)) for _ in range(nb)])
            yp[:, :, j] = Y[order, j].T
        correlations = (X.T @ yp.reshape(len(Y), -1)).reshape(X.shape[1], nb, Y.shape[1]).transpose(1, 0, 2)
        if a.matching == 'max':
            null[start:stop] = correlations.max(axis=1)
        else:
            for b in range(nb):
                null[start+b] = greedy(correlations[b])[0]
    rows = []
    for j, ct in enumerate(truth.columns):
        v = null[:, j]
        v = v[np.isfinite(v)]
        rows.append(dict(dataset=a.dataset, cell_type=ct, method=a.method,
                         assigned_LV=str(names[assigned[j]]) if assigned[j] >= 0 else '',
                         observed_r=observed[j], n_LVs_searched=X.shape[1],
                         null_mean_r=v.mean(), null_median_r=np.median(v), null_sd_r=v.std(ddof=1),
                         excess_r=observed[j]-np.median(v),
                         empirical_p=(1+np.count_nonzero(v >= observed[j]))/(len(v)+1) if np.isfinite(observed[j]) else np.nan,
                         n_null_assignments=len(v), n_donors=len(Y)))
    out = pd.DataFrame(rows)
    out['permutation_q'] = np.nan
    ok = out.empirical_p.notna()
    out.loc[ok, 'permutation_q'] = bh(out.loc[ok, 'empirical_p'])
    out.to_csv(a.output, sep='\t', index=False, float_format='%.12g')
    np.savez_compressed(str(a.output)+'.null.npz', null=null, cell_types=np.array(truth.columns, dtype=str))
    mean_null = np.nanmean(null, axis=1)
    pd.DataFrame({'permutation': np.arange(1, a.permutations+1), 'mean_r': mean_null}).to_csv(str(a.output)+'.mean_null.tsv.gz', sep='\t', index=False)
    provenance(a.output, [a.matrix, a.truth], dict(vars(a), func=None, n_targets=len(out), overall_observed=float(np.nanmean(observed)), overall_null_median=float(np.median(mean_null)), overall_excess=float(np.nanmean(observed)-np.median(mean_null)), overall_empirical_p=float((1+np.count_nonzero(mean_null >= np.nanmean(observed)))/(a.permutations+1))))


def embedding(a):
    from sklearn.preprocessing import StandardScaler
    fac = matrix(a.matrix).T
    fac.index = fac.index.astype(str).str.strip()
    meta = pd.read_csv(a.metadata, sep='\t', dtype=str, keep_default_na=False).set_index('SAMPID')
    shared = fac.index.intersection(meta.index, sort=False)
    labels = meta.loc[shared, 'SMTSD']
    counts = labels.value_counts()
    shared = shared[labels.isin(counts[counts >= a.min_samples].index)]
    labels = meta.loc[shared, 'SMTSD'].to_numpy(dtype=str)
    if len(np.unique(labels)) != a.n_clusters or len(shared) != len(fac):
        raise ValueError('GTEx projection cohort differs from the manuscript cohort')
    X = fac.loc[shared].to_numpy(dtype=np.float32)
    if not np.isfinite(X).all():
        raise ValueError('Nonfinite projected activities')
    X = StandardScaler().fit_transform(X)
    np.save(a.output, X)
    info = {'samples': shared.tolist(), 'labels': labels.tolist(), 'lv_names': fac.columns.astype(str).tolist(), 'n_clusters': a.n_clusters, 'matrix_sha256': digest(a.matrix)}
    Path(str(a.output)+'.metadata.json').write_text(json.dumps(info))
    provenance(a.output, [a.matrix, a.metadata], {'scaling': 'StandardScaler, float32', 'min_samples': a.min_samples, 'n_clusters': a.n_clusters, 'n_samples': len(shared)})


def subset(a):
    info = json.loads(Path(a.metadata).read_text())
    names = info['lv_names']
    if a.k > len(names):
        raise ValueError('Requested more LVs than available')
    if a.k == len(names):
        if a.rep != 0:
            raise ValueError('Full representation has a single subset (0)')
        ix = np.arange(len(names))
    else:
        rng = np.random.default_rng(np.random.SeedSequence([a.seed, a.k, a.rep]))
        ix = np.sort(rng.choice(len(names), a.k, replace=False))
    Path(a.output).write_text(json.dumps({'k': a.k, 'subset': a.rep, 'master_seed': a.seed, 'indices': ix.tolist(), 'lv_ids': [names[i] for i in ix]}, indent=2))


def cluster(a):
    import cupy as cp
    from cuml.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score
    info = json.loads(Path(a.metadata).read_text())
    X = np.load(a.embedding, mmap_mode='r')
    subset_info = json.loads(Path(a.subset).read_text()) if a.subset else {'k': X.shape[1], 'subset': 0, 'indices': list(range(X.shape[1]))}
    Xg = cp.asarray(np.asarray(X[:, subset_info['indices']], dtype=np.float32))
    rows = []
    for rep in range(a.repetitions):
        seed = a.base_seed + a.n_clusters * 1000 + rep
        km = KMeans(n_clusters=a.n_clusters, init='k-means++', random_state=seed, max_iter=300, tol=1e-4, verbose=0)
        labels = cp.asnumpy(km.fit_predict(Xg))
        rows.append({'k': subset_info['k'], 'subset': subset_info['subset'], 'rep': rep, 'seed': seed, 'ari': adjusted_rand_score(info['labels'], labels)})
    pd.DataFrame(rows).to_csv(a.output, index=False)
    provenance(a.output, [a.embedding, a.metadata] + ([a.subset] if a.subset else []), {'algorithm': 'cuml.KMeans', 'init': 'k-means++', 'max_iter': 300, 'tol': 1e-4, 'n_clusters': a.n_clusters, 'n_repetitions': a.repetitions, 'base_seed': a.base_seed})
    print(a.output, 'mean ARI', np.mean([r['ari'] for r in rows]), flush=True)


def aggregate(a):
    root = Path(a.root)
    cfg = json.loads(Path(a.settings).read_text())
    datasets, methods, seeds, ks = cfg['datasets'], cfg['methods'], cfg['null_seeds'], cfg['ks']
    pb = []
    for ds in datasets:
        for m in methods + ['projection']:
            t = pd.read_csv(root/'scores'/ds/m/'targets.tsv', sep='\t')
            pb.append(t[['dataset', 'cell_type', 'method', 'observed_r', 'null_median_r', 'excess_r']])
        null_rows = []
        for seed in seeds:
            t = pd.read_csv(root/'scores'/ds/f'null_{seed}'/'targets.tsv', sep='\t')
            null_rows.append(t)
        t = pd.concat(null_rows).groupby(['dataset','cell_type'], sort=False)[['observed_r','null_median_r','excess_r']].mean().reset_index()
        t['method'] = 'ARCHS4 projection null'
        pb.append(t)
    pb = pd.concat(pb, ignore_index=True)
    if pb.groupby('method').size().ne(47).any():
        raise ValueError('Expected 47 targets for each plotted condition')
    pb.to_csv(root/'pseudobulk_adjusted.tsv', sep='\t', index=False)
    published = pd.read_csv(a.local_ari)
    gt = []
    for m in ['CLAMPfull', 'CLAMPbase', 'RNA-Seq', 'PLIER']:
        d = published[published.method == m].copy()
        if len(d) != cfg['cluster_repetitions']:
            raise ValueError(f'{m} lacks the expected matched repetitions')
        d['rep'] = np.arange(len(d)); gt.append(d[['method','rep','ari']])
    full = pd.read_csv(root/'gtex/clusters/1728/0.csv')
    full['method'] = 'ARCHS4 projection'; gt.append(full[['method','rep','ari']])
    ns = []
    for seed in seeds:
        d = pd.read_csv(root/f'gtex/null_clusters/{seed}.csv');d['null_seed']=seed;ns.append(d)
    null_gt = pd.concat(ns).groupby('rep', sort=True).ari.mean().reset_index()
    null_gt['method'] = 'ARCHS4 projection null';gt.append(null_gt[['method','rep','ari']])
    pd.concat(gt,ignore_index=True).to_csv(root/'gtex_full_ari.tsv',sep='\t',index=False)
    sensitivity=[]
    for k in ks:
        for rep in ([0] if k==1728 else range(cfg['n_subsets'])):
            sensitivity.append(pd.read_csv(root/f'gtex/clusters/{k}/{rep}.csv'))
    all_reps=pd.concat(sensitivity,ignore_index=True)
    all_reps.to_csv(root/'gtex_sensitivity_repetitions.tsv',sep='\t',index=False)
    subset_means=all_reps.groupby(['k','subset']).ari.mean().reset_index()
    subset_means.to_csv(root/'gtex_sensitivity_subsets.tsv',sep='\t',index=False)
    summaries=[]
    for k,d in subset_means.groupby('k'):
        summaries.append({'k':k,'mean_ari':d.ari.mean(),'subset_p025':d.ari.quantile(.025),'subset_p975':d.ari.quantile(.975),'n_subsets':len(d)})
    pd.DataFrame(summaries).to_csv(root/'gtex_sensitivity_summary.tsv',sep='\t',index=False)
    inputs=[a.settings,a.local_ari]+list(root.glob('scores/*/*/targets.tsv'))+list(root.glob('gtex/clusters/*/*.csv'))+list(root.glob('gtex/null_clusters/*.csv'))
    provenance(root/'pseudobulk_adjusted.tsv',inputs,{'config':cfg,'gene_null_summary':'mean across gene-null replicates per target/cluster seed'})


def main():
    p=argparse.ArgumentParser(description='Score projection recovery and GTEx clustering benchmarks.'); sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('score');q.set_defaults(func=score)
    for name in ['matrix','truth','dataset','method','output']:q.add_argument('--'+name,required=True)
    q.add_argument('--matching',choices=['max','greedy'],required=True);q.add_argument('--exclude',nargs='*',default=[])
    q.add_argument('--permutations',type=int,default=10000);q.add_argument('--seed',type=int,default=123)
    q=sub.add_parser('embedding');q.set_defaults(func=embedding)
    for name in ['matrix','metadata','output']:q.add_argument('--'+name,required=True)
    q.add_argument('--min-samples',type=int,default=1);q.add_argument('--n-clusters',type=int,default=54)
    q=sub.add_parser('subset');q.set_defaults(func=subset)
    for name in ['metadata','output']:q.add_argument('--'+name,required=True)
    q.add_argument('--k',type=int,required=True);q.add_argument('--rep',type=int,required=True);q.add_argument('--seed',type=int,default=123)
    q=sub.add_parser('cluster');q.set_defaults(func=cluster)
    for name in ['embedding','metadata','output']:q.add_argument('--'+name,required=True)
    q.add_argument('--subset');q.add_argument('--n-clusters',type=int,default=54);q.add_argument('--repetitions',type=int,default=50);q.add_argument('--base-seed',type=int,default=10000)
    q=sub.add_parser('aggregate');q.set_defaults(func=aggregate)
    for name in ['root','settings','local-ari']:q.add_argument('--'+name,required=True)
    a=p.parse_args();Path(a.output).parent.mkdir(parents=True,exist_ok=True) if hasattr(a,'output') else None;a.func(a)

if __name__=='__main__':main()

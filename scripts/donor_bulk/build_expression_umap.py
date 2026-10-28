#!/usr/bin/env python3
# Build expression UMAPs for retained single cells.
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scanpy.preprocessing._highly_variable_genes import _get_disp_stats, _get_mean_bins
from sklearn.utils.extmath import randomized_svd
from umap.umap_ import find_ab_params, simplicial_set_embedding


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "pseudobulk"))
sys.path.insert(0, str(REPO_ROOT / "scripts" / "donor_bulk"))

import build_donor_bulk as db  # noqa: E402
from common import (  # noqa: E402
    h5_gene_names,
    iter_matrix_market,
    matrix_group,
    matrix_market_header,
    read_csr_rows,
    read_yaml,
    repo_path,
)


def path(value: str) -> Path:
    candidate = Path(value)
    return candidate if candidate.is_absolute() else REPO_ROOT / candidate


def cache_mtx(cfg: dict, keep: np.ndarray, n_genes: int, n_cells: int, cache: Path, chunk_nnz: int) -> None:
    partial = cache.with_suffix(".partial.h5")
    with h5py.File(partial, "w") as h5:
        data = h5.create_dataset("data", (0,), maxshape=(None,), dtype="float32", chunks=(1_000_000,))
        indices = h5.create_dataset("indices", (0,), maxshape=(None,), dtype="int32", chunks=(1_000_000,))
        pointers = [0]
        cursor = 0

        def append(row, col, value, end):
            nonlocal cursor
            block = sp.coo_matrix(
                (value, (col - cursor, row)), shape=(end - cursor, n_genes)
            ).tocsr()[keep[cursor:end]]
            old = len(data)
            new = old + block.nnz
            data.resize((new,))
            indices.resize((new,))
            data[old:new] = block.data
            indices[old:new] = block.indices
            pointers.extend((block.indptr[1:] + old).tolist())
            cursor = end

        carry = None
        for row, col, value in iter_matrix_market(repo_path(cfg["raw"]), chunk_nnz):
            if carry is not None:
                row, col, value = (np.concatenate([a, b]) for a, b in zip(carry, (row, col, value)))
            if np.any(col[1:] < col[:-1]):
                raise ValueError("Matrix Market input must be sorted by cell")
            end = int(col[-1])
            cut = int(np.searchsorted(col, end))
            append(row[:cut], col[:cut], value[:cut], end)
            carry = (row[cut:], col[cut:], value[cut:])
        append(*carry, len(keep))
        if len(pointers) != n_cells + 1:
            raise ValueError("retained-cell cache is incomplete")
        h5.create_dataset("indptr", data=np.asarray(pointers, dtype=np.int64))
    partial.rename(cache)


def retained_cells(args: argparse.Namespace, cfg: dict, work_dir: Path):
    manifest = db.load_manifest(cfg, args.manifest)
    if cfg["kind"] == "h5ad":
        with h5py.File(repo_path(cfg["raw"]), "r") as h5:
            ids, samples, labels, eligible = db.metadata_arrays_h5(h5, cfg)
            genes = h5_gene_names(h5, cfg.get("matrix", "X"))
    elif cfg["kind"] == "matrix_market":
        _, n_raw, _ = matrix_market_header(repo_path(cfg["raw"]))
        ids, samples, labels, eligible = db.metadata_arrays_mtx(cfg, n_raw)
        genes = pd.read_csv(repo_path(cfg["features"]), header=None, dtype=str)[0].to_numpy()
    else:
        raise ValueError(f"unsupported input kind: {cfg['kind']}")
    _, _, keep, _ = db.donor_codes(samples, manifest, eligible)
    index = np.flatnonzero(keep)
    expected = int(pd.read_csv(path(args.aggregation)).n_cells_retained.iloc[0])
    if len(index) != expected:
        raise ValueError(f"{args.dataset}: {len(index):,} retained cells, aggregation expects {expected:,}")
    cells = pd.DataFrame({
        "cell_index": index,
        "cell_id": ids[index],
        "sample": samples[index],
        "mapped_cell_type": labels[index],
    })
    if cfg["kind"] == "matrix_market":
        source = work_dir / "retained_counts.h5"
        if not source.exists():
            cache_mtx(cfg, keep, len(genes), len(index), source, args.chunk_nnz)
    else:
        source = repo_path(cfg["raw"])
    unique, aggregator = db.gene_aggregator(genes)
    aggregator = aggregator.astype(np.float32)

    def blocks():
        with h5py.File(source, "r") as h5:
            group = matrix_group(h5, cfg.get("matrix", "X")) if cfg["kind"] == "h5ad" else h5
            n_rows = len(keep) if cfg["kind"] == "h5ad" else len(index)
            for start in range(0, n_rows, args.chunk_cells):
                stop = min(start + args.chunk_cells, n_rows)
                block = read_csr_rows(group, start, stop, len(genes)).astype(np.float32)
                if cfg["kind"] == "h5ad":
                    block = block[keep[start:stop]]
                if block.shape[0]:
                    yield (block @ aggregator).tocsr()

    return cells, unique, blocks


def normalize(block: sp.csr_matrix, target_sum: float) -> sp.csr_matrix:
    totals = np.asarray(block.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise ValueError("retained cell with zero counts")
    return (sp.diags(target_sum / totals) @ block).tocsr()


def highly_variable(blocks, genes: np.ndarray, n_cells: int, args: argparse.Namespace):
    n_genes = len(genes)
    sums, squares, log_sums, log_squares = (np.zeros(n_genes) for _ in range(4))
    seen = 0
    for block in blocks():
        block = normalize(block, args.target_sum)
        sums += np.asarray(block.sum(axis=0)).ravel()
        squares += np.asarray(block.power(2).sum(axis=0)).ravel()
        np.log1p(block.data, out=block.data)
        log_sums += np.asarray(block.sum(axis=0)).ravel()
        log_squares += np.asarray(block.power(2).sum(axis=0)).ravel()
        seen += block.shape[0]
    if seen != n_cells:
        raise ValueError(f"{args.dataset}: streamed {seen:,} of {n_cells:,} retained cells")
    mean = sums / n_cells
    variance = (squares / n_cells - mean**2) * n_cells / (n_cells - 1)
    mean[mean == 0] = 1e-12
    dispersion = variance / mean
    dispersion[dispersion == 0] = np.nan
    stats = pd.DataFrame({"means": np.log1p(mean), "dispersions": np.log(dispersion)})
    stats["mean_bin"] = _get_mean_bins(stats.means, "seurat", 20)
    disp_stats = _get_disp_stats(stats, "seurat")
    normalized = ((stats.dispersions - disp_stats["avg"]) / disp_stats["dev"]).to_numpy()
    finite = np.isfinite(normalized)
    n_hvg = min(args.n_hvg, int(finite.sum()))
    if n_hvg < 2:
        raise ValueError(f"{args.dataset}: fewer than two highly variable genes")
    threshold = np.sort(normalized[finite])[::-1][n_hvg - 1]
    hvg = np.flatnonzero(np.nan_to_num(normalized, nan=-np.inf) >= threshold)
    log_mean = log_sums[hvg] / n_cells
    log_sd = np.sqrt(np.maximum((log_squares[hvg] / n_cells - log_mean**2) * n_cells / (n_cells - 1), 0))
    log_sd[log_sd == 0] = 1
    return hvg, log_mean, log_sd


def principal_components(blocks, n_cells: int, hvg, log_mean, log_sd, work_dir: Path, args) -> np.ndarray:
    scaled_path = work_dir / "scaled_expression.npy"
    scaled = np.lib.format.open_memmap(scaled_path, mode="w+", dtype=np.float32, shape=(n_cells, len(hvg)))
    offset = 0
    for block in blocks():
        block = normalize(block, args.target_sum)[:, hvg]
        np.log1p(block.data, out=block.data)
        dense = block.toarray()
        dense -= log_mean
        dense /= log_sd
        np.clip(dense, -10, 10, out=dense)
        scaled[offset:offset + len(dense)] = dense
        offset += len(dense)
    clipped_mean = np.zeros(len(hvg))
    for start in range(0, n_cells, args.chunk_cells):
        clipped_mean += scaled[start:start + args.chunk_cells].sum(axis=0, dtype=np.float64)
    clipped_mean /= n_cells
    for start in range(0, n_cells, args.chunk_cells):
        scaled[start:start + args.chunk_cells] -= clipped_mean
    scaled.flush()
    n_pcs = min(args.n_pcs, n_cells - 1, len(hvg) - 1)
    if n_pcs < 2:
        raise ValueError(f"{args.dataset}: expression matrix cannot support PCA")
    u, s, _ = randomized_svd(scaled, n_components=n_pcs, n_iter=args.pca_iterations, random_state=args.seed)
    del scaled
    scaled_path.unlink()
    return u * s


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--aggregation", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--chunk-cells", type=int, default=10_000)
    parser.add_argument("--chunk-nnz", type=int, default=5_000_000)
    parser.add_argument("--target-sum", type=float, default=10_000)
    parser.add_argument("--n-hvg", type=int, default=2_000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--pca-iterations", type=int, default=4)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.30)
    parser.add_argument("--n-epochs", type=int, default=200)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--threads", type=int, default=16)
    args = parser.parse_args()
    sc.settings.n_jobs = args.threads

    cfg = dict(read_yaml(path(args.config))["datasets"][args.dataset], dataset=args.dataset)
    work_dir = path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    cells, genes, blocks = retained_cells(args, cfg, work_dir)
    n_cells = len(cells)
    if cells["cell_index"].duplicated().any() or cells["cell_id"].duplicated().any():
        raise ValueError(f"{args.dataset}: retained cells are not unique")

    hvg, log_mean, log_sd = highly_variable(blocks, genes, n_cells, args)
    pcs = principal_components(blocks, n_cells, hvg, log_mean, log_sd, work_dir, args)
    n_pcs = pcs.shape[1]
    n_neighbors = min(args.n_neighbors, n_cells - 1)
    adata = ad.AnnData(X=sp.csr_matrix((n_cells, 0), dtype=np.float32))
    adata.obsm["X_pca"] = pcs
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep="X_pca", random_state=args.seed)
    a, b = find_ab_params(1.0, args.min_dist)
    coords, _ = simplicial_set_embedding(
        data=pcs, graph=adata.obsp["connectivities"].tocoo(), n_components=2,
        initial_alpha=1.0, a=a, b=b, gamma=1.0, negative_sample_rate=5,
        n_epochs=args.n_epochs, init="spectral", random_state=np.random.RandomState(args.seed),
        metric="euclidean", metric_kwds={}, densmap=False, densmap_kwds={},
        output_dens=False, parallel=True, verbose=True,
    )
    coords = np.asarray(coords, dtype=np.float64)
    if coords.shape != (n_cells, 2) or not np.isfinite(coords).all():
        raise ValueError(f"{args.dataset}: UMAP output is incomplete or non-finite")

    output = cells.copy()
    output["dataset"] = args.dataset
    output["umap1"] = coords[:, 0]
    output["umap2"] = coords[:, 1]
    output_path = path(args.output)
    summary_path = path(args.summary)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)
    pd.DataFrame(
        [
            {
                "dataset": args.dataset,
                "n_cells": n_cells,
                "n_genes_input": len(genes),
                "n_hvg": len(hvg),
                "n_pcs": n_pcs,
                "n_neighbors": n_neighbors,
                "target_sum": args.target_sum,
                "min_dist": args.min_dist,
                "n_epochs": args.n_epochs,
                "seed": args.seed,
                "pca_solver": f"randomized_svd_{args.pca_iterations}_power_iterations",
                "umap_optimizer": "umap_learn_parallel_sgd",
                "embedding_source": "all_retained_raw_expression_normalize_hvg_pca_neighbors_umap",
                "all_coordinates_finite": True,
            }
        ]
    ).to_csv(summary_path, index=False)
    shutil.rmtree(work_dir)
    print(f"[{args.dataset}] wrote expression UMAP for {n_cells:,} cells")


if __name__ == "__main__":
    main()

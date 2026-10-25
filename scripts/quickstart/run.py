#!/usr/bin/env python3
"""Exercise production pseudobulk processing and a notebook on synthetic counts."""
import argparse
import json
from pathlib import Path
import sys
import time

import h5py
import nbformat
from nbclient import NotebookClient
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'pseudobulk'))
import build_pseudobulk as build


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    started = time.monotonic()
    root = Path(args.output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    raw = root / 'synthetic.h5ad'
    matrix = sparse.csr_matrix([[1, 1, 2], [3, 1, 0], [5, 0, 5], [0, 2, 2]])
    with h5py.File(raw, 'w') as handle:
        x = handle.create_group('X')
        x.attrs['shape'] = matrix.shape
        x.attrs['encoding-type'] = 'csr_matrix'
        for name in ['data', 'indices', 'indptr']:
            x.create_dataset(name, data=getattr(matrix, name))
        obs = handle.create_group('obs')
        obs.create_dataset('donor', data=np.array(['b', 'a', 'b', 'a'], dtype='S'))
        obs.create_dataset('label', data=np.array(['T', 'B', 'excluded', 'T'], dtype='S'))
        handle.create_group('var').create_dataset('gene_name', data=np.array(['A', 'A', 'B'], dtype='S'))
    cfg = dict(raw=str(raw), matrix='X', sample_col='donor', truth_v0_col='label',
               exclude={'label': ['excluded']}, min_cells=1, chunk_cells=2)
    result = build.aggregate_h5ad('synthetic', cfg, pd.DataFrame({'sample': ['b', 'a']}))
    build.validate('synthetic', result)
    np.testing.assert_allclose(result.bulk.to_numpy(), [[250000, 500000], [375000, 250000]], rtol=0, atol=1e-8)
    assert result.info.nCells.tolist() == [1, 2]
    np.testing.assert_allclose(result.truth_v0.loc[['b', 'a'], ['B', 'T']], [[0, 1], [.5, .5]])
    result.bulk.to_csv(root / 'expression.csv')
    result.truth_v0.to_csv(root / 'fractions.csv')
    notebook = nbformat.v4.new_notebook(cells=[
        nbformat.v4.new_markdown_cell('# Synthetic pseudobulk example\nFour cells, two donors and two genes after duplicate-gene aggregation.'),
        nbformat.v4.new_code_cell("import pandas as pd\nimport numpy as np\nx = pd.read_csv('expression.csv', index_col=0)\nf = pd.read_csv('fractions.csv', index_col=0)\nassert x.index.tolist() == ['b', 'a']\nnp.testing.assert_allclose(x.to_numpy(), [[250000, 500000], [375000, 250000]])\nnp.testing.assert_allclose(f.sum(axis=1), 1)\nprint('PASS: expected expression and cell-type fractions')\nx")
    ])
    NotebookClient(notebook, timeout=60, kernel_name='python3', resources={'metadata': {'path': str(root)}}).execute()
    nbformat.write(notebook, root / 'quickstart.ipynb')
    report = {'status': 'PASS', 'cells': 4, 'retained_cells': 3, 'donors': 2, 'genes': 2,
              'elapsed_seconds': round(time.monotonic() - started, 3)}
    (root / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()

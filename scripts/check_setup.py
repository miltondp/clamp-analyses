#!/usr/bin/env python3
"""Check installed tools without downloading data or running analyses."""
import importlib.metadata
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    failures = []
    def check(label, command):
        result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        print(f"{'PASS' if result.returncode == 0 else 'FAIL'} {label}")
        if result.returncode:
            print((result.stderr or result.stdout).strip())
            failures.append(label)
    if not shutil.which('conda'):
        raise SystemExit('Conda is not on PATH. Source its etc/profile.d/conda.sh first.')
    try:
        version = importlib.metadata.version('snakemake')
        plugin = importlib.metadata.version('snakemake-executor-plugin-slurm')
        if not version.startswith('8.'):
            failures.append('Snakemake 8 required')
        print(f'Snakemake {version}; SLURM plugin {plugin}')
    except importlib.metadata.PackageNotFoundError:
        raise SystemExit('Activate the snakemake environment before running this check.')
    check('CLAMP revision', ['conda', 'run', '-n', 'clamp-analyses', 'Rscript', 'scripts/install_clamp.R', '--check'])
    check('comparison method revisions', ['conda', 'run', '-n', 'clamp-analyses', 'Rscript', 'scripts/install_analysis_packages.R', '--check'])
    check('analysis Python/notebook imports', ['conda', 'run', '-n', 'clamp-analyses', 'python', '-c',
          'import yaml,numpy,pandas,scipy,h5py,anndata,scanpy,mofaflex,nbformat,nbclient,nbconvert; from jupyter_client.kernelspec import KernelSpecManager; assert \"ir\" in KernelSpecManager().find_kernel_specs(), \"R kernel ir is missing\"'])
    check('R notebook and plotting packages', ['conda', 'run', '-n', 'clamp-analyses', 'Rscript', '-e',
          'p<-c("IRkernel","bigstatsr","CoGAPS","MOFA2","cowplot","patchwork","yaml"); stopifnot(all(vapply(p,requireNamespace,logical(1),quietly=TRUE)))'])
    check('GPU Python packages (driver checked at execution)', ['conda', 'run', '-n', 'gpu-kmeans', 'python', '-c',
          'import cuml,cupy,pyreadr,rpy2'])
    check('workflow parsing', [sys.executable, '-m', 'snakemake', '--snakefile', 'workflow/Snakefile', '--list-target-rules'])
    print('Inkscape: ' + (shutil.which('inkscape') or 'missing (required for Figure 4)'))
    print('Data availability and sufficient RAM/GPU capacity require a separate target dry-run and hardware check.')
    return bool(failures)


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""Adapt the pinned phenoplier CLI to the repository's summary/store targets."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def run(*args):
    subprocess.run([str(arg) for arg in args], check=True)


@contextmanager
def lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def publish(source, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + f'.tmp.{os.getpid()}')
    try:
        shutil.copyfile(source, tmp)
        tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['gls', 'store'])
    p.add_argument('--rds', '--model-input', dest='model', required=True, type=Path)
    p.add_argument('--name', required=True)
    p.add_argument('--summary-out', type=Path)
    p.add_argument('--store-out', type=Path)
    p.add_argument('--workspace', type=Path,
                   default=Path(os.environ.get('PHENOPLIER_HOME', '~/phenoplier')).expanduser())
    p.add_argument('--conda-env')  # activated by the shell entrypoint
    p.add_argument('--require-version', default='0.7.3')
    p.add_argument('--namespace', default='clamp')
    p.add_argument('--cohort', default='phenomexcan_rapid_gwas')
    p.add_argument('--lv-percentile', type=float, default=0.01)
    p.add_argument('--trait-filter', choices=['all', 'phenomexcan-curated'], default='phenomexcan-curated')
    p.add_argument('--fit-trait-filter', choices=['all', 'phenomexcan-curated'])
    p.add_argument('--expected-phenotypes', type=int, default=2824)
    p.add_argument('--expected-genes', type=int)
    p.add_argument('--expected-lvs', type=int)
    p.add_argument('--executor', choices=['local', 'slurm'], default='local')
    p.add_argument('--cluster', default='')
    p.add_argument('--n-jobs', type=int, default=14)
    p.add_argument('--mem-mb', type=int, default=100000)
    p.add_argument('--workers-step6', type=int, default=14)
    p.add_argument('--workers-step7', type=int, default=14)
    p.add_argument('--blas-step6', type=int, choices=[1], default=1)
    p.add_argument('--blas-step7', type=int, choices=[1], default=1)
    return p


def main():
    p = parser()
    a = p.parse_args()
    if a.mode == 'gls' and a.summary_out is None:
        p.error('--summary-out is required for gls')
    if a.mode == 'store' and a.store_out is None:
        p.error('--store-out is required for store')
    if min(a.n_jobs, a.workers_step6, a.workers_step7) < 1:
        p.error('CPU budgets must be positive')
    fit_filter = a.fit_trait_filter or a.trait_filter
    if fit_filter != 'all' and a.trait_filter == 'all':
        p.error('an all-trait summary cannot be derived from a curated fit')
    a.model = a.model.resolve(strict=True)
    a.workspace = a.workspace.resolve()
    a.summary_out = a.summary_out.resolve() if a.summary_out else None
    a.store_out = a.store_out.resolve() if a.store_out else None
    if not (a.workspace / 'config.toml').is_file():
        p.error(f'No workspace at {a.workspace}; see scripts/phenoplier/README.md')
    os.environ.update(PHENOPLIER_HOME=str(a.workspace), PHENOPLIER_ROOT_DIR=str(a.workspace))
    os.chdir(a.workspace)  # avoid importing a checkout from the working directory
    from phenoplier import __version__
    if __version__ != a.require_version:
        p.error(f'phenoplier {__version__} installed; require {a.require_version}')
    if a.expected_genes or a.expected_lvs:
        import pandas as pd
        if a.model.suffix in {'.pkl', '.pickle'}:
            matrix = pd.read_pickle(a.model)
        else:
            matrix = pd.read_csv(a.model, index_col=0,
                                 sep='\t' if a.model.suffix == '.tsv' else ',')
        if a.expected_genes and matrix.shape[0] != a.expected_genes:
            p.error(f'Model genes {matrix.shape[0]} != inventory {a.expected_genes}')
        if a.expected_lvs and matrix.shape[1] != a.expected_lvs:
            p.error(f'Model LVs {matrix.shape[1]} != inventory {a.expected_lvs}')
        del matrix
    digest = hashlib.sha256()
    with a.model.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    contract = dict(source_sha256=digest.hexdigest(), version=a.require_version,
                    cohort=a.cohort, lv_percentile=a.lv_percentile, trait_filter=fit_filter,
                    reference_panel='GTEX_V8', eqtl_model='MASHR')
    fingerprint = hashlib.sha256(json.dumps(contract, sort_keys=True).encode()).hexdigest()
    # New inputs/settings get independent projects. Never remove registrations,
    # overwrite someone else's config, or delete locks/sentinels to force a resume.
    key = f'{a.namespace}/{a.name}/{fingerprint[:20]}'
    project = a.workspace / 'projects' / key
    cfg = project / 'pipeline_config.yaml'
    results = project / 'results/gls/phenoplier'
    suffix = '' if a.trait_filter == 'all' else f'-{a.trait_filter}'
    summary = results / f'gls-summary-{a.cohort}{suffix}.tsv.gz'
    pickle = summary.with_name(summary.name.replace('.tsv.gz', '.pkl.gz'))
    pp = [sys.executable, '-m', 'phenoplier']
    with lock(a.workspace / '.clamp-locks' / f'{fingerprint}-{a.name}.lock'):
        if a.mode == 'gls':
            if not cfg.is_file():
                # Registration is serialized; expensive model runs are independent.
                with lock(a.workspace / '.clamp_registry.lock'):
                    run(*pp, 'model', 'register', a.model,
                        '--namespace', f'{a.namespace}/{a.name}', '--name', fingerprint[:20])
                args = [*pp, 'run', 'gls', '--model', key, '--project', key,
                        '--cohort', a.cohort, '--reference-panel', 'GTEX_V8', '--eqtl-model', 'MASHR',
                        '--lv-percentile', a.lv_percentile, '--trait-filter', fit_filter,
                        '--executor', a.executor, '--cores', a.n_jobs,
                        '--cores-step6', min(a.n_jobs, a.workers_step6),
                        '--cores-step7', min(a.n_jobs, a.workers_step7),
                        '--shell-prefix', 'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
                        '--init-only']
                if a.cluster:
                    args += ['--cluster', a.cluster]
                run(*args)
                (project / 'clamp_source.json').write_text(json.dumps(contract, indent=2) + '\n')
            # Bound the nested scheduler by the outer job's allocation. Preserve
            # per-chromosome requests (especially chr1); cap only total resources.
            import yaml
            profile = project / 'profiles' / f'gls_{a.executor}' / 'config.yaml'
            settings = yaml.safe_load(profile.read_text()) if profile.exists() else {}
            settings = settings or {}
            if a.executor == 'local':
                settings['cores'] = a.n_jobs
                resources = settings.setdefault('resources', {})
                if isinstance(resources, list):
                    settings['resources'] = [r for r in resources if not str(r).startswith('mem_mb=')] + [f'mem_mb={a.mem_mb}']
                else:
                    resources['mem_mb'] = a.mem_mb
            else:
                settings['jobs'] = max(1, a.n_jobs // 2)
            profile.parent.mkdir(parents=True, exist_ok=True)
            profile.write_text(yaml.safe_dump(settings, sort_keys=False))
            run(*pp, 'run', 'gls', '--config-file', cfg, '--executor', a.executor)
            for attempt in range(3):
                try:
                    run(*pp, 'run', 'status', '-p', project, '--workflow', 'gls', '--verify')
                    break
                except subprocess.CalledProcessError:
                    if attempt == 2:
                        raise
                    time.sleep(10)
            if fit_filter == 'all' and a.trait_filter != 'all':
                all_gate = [sys.executable, Path(__file__).with_name('verify_summary.py'),
                            results / f'gls-summary-{a.cohort}.pkl.gz', 4049]
                if a.expected_lvs:
                    all_gate.append(a.expected_lvs)
                run(*all_gate)
            if fit_filter != a.trait_filter:
                run(*pp, 'results', 'summarize', '--project-dir', project, '--cohort', a.cohort,
                    '--trait-filter', a.trait_filter, '--no-export-rds', '--no-dry-run')
            gate = [sys.executable, Path(__file__).with_name('verify_summary.py'),
                    pickle, a.expected_phenotypes]
            if a.expected_lvs:
                gate.append(a.expected_lvs)
            run(*gate)
            publish(pickle, a.summary_out.with_name(a.summary_out.name.replace('.tsv.gz', '.pkl.gz')))
            if fit_filter == 'all' and a.trait_filter != 'all':
                for ext in ('tsv.gz', 'pkl.gz'):
                    publish(results / f'gls-summary-{a.cohort}.{ext}',
                            a.summary_out.parent / f'all-traits.{ext}')
            excluded = results / f'trait-filter-excluded-{a.trait_filter}.tsv'
            if excluded.is_file():
                publish(excluded, a.summary_out.with_name(a.summary_out.name.replace('.tsv.gz', '_excluded.tsv')))
            publish(summary, a.summary_out)  # completion target is published last
        else:
            if not pickle.is_file():
                p.error(f'Missing exact summary {pickle}; run the GLS target first')
            a.store_out.parent.mkdir(parents=True, exist_ok=True)
            tmp = a.store_out.with_name(a.store_out.name + f'.tmp.{os.getpid()}')
            try:
                run(*pp, 'store', 'build', '--gls-dir', results / a.cohort,
                    '--cohort', a.cohort, '--trait-filter', a.trait_filter,
                    '--gls-summary', pickle, '--output', tmp, '--model', key, '--method-name', 'CLAMP')
                import h5py
                with h5py.File(tmp, 'r'):
                    pass
                tmp.replace(a.store_out)
            finally:
                tmp.unlink(missing_ok=True)


if __name__ == '__main__':
    main()

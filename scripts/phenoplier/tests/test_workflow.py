"""Contract tests: no reference data or GLS computation required."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('adapter', SCRIPT / 'workflow.py')
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


class WorkflowTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ws = self.root / 'workspace'
        self.ws.mkdir()
        (self.ws / 'config.toml').touch()
        self.model = self.root / 'Z.csv'
        pd.DataFrame({'LV1': [1., 2.]}, index=['A', 'B']).to_csv(self.model)
        self.out = self.root / 'out/curated.tsv.gz'
        self.calls = []
        self.project = None
        self.original_cwd = Path.cwd()
        self.addCleanup(os.chdir, self.original_cwd)
        self.df = pd.DataFrame({'phenotype': ['p1', 'p2'], 'lv': ['LV1', 'LV1'],
                                'pvalue': [.1, .2], 'fdr': [.2, .2]})

    def fake_run(self, *args):
        args = list(map(str, args))
        self.calls.append(args)
        if args[1] == str(SCRIPT / 'verify_summary.py'):
            return subprocess.run(args, check=True, capture_output=True)
        if '--init-only' in args:
            self.project = self.ws / 'projects' / args[args.index('--project') + 1]
            self.project.mkdir(parents=True)
            (self.project / 'pipeline_config.yaml').touch()
        if '--config-file' in args or 'summarize' in args:
            results = self.project / 'results/gls/phenoplier'
            results.mkdir(parents=True, exist_ok=True)
            for suffix in ['', '-phenomexcan-curated']:
                stem = results / f'gls-summary-phenomexcan_rapid_gwas{suffix}'
                frame = self.df if suffix else pd.DataFrame({
                    'phenotype': [f'p{i}' for i in range(4049)], 'lv': 'LV1',
                    'pvalue': .1, 'fdr': .2})
                frame.to_csv(str(stem) + '.tsv.gz', sep='\t', index=False)
                frame.to_pickle(str(stem) + '.pkl.gz')
        if 'build' in args:
            import h5py
            with h5py.File(args[args.index('--output') + 1], 'w'):
                pass

    def invoke(self, mode='gls', extra=()):
        argv = ['workflow.py', mode, '--rds', str(self.model), '--name', 'example',
                '--workspace', str(self.ws), '--expected-phenotypes', '2',
                '--expected-genes', '2', '--expected-lvs', '1',
                '--summary-out', str(self.out), '--store-out', str(self.root / 'store.h5'),
                *extra]
        with patch.object(sys, 'argv', argv), patch.object(adapter, 'run', self.fake_run), \
             patch.dict(sys.modules, phenoplier=types.SimpleNamespace(__version__='0.7.3')), \
             patch.dict(os.environ):
            adapter.main()

    def test_historical_summarizes_without_second_fit_and_store_uses_exact_family(self):
        self.invoke(extra=['--fit-trait-filter', 'all'])
        self.assertTrue(self.out.is_file())
        self.assertTrue((self.out.parent / 'all-traits.pkl.gz').is_file())
        self.assertEqual(sum('--config-file' in call for call in self.calls), 1)
        self.assertEqual(sum('summarize' in call for call in self.calls), 1)
        self.invoke('store', ['--fit-trait-filter', 'all'])
        store = next(call for call in self.calls if 'build' in call)
        self.assertIn('phenomexcan-curated', store)
        self.assertTrue(store[store.index('--gls-summary') + 1].endswith('-phenomexcan-curated.pkl.gz'))

    def test_changed_source_or_filter_gets_new_project_and_retries_keep_config(self):
        self.invoke()
        first = self.project
        self.invoke()
        self.assertEqual(sum('--init-only' in call for call in self.calls), 1)
        self.invoke(extra=['--fit-trait-filter', 'all'])
        self.assertNotEqual(first, self.project)
        second = self.project
        self.model.write_text(',LV1\nA,3\nB,4\n')
        self.invoke(extra=['--fit-trait-filter', 'all'])
        self.assertNotEqual(second, self.project)
        self.assertTrue(first.exists())

    def test_nested_scheduler_respects_outer_allocation(self):
        import yaml
        self.invoke(extra=['--n-jobs', '6', '--mem-mb', '64000'])
        profile = yaml.safe_load((self.project / 'profiles/gls_local/config.yaml').read_text())
        self.assertEqual(profile['cores'], 6)
        self.assertEqual(profile['resources']['mem_mb'], 64000)
        init = next(call for call in self.calls if '--init-only' in call)
        self.assertEqual(init[init.index('--cores-step7') + 1], '6')

    def test_bad_summary_is_not_published(self):
        self.df.loc[0, 'fdr'] = float('nan')
        with self.assertRaises(subprocess.CalledProcessError):
            self.invoke()
        self.assertFalse(self.out.exists())

    def test_wrong_model_dimensions_fail_before_registration(self):
        with self.assertRaises(SystemExit):
            self.invoke(extra=['--expected-lvs', '2'])
        self.assertFalse(self.calls)

    def test_gate_rejects_missing_grid_row_and_duplicate_pairs(self):
        frame = pd.concat([self.df, self.df.assign(lv='LV2')], ignore_index=True)
        for invalid in [frame.iloc[:-1], pd.concat([frame, frame.iloc[:1]])]:
            path = self.root / 'invalid.pkl.gz'
            invalid.to_pickle(path)
            result = subprocess.run([sys.executable, str(SCRIPT / 'verify_summary.py'), str(path), '2'], capture_output=True)
            self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()

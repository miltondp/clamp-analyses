"""Recorded timings are restored unchanged and never replace different results."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('recorded_timings', Path(__file__).resolve().parents[1] / 'scripts/inputs/restore_recorded_timings.py')
recorded = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(recorded)


class RecordedTimingTests(unittest.TestCase):
    def test_restore_both_reports_and_reuse(self):
        with tempfile.TemporaryDirectory() as temporary:
            for dataset in ['pseudobulk', 'gtex']:
                output = Path(temporary) / dataset
                recorded.restore(dataset, output)
                for source in (recorded.ROOT / 'references/runtime' / dataset).glob('*.csv'):
                    self.assertEqual((output / source.name).read_bytes(), source.read_bytes())
                recorded.restore(dataset, output)
                self.assertEqual(len(list(output.glob('*.csv'))), 3)

    def test_do_not_replace_new_measurements(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            local = output / 'runtime_long.csv'
            local.write_bytes(b'new measurements')
            with self.assertRaises(ValueError):
                recorded.restore('gtex', output)
            self.assertEqual(local.read_bytes(), b'new measurements')


if __name__ == '__main__':
    unittest.main()

import importlib.util
from pathlib import Path
import tempfile
import unittest
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/inputs"))

SPEC = importlib.util.spec_from_file_location('spredixcan_inputs', Path(__file__).resolve().parents[1] / 'scripts/inputs/link_spredixcan_inputs.py')
staging = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(staging)


class SpredixcanStagingTests(unittest.TestCase):
    def test_inventory_rejects_missing_or_changed_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            name = 'spredixcan-mashr-zscores-tissue.pkl'
            expected = {name: 3}
            with self.assertRaises(ValueError):
                staging.validate(root, expected)
            (root / name).write_bytes(b'abc')
            staging.validate(root, expected)
            (root / name).write_bytes(b'changed')
            with self.assertRaises(ValueError):
                staging.validate(root, expected)


if __name__ == '__main__':
    unittest.main()

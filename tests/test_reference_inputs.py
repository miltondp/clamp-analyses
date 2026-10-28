"""Small offline checks for pinned reference publication and integrity."""
import gzip
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('reference_inputs', Path(__file__).resolve().parents[1] / 'scripts/inputs/reference_inputs.py')
references = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(references)


class ReferenceInputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.payload = b'gene\tsample\nA\t7\n'
        self.source = self.root / 'source.txt'
        self.source.write_bytes(self.payload)
        self.asset = {'path': 'data/input.txt', 'bundled': False, 'url': self.source.as_uri(),
                      'sha256': hashlib.sha256(self.payload).hexdigest()}

    def test_download_and_reuse(self):
        references.fetch('tiny', self.asset, self.root)
        output = self.root / self.asset['path']
        before = output.stat().st_mtime_ns
        self.source.unlink()  # Reuse must not need a second request.
        references.fetch('tiny', self.asset, self.root)
        self.assertEqual(output.read_bytes(), self.payload)
        self.assertEqual(output.stat().st_mtime_ns, before)

    def test_gzip_is_verified_after_decompression(self):
        self.source.write_bytes(gzip.compress(self.payload))
        references.fetch('tiny', dict(self.asset, compression='gzip'), self.root)
        self.assertEqual((self.root / self.asset['path']).read_bytes(), self.payload)

    def test_changed_download_is_not_published(self):
        self.source.write_bytes(b'changed upstream')
        with self.assertRaises(ValueError):
            references.fetch('tiny', self.asset, self.root)
        self.assertFalse((self.root / self.asset['path']).exists())
        self.assertFalse(list((self.root / 'data').glob('.reference-*')))

    def test_existing_file_is_never_replaced(self):
        output = self.root / self.asset['path']
        output.parent.mkdir()
        output.write_bytes(b'local scientific input')
        with self.assertRaises(ValueError):
            references.fetch('tiny', self.asset, self.root)
        self.assertEqual(output.read_bytes(), b'local scientific input')

    def test_check_and_offline_do_not_download(self):
        with self.assertRaises(FileNotFoundError):
            references.fetch('tiny', self.asset, self.root, check=True)
        with self.assertRaises(RuntimeError):
            references.fetch('tiny', self.asset, self.root, offline=True)
        self.assertFalse((self.root / self.asset['path']).exists())

    def test_configured_destination_is_supported(self):
        output = self.root / 'custom' / 'reference.txt'
        references.fetch('tiny', self.asset, self.root, destination_override=output)
        self.assertEqual(output.read_bytes(), self.payload)
        self.assertFalse((self.root / self.asset['path']).exists())


class WorkspaceInputTests(unittest.TestCase):
    setUp = ReferenceInputTests.setUp
    def test_matching_workspace_input_is_linked(self):
        asset = dict(self.asset, source_candidates=['source.txt'])
        references.fetch('tiny', asset, self.root, source_root=self.root)
        output = self.root / asset['path']
        self.assertTrue(output.is_symlink())
        self.assertEqual(output.read_bytes(), self.payload)
        self.assertEqual(self.source.read_bytes(), self.payload)

    def test_wrong_workspace_input_is_rejected(self):
        asset = dict(self.asset, source_candidates=['source.txt'])
        self.source.write_bytes(b'other cohort')
        with self.assertRaises(ValueError):
            references.fetch('tiny', asset, self.root, source_root=self.root)
        self.assertFalse((self.root / asset['path']).exists())


if __name__ == '__main__':
    unittest.main()

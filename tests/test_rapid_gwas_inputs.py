"""Exercise workspace staging without downloading the 10 GB reference ZIP."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts/inputs'
sys.path.insert(0, str(SCRIPTS))
import rapid_gwas_inputs as inputs
from spredixcan_inputs import checksum


class RapidGwasTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.workspace = self.root / 'workspace'
        self.workspace.mkdir()
        (self.workspace / 'config.toml').write_text('')
        self.archive = self.root / 'rapid_gwas.zip'
        self.name = 'rapid_gwas/smultixcan_123_ccn30.tsv.gz'
        with zipfile.ZipFile(self.archive, 'w') as handle:
            handle.writestr(self.name, b'exact published bytes')
        self.asset = {'url': self.archive.as_uri(), 'bytes': self.archive.stat().st_size,
                      'md5': checksum(self.archive, 'md5'), 'files': 1}

    def test_download_and_resume_preserve_existing_input(self):
        cached = self.root / 'cache/rapid_gwas.zip'
        inputs.stage(self.workspace, cached, self.asset)
        output = self.workspace / 'data/smultixcan' / self.name
        before = output.stat().st_mtime_ns
        self.archive.unlink()
        inputs.stage(self.workspace, cached, self.asset)
        self.assertEqual(output.read_bytes(), b'exact published bytes')
        self.assertEqual(before, output.stat().st_mtime_ns)

    def test_changed_existing_input_is_rejected(self):
        inputs.stage(self.workspace, self.archive, self.asset)
        output = self.workspace / 'data/smultixcan' / self.name
        output.write_bytes(b'local edited bytes')
        with self.assertRaisesRegex(ValueError, 'differs from published'):
            inputs.stage(self.workspace, self.archive, self.asset)
        self.assertEqual(output.read_bytes(), b'local edited bytes')

    def test_unsafe_archive_member_is_rejected(self):
        with zipfile.ZipFile(self.archive, 'w') as handle:
            handle.writestr('../smultixcan_escape.tsv.gz', b'escape')
        asset = dict(self.asset, bytes=self.archive.stat().st_size, md5=checksum(self.archive, 'md5'))
        with self.assertRaisesRegex(ValueError, 'Unexpected archive path'):
            inputs.stage(self.workspace, self.archive, asset)
        self.assertFalse((self.root / 'smultixcan_escape.tsv.gz').exists())


if __name__ == '__main__':
    unittest.main()

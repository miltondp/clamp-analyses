"""Exercise Perturb-seq input staging without downloading the 7.7 GB Ota et al. ZIP."""
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

SCRIPTS = Path(__file__).resolve().parents[1] / 'scripts/inputs'
sys.path.insert(0, str(SCRIPTS))
import perturbseq_inputs as inputs
from spredixcan_inputs import checksum

MEMBER = 'data/Perturbseq/geneLevel/K562GW/limma_logFC_sum.txt'
PAYLOAD = b'GENE\tGATA1\nENSG00000000003\t0.5\n'
DESTINATION = 'data/perturbseq/K562GW/limma_logFC_sum.txt'
ARCHIVE_PATH = 'data/perturbseq/downloads/release.zip'


class PerturbseqInputTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.archive = self.root / 'release.zip'
        self.write_archive()
        self.output = self.root / 'repo'
        self.destination = self.output / DESTINATION

    def write_archive(self, name=MEMBER):
        with zipfile.ZipFile(self.archive, 'w') as handle:
            handle.writestr('README.md', b'release')
            handle.writestr(name, PAYLOAD)

    def manifest(self):
        return {
            'archive': {'url': self.archive.as_uri(), 'path': ARCHIVE_PATH,
                        'bytes': self.archive.stat().st_size, 'md5': checksum(self.archive, 'md5')},
            'members': {'k562gw_logfc': {'member': MEMBER, 'path': DESTINATION,
                                         'bytes': len(PAYLOAD), 'sha256': hashlib.sha256(PAYLOAD).hexdigest()}},
        }

    def test_download_extract_and_reuse(self):
        manifest = self.manifest()
        inputs.stage(manifest, self.output)
        self.assertEqual(self.destination.read_bytes(), PAYLOAD)
        before = self.destination.stat().st_mtime_ns
        self.archive.unlink()  # Reuse must not need the archive again.
        inputs.stage(manifest, self.output)
        self.assertEqual(self.destination.stat().st_mtime_ns, before)

    def test_member_mismatch_is_not_published(self):
        manifest = self.manifest()
        manifest['members']['k562gw_logfc']['sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'Checksum mismatch'):
            inputs.stage(manifest, self.output)
        self.assertFalse(self.destination.exists())
        self.assertFalse(list(self.destination.parent.glob('.perturbseq-*')))

    def test_existing_destination_is_never_replaced(self):
        self.destination.parent.mkdir(parents=True)
        self.destination.write_bytes(b'local scientific input')
        with self.assertRaises(ValueError):
            inputs.stage(self.manifest(), self.output)
        self.assertEqual(self.destination.read_bytes(), b'local scientific input')

    def test_source_root_links_verified_copy(self):
        source = self.root / 'extracted'
        (source / MEMBER).parent.mkdir(parents=True)
        (source / MEMBER).write_bytes(PAYLOAD)
        manifest = self.manifest()
        self.archive.unlink()  # Linking must not touch the archive.
        inputs.stage(manifest, self.output, source_root=source)
        self.assertTrue(self.destination.is_symlink())
        self.assertEqual(self.destination.read_bytes(), PAYLOAD)

    def test_source_root_rejects_wrong_copy(self):
        source = self.root / 'extracted'
        (source / MEMBER).parent.mkdir(parents=True)
        (source / MEMBER).write_bytes(b'other release')
        with self.assertRaisesRegex(ValueError, 'Checksum mismatch'):
            inputs.stage(self.manifest(), self.output, source_root=source)
        self.assertFalse(self.destination.exists())

    def test_check_does_not_download(self):
        with self.assertRaises(FileNotFoundError):
            inputs.stage(self.manifest(), self.output, check=True)
        self.assertFalse((self.output / ARCHIVE_PATH).exists())

    def test_missing_member_is_rejected(self):
        self.write_archive(name='data/elsewhere.txt')
        with self.assertRaisesRegex(ValueError, 'has no member'):
            inputs.stage(self.manifest(), self.output)
        self.assertFalse(self.destination.exists())

    def test_unsafe_member_is_rejected(self):
        manifest = self.manifest()
        manifest['members']['k562gw_logfc']['member'] = '../escape.txt'
        with self.assertRaisesRegex(ValueError, 'Unsafe archive member'):
            inputs.stage(manifest, self.output)


if __name__ == '__main__':
    unittest.main()

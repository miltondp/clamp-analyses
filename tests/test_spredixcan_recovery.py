"""Integrity and resume checks for the public S-PrediXcan recovery path."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

import pandas as pd

SPEC = importlib.util.spec_from_file_location('spredixcan_inputs', Path(__file__).resolve().parents[1] / 'scripts/inputs/spredixcan_inputs.py')
inputs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inputs)


class TissueRecoveryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.tsv = self.root / 'tissue.tsv.gz'
        self.frame = pd.DataFrame({'trait-A': [0.125, -2.0], 'trait-B': [3.5, 4.0]}, index=pd.Index(['ENSG1', 'ENSG2'], name='gene_name'))
        self.frame.to_csv(self.tsv, sep='\t')
        original = self.root / 'original.pkl'
        self.frame.to_pickle(original)
        self.asset = {'shape': [2, 2], 'index_name': 'gene_name', 'columns_name': None,
                      'semantic_sha256': inputs.fingerprint(self.frame),
                      'pickle_sha256': inputs.checksum(original),
                      'tsv_sha256': inputs.checksum(self.tsv), 'tsv_bytes': self.tsv.stat().st_size,
                      'archive_member': 'spredixcan/tissue.tsv.gz'}
        self.destination = self.root / 'pkl/model.pkl'

    def test_convert_checks_values_and_preserves_existing_output(self):
        inputs.convert(self.tsv, self.destination, self.asset)
        pd.testing.assert_frame_equal(pd.read_pickle(self.destination), self.frame)
        before = self.destination.stat().st_mtime_ns
        self.tsv.unlink()
        inputs.convert(self.tsv, self.destination, self.asset)
        self.assertEqual(before, self.destination.stat().st_mtime_ns)

    def test_rounded_values_are_rejected(self):
        self.frame.iloc[0, 0] += 0.00001
        self.frame.to_csv(self.tsv, sep='\t')
        asset = dict(self.asset, tsv_sha256=inputs.checksum(self.tsv))
        with self.assertRaisesRegex(ValueError, 'matrix mismatch'):
            inputs.convert(self.tsv, self.destination, asset)
        self.assertFalse(self.destination.exists())

    def test_existing_changed_matrix_is_never_replaced(self):
        self.destination.parent.mkdir()
        self.frame.iloc[0, 0] = 100
        self.frame.to_pickle(self.destination)
        before = self.destination.read_bytes()
        with self.assertRaises(ValueError):
            inputs.convert(self.tsv, self.destination, self.asset)
        self.assertEqual(before, self.destination.read_bytes())

    def test_alternate_pickle_serialization_requires_identical_values(self):
        self.destination.parent.mkdir()
        self.frame.to_pickle(self.destination, protocol=4)
        inputs.verify_pickle(self.destination, self.asset)

    def test_archive_download_extraction_and_resume(self):
        source = self.root / 'source.tar'
        with tarfile.open(source, 'w') as handle:
            handle.add(self.tsv, arcname=self.asset['archive_member'])
            # Unrequested members must never be extracted.
            unwanted = tarfile.TarInfo('../unrequested.txt')
            unwanted.size = 4
            handle.addfile(unwanted, io.BytesIO(b'nope'))
        manifest = {'archive': {'url': source.as_uri(), 'bytes': source.stat().st_size,
                               'md5': inputs.checksum(source, 'md5')},
                    'tissues': {'model.pkl': self.asset}}
        archive = self.root / 'cache/archive.tar'
        inputs.recover(self.destination.parent, manifest, archive=archive)
        pd.testing.assert_frame_equal(pd.read_pickle(self.destination), self.frame)
        self.assertFalse((self.root / 'unrequested.txt').exists())
        source.unlink()
        archive.unlink()
        inputs.recover(self.destination.parent, manifest, archive=archive)

    def test_cached_conversion_survives_workflow_output_cleanup(self):
        manifest = {'tissues': {'model.pkl': self.asset}}
        # The fixture name can differ from the public member basename.
        source = self.root / 'tissue.tsv.gz'
        inputs.recover(self.destination.parent, manifest, tsv_root=source.parent)
        inode = self.destination.stat().st_ino
        self.destination.unlink()
        self.destination.parent.rmdir()
        source.unlink()
        inputs.recover(self.destination.parent, manifest)
        self.assertEqual(inode, self.destination.stat().st_ino)
        pd.testing.assert_frame_equal(pd.read_pickle(self.destination), self.frame)

    def test_corrupt_archive_is_rejected(self):
        archive = self.root / 'corrupt.tar'
        archive.write_bytes(b'not the published file')
        with self.assertRaisesRegex(ValueError, 'archive checksum mismatch'):
            inputs.obtain_archive(archive, {'bytes': archive.stat().st_size, 'md5': '0' * 32})


if __name__ == '__main__':
    unittest.main()

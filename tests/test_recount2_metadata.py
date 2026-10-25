import importlib.util
import io
from pathlib import Path
import tarfile
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('metadata', Path(__file__).parents[1] / 'scripts/inputs/extract_recount2_metadata.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ExtractionTests(unittest.TestCase):
    def archive(self, root, name, symlink=False):
        path = root / 'source.tar.gz'
        with tarfile.open(path, 'w:gz') as tar:
            member = tarfile.TarInfo(name)
            if symlink:
                member.type = tarfile.SYMTYPE
                member.linkname = '/tmp/escape'
                tar.addfile(member)
            else:
                member.size = 3
                tar.addfile(member, io.BytesIO(b'abc'))
        return path

    def test_layout_and_existing_destination(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = self.archive(root, 'recount2/data/SRP001/rse_gene.Rdata')
            dst = root / 'metadata'
            module.extract(archive, dst)
            self.assertEqual((dst / 'SRP001/rse_gene.Rdata').read_bytes(), b'abc')
            with self.assertRaises(FileExistsError):
                module.extract(archive, dst)

    def test_unsafe_empty_and_wrong_layout(self):
        for name, link in [('recount2/data/../../escape', False), ('recount2/data/link', True), ('unexpected/file', False), ('recount2/data/SRP001/other.txt', False)]:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                archive = self.archive(root, name, link)
                with self.assertRaises(ValueError):
                    module.extract(archive, root / 'metadata')
                self.assertFalse((root / 'metadata').exists())

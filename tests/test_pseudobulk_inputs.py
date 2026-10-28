"""Pinned acquisition, original export recovery, and count aggregation checks."""
from pathlib import Path
import hashlib
import importlib.util
import io
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts/inputs'))
sys.path.insert(0, str(ROOT / 'scripts/pseudobulk'))
import pseudobulk_inputs as acquisition
import build_pseudobulk as build
from common import transform_ids, map_labels


class Response(io.BytesIO):
    def __init__(self, payload, start, stop, total, status=206):
        super().__init__(payload)
        self.status = status
        self.headers = {'Content-Range': f'bytes {start}-{stop}/{total}'}


class DownloadTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        payload = io.BytesIO()
        with zipfile.ZipFile(payload, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr('lab/bulk.csv', b'gene,donor\nA,1.25\n')
            archive.writestr('../unrequested.csv', b'never extracted')
        self.archive = payload.getvalue()
        self.asset = {'url': 'https://example.invalid/archive.zip', 'archive_bytes': len(self.archive),
                      'archive_member': 'lab/bulk.csv', 'bytes': 18,
                      'path': 'data/input.csv', 'sha256': hashlib.sha256(b'gene,donor\nA,1.25\n').hexdigest()}
        self.asset['bytes'] = len(b'gene,donor\nA,1.25\n')

    def request(self, request, **kwargs):
        start, stop = map(int, request.get_header('Range').removeprefix('bytes=').split('-'))
        return Response(self.archive[start:stop+1], start, stop, len(self.archive))

    def test_selected_member_integrity_and_existing_reuse(self):
        with patch('urllib.request.urlopen', side_effect=self.request):
            acquisition.fetch_member('tiny', self.asset, self.root)
        output = self.root / self.asset['path']
        self.assertEqual(output.read_bytes(), b'gene,donor\nA,1.25\n')
        self.assertFalse((self.root / 'unrequested.csv').exists())
        before = output.stat().st_mtime_ns
        with patch('urllib.request.urlopen', side_effect=AssertionError('No second download')):
            acquisition.fetch_member('tiny', self.asset, self.root)
        self.assertEqual(before, output.stat().st_mtime_ns)

    def test_unverified_member_and_existing_edits_are_rejected(self):
        bad = dict(self.asset, sha256='0' * 64)
        with patch('urllib.request.urlopen', side_effect=self.request), self.assertRaisesRegex(ValueError, 'Checksum mismatch'):
            acquisition.fetch_member('tiny', bad, self.root)
        output = self.root / self.asset['path']
        self.assertFalse(output.exists())
        output.write_text('author-owned data')
        with self.assertRaises(ValueError):
            acquisition.fetch_member('tiny', self.asset, self.root)
        self.assertEqual(output.read_text(), 'author-owned data')

    def test_server_ignoring_range_fails_without_whole_archive_fallback(self):
        with patch('urllib.request.urlopen', return_value=Response(self.archive, 0, len(self.archive)-1, len(self.archive), status=200)):
            with self.assertRaisesRegex(ValueError, 'byte range'):
                acquisition.RemoteZipFile(self.asset['url'],len(self.archive)).read(1)


class AggregationTests(unittest.TestCase):
    def test_fixed_selection_and_identifier_mapping(self):
        values=transform_ids(np.array(['1-1','2-2']), 'sample_prefix_hyphen_to_underscore')
        self.assertEqual(values.tolist(),['sample_1_1','sample_2_2'])
        manifest=pd.DataFrame({'sample':['sample_2_2','sample_1_1']})
        samples,mask,counts=build.cohort_selection(values,manifest,1)
        self.assertEqual(samples,manifest['sample'].tolist())
        np.testing.assert_array_equal(counts,[1,1])
        with self.assertRaisesRegex(ValueError,'absent'):
            build.cohort_selection(values,pd.DataFrame({'sample':['missing']}),None)
        with self.assertRaisesRegex(ValueError,'below'):
            build.cohort_selection(values,manifest,2)

    def test_count_normalization_exclusion_and_duplicate_genes(self):
        with tempfile.TemporaryDirectory() as temporary:
            p=Path(temporary)/'input.h5ad'
            matrix=sp.csr_matrix([[1,1,2],[3,1,0],[5,0,5],[0,2,2]])
            with h5py.File(p,'w') as f:
                g=f.create_group('X');g.attrs['shape']=matrix.shape;g.attrs['encoding-type']='csr_matrix'
                for key in ['data','indices','indptr']:g.create_dataset(key,data=getattr(matrix,key))
                o=f.create_group('obs');o.create_dataset('donor',data=np.array(['b','a','b','a'],dtype='S'))
                o.create_dataset('label',data=np.array(['T','B','excluded','T'],dtype='S'))
                v=f.create_group('var');v.create_dataset('gene_name',data=np.array(['A','A','B'],dtype='S'))
            cfg={'raw':str(p),'matrix':'X','sample_col':'donor','truth_v0_col':'label',
                 'exclude':{'label':['excluded']},'min_cells':1,'chunk_cells':2}
            manifest=pd.DataFrame({'sample':['b','a']})
            result=build.aggregate_h5ad('tiny',cfg,manifest)
            # Cell CPM before donor mean; duplicate A features averaged afterwards.
            np.testing.assert_allclose(result.bulk.loc['b'],[250000,500000])
            np.testing.assert_allclose(result.bulk.loc['a'],[375000,250000])
            self.assertEqual(result.info.nCells.tolist(),[1,2])
            self.assertEqual(result.truth_v0.loc['b','T'],1)
            np.testing.assert_allclose(result.truth_v0.loc['a',['B','T']],[.5,.5])

    def test_reference_gene_filtering_preserves_normalization_denominator(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            (root/'axes.json').write_text(json.dumps({'samples':['a'],'truthFrac_v0':{'samples':['a'],'labels':['T']}}))
            pd.DataFrame({'gene':['retained']}).to_csv(root/'genes.csv',index=False)
            pd.DataFrame({'sample':['a'],'nCells':[1]}).to_csv(root/'patient_info.csv',index=False)
            result=build.Result(pd.DataFrame({'retained':[250000.],'filtered':[750000.]},index=['a']),pd.DataFrame({'sample':['a'],'nCells':[1]}),pd.DataFrame({'T':[1.]},index=['a']),None,{'kind':'h5ad','count_like_fraction':1})
            filtered=build.reference_axes(result,{'reference_root':str(root)})
            self.assertEqual(filtered.bulk.columns.tolist(),['retained'])
            self.assertEqual(filtered.bulk.loc['a','retained'],250000.)
            build.validate('tiny',filtered)

    def test_reference_axes_reject_missing_genes_and_changed_cell_counts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            (root/'axes.json').write_text(json.dumps({'samples':['b','a'],'truthFrac_v0':{'samples':['b','a'],'labels':['T']}}))
            pd.DataFrame({'gene':['G']}).to_csv(root/'genes.csv',index=False)
            pd.DataFrame({'sample':['a','b'],'nCells':[2,1]}).to_csv(root/'patient_info.csv',index=False)
            result=build.Result(pd.DataFrame({'G':[1.,2.]},index=['a','b']),pd.DataFrame({'sample':['a','b'],'nCells':[2,1]}),pd.DataFrame({'T':[1.,1.]},index=['a','b']),None,{'kind':'prepared'})
            restored=build.reference_axes(result,{'reference_root':str(root)})
            self.assertEqual(restored.bulk.index.tolist(),['b','a'])
            self.assertEqual(restored.info['sample'].tolist(),['a','b'])
            build.validate('tiny',restored)
            restored.info.loc[0,'nCells']=999
            with self.assertRaisesRegex(ValueError,'Cell totals'):
                build.reference_axes(restored,{'reference_root':str(root)})
            restored.bulk.columns=['missing']
            with self.assertRaisesRegex(ValueError,'gene mismatch'):
                build.reference_axes(restored,{'reference_root':str(root)})


if __name__=='__main__':
    unittest.main()

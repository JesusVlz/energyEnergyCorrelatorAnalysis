from array import array
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import merge_response as merge

class ManifestGates(unittest.TestCase):
    def fixture(self,base,name,input_name,offset=0):
        folder=base/name;folder.mkdir()
        counts=[0,0]
        with gzip.open(folder/'event_identity.csv.gz','wt') as f:
            f.write('source_entry,run,lumi,event,event_key,fold,fold_entry\n')
            for i in range(4):
                key=merge.event_key(dict(sample='pp',sample_namespace='sample',id=[1,2,i+offset]));which=key%2
                f.write(','.join(map(str,[i,1,2,i+offset,key,which,counts[which]]))+'\n');counts[which]+=1
        record=dict(status='COMPLETE_ENGINEERING_ONLY',system='pp',source_input=input_name,
                    native_sha256='native',card_sha256='card',helper_sha256='helper',pair_counts_sha256='pair',
                    skim=dict(processed_entries=4,input_entries=4,input_size_bytes=7,sample_namespace='sample',fold_entries=counts,
                              ledger_sha256=merge.digest(folder/'event_identity.csv.gz')),folds={})
        for label in ('train','validation'):
            path=folder/(label+'_raw.root');path.write_bytes(b'fixture');record['folds'][label]=dict(sha256=merge.digest(path))
        (folder/'receipt.json').write_text(json.dumps(record))
        return folder
    def manifest(self,*names):
        return dict(tasks=[dict(system='pp',input=n,sample_namespace='sample',input_identity=dict(entries=4,size_bytes=7)) for n in names])
    def report(self):return dict(card={'WeightExponent':[1,2]},inclusive_activity_accounted=True)
    def test_full_manifest_missing_input_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);folder=self.fixture(p,'a','a.root')
            with patch.object(merge,'raw_report',return_value=self.report()),self.assertRaisesRegex(ValueError,'Incomplete frozen'):
                merge.validate_inputs([folder],'pp',self.manifest('a.root','b.root'),p/'db')
    def test_global_duplicate_event_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);a=self.fixture(p,'a','a.root');b=self.fixture(p,'b','b.root')
            with patch.object(merge,'raw_report',return_value=self.report()),self.assertRaisesRegex(ValueError,'Duplicate global'):
                merge.validate_inputs([a,b],'pp',self.manifest('a.root','b.root'),p/'db')
    def test_complete_disjoint_manifest_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);a=self.fixture(p,'a','a.root');b=self.fixture(p,'b','b.root',offset=4)
            with patch.object(merge,'raw_report',return_value=self.report()):
                report=merge.validate_inputs([a,b],'pp',self.manifest('a.root','b.root'),p/'db')
            self.assertEqual(report['events'],8);self.assertEqual(sum(report['fold_events']),8)
    def test_changed_output_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);a=self.fixture(p,'a','a.root');(a/'train_raw.root').write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'Changed native'):
                merge.validate_inputs([a],'pp',self.manifest('a.root'),p/'db')
try:
    import ROOT
except ImportError:
    ROOT = None


@unittest.skipIf(ROOT is None, 'PyROOT is required for native histogram merge tests')
class RecursiveHistogramMerge(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ROOT.gROOT.SetBatch(True)

    def fixture(self, path, weight, *, nested=True, extra=False, histogram_class='TH1D',
                edges=(0., 1., 2.), label='selected', sparse_bins=(2, 3)):
        handle = ROOT.TFile.Open(str(path), 'RECREATE')
        try:
            top = ROOT.TH1D('counts', '', 2, 0., 2.)
            top.Sumw2()
            top.Fill(.5, weight)
            top.Write()
            top.Fill(1.5, 2 * weight)
            top.Write()  # Only this latest cycle may contribute to the sum.
            ROOT.TObjString('per-input metadata').Write('ResponseProvenance')
            card = handle.mkdir('JCard')
            card.cd()
            ROOT.TObjString('fixture-card').Write('Configuration')
            qa = handle.mkdir('eventQA')
            qa.cd()
            ROOT.TObjString('not additive').Write('InputMetadata')
            if nested:
                histogram = getattr(ROOT, histogram_class)('counts', '', len(edges) - 1, array('d', edges))
                histogram.Sumw2()
                histogram.GetXaxis().SetBinLabel(1, label)
                for value, factor in ((-.5, 3), (.5, 1), (2.5, 4)):
                    histogram.Fill(value, factor * weight)
                histogram.Write()
            if extra:
                ROOT.TH1D('unexpected', '', 1, 0., 1.).Write()
            deep = qa.mkdir('deep')
            deep.cd()
            sparse = ROOT.THnSparseD('response', '', 2, array('i', sparse_bins),
                                     array('d', [0., 0.]), array('d', [2., 3.]))
            sparse.Sumw2()
            sparse.Fill(array('d', [.25, .5]), weight)
            sparse.Fill(array('d', [-.25, 1.5]), 2 * weight)
            sparse.Write()
        finally:
            handle.Close()

    def test_nested_histograms_errors_flows_and_latest_cycles_are_summed(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            inputs = [base/'a.root', base/'b.root']
            self.fixture(inputs[0], 2.)
            self.fixture(inputs[1], 3.)
            output = base/'merged.root'
            with patch.object(merge, 'raw_report', return_value={'fixture': True}):
                merge.merge_raw(inputs, output, {'fold': 'train'})
            result = ROOT.TFile.Open(str(output), 'READ')
            try:
                top = result.Get('counts')
                self.assertEqual(top.GetBinContent(1), 5.)
                self.assertEqual(top.GetBinContent(2), 10.)
                self.assertAlmostEqual(top.GetBinError(2) ** 2, 52.)
                nested = result.Get('eventQA/counts')
                for cell, factor in ((0, 3), (1, 1), (2, 0), (3, 4)):
                    self.assertEqual(nested.GetBinContent(cell), 5. * factor)
                    self.assertAlmostEqual(nested.GetBinError(cell) ** 2, 13. * factor ** 2)
                self.assertEqual(nested.GetXaxis().GetBinLabel(1), 'selected')
                sparse = result.Get('eventQA/deep/response')
                for coordinates, factor in (((1, 1), 1), ((0, 2), 2)):
                    cell = sparse.GetBin(array('i', coordinates))
                    self.assertEqual(sparse.GetBinContent(cell, array('i', [0, 0])), 5. * factor)
                    self.assertAlmostEqual(sparse.GetBinError2(cell), 13. * factor ** 2)
                self.assertEqual(merge.card_values(result.GetDirectory('JCard')), {'Configuration': 'fixture-card'})
                self.assertEqual(len(list(result.GetDirectory('JCard').GetListOfKeys())), 1)
                self.assertEqual(json.loads(str(result.Get('ResponseMergeProvenance').GetString())), {'fold': 'train'})
                self.assertFalse(result.Get('ResponseProvenance'))
                self.assertFalse(result.Get('eventQA/InputMetadata'))
            finally:
                result.Close()

    def test_missing_or_extra_nested_histograms_are_rejected(self):
        for change in ({'nested': False}, {'extra': True}):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                self.fixture(base/'a.root', 2.)
                self.fixture(base/'b.root', 3., **change)
                with self.assertRaisesRegex(ValueError, 'Different raw histogram paths'):
                    merge.merge_raw([base/'a.root', base/'b.root'], base/'merged.root', {})
                self.assertFalse((base/'merged.root').exists())

    def test_incompatible_classes_axes_and_labels_are_rejected(self):
        for change in ({'histogram_class': 'TH1F'}, {'edges': (0., .75, 2.)},
                       {'label': 'other-selection'}, {'sparse_bins': (3, 3)}):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                self.fixture(base/'a.root', 2.)
                self.fixture(base/'b.root', 3., **change)
                with self.assertRaisesRegex(ValueError, 'Incompatible raw histogram class/binning/labels'):
                    merge.merge_raw([base/'a.root', base/'b.root'], base/'merged.root', {})
                self.assertFalse((base/'merged.root').exists())

    def test_failed_root_add_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            self.fixture(base/'a.root', 2.)
            self.fixture(base/'b.root', 3.)
            with patch.object(ROOT.TH1D, 'Add', return_value=False), self.assertRaisesRegex(ValueError, 'ROOT failed to add'):
                merge.merge_raw([base/'a.root', base/'b.root'], base/'merged.root', {})
            self.assertFalse((base/'merged.root').exists())


if __name__ == '__main__':
    unittest.main()

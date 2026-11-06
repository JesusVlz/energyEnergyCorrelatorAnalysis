#!/usr/bin/env python3
"""ROOT integration regressions for campaign merging; run with PyROOT and hadd."""
from array import array
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import ROOT
from merge_root import inspect_root, merge_root


class RootMergeTest(unittest.TestCase):
    def setUp(self):
        ROOT.gROOT.SetBatch(True)
        self.temp = tempfile.TemporaryDirectory(prefix='eec-merge-test-')
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def write_input(self, name, weight, data_type=10, jet_bins=3, nonfinite=False):
        path = self.directory/name
        file = ROOT.TFile(str(path),'RECREATE')
        events = ROOT.TH1F('nEvents','events',9,-.5,8.5)
        for cut in range(9):
            events.Fill(cut)
        events.Write()
        # Weighted dense underflow/overflow and cancellation must survive hadd.
        dense = ROOT.TH1F('weighted','weighted',3,0,3)
        dense.Sumw2()
        for x,w in ((-1,2*weight),(.5,weight),(1.5,-weight),(4,3*weight)):
            dense.Fill(x,w)
        if nonfinite:
            dense.SetBinContent(1,float('nan'))
        dense.Write()
        for hist_name in ('inclusiveJet','energyEnergyCorrelator'):
            hist = ROOT.THnSparseF(hist_name,hist_name,2,array('i',[jet_bins,2]),
                                  array('d',[0,0]),array('d',[3,2]))
            hist.Sumw2()
            hist.Fill(array('d',[.5,.5]),weight)
            hist.Fill(array('d',[-1,3]),-weight)
            hist.Write()
        qa = file.mkdir('eventQA')
        qa.cd()
        histogram = ROOT.TH2D('coverage','coverage',2,0,2,2,0,2)
        histogram.Sumw2()
        histogram.Fill(.5,1.5,weight)
        histogram.Write()
        ROOT.TNamed('WeightState','Same immutable metadata').Write()
        file.cd()
        ROOT.TObjString('input-'+name).Write('InputFiles')
        card = file.mkdir('JCard')
        card.cd()
        for key,value in (('DataType',data_type),('McCorrelationType',3),('BackgroundMethods',6)):
            vector = ROOT.TVectorF(1,1)
            vector.GetMatrixArray()[0] = value
            vector.Write(key)
        ROOT.TObjString('/worker/'+name+'/pool.txt').Write('MixingFileListPath')
        ROOT.TObjString('same-git').Write('GitHash')
        file.Close()
        return path

    def test_weighted_dense_sparse_and_configuration_survive_real_hadd(self):
        first = self.write_input('a.root',2)
        second = self.write_input('b.root',-1)
        original = inspect_root(first)
        output = self.directory/'merged.root'
        result = merge_root([first,second],output,{'all_events':2,'inclusive_jet_entries':4,'eec_entries':4})
        self.assertEqual(result['status'],'PASS_ALL_HISTOGRAM_BINS_AND_CANONICAL_CARD')
        self.assertEqual(result['root_validation']['card'],original['card'])
        self.assertEqual(result['root_validation']['n_events_bins'],[0]+[2]*9+[0])
        file = ROOT.TFile.Open(str(output))
        dense = file.Get('weighted')
        self.assertAlmostEqual(dense.GetBinContent(0),2)
        self.assertAlmostEqual(dense.GetBinContent(2),-1)
        self.assertAlmostEqual(dense.GetBinContent(4),3)
        self.assertAlmostEqual(dense.GetBinError(1)**2,5)
        provenance = json.loads(str(file.Get('MergeProvenance').GetString()))
        self.assertEqual(provenance['input_count'],2)
        self.assertNotEqual(provenance['inputs'][0]['runtime_card'],provenance['inputs'][1]['runtime_card'])
        self.assertEqual(provenance['inputs'][1]['runtime_metadata']['InputFiles']['value'],'input-b.root')
        file.Close()

    def test_physics_card_and_axis_mismatches_rejected_before_hadd(self):
        first = self.write_input('a.root',1)
        for second in (self.write_input('wrong_card.root',1,data_type=12),
                       self.write_input('wrong_axis.root',1,jet_bins=4)):
            with self.subTest(second=second.name):
                with mock.patch('merge_root.subprocess.run') as process:
                    with self.assertRaises(ValueError):
                        merge_root([first,second],self.directory/'bad.root')
                    process.assert_not_called()

    def test_nonfinite_bins_duplicate_inputs_and_wrong_ledger_totals_rejected(self):
        good = self.write_input('a.root',1)
        bad = self.write_input('nonfinite.root',1,nonfinite=True)
        with self.assertRaises(Exception):
            inspect_root(bad)
        with self.assertRaises(ValueError):
            merge_root([good,good],self.directory/'duplicate.root')
        with mock.patch('merge_root.subprocess.run') as process:
            with self.assertRaisesRegex(ValueError,'production all_events'):
                merge_root([good],self.directory/'bad_total.root',{'all_events':99})
            process.assert_not_called()

    def test_bin_damage_with_conserved_entries_and_integral_is_detected(self):
        first = self.write_input('a.root',2)
        second = self.write_input('b.root',1)
        output = self.directory/'merged.root'
        real_restore = __import__('merge_root')._restore_card
        def damage(*args):
            real_restore(*args)
            file = ROOT.TFile.Open(str(output),'UPDATE')
            dense = file.Get('weighted')
            entries = dense.GetEntries()
            dense.SetBinContent(1,dense.GetBinContent(1)+1)
            dense.SetBinContent(2,dense.GetBinContent(2)-1)
            dense.SetEntries(entries)
            dense.Write('weighted',ROOT.TObject.kOverwrite)
            file.Close()
        with mock.patch('merge_root._restore_card',side_effect=damage):
            with self.assertRaisesRegex(Exception,'bin does not sum'):
                merge_root([first,second],output)


if __name__=='__main__':
    unittest.main()

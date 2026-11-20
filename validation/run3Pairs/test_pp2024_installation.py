"""Software installation/configuration tests; no detector closure assertion."""
import importlib.util
import json
from pathlib import Path
import unittest
R=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('pp_install',R/'cards/install_ppref_pair_tables.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class PpInstallation(unittest.TestCase):
    def setUp(self):
        self.manifest=m.manifest();self.payload=self.manifest['tables']['Nominal']['path']
        self.base=(R/'cardEECpp5v36TeVMCR04.input').read_text()
    def test_three_payloads(self):
        self.assertEqual(set(self.manifest['tables']),{'Nominal','Loose','Tight'})
        for row in self.manifest['tables'].values():
            self.assertEqual(len(row['sha256']),64);m.verify_single_track(row)
    def test_provenance(self):
        self.assertEqual(self.manifest['derivation_jet_id'],0)
        self.assertEqual(self.manifest['application_jet_id'],1)
    def test_six_cards(self):
        files=list((R/'cards/ppRef2024/R4').glob('*.input'));self.assertEqual(len(files),6)
        for path in files:
            c=m.card_values(path.read_text())
            for k,v in {'Run3JetID':'1','Run3UsePpL1FastJet':'0','DisableTrackPairEfficiencyCorrection':'0','MinPairTrackPt':'2','MinTrackPtCut':'0.5','WeightExponent':'1 2','MinJetPtCut':'15','MaxJetPtCut':'200','McCorrelationType':'0'}.items():
                self.assertEqual(c[k],v)
            self.assertEqual(int(c['FilledHistograms'])&(32|64|128|256),0)
            self.assertEqual(c['Run3TrackPairCorrectionFile'],self.manifest['tables'][c['Run3TrackSelection']]['path'])
    def test_idempotent_activation(self):
        new=m.activated_card(self.base,self.payload);self.assertEqual(new,m.activated_card(new,self.payload))
        before,after=m.card_values(self.base),m.card_values(new)
        for k in ('MinJetPtCut','MaxJetPtCut','MinJetPtUnfoldingReco','MinJetPtUnfoldingTruth','BackgroundMethods','JetPtBinEdgesEEC'):
            self.assertEqual(before[k],after[k])
    def test_bad_bytes(self):
        with self.assertRaises(ValueError):m.verify_bytes(b'not ROOT',self.manifest['tables']['Nominal'])
    def test_reject_oo(self):
        with self.assertRaises(ValueError):m.activated_card(self.base.replace('DataType 12','DataType 10'),self.payload)
    def test_reject_r8(self):
        with self.assertRaises(ValueError):m.activated_card(self.base.replace('JetRadius 0.4','JetRadius 0.8'),self.payload)
    def test_reject_response(self):
        with self.assertRaises(ValueError):m.activated_card(self.base.replace('FilledHistograms 31','FilledHistograms 279'),self.payload)
    def test_reject_gen_tracks(self):
        with self.assertRaises(ValueError):m.activated_card(self.base.replace('McCorrelationType 0','McCorrelationType 1'),self.payload)
    def test_reject_selection(self):
        with self.assertRaises(ValueError):m.activated_card(self.base.replace('Run3TrackSelection Nominal','Run3TrackSelection Tight'),self.payload)
    def test_reject_duplicate(self):
        with self.assertRaises(ValueError):m.activated_card(self.base+'\nRun3JetID 1\n',self.payload)
    def test_r8_bytes(self):
        records=json.loads((R/'cards/R8/moved.json').read_text());self.assertEqual(len(records),8)
        for old,row in records.items():
            self.assertFalse((R/old).exists())
            self.assertEqual(m.digest((R/row['path']).read_bytes()),row['sha256'])

if __name__=='__main__':unittest.main()

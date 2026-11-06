import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest
from array import array

import ROOT
ROOT.gROOT.SetBatch(True)
spec=importlib.util.spec_from_file_location('worker',Path(__file__).with_name('response_worker.py'))
worker=importlib.util.module_from_spec(spec)
spec.loader.exec_module(worker)

class SkimTest(unittest.TestCase):
    def fixture(self,path,duplicate=False,mismatch=False):
        f=ROOT.TFile.Open(str(path),'RECREATE')
        for dirname,treename in [('hiEvtAnalyzer','HiTree'),('ppTracks','trackTree')]:
            d=f.mkdir(dirname);d.cd();t=ROOT.TTree(treename,treename)
            run=array('I',[1]);lumi=array('I',[2]);event=array('L',[0]);value=array('i',[0])
            t.Branch('run',run,'run/i');t.Branch('lumi',lumi,'lumi/i');t.Branch('evt',event,'evt/l');t.Branch('value',value,'value/I')
            for i in range(19 if mismatch and dirname=='ppTracks' else 20):
                event[0]=0 if duplicate else i;value[0]=i*7;t.Fill()
                if i==9:t.Write()  # Older partial cycle must never be copied
            t.Write()
        f.mkdir('HiForestInfo').cd();t=ROOT.TTree('HiForest','metadata');x=array('i',[4]);t.Branch('x',x,'x/I');t.Fill();t.Write();f.Close()
    def test_synchronized_exact_partition(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);self.fixture(p/'input.root')
            record,files=worker.skim(p/'input.root',p/'split','pp','sample')
            self.assertEqual(sum(record['fold_entries']),20)
            seen=[]
            for which,path in enumerate(files):
                f=ROOT.TFile.Open(str(path));event=f.Get('hiEvtAnalyzer/HiTree');track=f.Get('ppTracks/trackTree')
                self.assertEqual(event.GetEntries(),track.GetEntries());self.assertEqual(f.Get('HiForestInfo/HiForest').GetEntries(),1)
                for i in range(event.GetEntries()):
                    event.GetEntry(i);track.GetEntry(i)
                    self.assertEqual(event.evt,track.evt);self.assertEqual(track.value,event.evt*7)
                    self.assertEqual(worker.fold(dict(sample='pp',sample_namespace='sample',id=[1,2,int(event.evt)])),which)
                    seen.append(int(event.evt))
                f.Close()
            self.assertEqual(sorted(seen),list(range(20)))
            with gzip.open(p/'split/event_identity.csv.gz','rt') as f:self.assertEqual(len(f.readlines()),21)
    def test_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);self.fixture(p/'input.root',mismatch=True)
            with self.assertRaisesRegex(ValueError,'Unsynchronized'):
                worker.skim(p/'input.root',p/'split','pp','sample')
    def test_duplicate_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);self.fixture(p/'input.root',duplicate=True)
            with self.assertRaisesRegex(ValueError,'Duplicate'):
                worker.skim(p/'input.root',p/'split','pp','sample')
if __name__=='__main__':unittest.main()

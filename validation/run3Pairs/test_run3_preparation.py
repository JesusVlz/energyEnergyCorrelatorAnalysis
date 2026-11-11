"""Software regression tests only; none of these assert detector/MC closure."""
import json
import math
import shlex
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
import preparation_contract as c
import prepare_inputs as prep
import pair_counts as pairs

CPP = r'''
#include "Run3TrackSelection.h"
#include "Run3TrackCharge.h"
#include "TrackingTableValidation.h"
#include "trackingEfficiency2024pp.h"
#include "trackingEfficiency2025OO.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2.h>
#include <TH3.h>
#include <cassert>
#include <limits>
#include <iostream>
using namespace EECRun3Tracking;
void writeCharges(const std::string& path, bool vector, int sign) {
  TFile file(path.c_str(),"RECREATE"); TTree tree("tracks","");
  int n=4, q[4]={-sign,sign,0,2}; std::vector<char> v={static_cast<char>(-sign),static_cast<char>(sign),0,2};
  tree.Branch("nTrk",&n,"nTrk/I");
  if(vector) tree.Branch("trkCharge",&v); else tree.Branch("trkCharge",q,"trkCharge[nTrk]/I");
  tree.Fill(); tree.Write();
}
int main(int argc,char** argv) {
  assert(argc==2); const std::string tmp=argv[1];
  assert(PassMomentum(2,1)); assert(PassMomentum(10,8));
  assert(PassMomentum(20,2)); assert(!PassMomentum(20,2.01));
  assert(!PassMomentum(2,-1)); assert(!PassMomentum(0,0));
  assert(!PassMomentum(std::numeric_limits<double>::quiet_NaN(),0));
  for(const char* name:{"Nominal","Loose","Tight"}) {
    double limit=VertexLimit(name);
    assert(PassVertex(limit,1,limit)); assert(PassVertex(-limit,1,limit));
    assert(!PassVertex(limit+.01,1,limit)); assert(!PassVertex(0,0,limit));
    assert(!PassVertex(0,-1,limit));
    assert(ResolveSelection(name,.1,limit)==name);
  }
  bool rejected=false; try { ResolveSelection("Tight",.05,2); } catch(const std::invalid_argument&) {rejected=true;} assert(rejected);
  rejected=false; try { ResolveSelection("Tight",.1,3); } catch(const std::invalid_argument&) {rejected=true;} assert(rejected);
  assert(Payload(true,"Loose").find("3D_Loose")!=std::string::npos);
  assert(Payload(false,"Tight").find("2D_Tight")!=std::string::npos);
  for(bool vector:{false,true}) {
    const std::string a=tmp+(vector?"/v1.root":"/a1.root"), b=tmp+(vector?"/v2.root":"/a2.root");
    writeCharges(a,vector,1); writeCharges(b,vector,-1);
    TChain chain("tracks"); chain.Add(a.c_str()); chain.Add(b.c_str());
    TrackCharge charge; charge.Bind(&chain);
    chain.GetEntry(0); assert(charge.Get(0)==-1); assert(charge.Get(1)==1);
    assert(charge.Get(2)==0); // Neutral candidates reach the analyzer charge cut.
    rejected=false; try { charge.Get(3); } catch(const std::runtime_error&) {rejected=true;} assert(rejected);
    chain.GetEntry(1); assert(charge.Get(0)==1); assert(charge.Get(1)==-1);
    assert(charge.Get(2)==0);
    rejected=false; try { charge.Get(3); } catch(const std::runtime_error&) {rejected=true;} assert(rejected);
  }
  TTree missing("missing",""); TrackCharge charge;
  rejected=false; try {charge.Bind(&missing);} catch(const std::runtime_error&){rejected=true;} assert(rejected);
  std::size_t checks=0;
  for(bool oxygen:{true,false}) for(const char* selection:{"Nominal","Loose","Tight"}) {
    const std::string path=Payload(oxygen,selection); TFile file(path.c_str(),"READ"); assert(!file.IsZombie());
    const std::string suffix=oxygen?"_3D":"_2D";
    auto* eff=dynamic_cast<TH1*>(file.Get(("hEff"+suffix).c_str()));
    auto* fak=dynamic_cast<TH1*>(file.Get(("hFak"+suffix).c_str()));
    auto* sec=dynamic_cast<TH1*>(file.Get(("hSec"+suffix).c_str()));
    auto* mul=dynamic_cast<TH1*>(file.Get(("hMul"+suffix).c_str()));
    ValidateTables({eff,fak,sec,mul},oxygen?3:2,path);
    std::unique_ptr<TrackingEfficiencyInterface> tool;
    if(oxygen) tool.reset(new TrkEff2025OO(true,path)); else tool.reset(new TrkEff2024ppRef(true,path));
    for(float eta:{-1.5f,0.f,1.5f}) for(float pt:{.6f,2.01f,6.1f,10.f,20.f,100.f}) for(int activity:{0,60,100,180,199}) {
      const float lookup=oxygen && pt>6.1f?6.1f:pt;
      const int bin=oxygen?dynamic_cast<TH3*>(eff)->FindBin(eta,lookup,activity):dynamic_cast<TH2*>(eff)->FindBin(eta,lookup);
      const double e=eff->GetBinContent(bin), f=fak->GetBinContent(bin), s=sec->GetBinContent(bin), m=mul->GetBinContent(bin);
      const double denominator=e*(1+m);
      const double expected=ValidContents(e,f,s,m) && (oxygen?e>.001:denominator>.001) ? (1-f)*(1-s)/denominator : 0.;
      const double actual=tool->getCorrection(pt,eta,activity);
      if(std::abs(actual-expected)>1e-5*std::max(1.,std::abs(expected))) {
        std::cerr<<path<<" eta="<<eta<<" pt="<<pt<<" hiBin="<<activity<<" got="<<actual<<" expected="<<expected<<std::endl; return 3;
      }
      ++checks;
    }
    assert(tool->getCorrection(std::numeric_limits<float>::quiet_NaN(),0,0)==0);
  }
  std::cout<<"PASS cut boundaries, AOD/MiniAOD physical charge across chained files, and "<<checks<<" actual-payload lookup comparisons"<<std::endl;
}
'''


class PreparationTests(unittest.TestCase):
    def test_cards_and_payload_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            cards=prep.write_cards(Path(tmp))
            self.assertEqual(len(cards),12)
            for (system,selection,population),path in cards.items():
                values=c.card_values(path)
                self.assertEqual(c.selection_from_card(path),selection)
                self.assertEqual(float(values['MinTrackPtCut']),.5)
                self.assertEqual(float(values['MinPairTrackPt']),2.)
                self.assertEqual(values['WeightExponent'],'1 2')
                self.assertEqual(float(values['HFEnergyShift']),0.)
                self.assertEqual(int(values['DisableTrackPairEfficiencyCorrection']),1)
                self.assertEqual(int(values['TrackEfficiencyVariation']),0)
                self.assertTrue(c.payload_path(system,selection).is_file())

    def test_bad_card_variations_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'bad.input'
            for content in ('Run3TrackSelection Tight\nMaxTrackPtRelativeError .1\nVertexMaxDistance 3\n',
                            'MaxTrackPtRelativeError .05\nVertexMaxDistance 2\n',
                            'MaxTrackPtRelativeError .15\nVertexMaxDistance 5\n'):
                path.write_text(content)
                with self.assertRaises(ValueError): c.selection_from_card(path)

    def test_input_contract_and_checksum(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'events.jsonl';path.write_text('{}\n')
            contract=c.export_contract('OO','Nominal',c.payload_path('OO','Nominal'))
            meta={'tracking_contract':contract,'output_sha256':c.digest(path)}
            Path(str(path)+'.metadata.json').write_text(json.dumps(meta))
            self.assertEqual(c.check_export_metadata([path]),{'OO':contract})
            path.write_text('{}\n{}\n')
            with self.assertRaises(ValueError):c.check_export_metadata([path])
            with self.assertRaises(ValueError):c.validate_contract({'version':'old'})

    def test_mixed_selections_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths=[]
            for selection in ('Nominal','Loose'):
                path=Path(tmp)/(selection+'.jsonl');path.write_text('{}\n');paths.append(path)
                meta={'tracking_contract':c.export_contract('OO',selection,c.payload_path('OO',selection)),
                      'output_sha256':c.digest(path)}
                Path(str(path)+'.metadata.json').write_text(json.dumps(meta))
            with self.assertRaises(ValueError):c.check_export_metadata(paths)

    def test_factor_once_for_both_powers(self):
        a={'pt':4.,'correction':2.};b={'pt':5.,'correction':3.}
        for n in (1,2):self.assertAlmostEqual(pairs.pair_weight(a,b,n,7.,1.4),7*20**n*2*3*1.4)

    def test_duplicate_track_identity_fails_before_counting(self):
        track={'index':1,'pt':3.,'eta':0.,'phi':0.,'correction':1.}
        with self.assertRaises(ValueError):pairs.validate_population([track,dict(track)])

    def test_source_wiring(self):
        source=(c.REPO/'src/EECAnalyzer.cxx').read_text()
        self.assertEqual(source.count('EECRun3Tracking::PassMomentum(trackPt,'),2)
        self.assertEqual(source.count('EECRun3Tracking::PassVertex('),4)
        self.assertIn('EECRun3Tracking::Payload(true, fTrackSelectionVariation)',source)
        self.assertIn('EECRun3Tracking::Payload(false, fTrackSelectionVariation)',source)
        self.assertIn('Run3TrackSelectionVersion',source)
        self.assertIn('if(trackReader->GetTrackCharge(iTrack) == 0) return false;',source)
        self.assertIn('selectedTrackCorrection[firstParticleType[iPairingType]]',source)

    @unittest.skipUnless(shutil.which('root-config'),'ROOT development environment required')
    def test_native_cut_charge_and_real_table_lookups(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'test.cpp';source.write_text(CPP);exe=Path(tmp)/'test'
            flags=shlex.split(subprocess.check_output(['root-config','--cflags','--libs'],text=True))
            command=['g++','-std=c++17','-I'+str(c.REPO/'src'),str(source),
                     str(c.REPO/'src/trackingEfficiency2025OO.cxx'),str(c.REPO/'src/trackingEfficiency2024pp.cxx'),
                     str(c.REPO/'src/TrackingEfficiencyInterface.cxx'),'-o',str(exe)]+flags
            subprocess.run(command,check=True,cwd=c.REPO)
            subprocess.run([str(exe),tmp],check=True,cwd=c.REPO)

    @unittest.skipUnless(shutil.which('root-config'),'ROOT required to open binary tables')
    def test_six_table_structure_audit(self):
        report=prep.audit_tables()
        self.assertEqual(len(report['tables']),6)
        self.assertEqual(report['physics_closure'],'PENDING_MC')
        for row in report['tables']:
            self.assertEqual(row['nonfinite_content'],0)
            self.assertEqual(len(row['sha256']),64)


if __name__=='__main__': unittest.main()

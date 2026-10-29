#!/usr/bin/env python3
"""Fresh-source test of ordinary EEC projection, background processing and reload.
Synthetic fixtures check software only, never detector or physics closure.
"""
import argparse
import math
import shlex
import subprocess
from pathlib import Path
import ROOT


def manager_library(folder):
    repo = Path(__file__).resolve().parents[2]
    library = folder / 'standardManager.so'
    sources = ['EECHistogramManager', 'EECCard', 'AlgorithmLibrary', 'EECBackgroundScale', 'EECSignalToBackgroundUnfoldingScale', 'EECUnfoldConfiguration']
    command = ['g++', '-shared', '-fPIC', '-O0', '-I'+str(repo),
               *[str(repo/'plotting'/(s+'.cxx')) for s in sources],
               str(repo/'src/ConfigurationCard.cxx'), str(repo/'src/EECHistograms.cxx'),
               *shlex.split(subprocess.check_output(['root-config', '--cflags', '--libs'], text=True)),
               '-lEG', '-o', str(library)]
    subprocess.run(command, check=True)
    if ROOT.gSystem.Load(str(library)) < 0:
        raise RuntimeError('Cannot load freshly compiled standard manager')
    ROOT.gInterpreter.AddIncludePath(str(repo))
    ROOT.gInterpreter.Declare('''
#include "src/ConfigurationCard.h"
#include "src/EECHistograms.h"
#include "plotting/EECHistogramManager.h"
#include "plotting/EECUnfoldConfiguration.h"
void standardResponse(TFile* source, const char* path){
  TDirectory::TContext context;
  EECCard card(source);
  const bool oo=card.GetDataType().Contains("OO");
  TFile out(path,"RECREATE");
  card.Write(&out);
  out.cd("JCard");
  TVectorT<float> type(1,1); type[1]=oo ? 10 : 12;
  type.Write("DataType",TObject::kOverwrite);
  TVectorT<float> mode(1,1); mode[1]=3;
  mode.Write("McCorrelationType",TObject::kOverwrite);
  out.cd();
  for(const char* name: {"jetPtUnfoldingMeasured","jetPtUnfoldingTruth","jetPtUnfoldingResponse"})
    out.mkdir(name);
}
void standardFixture(const char* cardName, const char* path){
  ConfigurationCard card(cardName);
  auto* setting=card.GetVector("BackgroundMethods"); (*setting)[1]=6;
  EECHistograms h(&card);h.CreateHistograms();
  // Real OO coordinate is retained; pp keeps its own inherited convention.
  const int type=card.Get("DataType"); const double centrality=(type==11||type==12)?-0.5:35.;
  const double radius=card.Get("JetRadius");
  for(int j=0;j<3;++j){
    for(int n=0;n<2;++n){
      for(int pairing: {0,3,4,5,6,7,8,9,10,11,12}){
        const double factors[]={20,0,0,2,1,4,2,1,2,1,4,2,1};
        double x[]={1.5*radius,60.+40*j,2.5,centrality,double(pairing),0.,double(n)};
        double factor=pairing==0 && j==2 ? 1 : factors[pairing];
        h.fhEnergyEnergyCorrelator->Fill(x,(j+1)*std::pow(15.,n+1)*factor);
      }
    }
  }
  TFile out(path,"RECREATE");h.Write();card.WriteCard(&out);out.Close();
}
''')
    return library


def configure(file, power, processed=False):
    m = ROOT.EECHistogramManager(file)
    m.SetLoadEnergyEnergyCorrelators(True)
    m.SetLoadedWeightExponent(power)
    m.SetCentralityBins(True, 0, ROOT.nullptr, not processed)
    m.SetJetPtBinsEEC(True, 0, ROOT.nullptr, not processed)
    m.SetTrackPtBinsEEC(True, 0, ROOT.nullptr, not processed)
    m.SetCentralityBinRange(0, 0)
    m.SetJetPtBinRangeEEC(0, 2)
    m.SetTrackPtBinRangeEEC(3, 3)
    return m


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch', type=Path)
    args = parser.parse_args()
    args.scratch.mkdir(parents=True, exist_ok=False)
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    manager_library(args.scratch)
    repo=Path(__file__).resolve().parents[2]
    tested=0
    for radius in (4,8):
        for system, prefix in [('OO','cardEECOOInclusive'),('OOMC','cardEECOOMCInclusive'),('pp','cardEECpp5v36TeV'),('ppMC','cardEECpp5v36TeVMC')]:
            raw=args.scratch/f'{system}_R0{radius}_standard.root'
            ROOT.standardFixture(str(repo/f'{prefix}R0{radius}.input'),str(raw))
            response=args.scratch/f'{system}_R0{radius}_response.root'
            response_source=ROOT.TFile(str(raw))
            ROOT.standardResponse(response_source,str(response))
            response_source.Close()
            for power in (1,2):
                source=ROOT.TFile(str(raw))
                # Standard unfolding configuration rejects unqualified legacy fallback.
                card=ROOT.EECCard(source)
                rejected=False
                try:
                    ROOT.EECUnfoldConfiguration(card,0,0,power)
                except Exception:
                    rejected=True
                assert rejected
                configured=ROOT.EECUnfoldConfiguration(card,0,0,power,str(response),7)
                bins=ROOT.EECUnfoldConfiguration.SelectUnfoldedJetPtBins(card)
                assert [(b.first,b.second) for b in bins]==[(40.,80.),(80.,120.),(120.,160.)]
                assert configured.GetNumberOfIterations(ROOT.std.pair('double,double')(-1.,100.),2.)==7
                assert card.GetFirstUnfoldedJetPtBin()==0 and card.GetLastUnfoldedJetPtBin()==2
                m=configure(source,power)
                m.LoadHistograms()
                for j in range(3):
                    h=m.GetHistogramEnergyEnergyCorrelator(0,0,j,3,0,4)
                    factor=1 if j==2 else 20
                    expected=(j+1)*15**power*factor
                    assert h and math.isclose(h.Integral('width'),expected,rel_tol=1e-7)
                    assert h.GetXaxis().GetXmax() >= 2*radius/10-1e-6
                    assert h.GetBinContent(h.FindBin(1.5*radius/10))>0
                projected=args.scratch/f'{system}_R0{radius}_N{power}_projected.root'
                m.Write(str(projected),'RECREATE')
                for method in (0,1):
                    m.SubtractBackground(method,0)
                    for j in range(3):
                        background=m.GetHistogramEnergyEnergyCorrelatorProcessed(0,0,j,3,1)
                        signal=m.GetHistogramEnergyEnergyCorrelatorProcessed(0,0,j,3,2)
                        scale=(j+1)*15**power
                        assert math.isclose(background.Integral('width'),2*scale,rel_tol=1e-7)
                        assert math.isclose(signal.Integral('width'),(-1 if j==2 else 18)*scale,rel_tol=1e-7)
                # The ordinary processed file is exactly the unfolding interface input.
                loaded=ROOT.TFile(str(projected))
                restored=configure(loaded,power,True)
                restored.LoadProcessedHistograms()
                for j in range(3):
                    h=restored.GetHistogramEnergyEnergyCorrelator(0,0,j,3,0,4)
                    expected=(j+1)*15**power*(1 if j==2 else 20)
                    assert h and math.isclose(h.Integral('width'),expected,rel_tol=1e-7)
                del restored
                loaded.Close()
                del m
                source.Close()
                tested+=1
    print(f'PASS: {tested} standard manager OO/pp data/MC radius/power fixtures; three jet bins, mixed/perpendicular categories, negative signal retained, wide support, save/load and unfolding-input schema; software only')


if __name__=='__main__':
    main()

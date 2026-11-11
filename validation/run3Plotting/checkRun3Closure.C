#include "../../plotting/Run3ClosureTools.h"
R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include <iostream>
#include <memory>
// Regression test for the supplied OO GenGen split-centrality fixture.
// No input data are stored in this repository. Supply the local input and plots.
void checkRun3Closure(TString input,TString directory,double low=80.,double high=120.) {
  using namespace Run3Closure;
  std::unique_ptr<TFile> in(Open(input));
  auto manager=std::make_unique<EECHistogramManager>(in.get());
  EECCard& card=*manager->GetCard();CheckCard(card);
  const int c=card.FindBinIndexCentrality(0.,30.),j=card.FindBinIndexJetPtEEC(low,high),t=Track(card,2.);
  Require(c==1 && t==3,"Unexpected indexing for this supplied split fixture");
  Preflight(*manager,c,j,t);
  std::unique_ptr<TH1D> total(Read(*manager,c,j,t,EECHistograms::kSameJetPair,EECHistograms::knSubeventCombinations,"check_total"));
  std::unique_ptr<TH1D> truth(Read(*manager,c,j,t,EECHistograms::kSameJetPair,EECHistograms::kPythiaPythia,"check_signal"));
  std::unique_ptr<TH1D> bt(TrueBackground(*manager,c,j,t));
  TString stem=directory+"/"+Form("C0-30_J%.0f-%.0f_T2.0_n1_signal_legacyC1_raw_",low,high);
  TFile independent(stem+"uncorrelated.root"), shared(stem+"sharedSignalApprox.root");
  Require(!independent.IsZombie() && !shared.IsZombie(),"Missing comparison outputs");
  int bins=0;
  for(int k=0;k<3;++k) {
    auto* a=dynamic_cast<TH1D*>(independent.Get(Form("curve%d",k)));
    auto* b=dynamic_cast<TH1D*>(shared.Get(Form("curve%d",k)));
    auto* ra=dynamic_cast<TH1D*>(independent.Get(Form("ratio%d",k)));
    auto* rb=dynamic_cast<TH1D*>(shared.Get(Form("ratio%d",k)));
    Compatible(a,b);Compatible(ra,rb);
    std::unique_ptr<TH1D> estimate(k ? Background(*manager,c,j,t,k-1,false) : nullptr);
    for(int q=1;q<=a->GetNbinsX();++q) {
      Require(a->GetBinContent(q)==b->GetBinContent(q),"Upper central value changed");
      Require(a->GetBinError(q)==b->GetBinError(q),"Upper uncertainty changed");
      Require(ra->GetBinContent(q)==rb->GetBinContent(q),"Ratio central value changed");
      if(!k || !(truth->GetBinContent(q)>0.)) continue;
      const double s=truth->GetBinContent(q),x=total->GetBinContent(q)-estimate->GetBinContent(q),r=x/s;
      const double vs=std::pow(truth->GetBinError(q),2),vt=std::pow(total->GetBinError(q),2),ve=std::pow(estimate->GetBinError(q),2),vb=std::pow(bt->GetBinError(q),2);
      Require(Close(a->GetBinContent(q),x),"Legacy arithmetic changed");
      Require(Close(ra->GetBinError(q)*ra->GetBinError(q),(vt+ve+r*r*vs)/(s*s),1.e-12),"Independent ratio formula mismatch");
      Require(Close(rb->GetBinError(q)*rb->GetBinError(q),(vb+ve+(r-1)*(r-1)*vs)/(s*s),1.e-12),"Shared-signal formula mismatch");
      ++bins;
    }
  }
  bool caught=false;try {Track(card,2.1);} catch(const std::exception&) {caught=true;}
  Require(caught,"Unavailable threshold silently accepted");
  Require(card.FindBinIndexJetPtEEC(80.,121.)<0,"Missing jet boundary accepted");
  std::cout<<"RUN3_NUMERICAL_PASS C="<<c<<" J="<<j<<" T="<<t<<" checkedRatioBins="<<bins
    <<" identicalCentralValues=true identicalUpperErrors=true legacyArithmetic=true exactThresholdGuard=true"<<std::endl;
}

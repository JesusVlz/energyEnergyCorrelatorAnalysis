// Synthetic software tests only. These are NOT collision data or physics plots.
#include "../../plotting/trackClosurePlotterRun3.C"
#include <TObjString.h>
#include <map>
#include <string>

namespace TrackingFixture {
void Vector(TDirectory& dir,const char* name,const std::vector<double>& x) {
  dir.cd();TVectorT<float> v(1,x.size());for(size_t i=0;i<x.size();++i)v[i+1]=x[i];v.Write(name);
}
void Make(const TString& name,int mode,int power,bool badJets=false) {
  TFile f(name,"RECREATE");auto* card=f.mkdir("JCard");
  std::map<std::string,double> vals={
    {"DataType",10},{"McCorrelationType",double(mode)},{"MatchJets",0},{"TriggerSelection",0},
    {"JetType",1},{"JetAxis",1},{"JetEtaCut",1.6},{"MinJetPtCut",40},{"MaxJetPtCut",160},
    {"CutBadPhi",0},{"MinMaxTrackPtFraction",0.01},{"MaxMaxTrackPtFraction",0.98},
    {"JetUncertainty",0},{"TrackEtaCut",2.4},{"MinTrackPtCut",0.5},{"MaxTrackPtCut",500},
    {"MaxTrackPtRelativeError",0.1},{"VertexMaxDistance",3},{"SubeventCut",2},
    {"TrackEfficiencyVariation",0},{"JetPtWeight",0},{"DisableTrackPairEfficiencyCorrection",1},
    {"ZVertexCut",15},{"LowPtHatCut",-15},{"HighPtHatCut",5020},{"MultiplicityMode",0},
    {"JetRadius",0.4},{"WeightExponent",double(power)},{"SkipCovarianceMatrix",1},
    {"MinJetPtUnfoldingReco",20},{"MinJetPtUnfoldingTruth",10},{"BackgroundMethods",0},
    {"MegaSkimMixing",0},{"FirstUnfoldedJetPtBin",0},{"LastUnfoldedJetPtBin",2}};
  for(auto& v:vals)Vector(*card,v.first.c_str(),{v.second});
  Vector(*card,"CentralityBinEdges",{0,30,50,90,100});
  Vector(*card,"JetPtBinEdgesEEC",{40,80,120,160});
  Vector(*card,"TrackPtBinEdgesEEC",{0.5,1,1.5,2,2.5,3,3.5,4,300});
  Vector(*card,"TrackPtBinEdges",{0.7,1,2,3,4,8,12,300});
  Vector(*card,"PtHatBinEdges",{0,30,50,80,120,170,220,280,370,460});
  card->cd();TObjString("SYNTHETIC_SOFTWARE_TEST").Write("GitHash");
  auto* mult=f.mkdir("multiplicity");auto* jets=f.mkdir("inclusiveJet");
  auto* track=f.mkdir("track");auto* raw=f.mkdir("trackUncorrected");
  auto* eec=f.mkdir("energyEnergyCorrelator")->mkdir("SignalConePair");
  for(int c=0;c<3;++c) {
    mult->cd();TH1D m(Form("multiplicity_C%d",c),"",20,0,200);
    m.SetBinContent(2,1000./10.);m.SetBinError(2,std::sqrt(1000.)/10.);m.Write();
    jets->cd();TH1D j(Form("inclusiveJetPt_C%d",c),"",200,0,200);
    for(int b=41;b<=160;++b){j.SetBinContent(b,100.+b+(badJets?2.:0.));j.SetBinError(b,1.);}j.Write();
    for(int t=0;t<2;++t) {
      auto* dir=t==0?track:raw;const char* prefix=t==0?"track":"trackUncorrected";
      dir->cd();TH1D p(Form("%sPt_C%d",prefix,c),"",1000,0,100);
      const double scale=mode==1?1.:(t==0?0.96:0.60);
      for(int b=1;b<=p.GetNbinsX();++b){double v=5000.*std::exp(-p.GetBinCenter(b)/7.);p.SetBinContent(b,v*scale);p.SetBinError(b,std::sqrt(v)*scale);}p.Write();
      for(int a=0;a<2;++a){TH1D h(Form("%s%s_C%dT2",prefix,a==0?"Eta":"Phi",c),"",40,a==0?-2.4:-3.14159265,a==0?2.4:3.14159265);
        for(int b=1;b<=40;++b){h.SetBinContent(b,2000.*scale);h.SetBinError(b,20.*scale);}h.Write();}
    }
    eec->cd();double edges[31];for(int b=0;b<=30;++b)edges[b]=0.008*std::pow(0.4/0.008,b/30.);
    for(int jb=0;jb<3;++jb){TH1D h(Form("energyEnergyCorrelatorSignalConePair_C%dT3J%d",c,jb),"",30,edges);
      for(int b=1;b<=30;++b){double v=1.e5*(jb+1)*std::exp(-h.GetBinCenter(b)*4.);double s=mode==1?1.:0.95;
        h.SetBinContent(b,v*s);h.SetBinError(b,std::sqrt(v)*s);}h.Write();}
  }
  f.Close();
}
}
void testTrackClosurePlotterRun3(const char* output="tracking_fixture") {
  using namespace Run3TrackingPlot;
  gSystem->mkdir(output,true);
  TString base(output);TrackingFixture::Make(base+"/rr1.root",0,1);TrackingFixture::Make(base+"/rg1.root",1,1);
  TrackingFixture::Make(base+"/rr2.root",0,2);TrackingFixture::Make(base+"/rg2.root",1,2);
  TrackingFixture::Make(base+"/badJets.root",0,1,true);
  trackClosurePlotterRun3(base+"/rr1.root",base+"/rg1.root",base+"/plots1",0,30,2,"Nominal",true,true);
  trackClosurePlotterRun3(base+"/rr2.root",base+"/rg2.root",base+"/plots2",0,30,2,"Nominal",false,true);
  TFile check(base+"/plots1/Nominal_C0-30_n1_EEC_J40-80_T2.root","READ");
  auto* graph=dynamic_cast<TGraphErrors*>(check.Get("corrected_over_gen"));Require(graph,"Missing ratio output");
  Require(graph->GetN()==30,"Wrong number of EEC points");
  for(int p=0;p<graph->GetN();++p){double x,y;graph->GetPoint(p,x,y);Require(Near(y,0.95,1.e-12),"Central values changed");}
  TH1D n("n","",2,0,2),d("d","",2,0,2);n.SetBinContent(1,6);n.SetBinError(1,2);d.SetBinContent(1,3);d.SetBinError(1,1);
  auto ratio=Ratio(n,d,"testRatio");Require(ratio.graph->GetN()==1,"Zero denominator was not omitted");
  Require(Near(ratio.value[1],2.) && Near(ratio.error[1],std::sqrt(8.)/3.),"Ratio error formula failed");
  bool caught=false;try{TFile a(base+"/rr1.root"),b(base+"/rg2.root");CheckInputs(a,b);}catch(const std::exception&){caught=true;}
  Require(caught,"Different powers were accepted");
  caught=false;try{TFile a(base+"/rg1.root"),b(base+"/rr1.root");CheckInputs(a,b);}catch(const std::exception&){caught=true;}
  Require(caught,"Wrong mode order was accepted");
  caught=false;try{TFile a(base+"/badJets.root"),b(base+"/rg1.root");auto ma=std::make_unique<EECHistogramManager>(&a);auto mb=std::make_unique<EECHistogramManager>(&b);CheckJets(*ma,*mb,0,0);}catch(const std::exception&){caught=true;}
  Require(caught,"Jet-population mismatch was accepted");
  std::cout<<"SYNTHETIC_TRACKING_TESTS_PASS\n";
}

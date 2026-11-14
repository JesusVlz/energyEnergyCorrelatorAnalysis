// Synthetic fixtures only, never collision-data evidence.
#include "../../plotting/Run3ResponsePlotTools.h"
R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include <TObjString.h>
#include <map>
namespace Run3ResponseFixture {
void Vector(TDirectory& d,const char* name,std::vector<double> x){d.cd();TVectorT<float> v(1,x.size());for(size_t i=0;i<x.size();++i)v[i+1]=x[i];v.Write(name);}
void Make(const TString& path,int mode,int power,bool integrated=false){
  TFile f(path,"RECREATE");auto* c=f.mkdir("JCard");
  std::map<std::string,double> scalars={{"DataType",10},{"McCorrelationType",double(mode)},{"MatchJets",0},{"TriggerSelection",0},{"JetType",1},{"JetAxis",1},{"JetEtaCut",1.6},{"MinJetPtCut",40},{"MaxJetPtCut",5020},{"JetRadius",.4},{"WeightExponent",double(power)},{"SubeventCut",2},{"BackgroundMethods",0},{"SkipCovarianceMatrix",0},{"MinJetPtUnfoldingReco",40},{"MinJetPtUnfoldingTruth",5},{"JetPtWeight",0},{"DisableTrackPairEfficiencyCorrection",1},{"MultiplicityMode",0},{"TrackEtaCut",2.4},{"MinTrackPtCut",mode==1?2.:.5},{"MaxTrackPtCut",500}};
  for(auto& kv:scalars)Vector(*c,kv.first.c_str(),{kv.second});
  Vector(*c,"CentralityBinEdges",integrated?std::vector<double>{-1,100}:std::vector<double>{-1,0,30,50,90,100});
  Vector(*c,"JetPtBinEdgesEEC",{5,10,20,40,80,120,160,200,300,500,5020});
  Vector(*c,"TrackPtBinEdgesEEC",{.5,1,1.5,2,2.5,3,3.5,4,300});Vector(*c,"TrackPtBinEdges",{.7,1,2,3,4,8,12,300});Vector(*c,"PtHatBinEdges",{0,30,50,80,120,170,220,280,370,460});c->cd();TObjString("SYNTHETIC_ONLY").Write("GitHash");
  const int cent=integrated?0:1;
  auto write=[&](TH1& h,const TString& dir,const TString& name){TDirectory* d=f.GetDirectory(dir);if(!d)d=f.mkdir(dir);d->cd();h.Write(name);};
  double a[]={0,.008,.016,.04,.1,.2,.4,.8},p[]={4,6,9,16,25,100,1000,90000};
  for(int j=3;j<=5;++j){
    TH2D angle("angle","",7,a,7,a),mom("mom","",7,p,7,p);TH1D ratio("ratio","",40,.025,2.025);
    for(int x=1;x<=7;++x)for(int y=1;y<=7;++y){const double v=(x==7)?0.:(x==y?100.:2.);angle.SetBinContent(x,y,v);angle.SetBinError(x,y,std::sqrt(v));mom.SetBinContent(x,y,v*power);mom.SetBinError(x,y,std::sqrt(v*power));}
    angle.SetBinContent(0,1,3);angle.SetBinError(0,1,std::sqrt(3.));
    for(int b=1;b<=40;++b){double v=100.*std::exp(-std::pow((ratio.GetBinCenter(b)-1.)/.12,2));ratio.SetBinContent(b,v);ratio.SetBinError(b,std::sqrt(v));}
    write(angle,"particleDeltaRResponseMatrix",Form("particleDeltaRResponseMatrix_C%dJ%dT3",cent,j));
    write(mom,"particlePtResponseMatrix",Form("particlePtResponseMatrix_C%dJ%dT3",cent,j));
    write(ratio,"particlePairPtClosure",Form("particlePairPtClosure_C%dJ%dT3",cent,j));
  }
  double legacy[]={50,80,100,120,140,160,200,500},native[]={5,10,20,40,80,120,160,200,300,500,5020};
  TH2D jet("jet","",7,legacy,7,legacy),nat("native","",10,native,10,native),eec("eec","",20,0,8,20,0,8),cov("cov","",20,0,8,20,0,8);
  for(TH2D* h:{&jet,&nat,&eec,&cov})for(int x=1;x<=h->GetNbinsX();++x)for(int y=1;y<=h->GetNbinsY();++y){const double v=x==y?100.:1.;h->SetBinContent(x,y,v);h->SetBinError(x,y,std::sqrt(v));}
  write(jet,"jetPtResponseMatrix",Form("jetPtResponseMatrix_C%d",cent));
  write(nat,"oneDimensionalJetPtUnfoldingResponse",Form("oneDimensionalJetPtUnfoldingResponse_C%d",cent));
  write(eec,"jetPtUnfoldingResponse",Form("jetPtUnfoldingResponse_C%dT3",cent));
  write(cov,"jetPtUnfoldingCovariance",Form("jetPtUnfoldingCovariance_C%dT3",cent));
  for(const char* name:{"jetPtUnfoldingMeasured","jetPtUnfoldingTruth"}){TH1D h("marginal","",20,0,8);for(int b=1;b<=20;++b){h.SetBinContent(b,125.);h.SetBinError(b,std::sqrt(125.));}write(h,name,Form("%s_C%dT3",name,cent));}
  f.Close();
}
}
void testRun3ResponsePlotters(const char* output="/tmp/run3_response_fixture"){
  using namespace Run3ResponsePlot;
  try {
    Require(gSystem->AccessPathName(output),"Fixture output already exists");gSystem->mkdir(output,true);gROOT->SetBatch(true);gStyle->SetOptStat(0);
    TString base(output);Run3ResponseFixture::Make(base+"/pair1.root",0,1);Run3ResponseFixture::Make(base+"/pair2.root",0,2,true);Run3ResponseFixture::Make(base+"/response1.root",1,1);Run3ResponseFixture::Make(base+"/response2.root",1,2,true);
    TH2D raw("rawtest","",3,0,3,2,0,2);raw.SetBinContent(1,1,2);raw.SetBinContent(1,2,6);raw.SetBinError(1,1,2);auto valid=Copy(&raw);auto conditional=ConditionalX(raw,*valid);
    Require(Near(conditional->GetBinContent(1,1),.25)&&Near(conditional->GetBinContent(1,2),.75)&&Near(conditional->GetBinError(1,1),.25),"Wrong conditional direction or scale");Require(valid->GetBinContent(2,1)==0.,"Empty conditional slice not masked");Require(raw.GetBinContent(1,1)==2.,"Input mutated");
    raw.SetBinContent(3,1,-1.);bool caught=false;try{ConditionalX(raw,*valid);}catch(...){caught=true;}Require(caught,"Signed conditional input accepted");
    PairPlots(base+"/pair1.root",base+"/pairs1",0,30,2);PairPlots(base+"/pair2.root",base+"/pairs2",-1,100,2);
    JetPlots(base+"/pair1.root",base+"/jets",0,30);
    UnfoldingPlots(base+"/response1.root",base+"/unfold1",0,30,2);UnfoldingPlots(base+"/response2.root",base+"/unfold2",-1,100,2);
    caught=false;try{Input in(base+"/pair1.root",base+"/invalid",0,30,2,1);}catch(...){caught=true;}Require(caught,"Wrong reconstruction mode accepted");
    caught=false;try{Input in(base+"/pair1.root",base+"/invalid",0,10,2,0);}catch(...){caught=true;}Require(caught,"Unavailable centrality accepted");
    caught=false;try{Input in(base+"/pair1.root",base+"/invalid",0,30,2.1,0);}catch(...){caught=true;}Require(caught,"Unavailable threshold accepted");
    TFile check(base+"/jets/legacyJetResponse.root","READ");auto* h=dynamic_cast<TH2D*>(check.Get("raw"));Require(h&&Exact(*h->GetXaxis(),40,80).first<0,"Legacy 50 GeV floor not preserved");Require(Exact(*h->GetXaxis(),80,120).first>0,"Valid legacy interval rejected");
    Require(gSystem->AccessPathName(base+"/jets/legacyJetSlice_fixedGen40-80.root"),"Invented legacy 40-80 plot");
    Require(!gSystem->AccessPathName(base+"/unfold1/nativeJetSlice_fixedGen40-80.root"),"Native 40-80 plot missing");
    std::cout<<"SYNTHETIC_RUN3_RESPONSE_TESTS_PASS\n";
  }catch(const std::exception& e){std::cerr<<"SYNTHETIC_RUN3_RESPONSE_TESTS_FAIL "<<e.what()<<'\n';gSystem->Exit(1);}
}

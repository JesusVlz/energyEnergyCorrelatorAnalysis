// Standard EEC pipeline regression. Private access is confined to this fixture.
#include "../../src/EECAnalyzer.h"
#include "../../src/OOEventIdentity.h"
#include <TMemFile.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
void Require(bool ok, const char* message){ if(!ok) throw std::runtime_error(message); }
void Set(ConfigurationCard& c,const char* key,float v){
  if(!c.HasKey(key)) c.ReadInputLine((std::string(key)+" "+std::to_string(v)).c_str());
  auto* x=c.GetVector(key);(*x)[x->GetLwb()]=v;
}
void SetString(ConfigurationCard& c,const char* key,const char* value){
  if(!c.HasKey(key)) c.ReadInputLine((std::string(key)+" "+value).c_str());
  else c.fValueString[c.GetTVectorIndex(key)]=value;
}
std::unique_ptr<EECAnalyzer> Analyzer(ConfigurationCard& card){
  // Pairing-only fixtures provide no native jet payload and explicitly turn
  // off the independently tested PF Jet ID/event-veto integration.
  Set(card,"Run3JetID",0);Set(card,"ApplyJetVeto",0);
  auto a=std::make_unique<EECAnalyzer>(std::vector<TString>{},&card,1,-1);
  a->fJetReader=new HighForestReader();
  a->fTrackReader=new HighForestReader();
  a->fMcCorrelationType=EECAnalyzer::kRecoGen;
  return a;
}
void Pair(EECAnalyzer& a,int hiBin,double weight,double pt=2.5){
  a.fTrackReader->fHiBin=hiBin;
  a.fTotalEventWeight=weight;
  std::vector<double> p[EECHistograms::knJetConeTypes],e[EECHistograms::knJetConeTypes],f[EECHistograms::knJetConeTypes];
  std::vector<int> sub[EECHistograms::knJetConeTypes];
  // Both mixed and perpendicular controls exercise the ordinary pairing categories.
  for(int cone: {0,2,3,4,5}){p[cone]={pt,3};e[cone]={0,.1};f[cone]={0,0};sub[cone]={0,0};}
  a.CalculateEnergyEnergyCorrelator(p,e,f,sub,60);
}
double Yield(THnSparse* h,int power,double low,double high){
  h->GetAxis(6)->SetRange(power,power);
  h->GetAxis(4)->SetRange(1,1);
  h->GetAxis(3)->SetRange(h->GetAxis(3)->FindBin(low),h->GetAxis(3)->FindBin(high)-1);
  std::unique_ptr<TH1D> p(h->Projection(0));
  double result=p->Integral(0,p->GetNbinsX()+1);
  for(int i=0;i<7;++i)h->GetAxis(i)->SetRange(0,0);
  return result;
}
}
int main(){
 try{
  TH1::AddDirectory(false);
  ConfigurationCard broad("cardEECOOMCInclusiveR04.input"),split("cardEECOOMC.input");
  Set(broad,"BackgroundMethods",6);Set(split,"BackgroundMethods",6);
  Set(split,"SkipCovarianceMatrix",1);Set(split,"MinPairTrackPt",2);
  auto a=Analyzer(broad),b=Analyzer(split);
  // RecoGen keeps synthetic boundary coordinates independent of payload-domain behavior.
  a->fMcCorrelationType=b->fMcCorrelationType=EECAnalyzer::kRecoGen;
  const int hibins[]={-3,-2,-1,0,59,60,99,100,179,180,199,200,202};
  for(int h:hibins){
   const double weight=(h==100?-2.:1.5)*a->GetVzWeight(.5)*a->GetCentralityWeight(h);
   Require(weight==((h==100?-2.:1.5)*b->GetVzWeight(.5)*b->GetCentralityWeight(h)),"Binning changed event correction weight");
   if(h>=0 && h<180){Pair(*a,h,weight);Pair(*b,h,weight);}
   else {a->fTrackReader->fHiBin=h;b->fTrackReader->fHiBin=h;}
   Require(a->fTrackReader->GetStoredHiBin()==h,"Stored hiBin changed");
   Require(a->fTrackReader->GetCentrality()==h/2.,"OO acquired a pp centrality placeholder");
   for(int stage=0;stage<4;++stage)a->GetHistograms()->FillEventCoverage(h,stage,true);
   a->GetHistograms()->FillSelectedEvent(h,100,.5,weight);
  }
  auto* raw=a->GetHistograms()->fhEnergyEnergyCorrelator;
  auto* other=b->GetHistograms()->fhEnergyEnergyCorrelator;
  Require(raw->GetNdimensions()==7,"Standard EEC schema changed");
  for(int n: {1,2}){
   double x=Yield(raw,n,-1,100),y=Yield(other,n,0,90);
   Require(std::abs(x-y)<1e-5*std::max(1.,std::abs(y)),"Broad and split raw yields differ over common [0,90) range");
  }
  auto* qa=a->GetHistograms()->fhEventCoverage;
  Require(qa->GetBinContent(3,4)==3 && qa->GetBinContent(4,4)==2,"Invalid stored coordinates were not reported separately");
  auto* axis=raw->GetAxis(3);
  Require(axis->FindBin(-1)==1 && axis->FindBin(99.5)==1 && axis->FindBin(100)==2 && axis->FindBin(-1.5)==0,"Ordinary centrality boundaries/flows changed");
  a->fDoMixedCone=a->fDoPerpendicularCone=false;
  raw->Reset();Pair(*a,20,1,2.0);Require(raw->GetEntries()==0,"Exact pair threshold accepted");
  Pair(*a,20,1,2.000001);Require(raw->GetEntries()>0,"Above-threshold pair rejected");
  // Reconstructed-track lookup factors are unchanged by output binning, too.
  raw->Reset();other->Reset();double expected[2]={0,0};
  for(int h: {20,60,100}){
    const double w=2*a->GetVzWeight(.5)*a->GetCentralityWeight(h);
    a->fMcCorrelationType=b->fMcCorrelationType=EECAnalyzer::kRecoReco;
    const double correction=a->GetTrackEfficiencyCorrection(2.5,0,h)*a->GetTrackEfficiencyCorrection(3,.1,h);
    a->fMcCorrelationType=b->fMcCorrelationType=EECAnalyzer::kRecoReco;
    Pair(*a,h,w);Pair(*b,h,w);
    a->fMcCorrelationType=b->fMcCorrelationType=EECAnalyzer::kRecoGen;
    for(int n=1;n<=2;++n)expected[n-1]+=w*correction*std::pow(7.5,n);
  }
  for(int n=1;n<=2;++n){
    const double yield=Yield(raw,n,-1,100);
    Require(std::abs(yield-expected[n-1])<1e-5*std::max(1.,std::abs(yield)),"Standard pair weight changed or factors repeated");
    Require(std::abs(yield-Yield(other,n,0,90))<1e-5*std::max(1.,std::abs(yield)),"Reco pair corrections depend on output binning");
  }
  TMemFile file("standardOOFixture.root","RECREATE");a->GetHistograms()->Write();
  Require(file.Get("energyEnergyCorrelator") && file.Get("eventQA/eventCoverage"),"Standard spectrum or coverage QA not saved");
  Require(!file.Get("inclusive"),"Retired raw-only namespace recreated");
  // pp keeps its existing centrality convention, standard corrections and ordinary histogram.
  ConfigurationCard pp("cardEECpp5v36TeVMCR04.input");Set(pp,"BackgroundMethods",0);
  auto p=Analyzer(pp);p->fMcCorrelationType=EECAnalyzer::kRecoGen;Pair(*p,-1,2);
  Require(p->GetHistograms()->fhEnergyEnergyCorrelator->GetEntries()==2,"pp standard pairing regression");
  // A broad output interval still uses actual event hiBin/vz for pool matching.
  a->fnEventsInMixingFile=6;a->fMixingStartIndex=0;
  a->fMixedEventVz={0,0,0,0,0,0};a->fMixedEventHiBin={20,99,20,60,20,20};
  a->fMixedEventEventNumber={123,999,123,44,555,666};
  if(!a->fRng)a->fRng=new TRandom3(1);
  std::vector<int> partners;
  a->FindHiBinMatchedEvents(partners,2,0,20,123,0,false);
  Require(partners==std::vector<int>({4,5}),"Integrated output disabled event-by-event activity matching or allowed self mixing");
  // Ordinary Run-3 analysis requires an explicit usable pool, even for a broad interval.
  ConfigurationCard noPool("validation/run3Workflow/cards/OOMC_R04_diagnostic.input");
  Set(noPool,"BackgroundMethods",6);
  Set(noPool,"Run3JetID",0);Set(noPool,"ApplyJetVeto",0);
  SetString(noPool,"MixingFileList","/tmp/oo-eec-deliberately-unavailable-pool-list.txt");
  bool missingRejected=false;
  try{EECAnalyzer missing({},&noPool,1,-1);missing.RunAnalysis();}
  catch(const std::runtime_error& error){missingRejected=std::string(error.what()).find("MixingFileList")!=std::string::npos;}
  Require(missingRejected,"Missing configured pool silently fell back to signal input");
  const auto key=oo_eec::EventKey("OO2025_MC",1,2,3);
  Require(key==UINT64_C(12337214294153387831) && oo_eec::EventFold(key)==1,"Physical-event identity/fold changed");
  const unsigned replicas[]={1,1,0,2,0,0};
  for(unsigned i=0;i<6;++i)Require(oo_eec::ReplicaWeight(key,i)==replicas[i],"Replica parity changed");
  std::cout<<"PASS: standard OO/pp EEC; binning-invariant weights/yields; centrality boundaries/invalid QA; strict pair threshold; ordinary pairing categories; physical identity\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

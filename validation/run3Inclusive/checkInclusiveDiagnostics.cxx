// Built by run_tests.py with -fno-access-control for this fixture only.
#include "../../src/EECAnalyzer.h"
#include "../../src/OOEventIdentity.h"
#include <TMemFile.h>
#include <TNamed.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
void Require(bool ok, const char* message){ if(!ok) throw std::runtime_error(message); }
void Set(ConfigurationCard& card, const char* key, float value){
  TVector* entry = card.GetVector(key); (*entry)[entry->GetLwb()] = value;
}
void Configure(ConfigurationCard& card){
  card.ReadInputLine("InclusiveMode 1");
  Set(card,"BackgroundMethods",0);
  Set(card,"SkipCovarianceMatrix",1);
}
std::unique_ptr<EECAnalyzer> Fill(ConfigurationCard& card){
  auto analyzer = std::make_unique<EECAnalyzer>(std::vector<TString>{}, &card, 1, -1);
  analyzer->fTrackReader = new HighForestReader();
  // Own through fJetReader, as in the normal one-reader production path.
  analyzer->fJetReader = analyzer->fTrackReader;
  std::vector<double> pt[EECHistograms::knJetConeTypes], eta[EECHistograms::knJetConeTypes], phi[EECHistograms::knJetConeTypes];
  std::vector<int> subevent[EECHistograms::knJetConeTypes];
  pt[0]={2.5,3}; eta[0]={0,0.1}; phi[0]={0,0}; subevent[0]={0,0};
  const int hiBin[]={20,190,-1,200};
  const double weights[]={2,3,5,7};
  for(int i=0;i<4;++i){
    analyzer->fTrackReader->fHiBin=hiBin[i];
    Require(analyzer->fTrackReader->GetStoredHiBin()==hiBin[i],"Stored hiBin was remapped");
    analyzer->fInclusiveEventWeight=weights[i];
    analyzer->fTotalEventWeight=12345; // Must never enter the dedicated raw sum.
    auto* hist=analyzer->GetHistograms();
    for(int stage=0;stage<4;++stage) hist->FillInclusiveCoverage(hiBin[i],stage,true);
    hist->FillInclusiveSelectedEvent(hiBin[i],100+i,0.5,weights[i]);
    hist->fhInclusiveJetExposure->Fill(60,weights[i]);
    analyzer->CalculateEnergyEnergyCorrelator(pt,eta,phi,subevent,60);
  }
  return analyzer;
}
}
int main(){
  try {
    TH1::AddDirectory(false);
    ConfigurationCard card("cardEECOO.input"); Configure(card);
    auto analyzer=Fill(card);
    auto* hist=analyzer->GetHistograms();
    analyzer->fIsRealData=false;
    analyzer->fVzWeight=2;
    analyzer->fJetReader->fEventWeight=3;
    Require(analyzer->GetInclusiveDiagnosticEventWeight()==6,"Actual stored MC generator weight was lost or applied twice");
    analyzer->fIsRealData=true;
    Require(analyzer->GetInclusiveDiagnosticEventWeight()==2,"Data diagnostic acquired an MC generator weight");
    Require(hist->fhInclusiveRawEEC->GetNdimensions()==6,"Inclusive raw schema acquired an overlapping centrality axis");
    Require(hist->fhInclusiveSelectedActivityWeighted->GetSumOfWeights()==17,"Accepted weighted events did not add");
    Require(std::abs(std::pow(hist->fhInclusiveJetExposure->GetBinError(1),2)-87)<1e-10,"Jet sum of weights squared lost");
    for(int category=0;category<4;++category) Require(hist->fhInclusiveEventCoverage->GetBinContent(category+1,4)==1,"Peripheral or invalid event coverage lost");
    const double coordinates1[6]={0.1,60,2.5,0,0,0}, coordinates2[6]={0.1,60,2.5,0,0,1};
    Require(std::abs(hist->fhInclusiveRawEEC->GetBinContent(hist->fhInclusiveRawEEC->GetBin(coordinates1))-127.5)<1e-10,"n=1 raw pair yield did not sum event weights once");
    Require(std::abs(hist->fhInclusiveRawEEC->GetBinContent(hist->fhInclusiveRawEEC->GetBin(coordinates2))-956.25)<1e-10,"n=2 raw pair yield was not independently constructed");
    Require(std::abs(std::pow(hist->fhInclusiveRawEEC->GetBinError(hist->fhInclusiveRawEEC->GetBin(coordinates1)),2)-4893.75)<1e-9,"Raw pair Sumw2 lost");
    // Changing disjoint diagnostic edges leaves the accepted raw event sum invariant.
    ConfigurationCard alternate("cardEECOO.input"); Configure(alternate);
    TVector* edges=alternate.GetVector("CentralityBinEdges");
    (*edges)[edges->GetLwb()+1]=10; (*edges)[edges->GetLwb()+2]=80; (*edges)[edges->GetLwb()+3]=100;
    auto other=Fill(alternate);
    Require(other->GetHistograms()->fhInclusiveRawEEC->GetBinContent(other->GetHistograms()->fhInclusiveRawEEC->GetBin(coordinates1))==127.5,"Diagnostic split bins changed the inclusive weighted sum");
    // Exact boundary control: a bin beginning at 2 GeV alone implements >=2,
    // so the measurement threshold must be applied to the pair before binning.
    other->GetHistograms()->fhInclusiveRawEEC->Reset();
    other->fInclusiveEventWeight=1;
    std::vector<double> boundaryPt[EECHistograms::knJetConeTypes], boundaryEta[EECHistograms::knJetConeTypes], boundaryPhi[EECHistograms::knJetConeTypes];
    std::vector<int> boundarySubevent[EECHistograms::knJetConeTypes];
    boundaryPt[0]={2.0,3.0}; boundaryEta[0]={0,0.1}; boundaryPhi[0]={0,0}; boundarySubevent[0]={0,0};
    other->CalculateEnergyEnergyCorrelator(boundaryPt,boundaryEta,boundaryPhi,boundarySubevent,60);
    Require(other->GetHistograms()->fhInclusiveRawEEC->GetEntries()==0,"Exact-2 GeV track entered the strict inclusive measurement");
    boundaryPt[0][0]=2.000001;
    other->CalculateEnergyEnergyCorrelator(boundaryPt,boundaryEta,boundaryPhi,boundarySubevent,60);
    Require(other->GetHistograms()->fhInclusiveRawEEC->GetEntries()==2,"Track strictly above 2 GeV did not fill both powers");
    const double boundaryCoordinates[6]={0.1,60,2.000001,0,0,0};
    Require(std::abs(other->GetHistograms()->fhInclusiveRawEEC->GetBinContent(other->GetHistograms()->fhInclusiveRawEEC->GetBin(boundaryCoordinates))-6.000003)<1e-10,"Above-threshold pair weight changed");
    TMemFile output("inclusiveDiagnosticFixture.root","RECREATE"); hist->Write();
    Require(output.Get("inclusive/inclusiveRawEEC")!=nullptr,"Dedicated inclusive output namespace missing");
    auto* threshold=dynamic_cast<TNamed*>(output.Get("inclusive/TrackSelection"));
    Require(threshold && std::string(threshold->GetTitle())=="strict pT > 2 GeV; input threshold remains configured separately","Strict measurement threshold provenance missing");
    auto* state=dynamic_cast<TNamed*>(output.Get("inclusive/CorrectionState"));
    Require(state && std::string(state->GetTitle()).find("BLOCKED")!=std::string::npos,"Uncalibrated diagnostic lost its blocked correction metadata");
    ConfigurationCard corrected("cardEECOO.input"); corrected.ReadInputLine("InclusiveMode 2");
    bool rejected=false;
    try { EECAnalyzer forbidden({},&corrected,1,-1); }
    catch(const std::invalid_argument& e){ rejected=std::string(e.what()).find("BLOCKED")!=std::string::npos; }
    Require(rejected,"Unsupported corrected inclusive request was accepted");
    const auto key=oo_eec::EventKey("OO2025_MC",1,2,3);
    Require(key==UINT64_C(12337214294153387831),"Physical-event key changed");
    Require(oo_eec::EventFold(key)==1,"Common physical-event fold changed");
    const unsigned int expected[]={1,1,0,2,0,0};
    for(unsigned int i=0;i<6;++i) Require(oo_eec::ReplicaWeight(key,i)==expected[i],"Common physical-event replica changed");
    Require(oo_eec::EventKey("ppRef2024_MC",1,2,3)!=key,"Independent collision samples share identity");
    std::cout<<"PASS: weighted inclusive raw sums, Sumw2, peripheral/invalid audit, split independence, strict pT>2 boundary, correction gate, physical identity\n";
  } catch(const std::exception& error){ std::cerr<<"FAIL: "<<error.what()<<'\n'; return 1; }
}

// Build with -fno-access-control for this synthetic regression only.
#include "../../src/EECAnalyzer.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH3D.h>
#include <TSystem.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
void Need(bool pass,const char* why) { if(!pass) throw std::runtime_error(why); }
template<class F> void Reject(F f,const char* why) { try {f();} catch(const std::exception&) {return;} throw std::runtime_error(why); }
void Set(ConfigurationCard& c,const char* key,double v) {
  if(!c.HasKey(key)) c.ReadInputLine((std::string(key)+" "+std::to_string(v)).c_str());
  (*c.GetVector(key))[1]=v;
}
void Str(ConfigurationCard& c,const char* key,const std::string& v) {
  if(!c.HasKey(key)) c.ReadInputLine((std::string(key)+" "+v).c_str());
  else c.fValueString[c.GetTVectorIndex(key)]=v.c_str();
}
EECJetSelection::JetInputs Jet(double eta=.3,double phi=.1) {
  return {eta,phi,50,80,.2,.1,.5,.1,.1,5,20,25,false};
}
void IDAndTrackMax() {
  using namespace EECJetSelection;
  auto j=Jet(); j.neutralHadronFraction=.95;
  Need(PassJetID(j,JetIDPrescription::MiniPOGpp),"New central pp NHF<.99 not implemented");
  Need(!PassJetID(j,JetIDPrescription::HINpp),"Historical table was silently overwritten");
  j=Jet(2.8);j.neutralHadronFraction=.995;
  Need(!PassJetID(j,JetIDPrescription::MiniPOGpp),"Transition NHF requirement missing");
  Need(PassJetID(j,JetIDPrescription::HINpp),"Historical transition table changed");
  j=Jet(3.2);j.neutralHadronFraction=.5;j.neutralEmFraction=.5;
  Need(!PassJetID(j,JetIDPrescription::MiniPOGpp),"Forward pp NEF<.4 missing");
  Need(!PassJetID(j,JetIDPrescription::MiniPOGOO),"Forward OO NEF<.4 missing");
  Need(PassJetID(j,JetIDPrescription::HINpp),"Historical forward pp changed");
  Need(PassJetID(j,JetIDPrescription::HINPbPb),"Historical forward PbPb changed");
  j.neutralEmFraction=.4;Need(!PassJetID(j,JetIDPrescription::MiniPOGOO),"NEF boundary accepted");
  j.neutralEmFraction=std::nextafter(.4,0.);Need(PassJetID(j,JetIDPrescription::MiniPOGOO),"Allowed forward OO jet rejected");
  j.neutralMultiplicity=10;j.multiplicitiesAreLowerBounds=true;
  Reject([&]{PassJetID(j,JetIDPrescription::MiniPOGpp);},"Unresolved HF multiplicity was guessed");
  Need(PassTrackMaxFraction(2,100,.01,.98),"Raw-pT cut should pass");
  Need(!PassTrackMaxFraction(2,400,.01,.98),"Fixture does not distinguish raw from corrected denominator");
  Need(!PassTrackMaxFraction(1,100,.01,.98),"Lower strict boundary accepted");
  Need(!PassTrackMaxFraction(98,100,.01,.98),"Upper strict boundary accepted");
  Need(!PassTrackMaxFraction(0,0,.01,.98),"Invalid rawpt accepted");
}
void WritePpL1Fixture(const std::string& path,float rhoValue,float areaValue) {
  TFile f(path.c_str(),"RECREATE");
  f.mkdir("hiEvtAnalyzer")->cd(); TTree events("HiTree",""); float vz=0; events.Branch("vz",&vz,"vz/F"); events.Fill(); events.Write();
  f.mkdir("ak4PFJetAnalyzer")->cd(); TTree jets("t",""); int nref=1; float area[250]{}; area[0]=areaValue;
  jets.Branch("nref",&nref,"nref/I"); jets.Branch("jtarea",area,"jtarea[nref]/F"); jets.Fill(); jets.Write();
  f.mkdir("ggHiNtuplizer")->cd(); TTree rho("EventTree",""); float r=rhoValue; rho.Branch("rho",&r,"rho/F"); rho.Fill(); rho.Write();
}
void Calibration(const std::string& dir) {
  const std::string ppBase="jetEnergyCorrections/pp2024/";
  {
    ConfigurationCard card("cardEECpp5v36TeVMCR04.input");
    Need(card.HasKey("Run3UsePpL1FastJet") && card.Get("Run3UsePpL1FastJet")==0,
         "Current ppRef MC production card must keep L1FastJet disabled");
    JetCorrector nominal(ppBase+"Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt");
    nominal.SetJetPT(50); nominal.SetJetEta(.3); nominal.SetJetPhi(.1);
    auto stages=nominal.GetCorrectedPTStages();
    Need(stages.size()==2 && stages[0]==50 && std::isfinite(stages[1]),
         "L1-disabled pp MC chain must be raw -> L2Relative");
  }
  const std::string fixture=dir+"/pp_l1_fixture.root";
  WritePpL1Fixture(fixture,2.5,.5);
  Run3JetL1Inputs inputs("ggHiNtuplizer/EventTree","rho","jtarea");
  inputs.Bind({fixture.c_str()},"ak4PFJetAnalyzer/t");
  inputs.Read(0,1); auto l1=inputs.Get(0);
  Need(std::abs(l1.first-2.5)<1e-6 && std::abs(l1.second-.5)<1e-6,
       "ppRef rho/jtarea branches were not read as calibration inputs");

  {
    JetCorrector mc({ppBase+"Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt",
                     ppBase+"Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt"});
    mc.SetJetPT(50); mc.SetJetEta(.3); mc.SetJetPhi(.1); mc.SetRho(l1.first); mc.SetJetArea(l1.second);
    auto stages=mc.GetCorrectedPTStages();
    Need(stages.size()==3 && stages[0]==50 && std::isfinite(stages[1]) && std::isfinite(stages[2]),
         "pp MC must be raw -> L1FastJet -> L2Relative");
  }
  {
    JetCorrector data({ppBase+"Prompt24HIpp_V1_MC_L1FastJet_AK4PF.txt",
                       ppBase+"Prompt24HIpp_V2_MC_L2Relative_AK4PF.txt",
                       ppBase+"Prompt24HIpp_V1_DATA_L2Residual_AK4PF.txt"});
    data.SetJetPT(50); data.SetJetEta(.3); data.SetJetPhi(.1); data.SetRho(l1.first); data.SetJetArea(l1.second);
    auto stages=data.GetCorrectedPTStages();
    Need(stages.size()==4 && stages[0]==50 && std::isfinite(stages[1]) && std::isfinite(stages[2]) && std::isfinite(stages[3]),
         "pp data must be raw -> L1FastJet -> L2Relative -> residual");
  }
}
class Reco : public HighForestReader {
public:
  std::vector<EECJetSelection::JetInputs> jets;
  std::vector<float> trackMax;
  Reco() {
    fVertexZ=0;fHiBin=120;
    fPrimaryVertexFilterBit=fHfCoincidenceFilterBit=fClusterCompatibilityFilterBit=1;
    fHBHENoiseFilterBit=fBeamScrapingFilterBit=fPileupFilterBit=1;
  }
  Int_t GetNJets() const override{return jets.size();}
  Float_t GetJetRawPt(Int_t i) const override{return jets.at(i).rawPt;}
  Float_t GetJetPt(Int_t i) const override{return jets.at(i).storedPt;}
  Float_t GetJetEta(Int_t i) const override{return jets.at(i).eta;}
  Float_t GetJetPhi(Int_t i) const override{return jets.at(i).phi;}
  Float_t GetJetMaxTrackPt(Int_t i) const override{return trackMax.at(i);}
  EECJetSelection::JetInputs GetJetSelectionInputs(Int_t i) const override{return jets.at(i);}
};
void StagesBeforeVeto(const std::string& dir) {
  const auto mapPath=dir+"/mask.root",correction=dir+"/factor2.txt";
  {
    TFile f(mapPath.c_str(),"RECREATE");TH2D h("jetvetomap_all","",104,-5.2,5.2,72,-TMath::Pi(),TMath::Pi());
    h.SetBinContent(h.FindBin(3.2,2.3),-100);h.Write();
  }
  {std::ofstream f(correction);f<<"{1 JetEta 1 JetPt [0] Correction L2Relative}\n-5 5 3 1 1000 2\n";}
  ConfigurationCard c("cardEECOOInclusiveR04.input");Set(c,"Run3JEC",1);Set(c,"Run3JetQA",1);Set(c,"BackgroundMethods",0);
  Str(c,"JetVetoMapFile",mapPath);
  EECAnalyzer a({},&c,1,-1);a.fJetCorrector=new JetCorrector(correction);a.fTotalEventWeight=2;
  auto* reco=new Reco();a.fJetReader=reco;
  reco->jets={Jet(),Jet(),Jet(),Jet(3.2,2.3)};reco->jets[1].neutralEmFraction=.95;
  reco->trackMax={10,10,0,0}; // jet 3 fails trackMax but must still reject the event via veto.
  Need(!a.PassEventCuts(reco,false),"Event-veto eligibility incorrectly included trackMax");
  reco->jets.pop_back();reco->trackMax.pop_back();
  Need(a.PassEventCuts(reco,false),"Unmasked event failed");
  const auto path=dir+"/synthetic_stages.root";
  {TFile f(path.c_str(),"RECREATE");a.fHistograms->Write();}
  TFile f(path.c_str(),"READ");auto* h=dynamic_cast<TH3D*>(f.Get("run3JetQA/selectionPtCounts"));
  Need(h,"No saved stage spectra");
  const int x=h->GetXaxis()->FindBin(100),z=h->GetZaxis()->FindBin(60);
  Need(h->GetBinContent(x,1,z)==6,"No-veto denominator lost rejected-event jets");
  Need(h->GetBinContent(x,2,z)==4 && h->GetBinContent(x,3,z)==4,"Independent cut spectra incorrect");
  Need(h->GetBinContent(x,4,z)==2 && h->GetBinContent(x,6,z)==1,"Additive selection or final event veto incorrect");
  auto* calibration=dynamic_cast<TH3D*>(f.Get("run3JetQA/calibrationPtWeighted"));
  Need(calibration && calibration->GetBinContent(calibration->GetXaxis()->FindBin(50),1,z)==2,"Raw stage was not filled for same final-selected jet");
  Need(f.Get("run3JetQA/configuredJetVetoMap"),"Actual veto payload not persisted");
}
}
int main(int argc,char** argv) {
  try {
    if(argc!=2) throw std::invalid_argument("Pass a new test directory");
    std::filesystem::create_directories(argv[1]);
    IDAndTrackMax();Calibration(argv[1]);StagesBeforeVeto(argv[1]);
    std::cout<<"PASS: updated ID table, strict additive trackMax/rawPt, card-disabled pp L1FastJet with optional L1 regression, data-only residual, pre-veto denominators, fixed-jet calibration stages and saved mask\n";
  } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}

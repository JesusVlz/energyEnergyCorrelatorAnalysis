// Build with -fno-access-control only for the regression executable.
#include "../../src/EECAnalyzer.h"
#include <TDirectory.h>
#include <TObjString.h>
#include <TSystem.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include <stdexcept>

namespace {
void Require(bool value, const char* why){ if(!value) throw std::runtime_error(why); }
class RecoFixture : public HighForestReader {
public:
  std::vector<EECJetSelection::JetInputs> jets;
  RecoFixture(){
    fVertexZ=0;
    fPrimaryVertexFilterBit=fHfCoincidenceFilterBit=fClusterCompatibilityFilterBit=1;
    fHBHENoiseFilterBit=fBeamScrapingFilterBit=fPileupFilterBit=1;
  }
  Int_t GetNJets() const override { return jets.size(); }
  EECJetSelection::JetInputs GetJetSelectionInputs(Int_t index) const override { return jets.at(index); }
};
EECJetSelection::JetInputs Jet(double pt=20,double eta=2.8,double phi=.2){
  return {eta,phi,pt,pt,.1,.1,.5,.1,.1,3,3,6,false};
}
void Add(ConfigurationCard& card, const std::string& line){ card.ReadInputLine(line.c_str()); }
void Set(ConfigurationCard& card, const char* key, float value){
  if(!card.HasKey(key)) Add(card,std::string(key)+" "+std::to_string(value));
  auto* entry=card.GetVector(key); (*entry)[entry->GetLwb()]=value;
}
void SetString(ConfigurationCard& card,const char* key,const std::string& value){
  if(!card.HasKey(key)) Add(card,std::string(key)+" "+value);
  else card.fValueString[card.GetTVectorIndex(key)]=value.c_str();
}
void LoadWithoutVetoSettings(ConfigurationCard& card,const char* path){
  std::ifstream input(path);
  if(!input) throw std::runtime_error(std::string("Cannot open fixture base card: ")+path);
  std::string line;
  while(std::getline(input,line)){
    std::istringstream tokens(line); std::string key; tokens>>key;
    if(key=="Run3JetID" || key=="ApplyJetVeto" || key=="JetVetoMapFile" ||
       key=="JetVetoMapHistogram" || key=="JetVetoPtConvention" || key=="JetVetoInputMinPt") continue;
    card.ReadInputLine(line.c_str());
  }
}
void WriteMap(const std::string& path){
  TDirectory::TContext context;
  TFile file(path.c_str(),"RECREATE");
  TH2D map("explicit_test_key",";eta;phi",104,-5.2,5.2,72,-TMath::Pi(),TMath::Pi());
  map.SetBinContent(map.FindBin(2.8,.2),-1);
  map.SetBinContent(map.FindBin(.3,.2),-1);
  map.Write();
}
void Configure(ConfigurationCard& card,const std::string& map){
  Set(card,"Run3JetID",1); Set(card,"ApplyJetVeto",1);
  SetString(card,"JetVetoMapFile",map); SetString(card,"JetVetoMapHistogram","explicit_test_key");
  SetString(card,"JetVetoPtConvention","stored"); Set(card,"JetVetoInputMinPt",15);
}
void CheckSystem(const char* name,const std::string& map,const std::string& output){
  ConfigurationCard card(name); Configure(card,map);
  auto analyzer=std::make_unique<EECAnalyzer>(std::vector<TString>{},&card,1,-1);
  Require(analyzer->fRun3JetID==1 && analyzer->fApplyJetVeto,"Explicit configuration not enabled");
  auto* signal=new RecoFixture();
  auto mixing=std::make_unique<RecoFixture>();
  analyzer->fJetReader=signal; analyzer->fMixedEventReader=mixing.get();
  analyzer->fJetMinimumPtCut=40; // explicit EEC study threshold, independent of veto
  analyzer->fJetEtaCut=1.6;
  const auto outside=Jet();
  Require(outside.storedPt<analyzer->fJetMinimumPtCut && std::abs(outside.eta)>analyzer->fJetEtaCut,
          "Fixture must sit below EEC pT and outside its eta acceptance");
  signal->jets={Jet(80,.3,-1.1),outside};
  Require(!analyzer->PassEventCuts(signal,false),"Vetoed low-pT/outside-EEC jet failed to reject whole event");
  mixing->jets=signal->jets;
  Require(!analyzer->PassEventCuts(mixing.get(),false),"Signal/mixing masks differ");
  signal->jets[1].storedPt=15;
  Require(analyzer->PassEventCuts(signal,false),"Exactly 15 GeV must not veto");
  signal->jets[1]=outside; signal->jets[1].neutralEmFraction=.8;
  Require(analyzer->PassEventCuts(signal,false),"EM sum exactly .9 must not veto");
  signal->jets[1]=Jet(20,.3,.2); signal->jets[1].muonFraction=.8;
  Require(analyzer->PassEventCuts(signal,false),"Jet failing ID must not veto");
  signal->jets={Jet(80,.3,-1.1)};
  Require(analyzer->PassEventCuts(signal,false),"Unmasked event must pass");
  // Gen EEC event selection must still be based on the reconstructed companion.
  auto* companion=new RecoFixture(); companion->jets={outside};
  analyzer->fRecoJetReader=companion;
  Require(!analyzer->PassEventCuts(signal,false),"Generator-mode event lost reconstructed veto");
  analyzer->fRecoJetReader=nullptr; delete companion;
  auto* counts=analyzer->fHistograms->fhJetVetoEventCounts;
  Require(counts && counts->GetBinContent(1)==6 && counts->GetBinContent(2)==4 && counts->GetBinContent(3)==2,
          "Signal event QA does not account for full event rejection");
  Require(counts->GetBinContent(4)==1 && counts->GetBinContent(6)==1,"Mixing event QA wrong");
  Require(analyzer->fHistograms->fhJetVetoEtaPhiBefore[0]->GetEntries()>
          analyzer->fHistograms->fhJetVetoEtaPhiAfter[0]->GetEntries(),"After QA must exclude all jets from vetoed events");
  {
    TDirectory::TContext context;
    TFile file(output.c_str(),"RECREATE"); analyzer->fHistograms->Write(); card.WriteCard(&file);
  }
  TFile saved(output.c_str(),"READ");
  auto* key=dynamic_cast<TObjString*>(saved.Get("JCard/JetVetoMapHistogramString"));
  Require(key && key->String()=="explicit_test_key","Map histogram string was not persisted");
  Require(saved.Get("eventQA/jetVetoEventCounts"),"Veto QA not persisted");
}
void ExpectBadConfig(const std::string& map,bool omitPt){
  // Explicitly omit the tested settings instead of relying on the current
  // production defaults or appending duplicate configuration keys.
  ConfigurationCard card(""); LoadWithoutVetoSettings(card,"cardEECOO.input");
  Set(card,"Run3JetID",1); Set(card,"ApplyJetVeto",1);
  SetString(card,"JetVetoMapFile",map); SetString(card,"JetVetoMapHistogram","explicit_test_key");
  Set(card,"JetVetoInputMinPt",20);
  if(!omitPt) SetString(card,"JetVetoPtConvention","stored");
  try { EECAnalyzer analyzer({},&card,1,-1); }
  catch(const std::invalid_argument&){return;}
  throw std::runtime_error("Incomplete or >15 GeV input coverage was accepted");
}
void CheckMegaSkimScope(const std::string& map){
  ConfigurationCard noMixing("cardEECOOMC.input"); Configure(noMixing,map);
  Set(noMixing,"BackgroundMethods",0); Set(noMixing,"MegaSkimMixing",1);
  auto accepted=std::make_unique<EECAnalyzer>(std::vector<TString>{},&noMixing,1,-1);
  Require(accepted->fApplyJetVeto && !accepted->fDoMixedCone && accepted->fMegaSkimMode,
          "Mega-skim scope fixture was not configured explicitly");
  accepted->RunAnalysis(); // Empty signal input still initializes production machinery.
  ConfigurationCard mixing("cardEECOOMC.input"); Configure(mixing,map);
  Set(mixing,"BackgroundMethods",4); Set(mixing,"MegaSkimMixing",1);
  try {
    EECAnalyzer rejected({},&mixing,1,-1); rejected.RunAnalysis();
  } catch(const std::invalid_argument& error){
    Require(std::string(error.what()).find("mega skims")!=std::string::npos,
            "Mixed mega-skim rejection lost its explicit diagnostic");
    return;
  }
  throw std::runtime_error("Actual mega-skim mixing accepted an unverifiable veto mask");
}
}
int main(int argc,char** argv){
  try {
    const std::string directory=argc>1?argv[1]:"/tmp/eecRun3JetVetoAnalyzer";
    gSystem->mkdir(directory.c_str(),true);
    const auto map=directory+"/synthetic_veto.root"; WriteMap(map);
    int index=0;
    for(const char* name:{"cardEECOO.input","cardEECOOMC.input","cardEECpp5v36TeV.input","cardEECpp5v36TeVMC.input"})
      CheckSystem(name,map,directory+"/qa"+std::to_string(index++)+".root");
    ExpectBadConfig(map,true); ExpectBadConfig(map,false);
    CheckMegaSkimScope(map);
    ConfigurationCard disabled("cardEECOO.input");
    Set(disabled,"Run3JetID",0); Set(disabled,"ApplyJetVeto",0);
    EECAnalyzer baseline({},&disabled,1,-1);
    Require(!baseline.fApplyJetVeto && !baseline.fHistograms->fhJetVetoEventCounts,"Disabled veto fabricated measured loss QA");
    std::cout<<"PASS: all four systems, whole-event/15 GeV/EM/ID boundaries, gen companion, same mixing mask, QA, explicit config and persisted strings\n";
  } catch(const std::exception& error){ std::cerr<<"FAIL: "<<error.what()<<'\n';return 1; }
}

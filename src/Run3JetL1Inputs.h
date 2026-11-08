#ifndef EEC_RUN3_JET_L1_INPUTS_H
#define EEC_RUN3_JET_L1_INPUTS_H

#include "Run3JetSelectionInputs.h"
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Independent branch-only chains avoid changing any existing reader address.
// The caller MUST identify the applicable event rho and measured per-jet area.
// No rho=0, area=pi*R^2, alternative branch, or stored-jtpt fallback exists.
class Run3JetL1Inputs {
public:
  Run3JetL1Inputs(std::string rhoTree, std::string rhoBranch, std::string areaBranch)
    : fRhoTreeName(std::move(rhoTree)), fRhoName(std::move(rhoBranch)), fAreaName(std::move(areaBranch)) {
    if(fRhoTreeName.empty() || fRhoName.empty() || fAreaName.empty())
      throw std::invalid_argument("Run-3 L1 requires explicit rho tree, scalar rho branch and per-jet area branch");
  }

  void Bind(const std::vector<TString>& files, const char* jetTreeName) {
    if(files.empty()) throw std::invalid_argument("Run-3 L1 input file list is empty");
    fAreas.reset(new TChain(jetTreeName));
    fRhos.reset(new TChain(fRhoTreeName.c_str()));
    bool first = true;
    for(const auto& path : files) {
      std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
      if(!file || file->IsZombie()) throw std::runtime_error("Run-3 L1 cannot open input forest");
      auto* jets = dynamic_cast<TTree*>(file->Get(jetTreeName));
      auto* rho = dynamic_cast<TTree*>(file->Get(fRhoTreeName.c_str()));
      auto* events = dynamic_cast<TTree*>(file->Get("hiEvtAnalyzer/HiTree"));
      if(!jets || !rho || !events || jets->GetEntries()!=events->GetEntries() || rho->GetEntries()!=events->GetEntries())
        throw std::runtime_error("Run-3 L1 needs aligned jet, rho and event trees with identical entry counts in every file");
      const std::string areaType = Type(jets, fAreaName, "nref");
      const std::string rhoType = Type(rho, fRhoName);
      if(first) { fAreaDouble = areaType=="Double_t"; fRhoDouble = rhoType=="Double_t"; first = false; }
      if(areaType!=(fAreaDouble?"Double_t":"Float_t") || rhoType!=(fRhoDouble?"Double_t":"Float_t"))
        throw std::runtime_error("Run-3 L1 branch representation changes between input files");
      EECJetSelection::Detail::RequireJetLeaf(jets, "nref", "Int_t");
      if(fAreas->Add(path)==0 || fRhos->Add(path)==0) throw std::runtime_error("Run-3 L1 cannot add input tree");
    }
    fAreas->SetBranchStatus("*",0); fRhos->SetBranchStatus("*",0);
    fAreas->SetBranchStatus("nref",1); fAreas->SetBranchStatus(fAreaName.c_str(),1);
    fRhos->SetBranchStatus(fRhoName.c_str(),1);
    if(fAreas->SetBranchAddress("nref",&fCount)<0 ||
       (fAreaDouble ? fAreas->SetBranchAddress(fAreaName.c_str(),fAreaD.data()) : fAreas->SetBranchAddress(fAreaName.c_str(),fAreaF.data()))<0 ||
       (fRhoDouble ? fRhos->SetBranchAddress(fRhoName.c_str(),&fRhoD) : fRhos->SetBranchAddress(fRhoName.c_str(),&fRhoF))<0)
      throw std::runtime_error("Run-3 L1 cannot bind configured rho/area inputs");
    fRead = false;
  }

  void Read(Long64_t entry, int expectedJetCount) {
    fRead = false;
    if(!fAreas || !fRhos || entry<0 || entry>=fAreas->GetEntries()) throw std::out_of_range("Run-3 L1 invalid event entry");
    const auto areaEntry = fAreas->LoadTree(entry), rhoEntry = fRhos->LoadTree(entry);
    auto* areaTree = fAreas->GetTree(); auto* rhoTree = fRhos->GetTree();
    if(areaEntry<0 || rhoEntry<0 || !areaTree || !rhoTree) throw std::runtime_error("Run-3 L1 cannot load entry");
    // Revalidate each transition before reading any variable-length payload.
    EECJetSelection::Detail::RequireJetLeaf(areaTree,"nref","Int_t");
    if(Type(areaTree,fAreaName,"nref")!=(fAreaDouble?"Double_t":"Float_t") ||
       Type(rhoTree,fRhoName)!=(fRhoDouble?"Double_t":"Float_t"))
      throw std::runtime_error("Run-3 L1 input schema changed at a chain transition");
    if(areaTree->GetBranch("nref")->GetEntry(areaEntry)<=0 || fCount<0 || fCount>int(fAreaD.size()) || fCount!=expectedJetCount)
      throw std::runtime_error("Run-3 L1 native jet count is invalid or disagrees with the analysis reader");
    if((fCount && areaTree->GetBranch(fAreaName.c_str())->GetEntry(areaEntry)<=0) ||
       rhoTree->GetBranch(fRhoName.c_str())->GetEntry(rhoEntry)<=0)
      throw std::runtime_error("Run-3 L1 could not read configured rho/area payload");
    const double value = fRhoDouble ? fRhoD : fRhoF;
    if(!std::isfinite(value) || value<0) throw std::runtime_error("Run-3 L1 rho must be finite and nonnegative");
    fRead = true;
  }

  std::pair<double,double> Get(int index) const {
    if(!fRead || index<0 || index>=fCount) throw std::out_of_range("Run-3 L1 inputs not loaded for this jet");
    const double area = fAreaDouble ? fAreaD[index] : fAreaF[index];
    if(!std::isfinite(area) || area<=0) throw std::runtime_error("Run-3 L1 requires positive measured jet area; no geometric fallback");
    return {fRhoDouble ? fRhoD : fRhoF,area};
  }
private:
  static std::string Type(TTree* tree,const std::string& name,const char* count=nullptr) {
    auto* leaf=tree?tree->GetLeaf(name.c_str()):nullptr;
    if(!leaf) throw std::runtime_error("Run-3 L1 missing configured branch "+name);
    const std::string type=leaf->GetTypeName();
    if(type!="Float_t" && type!="Double_t") throw std::runtime_error("Run-3 L1 branch must be Float_t or Double_t: "+name);
    EECJetSelection::Detail::RequireJetLeaf(tree,name.c_str(),type.c_str(),count);
    return type;
  }
  std::string fRhoTreeName,fRhoName,fAreaName;
  std::unique_ptr<TChain> fAreas,fRhos;
  std::array<Float_t,250> fAreaF{};
  std::array<Double_t,250> fAreaD{};
  Float_t fRhoF=0;
  Double_t fRhoD=0;
  Int_t fCount=0;
  bool fAreaDouble=false,fRhoDouble=false,fRead=false;
};
#endif

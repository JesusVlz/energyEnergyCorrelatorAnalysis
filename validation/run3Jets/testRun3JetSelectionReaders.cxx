#include "HighForestReader.h"
#include <TDirectory.h>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

// Readers contain large legacy track buffers; allocate test readers on the heap.
// All fixtures are synthetic ROOT files under a unique temporary directory.
struct FixtureDirectory {
  std::filesystem::path original = std::filesystem::current_path();
  std::filesystem::path directory;
  FixtureDirectory(){
    char pattern[] = "/tmp/eec-run3-reader-inputs-XXXXXX";
    const char* created = mkdtemp(pattern);
    if(!created) throw std::runtime_error("Cannot create reader fixture directory");
    directory = created;
    std::filesystem::current_path(directory);
  }
  ~FixtureDirectory(){
    std::error_code error;
    std::filesystem::current_path(original, error);
    std::filesystem::remove_all(directory, error);
  }
};

void fixture(const char* path, float eta, int count=2, bool badType=false, bool badShape=false, bool omit=false){
  TFile file(path,"RECREATE");
  auto hi=file.mkdir("hiEvtAnalyzer"); hi->cd();
  TTree events("HiTree","events");
  Float_t vz=0, hf=3; ULong64_t evt=17;
  events.Branch("vz",&vz,"vz/F"); events.Branch("hiHFPlus_pf",&hf,"hiHFPlus_pf/F");
  events.Branch("hiHFMinus_pf",&hf,"hiHFMinus_pf/F"); events.Branch("evt",&evt,"evt/l");
  events.Fill(); events.Write();
  file.cd(); auto skim=file.mkdir("skimanalysis"); skim->cd(); TTree filters("HltTree","filters");
  Int_t primary=1; filters.Branch("pprimaryVertexFilter",&primary,"pprimaryVertexFilter/I"); filters.Fill(); filters.Write();
  file.cd(); auto jetdir=file.mkdir("ak4PFJetAnalyzer"); jetdir->cd(); TTree jets("t","jets");
  Int_t nref=count; Float_t values[600]={}; Int_t counts[600]={};
  for(int i=0;i<600;++i){values[i]=eta;counts[i]=3;}
  jets.Branch("nref",&nref,"nref/I");
  for(const char* name:{"jtpt","rawpt","jteta","jtphi","trackMax","jtPfCHF","jtPfNHF","jtPfCEF","jtPfNEF","jtPfMUF"}){
    if(omit && std::string(name)=="jtPfMUF") continue;
    const std::string shape=std::string(name)+((badShape&&std::string(name)=="jtPfCHF")?"[nref][2]/F":"[nref]/F");
    jets.Branch(name,values,shape.c_str());
  }
  Float_t wta[300]={}; for(float& value:wta)value=.25;
  jets.Branch("WTAeta",wta,"WTAeta[nref]/F"); jets.Branch("WTAphi",wta,"WTAphi[nref]/F");
  for(const char* name:{"jtPfCHM","jtPfNHM","jtPfCEM","jtPfNEM","jtPfMUM"}){
    const std::string shape=std::string(name)+((badType&&std::string(name)=="jtPfCHM")?"[nref]/F":"[nref]/I");
    jets.Branch(name,counts,shape.c_str());
  }
  jets.Fill(); jets.Write(); file.Close();
}

template<class F> void mustThrow(F operation){bool threw=false;try{operation();}catch(const std::exception& error){threw=true;std::cout<<"Expected: "<<error.what()<<'\n';}assert(threw);}
int main(){
  FixtureDirectory fixtureDirectory;
  fixture("run3_reader_a.root",1.25);fixture("run3_reader_b.root",1.75);
  auto normal=std::make_unique<HighForestReader>(ForestReader::kPpRef5p36TeV,0,1,1,0,false,false,false);
  normal->EnableRun3JetSelection(); normal->ReadForestFromFileList({"run3_reader_a.root","run3_reader_b.root"});
  assert(normal->GetNJets()==2); assert(normal->GetJetEta(0)==.25);
  auto inputs=normal->GetJetSelectionInputs(0);assert(inputs.eta==1.25);assert(inputs.chargedMultiplicity==9);
  assert(inputs.neutralMultiplicity==6);assert(inputs.numberOfConstituents==15);assert(inputs.multiplicitiesAreLowerBounds);
  normal->GetEvent(1); assert(normal->GetJetSelectionInputs(0).eta==1.75);assert(normal->GetJetEta(0)==.25);
  auto standard=std::make_unique<HighForestReader>(ForestReader::kPpRef5p36TeV,0,1,0,0,false,false,false);
  standard->EnableRun3JetSelection();standard->ReadForestFromFileList({"run3_reader_a.root"});
  assert(standard->GetJetEta(0)==1.25);assert(standard->GetJetSelectionInputs(0).eta==1.25);
  mustThrow([&]{normal->GetJetSelectionInputs(2);});mustThrow([&]{normal->EnableRun3JetSelection(false);});
  auto mixing=std::make_unique<HighForestReader>(ForestReader::kPpRef5p36TeV,0,1,1,0,false,true,false);
  mixing->EnableRun3JetSelection();mixing->ReadForestFromFileList({"run3_reader_a.root","run3_reader_b.root"});
  mixing->GetEvent(1);assert(mixing->GetNJets()==2);assert(mixing->GetJetSelectionInputs(1).eta==1.75);
  auto disabled=std::make_unique<HighForestReader>(ForestReader::kPpRef5p36TeV,0,1,1,0,false,true,false);
  disabled->ReadForestFromFileList({"run3_reader_a.root"});assert(disabled->GetNJets()==0);
  mustThrow([&]{disabled->GetJetSelectionInputs(0);});
  fixture("run3_reader_big.root",1.,251);
  auto growing=std::make_unique<HighForestReader>(ForestReader::kPpRef5p36TeV,0,1,1,0,false,true,false);
  growing->EnableRun3JetSelection();growing->ReadForestFromFileList({"run3_reader_a.root","run3_reader_big.root"});
  mustThrow([&]{growing->GetEvent(1);});
  fixture("run3_reader_type.root",1.,2,true);
  fixture("run3_reader_shape.root",1.,2,false,true);
  fixture("run3_reader_missing.root",1.,2,false,false,true);
  for(const char* file:{"run3_reader_big.root","run3_reader_type.root","run3_reader_shape.root","run3_reader_missing.root"}){
    mustThrow([&]{auto invalid=std::make_unique<HighForestReader>(ForestReader::kPpRef5p36TeV,0,1,1,0,false,true,false);invalid->EnableRun3JetSelection();invalid->ReadForestFromFileList({file});});
  }
  std::cout<<"PASS: Run-3 native reader inputs, WTA separation, chain transitions, mixing, disabled behavior, missing/type/shape/capacity failures\n";
}

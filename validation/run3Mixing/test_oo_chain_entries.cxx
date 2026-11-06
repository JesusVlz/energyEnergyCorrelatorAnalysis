#include "GeneratorLevelForestReader.h"
#include "HighForestReader.h"
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cassert>
#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>

static void fixture(const std::string& path,int count,ULong64_t first){
  TFile file(path.c_str(),"RECREATE");
  file.mkdir("hiEvtAnalyzer")->cd();
  TTree event("HiTree","physical events");
  Float_t vz=0,hf=.125,pthat=0,weight=1;Int_t hiBin=20;UInt_t run=1,lumi=1;ULong64_t evt=first;
  event.Branch("run",&run,"run/i");event.Branch("lumi",&lumi,"lumi/i");event.Branch("evt",&evt,"evt/l");
  event.Branch("vz",&vz,"vz/F");event.Branch("hiBin",&hiBin,"hiBin/I");
  event.Branch("hiHF_pf",&hf,"hiHF_pf/F");event.Branch("hiHFPlus_pf",&hf,"hiHFPlus_pf/F");event.Branch("hiHFMinus_pf",&hf,"hiHFMinus_pf/F");
  event.Branch("pthat",&pthat,"pthat/F");event.Branch("weight",&weight,"weight/F");
  for(int i=0;i<count;++i){evt=first+i;event.Fill();}event.Write();
  file.cd();file.mkdir("skimanalysis")->cd();TTree skim("HltTree","filters");Int_t pass=1;
  for(const char* name:{"pprimaryVertexFilter","OOphfCoincFilterPF2Th4","pileupVertexFilter"}) skim.Branch(name,&pass,(std::string(name)+"/I").c_str());
  for(int i=0;i<count;++i)skim.Fill();skim.Write();
  file.cd();file.mkdir("HiGenParticleAna")->cd();TTree gen("hi","particles");
  std::vector<float> pt{1},eta{0},phi{0};std::vector<int> chg{1},sube{0};
  gen.Branch("pt",&pt);gen.Branch("eta",&eta);gen.Branch("phi",&phi);gen.Branch("chg",&chg);gen.Branch("sube",&sube);
  for(int i=0;i<count;++i)gen.Fill();gen.Write();file.Close();
}
int main(int argc,char** argv){
  assert(argc==2);const std::filesystem::path directory=argv[1];std::filesystem::create_directories(directory);
  const auto a=(directory/"first.root").string(),b=(directory/"second.root").string();
  const ULong64_t large=ULong64_t(1)<<42;
  fixture(a,2,large+10);fixture(b,3,large+20);
  auto gen=std::make_unique<GeneratorLevelForestReader>(ForestReader::kOOMC,0,1,1,0,true,true,false);
  auto reco=std::make_unique<HighForestReader>(ForestReader::kOOMC,0,1,1,0,false,true,false);
  gen->ReadForestFromFileList({a,b});reco->ReadForestFromFileList({a,b});
  assert(gen->GetNEvents()==2 && reco->GetNEvents()==2); // legacy first-source behavior retained
  assert(gen->GetTotalChainEntries()==5 && reco->GetTotalChainEntries()==5);
  const int globalCap=4;
  const int limit=static_cast<int>(std::min<Long64_t>(gen->GetTotalChainEntries(),globalCap));
  assert(limit==4 && limit>gen->GetNEvents());
  std::vector<ULong64_t> seen;
  for(int i=0;i<limit;++i){gen->GetEvent(i);reco->GetEvent(i);assert(gen->GetEventNumber()==reco->GetEventNumber());seen.push_back(gen->GetEventNumber());}
  assert((seen==std::vector<ULong64_t>{large+10,large+11,large+20,large+21}));
  assert(gen->GetTotalChainEntries()==5 && reco->GetTotalChainEntries()==5);
  assert(gen->GetRawHFSum()==.125f);
  std::cout<<"PASS: explicit GenGen/reco total-chain count spans 2+3 entries, global cap 4 crosses files, 64-bit event IDs/read alignment retained, legacy count unchanged\n";
}

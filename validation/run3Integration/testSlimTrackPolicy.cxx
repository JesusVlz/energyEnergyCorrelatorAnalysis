// Synthetic reader policy test. Compile with -fno-access-control for tree status.
#include "../../src/HighForestReader.h"
#define EEC_READER_FIXTURES_ONLY
#include "../run3Response/validateUnfoldingForestReader.cxx"
#include <TError.h>
#include <TBranch.h>
#include <functional>
#include <vector>
#include <iostream>

namespace {
const std::vector<std::string> unused = {"trkNormChi2","trkNLayers","trkNHits",
  "pfEcal","pfHcal","trkChi2","trkNdof","trkNlayer","trkNHit"};
int unexpectedWarnings = 0;
void Handler(int level, Bool_t abort, const char* location, const char* message){
  if(level >= kError) ++unexpectedWarnings;
  else if(level >= kWarning) for(const auto& name : unused)
    if(std::string(message).find(name) != std::string::npos) ++unexpectedWarnings;
  DefaultErrorHandler(level,abort,location,message);
}
void Require(bool ok,const char* what){ if(!ok) throw std::runtime_error(what); }
void Throws(const std::function<void()>& operation){
  bool rejected=false;
  try {operation();} catch(const std::exception&) {rejected=true;}
  Require(rejected,"Run-3 legacy getter was not rejected");
}
template<class Reader> void Check(Reader& reader,TTree* tracks){
  for(const auto& name : unused) if(tracks->GetBranch(name.c_str()))
    Require(!tracks->GetBranchStatus(name.c_str()),"Unused track branch was enabled");
  Require(reader.GetNTracks()==2,"Lost reconstructed tracks");
  Require(reader.GetTrackCharge(1)==-1,"Lost physical charge");
  Require(std::abs(reader.GetTrackPt(0)-4)<1e-5,"Lost track momentum");
  Throws([&]{reader.GetTrackNormalizedChi2(0);});
  Throws([&]{reader.GetTrackChi2(0);});
  Throws([&]{reader.GetNTrackDegreesOfFreedom(0);});
  Throws([&]{reader.GetNHitsTrackerLayer(0);});
  Throws([&]{reader.GetNHitsTrack(0);});
  Throws([&]{reader.GetTrackEnergyEcal(0);});
  Throws([&]{reader.GetTrackEnergyHcal(0);});
}
// Complete the event/skim inputs needed by the ordinary reader, independently
// of the response-only fixture. Append branch entries, never duplicate events.
void CompleteFixture(const std::string& path,int type,bool aod){
  TFile file(path.c_str(),"UPDATE");
  auto* events=dynamic_cast<TTree*>(file.Get("hiEvtAnalyzer/HiTree"));
  ULong64_t eventId=101; float hf=4; unsigned int run=1,lumi=1;
  std::vector<TBranch*> branches;
  branches.push_back(events->Branch("evt",&eventId,"evt/l"));
  branches.push_back(events->Branch("run",&run,"run/i"));
  branches.push_back(events->Branch("lumi",&lumi,"lumi/i"));
  branches.push_back(events->Branch("hiHFPlus_pf",&hf,"hiHFPlus_pf/F"));
  branches.push_back(events->Branch("hiHFMinus_pf",&hf,"hiHFMinus_pf/F"));
  for(int i=0;i<2;++i){eventId=101+i;for(auto* branch:branches)branch->Fill();}
  events->GetDirectory()->cd();events->Write("",TObject::kOverwrite);
  if(!aod){
    auto* tracks=dynamic_cast<TTree*>(file.Get(type<=10 ? "PbPbTracks/trackTree":"ppTracks/trackTree"));
    int n=2;auto* branch=tracks->Branch("nTrk",&n,"nTrk/I");
    for(int i=0;i<2;++i)branch->Fill();
    tracks->GetDirectory()->cd();tracks->Write("",TObject::kOverwrite);
  }
  if(type==10 || type==12){
    auto* jets=dynamic_cast<TTree*>(file.Get(type==10?"akCs4PFJetAnalyzer/t":"ak4PFJetAnalyzer/t"));
    int flavor[2]={1,21};auto* branch=jets->Branch("refparton_flavorForB",flavor,"refparton_flavorForB[nref]/I");
    for(int i=0;i<2;++i){jets->GetEntry(i);branch->Fill();}
    jets->GetDirectory()->cd();jets->Write("",TObject::kOverwrite);
  }
  file.mkdir("skimanalysis")->cd();TTree skim("HltTree","");int pass=1;
  for(const char* name:{"pprimaryVertexFilter","phfCoincFilter2Th4","pOOprimaryVertexFilter",
                       "OOpfCoincFilterPF2Th4","OOphfCoincFilterPF2Th4","pileupVertexFilter"})
    skim.Branch(name,&pass,(std::string(name)+"/I").c_str());
  for(int i=0;i<2;++i)skim.Fill();skim.Write();
}
}
int main(int argc,char** argv){
  try {
    const std::string directory=argc>1?argv[1]:"/tmp/slim-track-policy";
    gSystem->mkdir(directory.c_str(),true);
    auto oldHandler=SetErrorHandler(Handler);
    int checked=0;
    for(int type:{9,10,11,12}) for(bool aod:{false,true}) for(bool present:{false,true}){
      const std::string stem=directory+"/"+std::to_string(type)+(aod?"AOD":"Mini")+(present?"Full":"Slim");
      const std::string a=stem+"A.root",b=stem+"B.root";
      for(const auto& path:{a,b}){
        WriteRun3ReaderFixture(path,type,1,"",aod,present);CompleteFixture(path,type,aod);
      }
      for(bool mixing:{false,true}){
        auto reader=std::make_unique<HighForestReader>(type,0,1,1,0,true,mixing,false);
        reader->ReadForestFromFileList({a.c_str(),b.c_str()});
        for(int entry:{0,2}){reader->GetEvent(entry);Check(*reader,reader->fTrackTree);++checked;}
      }
      TFile input(a.c_str(),"READ");auto reader=std::make_unique<UnfoldingForestReader>(type,1,1);
      reader->ReadForestFromFile(&input);reader->GetEvent(0);Check(*reader,reader->fTrackTree);++checked;
    }
    SetErrorHandler(oldHandler);
    Require(unexpectedWarnings==0,"Reader reported an unexpected error or requested disabled legacy track branches");
    std::cout<<"PASS: "<<checked<<" reader cases; all four Run-3 types, AOD/MiniAOD, fields present/absent, ordinary/mixed chains and unfolding; no reader errors or unused-branch warnings.\n";
  } catch(const std::exception& error){std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;}
}

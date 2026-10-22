// Synthetic ROOT schema regression, not a validation of a production MC forest.
/* Build from the repository root:
g++ $(root-config --cflags) validation/run3Response/validateUnfoldingForestReader.cxx \
  src/UnfoldingForestReader.cxx $(root-config --libs) -o /tmp/validateUnfoldingForestReader
/tmp/validateUnfoldingForestReader /tmp/eecRun3ReaderValidation
*/

#include "../../src/UnfoldingForestReader.h"

#include <TSystem.h>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

// Kept available when this source is included with EEC_READER_FIXTURES_ONLY, so
// analysis selection tests can reuse exactly the same small input forest.
void WriteRun3ReaderFixture(const std::string& path, int dataType, int jetType = 1,
                           const std::string& issue = "", bool aod = false){
  const bool isMC = dataType == 2 || dataType == 3 || dataType == 10 || dataType == 12;
  const bool heavyIon = dataType == 1 || dataType == 3 || dataType == 9 || dataType == 10;
  const bool run3 = dataType >= 9;
  const char* ooJets[] = {"akCs2PFJetAnalyzer", "akCs4PFJetAnalyzer", "akCs6PFJetAnalyzer", "akCs8PFJetAnalyzer"};
  const char* pbpbJets[] = {"akPu4CaloJetAnalyzer", "akCs4PFJetAnalyzer", "akPu4PFJetAnalyzer", "akFlowPuCs4PFJetAnalyzer"};
  const char* ppJets[] = {"ak4CaloJetAnalyzer", "ak4PFJetAnalyzer", "akCs4PFJetAnalyzer"};
  const char* jetDirectory = (dataType == 9 || dataType == 10) ? ooJets[jetType] :
                             (dataType == 1 || dataType == 3) ? pbpbJets[jetType] : ppJets[jetType];

  TFile file(path.c_str(), "RECREATE");
  if(file.IsZombie()) throw std::runtime_error("Cannot create fixture " + path);
  file.mkdir("hiEvtAnalyzer")->cd();
  TTree event("HiTree", "Synthetic event information");
  float vz = 2.5, pthat = 80, weight = 0.25;
  int hiBin = 60;
  event.Branch("vz", &vz, "vz/F");
  if(heavyIon || !run3) event.Branch("hiBin", &hiBin, "hiBin/I");
  if(isMC){
    event.Branch("pthat", &pthat, "pthat/F");
    if(issue != "weight") event.Branch("weight", &weight, "weight/F");
  }

  file.mkdir(jetDirectory)->cd();
  TTree jets("t", "Synthetic reconstructed and generated jets");
  int nref = issue == "overflowJets" ? 251 : 2, ngen = 1;
  // Larger fixture buffers let us create an input which must be rejected before
  // it is copied into the reader's 250-element buffers.
  float jtpt[300] = {128, 80}, rawpt[300] = {120, 75}, trackMax[300] = {12, 8};
  float jtphi[300] = {0.2, 2.0}, jteta[300] = {0.15, 1.2};
  float wtaPhi[300] = {0.25, 2.05}, wtaEta[300] = {0.18, 1.25};
  float refpt[300] = {130, -999}, refphi[300] = {0.2, -999}, refeta[300] = {0.15, -999};
  float genpt[300] = {130}, genphi[300] = {0.2}, geneta[300] = {0.15};
  float wtaGenPhi[300] = {0.22}, wtaGenEta[300] = {0.17};
  jets.Branch("nref", &nref, "nref/I");
  jets.Branch("jtpt", jtpt, "jtpt[nref]/F");
  jets.Branch("rawpt", rawpt, "rawpt[nref]/F");
  jets.Branch("trackMax", trackMax, "trackMax[nref]/F");
  jets.Branch("jtphi", jtphi, "jtphi[nref]/F");
  jets.Branch("jteta", jteta, "jteta[nref]/F");
  if(issue != "WTA"){
    jets.Branch("WTAphi", wtaPhi, "WTAphi[nref]/F");
    jets.Branch("WTAeta", wtaEta, "WTAeta[nref]/F");
  }
  if(isMC){
    if(issue != "refpt") jets.Branch("refpt", refpt, "refpt[nref]/F");
    jets.Branch("refphi", refphi, "refphi[nref]/F");
    jets.Branch("refeta", refeta, "refeta[nref]/F");
    jets.Branch("ngen", &ngen, "ngen/I");
    jets.Branch("genpt", genpt, "genpt[ngen]/F");
    jets.Branch("genphi", genphi, "genphi[ngen]/F");
    jets.Branch("geneta", geneta, "geneta[ngen]/F");
    if(issue != "WTA"){
      jets.Branch("WTAgenphi", wtaGenPhi, "WTAgenphi[ngen]/F");
      jets.Branch("WTAgeneta", wtaGenEta, "WTAgeneta[ngen]/F");
    }
  }

  file.mkdir(aod ? "ppTrack" : heavyIon ? "PbPbTracks" : "ppTracks")->cd();
  TTree tracks("trackTree", "Synthetic reconstructed tracks");
  int nTrk = 2;
  std::vector<float> pt = {4, 1}, ptError = {0.04, 0.02}, phi = {0.2, 0.25}, eta = {0.15, 0.18};
  std::vector<float> dz = {0.01, 0.01}, dzError = {0.1, 0.1}, dxy = {0.01, 0.01}, dxyError = {0.1, 0.1};
  std::vector<float> chi2 = {1, 1}, ecal = {2, 1}, hcal = {3, 1};
  std::vector<bool> highPurity = {true, true};
  std::vector<char> charge = {1, -1}, layers = {10, 10}, hits = {12, 12};
  float arrayPt[] = {4, 1}, arrayPtError[] = {0.04, 0.02}, arrayPhi[] = {0.2, 0.25}, arrayEta[] = {0.15, 0.18};
  float arrayDz[] = {0.01, 0.01}, arrayDzError[] = {0.1, 0.1}, arrayDxy[] = {0.01, 0.01}, arrayDxyError[] = {0.1, 0.1};
  float arrayChi2[] = {10, 10}, arrayEcal[] = {2, 1}, arrayHcal[] = {3, 1};
  bool arrayPurity[] = {true, true};
  unsigned char arrayNdof[] = {10, 10}, arrayLayers[] = {10, 10}, arrayHits[] = {12, 12};
  int arrayCharge[] = {1, -1};
  if(aod){
    tracks.Branch("nTrk", &nTrk, "nTrk/I");
    tracks.Branch("trkPt", arrayPt, "trkPt[nTrk]/F");
    tracks.Branch("trkPtError", arrayPtError, "trkPtError[nTrk]/F");
    tracks.Branch("trkPhi", arrayPhi, "trkPhi[nTrk]/F");
    tracks.Branch("trkEta", arrayEta, "trkEta[nTrk]/F");
    tracks.Branch("highPurity", arrayPurity, "highPurity[nTrk]/O");
    tracks.Branch("trkDz1", arrayDz, "trkDz1[nTrk]/F");
    tracks.Branch("trkDzError1", arrayDzError, "trkDzError1[nTrk]/F");
    tracks.Branch("trkDxy1", arrayDxy, "trkDxy1[nTrk]/F");
    tracks.Branch("trkDxyError1", arrayDxyError, "trkDxyError1[nTrk]/F");
    tracks.Branch("trkChi2", arrayChi2, "trkChi2[nTrk]/F");
    tracks.Branch("trkNdof", arrayNdof, "trkNdof[nTrk]/b");
    tracks.Branch("trkNlayer", arrayLayers, "trkNlayer[nTrk]/b");
    tracks.Branch("trkNHit", arrayHits, "trkNHit[nTrk]/b");
    tracks.Branch("pfEcal", arrayEcal, "pfEcal[nTrk]/F");
    tracks.Branch("pfHcal", arrayHcal, "pfHcal[nTrk]/F");
    tracks.Branch("trkCharge", arrayCharge, "trkCharge[nTrk]/I");
  } else {
    // Deliberately omit nTrk in Run-3 fixtures: the vectors define their size.
    if(!run3) tracks.Branch("nTrk", &nTrk, "nTrk/I");
    tracks.Branch("trkPt", &pt);
    tracks.Branch("trkPtError", &ptError);
    tracks.Branch("trkPhi", &phi);
    tracks.Branch("trkEta", &eta);
    tracks.Branch("highPurity", &highPurity);
    tracks.Branch("trkDzFirstVtx", &dz);
    tracks.Branch("trkDzErrFirstVtx", &dzError);
    tracks.Branch("trkDxyFirstVtx", &dxy);
    tracks.Branch("trkDxyErrFirstVtx", &dxyError);
    if(issue != "trkCharge") tracks.Branch("trkCharge", &charge);
    if(!run3){
      tracks.Branch("trkNormChi2", &chi2);
      tracks.Branch("trkNLayers", &layers);
      tracks.Branch("trkNHits", &hits);
      tracks.Branch("pfEcal", &ecal);
      tracks.Branch("pfHcal", &hcal);
    }
  }

  file.cd();
  std::unique_ptr<TTree> particles;
  std::vector<float> particlePt = {4.2, 1.1, 0.5}, particlePhi = {0.2, 0.25, 2.0}, particleEta = {0.15, 0.18, 1.2};
  std::vector<int> particleCharge = {1, -1, 0}, subevent = {0, 1, 0};
  std::string wrongChargeType = "invalid generator charge schema";
  if(isMC && issue != "genTree"){
    file.mkdir("HiGenParticleAna")->cd();
    particles.reset(new TTree("hi", "Synthetic generator particles"));
    particles->SetDirectory(nullptr);
    particles->Branch("pt", &particlePt);
    particles->Branch("phi", &particlePhi);
    particles->Branch("eta", &particleEta);
    if(issue == "genChargeType") particles->Branch("chg", &wrongChargeType);
    else particles->Branch("chg", &particleCharge);
    if(issue != "sube") particles->Branch("sube", &subevent);
  }
  for(int iEvent = 0; iEvent < 2; ++iEvent){
    if(iEvent == 1){ vz = 3.5; weight = 0.5; hiBin = 140; }
    if(issue == "vectorSize") charge.resize(1);
    event.Fill();
    jets.Fill();
    if(issue != "eventCount" || iEvent == 0) tracks.Fill();
    if(particles) particles->Fill();
  }
  file.GetDirectory("hiEvtAnalyzer")->cd(); event.Write();
  file.GetDirectory(jetDirectory)->cd(); jets.Write();
  file.GetDirectory(aod ? "ppTrack" : heavyIon ? "PbPbTracks" : "ppTracks")->cd(); tracks.Write();
  if(particles){ file.GetDirectory("HiGenParticleAna")->cd(); particles->Write(); }
}

#ifndef EEC_READER_FIXTURES_ONLY
namespace {
int checks = 0;
void Check(bool condition, const std::string& context){
  ++checks;
  if(!condition) throw std::runtime_error("Check failed: " + context);
}
bool Near(double a, double b){ return std::abs(a - b) < 1e-5; }
void ExpectFailure(const std::function<void()>& action, const std::string& fragment){
  try { action(); }
  catch(const std::exception& error){
    Check(std::string(error.what()).find(fragment) != std::string::npos, error.what());
    return;
  }
  throw std::runtime_error("Expected failure containing: " + fragment);
}

void ValidateFixture(const std::string& path, int dataType, int jetType, int axis, bool aod = false){
  const bool isMC = dataType == 2 || dataType == 3 || dataType == 10 || dataType == 12;
  const bool ppRef = dataType == 11 || dataType == 12;
  TFile file(path.c_str(), "READ");
  // The reader contains sizable legacy arrays, so keep it off the stack.
  std::unique_ptr<UnfoldingForestReader> reader(new UnfoldingForestReader(dataType, jetType, axis));
  reader->ReadForestFromFile(&file);
  Check(reader->GetNEvents() == 2, "event count");
  reader->GetEvent(0);
  Check(Near(reader->GetVz(), 2.5), "vertex z");
  Check(Near(reader->GetCentrality(), ppRef ? -0.5 : 30), "centrality");
  Check(Near(reader->GetEventWeight(), isMC ? 0.25 : 1), "MC weight / data unity");
  Check(Near(reader->GetPtHat(), isMC ? 80 : 0), "pthat");
  Check(reader->GetNJets() == 2, "reconstructed jet count");
  Check(Near(reader->GetJetPt(0), 128) && Near(reader->GetJetRawPt(0), 120), "raw and corrected pT remain distinct");
  Check(Near(reader->GetJetPhi(0), axis ? 0.25 : 0.2), "reconstructed jet axis phi");
  Check(Near(reader->GetJetEta(0), axis ? 0.18 : 0.15), "reconstructed jet axis eta");
  Check(reader->GetNTracks() == 2, "track count from vectors or AOD arrays");
  Check(Near(reader->GetTrackPt(0), 4) && Near(reader->GetTrackPtError(0), 0.04), "track pT and uncertainty");
  Check(reader->GetTrackHighPurity(0), "track quality");
  Check(reader->GetTrackCharge(1) == -1, "negative reconstructed charge");
  Check(Near(reader->GetTrackVertexDistanceZ(0), 0.01), "track vertex association");
  if(dataType >= 9 && !aod){
    ExpectFailure([&]{ reader->GetTrackEnergyEcal(0); }, "pfEcal is unavailable");
    ExpectFailure([&]{ reader->GetTrackNormalizedChi2(0); }, "trkNormChi2 is unavailable");
  } else {
    Check(Near(reader->GetTrackNormalizedChi2(0), 1), "legacy track normalized chi2");
    Check(reader->GetNHitsTrackerLayer(0) == 10, "legacy track layers");
    Check(Near(reader->GetTrackEnergyEcal(0), 2), "legacy track calorimeter energy");
  }
  if(isMC){
    Check(reader->GetNGeneratorJets() == 1, "generated jet count");
    Check(reader->HasMatchingGenJet(0) && !reader->HasMatchingGenJet(1), "matched and unmatched reconstructed jets");
    Check(reader->GetMatchingGenIndex(0) == 0 && reader->GetMatchingGenIndex(1) == -1, "reco to gen indices");
    Check(reader->HasMatchingRecoJet(0) && reader->GetMatchingRecoIndex(0) == 0, "gen to reco matching");
    Check(Near(reader->GetMatchedGenPt(0), 130) && Near(reader->GetMatchedRecoPt(0), 120), "matched raw reco and gen pT");
    Check(Near(reader->GetMatchedGenPhi(0), axis ? 0.22 : 0.2), "matched gen jet axis");
    Check(Near(reader->GetMatchedRecoEta(0), axis ? 0.18 : 0.15), "matched reco jet axis");
    Check(reader->GetNGenParticles() == 3 && Near(reader->GetGenParticlePt(0), 4.2), "generator particle vectors");
    Check(reader->GetGenParticleCharge(1) == -1 && reader->GetGenParticleSubevent(1) == 1, "generator charge and subevent");
  } else {
    Check(reader->GetNGeneratorJets() == 0 && reader->GetNGenParticles() == 0, "data has no generator entries");
    Check(reader->GetMatchingGenIndex(0) == -1 && Near(reader->GetMatchedGenPt(0), 0), "data matching guards");
  }
  reader->GetEvent(1);
  Check(Near(reader->GetVz(), 3.5) && Near(reader->GetEventWeight(), isMC ? 0.5 : 1), "second event values");
  ExpectFailure([&]{ reader->GetEvent(2); }, "outside the open forest");
  reader->BurnForest();
  reader->BurnForest();
  Check(reader->GetNEvents() == 0 && reader->GetNTracks() == 0, "repeated cleanup");
}
} // namespace

int main(int argc, char** argv){
  try {
    const std::string directory = argc > 1 ? argv[1] : "/tmp/eecRun3ReaderValidation";
    gSystem->mkdir(directory.c_str(), true);
    for(int dataType : {10, 12}){
      const int nJetTypes = dataType == 10 ? 4 : 3;
      for(int jetType = 0; jetType < nJetTypes; ++jetType){
        const std::string path = directory + "/mc" + std::to_string(dataType) + "_jet" + std::to_string(jetType) + ".root";
        WriteRun3ReaderFixture(path, dataType, jetType);
        for(int axis : {0, 1}) ValidateFixture(path, dataType, jetType, axis);
      }
    }
    for(int dataType : {9, 11, 2, 3}){
      const std::string path = directory + "/type" + std::to_string(dataType) + ".root";
      WriteRun3ReaderFixture(path, dataType);
      ValidateFixture(path, dataType, 1, 1);
    }
    const std::string legacyAOD = directory + "/legacyAOD.root";
    WriteRun3ReaderFixture(legacyAOD, 2, 1, "", true);
    ValidateFixture(legacyAOD, 2, 1, 0, true);
    for(const std::string issue : {"genTree", "weight", "refpt", "trkCharge", "sube", "genChargeType", "vectorSize", "eventCount", "overflowJets", "WTA"}){
      const std::string path = directory + "/bad_" + issue + ".root";
      WriteRun3ReaderFixture(path, 10, 1, issue);
      const std::string fragment = issue == "genTree" ? "HiGenParticleAna/hi" :
                                   issue == "genChargeType" ? "incompatible branch hi/chg" :
                                   issue == "vectorSize" ? "inconsistent particle vector trkCharge" :
                                   issue == "eventCount" ? "event counts differ" :
                                   issue == "overflowJets" ? "jet multiplicity exceeds" : issue;
      ExpectFailure([&]{ ValidateFixture(path, 10, 1, 1); }, fragment);
      if(issue == "WTA") ValidateFixture(path, 10, 1, 0);
    }
    for(int dataType : {-1, 4, 7, 8, 13}){
      ExpectFailure([&]{ std::unique_ptr<UnfoldingForestReader> reader(new UnfoldingForestReader(dataType, 1, 0)); }, "unsupported data type");
    }
    std::unique_ptr<UnfoldingForestReader> empty(new UnfoldingForestReader(10, 1, 0));
    empty->BurnForest();
    Check(empty->GetNEvents() == 0, "unopened forest");
    ExpectFailure([&]{ empty->ReadForestFromFile(nullptr); }, "cannot open input forest");
    ExpectFailure([&]{ empty->ReadForestFromFileList({}); }, "empty input file list");
    std::cout << "PASS: " << checks << " synthetic reader checks; fixtures in " << directory << std::endl;
    return 0;
  } catch(const std::exception& error){
    std::cerr << "FAIL: " << error.what() << std::endl;
    return 1;
  }
}
#endif

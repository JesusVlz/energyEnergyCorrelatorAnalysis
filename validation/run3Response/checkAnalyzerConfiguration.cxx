// Run after make, from the repository root:
// g++ $(root-config --cflags) -fno-access-control -I. \
//   validation/run3Response/checkAnalyzerConfiguration.cxx src/*.o \
//   $(root-config --libs) -lEG -o /tmp/checkRun3Analyzer
// /tmp/checkRun3Analyzer
// Private access is enabled only for this regression executable, so production
// interfaces need not expose configuration internals or track-selection helpers.

#include "../../src/EECAnalyzer.h"

#define EEC_READER_FIXTURES_ONLY
#include "validateUnfoldingForestReader.cxx"

#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

void Set(ConfigurationCard& card, const char* key, float value){
  if(!card.HasKey(key)) card.ReadInputLine((std::string(key)+" "+std::to_string(value)).c_str());
  TVector* entry = card.GetVector(key);
  (*entry)[entry->GetLwb()] = value;
}

void DisableJetSelection(ConfigurationCard& card){
  // This fixture isolates response and track behavior; native PF composition
  // and event-veto coverage are exercised by the separate jet-selection suite.
  Set(card,"Run3JetID",0); Set(card,"ApplyJetVeto",0);
}

void Require(bool condition, const char* message){
  if(!condition) throw std::runtime_error(message);
}

void ExpectRejected(const char* cardName, const char* key, float value, const char* diagnostic){
  ConfigurationCard card(cardName);
  DisableJetSelection(card);
  Set(card, key, value);
  try {
    EECAnalyzer analyzer({}, &card, 1, -1);
  } catch(const std::invalid_argument& error){
    Require(std::string(error.what()).find(diagnostic) != std::string::npos,
            "Configuration failed without the expected diagnostic");
    return;
  }
  throw std::runtime_error("Unsupported Run-3 configuration was accepted");
}

void CheckTrackSelection(const char* cardName){
  ConfigurationCard card(cardName);
  DisableJetSelection(card);
  auto analyzer = std::make_unique<EECAnalyzer>(std::vector<TString>{}, &card, 1, -1);
  auto ordinary = std::make_unique<HighForestReader>();
  auto response = std::make_unique<UnfoldingForestReader>();
  ordinary->fIsMiniAOD = true;
  response->fIsMiniAOD = true;

  // Omit all calorimeter, normalized-chi2, and hit/layer fields. Run-3
  // selection must accept a good high-pT track without dereferencing them.
  std::vector<float> pt{50}, eta{0.2}, error{0.5}, dz{0.01}, dzError{0.1}, dxy{0.01}, dxyError{0.1};
  std::vector<bool> purity{true};
  std::vector<char> charge{1};
  ordinary->fTrackPtVector = response->fTrackPtVector = &pt;
  ordinary->fTrackEtaVector = response->fTrackEtaVector = &eta;
  ordinary->fTrackPtErrorVector = response->fTrackPtErrorVector = &error;
  ordinary->fHighPurityTrackVector = response->fHighPurityTrackVector = &purity;
  ordinary->fTrackVertexDistanceZVector = response->fTrackVertexDistanceZVector = &dz;
  ordinary->fTrackVertexDistanceZErrorVector = response->fTrackVertexDistanceZErrorVector = &dzError;
  ordinary->fTrackVertexDistanceXYVector = response->fTrackVertexDistanceXYVector = &dxy;
  ordinary->fTrackVertexDistanceXYErrorVector = response->fTrackVertexDistanceXYErrorVector = &dxyError;
  response->fTrackChargeVector = &charge;

  auto check = [&](bool expected){
    Require(analyzer->PassTrackCuts(static_cast<ForestReader*>(ordinary.get()), 0, nullptr, true) == expected,
            "Unexpected ordinary Run-3 track selection");
    Require(analyzer->PassTrackCuts(response.get(), 0) == expected,
            "Response track selection differs from the ordinary Run-3 path");
  };
  check(true);
  purity[0] = false; check(false); purity[0] = true;
  error[0] = 10; check(false); error[0] = 0.5;
  dz[0] = 1; check(false); dz[0] = 0.01;
  dxy[0] = 1; check(false); dxy[0] = 0.01;
  eta[0] = 2.5; check(false); eta[0] = 0.2;
  pt[0] = 0.1; check(false); pt[0] = 50;
  charge[0] = 0;
  Require(!analyzer->PassTrackCuts(response.get(), 0), "Neutral response track was accepted");
}

void CheckResponseFilling(const char* cardName, int dataType, const std::string& directory){
  const std::string path = directory + "/responseFilling" + std::to_string(dataType) + ".root";
  WriteRun3ReaderFixture(path, dataType);
  TFile input(path.c_str(), "READ");
  ConfigurationCard card(cardName);
  DisableJetSelection(card);
  auto analyzer = std::make_unique<EECAnalyzer>(std::vector<TString>{}, &card, 1, -1);
  // Empty normal input initializes the same JEC and histogram machinery used
  // in production, without requiring unrelated HLT/skim fixture branches.
  analyzer->RunAnalysis();
  auto reader = std::make_unique<UnfoldingForestReader>(dataType, 1, 1);
  reader->ReadForestFromFile(&input);
  analyzer->fUnfoldingForestReader = reader.get();
  for(int event = 0; event < reader->GetNEvents(); ++event){
    reader->GetEvent(event);
    analyzer->fTotalEventWeight = analyzer->GetVzWeight(reader->GetVz()) * analyzer->GetCentralityWeight(reader->GetHiBin());
    analyzer->FillUnfoldingResponse();
    analyzer->ConstructParticleResponses();
  }
  Require(analyzer->fHistograms->fhUnfoldingResponse->GetEntries() > 0, "Empty EEC response from a matched synthetic jet");
  Require(analyzer->fHistograms->fhJetPtUnfoldingResponse->GetEntries() > 0, "Empty jet-pT response from a matched synthetic jet");
  Require(analyzer->fHistograms->fhParticleDeltaRResponse->GetEntries() > 0, "Empty particle angular response from matched synthetic tracks");
  Require(analyzer->fHistograms->fhParticlePtResponse->GetEntries() > 0, "Empty particle-pT response from matched synthetic tracks");
  analyzer->fUnfoldingForestReader = nullptr;
}

} // namespace

int main(int argc, char** argv){
  try {
    const std::string directory = argc > 1 ? argv[1] : "/tmp/eecRun3AnalyzerFixtures";
    gSystem->mkdir(directory.c_str(), true);
    for(const char* card : {"cardEECOO.input", "cardEECOOMC.input", "cardEECpp5v36TeV.input", "cardEECpp5v36TeVMC.input"}){
      ExpectRejected(card, "DisableTrackPairEfficiencyCorrection", 0, "DisableTrackPairEfficiencyCorrection 1");
      ExpectRejected(card, "JetUncertainty", 1, "JetUncertainty 0");
      ExpectRejected(card, "JetUncertainty", 3, "JetUncertainty 0");
      ExpectRejected(card, "JetPtWeight", 1, "JetPtWeight 0");
      ExpectRejected(card, "SmearDeltaR", 1, "SmearDeltaR 0");
      ExpectRejected(card, "SmearEnergyWeight", 1, "SmearEnergyWeight 0");
    }
    for(const char* card : {"cardEECOO.input", "cardEECpp5v36TeV.input"}){
      ExpectRejected(card, "FilledHistograms", 128, "require MC");
      ExpectRejected(card, "FilledHistograms", 256, "require MC");
    }
    CheckTrackSelection("cardEECOOMC.input");
    CheckTrackSelection("cardEECpp5v36TeVMC.input");
    CheckResponseFilling("cardEECOOMC.input", 10, directory);
    CheckResponseFilling("cardEECpp5v36TeVMC.input", 12, directory);
    std::cout << "PASS: Run-3 guards, consistent track selections, and nonzero synthetic OO/ppRef EEC/jet/particle responses\n";
  } catch(const std::exception& error){
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}

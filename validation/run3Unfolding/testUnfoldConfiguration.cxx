// Synthetic ROOT metadata tests only: no detector calibration or physical response is supplied here.
#include "EECUnfoldConfiguration.h"
#include <TFile.h>
#include <TVectorT.h>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void WriteVector(const char* name, const std::vector<float>& values){
  TVectorT<float> vector(1, values.size());
  for(unsigned i = 0; i < values.size(); ++i) vector[i+1] = values[i];
  vector.Write(name);
}

void WriteFixture(const std::string& path, int dataType, int correlation = 3, int axis = 1, int exponent = 1, float centralityEdge = 30, bool response = true, bool weightMetadata = true, float radius = 0.4, float trackEdge = 1, float jetEdge = 120){
  TFile file(path.c_str(), "RECREATE");
  file.mkdir("JCard")->cd();
  WriteVector("DataType", {float(dataType)});
  WriteVector("McCorrelationType", {float(correlation)});
  WriteVector("JetType", {1});
  WriteVector("JetAxis", {float(axis)});
  WriteVector("JetRadius", {radius});
  if(weightMetadata) WriteVector("WeightExponent", {float(exponent)});
  WriteVector("CentralityBinEdges", {0, 10, centralityEdge, 90});
  WriteVector("JetPtBinEdgesEEC", {100, jetEdge, 140});
  WriteVector("TrackPtBinEdgesEEC", {0.7, trackEdge, 2});
  file.cd();
  if(response){
    // Configuration checks directories; unfolding itself checks each loaded histogram.
    file.mkdir("jetPtUnfoldingMeasured");
    file.mkdir("jetPtUnfoldingTruth");
    file.mkdir("jetPtUnfoldingResponse");
  }
}

void Require(bool condition, const char* message){
  if(!condition) throw std::runtime_error(message);
}

void Reject(const std::function<void()>& action, const char* message){
  try { action(); } catch(const std::invalid_argument&) { return; }
  throw std::runtime_error(message);
}
}

int main(int argc, char** argv){
  if(argc != 2){
    std::cerr << "Usage: testUnfoldConfiguration EXISTING_SCRATCH_DIRECTORY\n";
    return 2;
  }
  const std::string prefix = std::string(argv[1]) + "/unfoldConfig_";
  WriteFixture(prefix + "OOdata.root", 9);
  WriteFixture(prefix + "ppRefdata.root", 11);
  WriteFixture(prefix + "ppdata.root", 0);
  WriteFixture(prefix + "PbPbdata.root", 1);
  WriteFixture(prefix + "OOresponse.root", 10);
  WriteFixture(prefix + "ppRefresponse.root", 12);
  WriteFixture(prefix + "wrongMode.root", 10, 4);
  WriteFixture(prefix + "RecoGen.root", 10, 1);
  WriteFixture(prefix + "wrongAxis.root", 10, 3, 0);
  WriteFixture(prefix + "wrongExponent.root", 10, 3, 1, 2);
  WriteFixture(prefix + "wrongCentrality.root", 10, 3, 1, 1, 40);
  WriteFixture(prefix + "missingResponse.root", 10, 3, 1, 1, 30, false);
  WriteFixture(prefix + "missingMetadata.root", 10, 3, 1, 1, 30, true, false);
  WriteFixture(prefix + "wrongRadius.root", 10, 3, 1, 1, 30, true, true, 0.3);
  WriteFixture(prefix + "wrongTrack.root", 10, 3, 1, 1, 30, true, true, 0.4, 1.5);
  WriteFixture(prefix + "wrongJet.root", 10, 3, 1, 1, 30, true, true, 0.4, 1, 125);
  WriteFixture(prefix + "nonfiniteRadius.root", 10, 3, 1, 1, 30, true, true, std::numeric_limits<float>::quiet_NaN());
  WriteFixture(prefix + "RecoReco.root", 10, 0);
  WriteFixture(prefix + "GenReco.root", 10, 2);
  TFile ooFile((prefix + "OOdata.root").c_str()); EECCard oo(&ooFile);
  TFile ppRefFile((prefix + "ppRefdata.root").c_str()); EECCard ppRef(&ppRefFile);
  TFile ppFile((prefix + "ppdata.root").c_str()); EECCard pp(&ppFile);
  TFile pbpbFile((prefix + "PbPbdata.root").c_str()); EECCard pbpb(&pbpbFile);
  const TString ooResponse = (prefix + "OOresponse.root").c_str();

  Reject([&]{ EECUnfoldConfiguration c(&oo); }, "OO selected a legacy response");
  Reject([&]{ EECUnfoldConfiguration c(&ppRef); }, "ppRef selected a legacy response");
  Reject([&]{ EECUnfoldConfiguration c(&oo, 0, 0, 1, ooResponse); }, "Missing explicit iterations accepted");
  Reject([&]{ EECUnfoldConfiguration c(&oo, 0, 0, 1, ooResponse, -2); }, "Negative iterations accepted");
  Reject([&]{ EECUnfoldConfiguration c(&oo, 0, 0, 1, (prefix + "ppRefresponse.root").c_str(), 6); }, "Wrong MC species accepted");
  Reject([&]{ EECUnfoldConfiguration c(&oo, 0, 0, 1, (prefix + "OOdata.root").c_str(), 6); }, "Data response accepted");
  for(const char* name : {"wrongMode", "wrongAxis", "wrongExponent", "wrongCentrality", "missingResponse", "missingMetadata", "wrongRadius", "wrongTrack", "wrongJet", "nonfiniteRadius"}){
    Reject([&]{ EECUnfoldConfiguration c(&oo, 0, 0, 1, (prefix + name + ".root").c_str(), 6); }, name);
  }
  Reject([&]{ EECUnfoldConfiguration c(&oo, 0, 0, 1, (prefix + "not_present.root").c_str(), 6); }, "Nonexistent response accepted");
  EECUnfoldConfiguration explicitOO(&oo, 1, 0, 1, ooResponse, 6);
  Require(explicitOO.GetResponseFileName() == ooResponse, "Explicit response path changed");
  Require(explicitOO.GetNumberOfIterations({30,90}, 0.7) == 6, "Explicit iterations replaced by a legacy count");
  EECUnfoldConfiguration explicitPP(&ppRef, 0, 0, 1, (prefix + "ppRefresponse.root").c_str(), 9);
  Require(explicitPP.GetNumberOfIterations({0,10}, 0.7) == 9, "ppRef explicit iterations lost");
  EECUnfoldConfiguration recoGen(&oo, 0, 0, 1, (prefix + "RecoGen.root").c_str(), 6);
  Require(recoGen.GetNumberOfIterations({0,10}, 0.7) == 6, "Valid RecoGen-produced particle response rejected");
  for(const char* name : {"RecoReco", "GenReco"}){
    EECUnfoldConfiguration validMode(&oo, 0, 0, 1, (prefix + name + ".root").c_str(), 6);
    Require(validMode.GetNumberOfIterations({0,10}, 0.7) == 6, "Valid accompanying MC mode rejected");
  }
  EECUnfoldConfiguration legacyPP(&pp);
  EECUnfoldConfiguration legacyPbPb(&pbpb);
  Require(legacyPP.GetResponseFileName().Contains("ppMC2017"), "Legacy pp response changed");
  Require(legacyPbPb.GetResponseFileName().Contains("PbPbMC2018"), "Legacy PbPb response changed");
  Require(legacyPP.GetNumberOfIterations({0,10}, 2) == 4, "Legacy default iterations changed");

  const std::vector<std::pair<double,double>> legacyBins = {{120,140}, {140,160}, {160,180}, {180,200}};
  Require(EECUnfoldConfiguration::SelectUnfoldedJetPtBins(pp) == legacyBins, "Legacy selected jet pT bins changed");
  Require(EECUnfoldConfiguration::SelectUnfoldedJetPtBins(pbpb) == legacyBins, "Legacy PbPb selected jet pT bins changed");
  double nativeEdges[] = {60,100,120,160,200,240};
  oo.AddVector(EECCard::kJetPtBinEdgesEEC, 6, nativeEdges);
  ppRef.AddVector(EECCard::kJetPtBinEdgesEEC, 6, nativeEdges);
  const std::vector<std::pair<double,double>> nativeBins = {{120,160}, {160,200}};
  Require(EECUnfoldConfiguration::SelectUnfoldedJetPtBins(oo) == nativeBins, "OO native jet pT bins were not selected");
  Require(EECUnfoldConfiguration::SelectUnfoldedJetPtBins(ppRef) == nativeBins, "ppRef native jet pT bins were not selected");
  double partialEdges[] = {100,110,130,170,210,250};
  oo.AddVector(EECCard::kJetPtBinEdgesEEC, 6, partialEdges);
  Require(EECUnfoldConfiguration::SelectUnfoldedJetPtBins(oo) == std::vector<std::pair<double,double>>{{130,170}}, "Partial or outside jet pT bins selected");
  double emptyEdges[] = {100,250};
  oo.AddVector(EECCard::kJetPtBinEdgesEEC, 2, emptyEdges);
  Reject([&]{ EECUnfoldConfiguration::SelectUnfoldedJetPtBins(oo); }, "Empty Run-3 jet pT range accepted");
  double integratedEdges[] = {40,80,120,160,200};
  oo.AddVector(EECCard::kJetPtBinEdgesEEC, 5, integratedEdges);
  oo.AddOneDimensionalVector(EECCard::kFirstUnfoldedJetPtBin, 0);
  oo.AddOneDimensionalVector(EECCard::kLastUnfoldedJetPtBin, 2);
  Require(EECUnfoldConfiguration::SelectUnfoldedJetPtBins(oo) == std::vector<std::pair<double,double>>{{40,80},{80,120},{120,160}}, "Explicit integrated bins lost precedence over native fallback");
  oo.AddOneDimensionalVector(EECCard::kFirstUnfoldedJetPtBin, -1);
  Reject([&]{ EECUnfoldConfiguration::SelectUnfoldedJetPtBins(oo); }, "Partially configured range accepted");
  oo.AddOneDimensionalVector(EECCard::kFirstUnfoldedJetPtBin, 3);
  Reject([&]{ EECUnfoldConfiguration::SelectUnfoldedJetPtBins(oo); }, "Reversed range accepted");
  oo.AddOneDimensionalVector(EECCard::kFirstUnfoldedJetPtBin, 0);
  oo.AddOneDimensionalVector(EECCard::kLastUnfoldedJetPtBin, 4);
  Reject([&]{ EECUnfoldConfiguration::SelectUnfoldedJetPtBins(oo); }, "Out-of-bounds range accepted");
  oo.AddOneDimensionalVector(EECCard::kLastUnfoldedJetPtBin, 1.5);
  Reject([&]{ EECUnfoldConfiguration::SelectUnfoldedJetPtBins(oo); }, "Fractional range accepted");
  oo.AddOneDimensionalVector(EECCard::kLastUnfoldedJetPtBin, std::numeric_limits<float>::quiet_NaN());
  Reject([&]{ EECUnfoldConfiguration::SelectUnfoldedJetPtBins(oo); }, "Nonfinite range accepted");
  std::cout << "PASS: explicit OO/ppRef response selection, metadata rejection, native jet bins and legacy defaults.\n";
}

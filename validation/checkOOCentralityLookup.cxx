// Compile and run from the repository root:
// g++ $(root-config --cflags) -O1 -ffunction-sections -fdata-sections \
//   validation/checkOOCentralityLookup.cxx plotting/EECHistogramManager.cxx \
//   plotting/EECCard.cxx plotting/AlgorithmLibrary.cxx \
//   plotting/EECBackgroundScale.cxx plotting/EECSignalToBackgroundUnfoldingScale.cxx \
//   -Wl,--gc-sections $(root-config --libs) \
//   -o /tmp/checkOOCentralityLookup && /tmp/checkOOCentralityLookup

#include "../plotting/EECHistogramManager.h"

#include <TMemFile.h>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

void WriteVector(const char* name, std::initializer_list<float> values){
  TVectorF entry(1, static_cast<int>(values.size()));
  int index = 1;
  for(float value : values) entry[index++] = value;
  entry.Write(name);
}

void Require(bool condition, const char* message){
  if(!condition) throw std::runtime_error(message);
}

void RequireContent(TH1D* histogram, double expected){
  Require(histogram != nullptr, "EEC histogram was not found");
  Require(histogram->GetBinContent(1) == expected, "EEC lookup returned the wrong centrality bin");
}

void CheckDataType(int dataType, bool hasCentrality, bool legacyNames){
  TMemFile input("centralityLookup.root", "RECREATE");
  input.mkdir("JCard")->cd();
  WriteVector("DataType", {static_cast<float>(dataType)});
  WriteVector("McCorrelationType", {0});
  WriteVector("BackgroundMethods", {0});
  WriteVector("CentralityBinEdges", {0, 10, 30, 50, 90});
  WriteVector("TrackPtBinEdges", {1, 2});
  WriteVector("JetPtBinEdgesEEC", {120, 140});
  WriteVector("TrackPtBinEdgesEEC", {1, 2});

  TDirectory* raw = input.mkdir("energyEnergyCorrelator");
  if(!legacyNames) raw = raw->mkdir("SignalConePair");
  TDirectory* processed = input.mkdir("energyEnergyCorrelatorProcessed");
  for(int centrality = 0; centrality < 4; ++centrality){
    raw->cd();
    const char* rawName = legacyNames
      ? "energyEnergyCorrelator_C%dT0J0"
      : "energyEnergyCorrelatorSignalConePair_C%dT0J0";
    TH1D rawHistogram(Form(rawName, centrality), "", 1, 0, 1);
    rawHistogram.SetBinContent(1, 100 + centrality);
    rawHistogram.Write();

    processed->cd();
    TH1D processedHistogram(Form("energyEnergyCorrelatorUnfolded_C%dT0J0", centrality), "", 1, 0, 1);
    processedHistogram.SetBinContent(1, 200 + centrality);
    processedHistogram.Write();
  }

  // The manager contains large histogram-pointer arrays and must live on the heap.
  auto manager = std::make_unique<EECHistogramManager>(&input);
  const int eec = EECHistogramManager::kEnergyEnergyCorrelator;
  const std::pair<int,int> jetPt{120, 140};
  Require(manager->GetLastCentralityBin() == (hasCentrality ? 3 : 0), "Wrong loaded centrality range");
  RequireContent(manager->GetHistogramEnergyEnergyCorrelator(eec, {30, 50}, jetPt, 1.0), hasCentrality ? 102 : 100);
  RequireContent(manager->GetHistogramEnergyEnergyCorrelator(eec, {50, 90}, jetPt, 1.0), hasCentrality ? 103 : 100);

  // Processed EECs have an index getter, without a boundary-based overload.
  // Confirm that its existing OO centrality handling retains distinct unfolded bins.
  const int unfolded = EECHistogramManager::kEnergyEnergyCorrelatorUnfolded;
  RequireContent(manager->GetHistogramEnergyEnergyCorrelatorProcessed(eec, 2, 0, 0, unfolded), 202);
  RequireContent(manager->GetHistogramEnergyEnergyCorrelatorProcessed(eec, 3, 0, 0, unfolded), 203);

  if(hasCentrality){
    Require(manager->GetHistogramEnergyEnergyCorrelator(eec, {25, 50}, jetPt, 1.0) == nullptr,
            "Unknown centrality boundaries must not index histogram arrays");
  }
  Require(manager->GetHistogramEnergyEnergyCorrelator(eec, {30, 50}, {125, 140}, 1.0) == nullptr,
          "Unknown jet boundaries must not index histogram arrays");
  Require(manager->GetHistogramEnergyEnergyCorrelator(eec, {30, 50}, jetPt, 0.5) == nullptr,
          "Track thresholds outside the card must not index histogram arrays");
}

} // namespace

int main(){
  try {
    for(bool legacyNames : {false, true}){
      for(int dataType : {1, 3, 9, 10}) CheckDataType(dataType, true, legacyNames);
      for(int dataType : {0, 2, 11, 12}) CheckDataType(dataType, false, legacyNames);
    }
    std::cout << "PASS: OO data/MC centrality lookup, processed bins, PbPb/pp behavior, both naming conventions, and invalid bounds\n";
  } catch(const std::exception& error){
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}

#include <TFile.h>
#include <TH1.h>
#include <THnBase.h>
#include <TSystem.h>
#include <TVectorT.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace {

void OOMCCheck(bool condition, const std::string& message, int& failures) {
  if (!condition) {
    std::cerr << "OOMC validation ERROR: " << message << std::endl;
    ++failures;
  }
}

void OOMCCheckCard(TFile& file, const char* key, int expected, int& failures) {
  const std::string path = std::string("JCard/") + key;
  auto* value = dynamic_cast<TVectorT<float>*>(file.Get(path.c_str()));
  OOMCCheck(value && value->GetNrows() == 1, "missing scalar card entry " + path, failures);
  if (!value || value->GetNrows() != 1) return;
  const double actual = (*value)[value->GetLwb()];
  OOMCCheck(std::isfinite(actual) && actual == expected,
            path + " must equal " + std::to_string(expected) + ", got " + std::to_string(actual),
            failures);
}

TH1* OOMCCheckHistogram(TFile& file, const char* name, bool requireEntries, int& failures) {
  auto* histogram = dynamic_cast<TH1*>(file.Get(name));
  OOMCCheck(histogram && histogram->GetDimension() == 1,
            std::string("missing one-dimensional histogram ") + name, failures);
  if (!histogram || histogram->GetDimension() != 1) return nullptr;
  const double entries = histogram->GetEntries();
  OOMCCheck(std::isfinite(entries) && entries >= 0,
            std::string(name) + " has invalid entries", failures);
  if (requireEntries) OOMCCheck(entries > 0, std::string(name) + " has no entries", failures);
  for (int bin = 0; bin < histogram->GetNcells(); ++bin) {
    if (!std::isfinite(histogram->GetBinContent(bin)) ||
        !std::isfinite(histogram->GetBinError(bin))) {
      OOMCCheck(false, std::string(name) + " has nonfinite content/error in bin " +
                        std::to_string(bin), failures);
      break;
    }
  }
  std::cout << "OOMC histogram " << name << " entries=" << entries << std::endl;
  return histogram;
}

void OOMCCheckSparse(TFile& file, const char* name, bool requireEntries, int& failures) {
  auto* histogram = dynamic_cast<THnBase*>(file.Get(name));
  OOMCCheck(histogram != nullptr, std::string("missing THnBase histogram ") + name, failures);
  if (!histogram) return;
  const double entries = histogram->GetEntries();
  OOMCCheck(std::isfinite(entries) && entries >= 0 &&
              std::isfinite(histogram->GetSumw()) && std::isfinite(histogram->GetSumw2()) &&
              histogram->GetNdimensions() > 0 && histogram->GetNbins() >= 0,
            std::string(name) + " has invalid entries, weight sums, or dimensions", failures);
  if (requireEntries) OOMCCheck(entries > 0, std::string(name) + " has no entries", failures);
  std::cout << "OOMC histogram " << name << " entries=" << entries
            << " bins=" << histogram->GetNbins() << std::endl;
}

}  // namespace

// Run with: root -l -b -q 'condor/checkOOMCOutput.C("output.root","RecoReco")'
// Response entries may be zero in an individual minimum-bias input file.
void checkOOMCOutput(const char* filename, const char* mode, int expectedDataType = 10) {
  int failures = 0;
  const std::string selectedMode = mode ? mode : "";
  int correlationType = -1;
  if (selectedMode == "RecoReco") correlationType = 0;
  else if (selectedMode == "RecoGen") correlationType = 1;
  else if (selectedMode == "GenReco") correlationType = 2;
  else if (selectedMode == "GenGen" || selectedMode == "Response") correlationType = 3;
  if (correlationType < 0 || !filename || !*filename ||
      (expectedDataType != 10 && expectedDataType != 12)) {
    std::cerr << "OOMC validation ERROR: provide an output file and mode "
              << "RecoReco, RecoGen, GenReco, GenGen, or Response, and MC data type 10 or 12" << std::endl;
    gSystem->Exit(2);
    return;
  }

  std::unique_ptr<TFile> file(TFile::Open(filename, "READ"));
  if (!file || file->IsZombie() || !file->IsOpen() || file->TestBit(TFile::kRecovered)) {
    std::cerr << "OOMC validation ERROR: unreadable, zombie, or recovered file "
              << filename << std::endl;
    gSystem->Exit(1);
    return;
  }
  const bool response = selectedMode == "Response";
  OOMCCheckCard(*file, "DataType", expectedDataType, failures);
  OOMCCheckCard(*file, "TriggerSelection", 0, failures);
  OOMCCheckCard(*file, "McCorrelationType", correlationType, failures);
  if (response) {
    OOMCCheckCard(*file, "FilledHistograms", 479, failures);
    OOMCCheckCard(*file, "MatchJets", 1, failures);
  }

  auto* events = OOMCCheckHistogram(*file, "nEvents", true, failures);
  if (events) {
    // The event counter has bin centers 0=All, 1=Triggered, ..., 8=vz cut.
    OOMCCheck(events->GetNbinsX() >= 9, "nEvents is missing event-cut bins", failures);
    if (events->GetNbinsX() >= 9) {
      const double all = events->GetBinContent(1);
      const double triggered = events->GetBinContent(2);
      const double selected = events->GetBinContent(9);
      OOMCCheck(all > 0 && triggered > 0 && selected > 0,
                "no events reached All, Triggered, or the vertex-z cut", failures);
      OOMCCheck(all == triggered, "TriggerSelection 0 must keep every event at the trigger step", failures);
      std::cout << "OOMC events All=" << all << " Triggered=" << triggered
                << " VzCut=" << selected << std::endl;
    }
  }
  OOMCCheckHistogram(*file, "vertexZ", true, failures);
  OOMCCheckHistogram(*file, "pthat", true, failures);
  for (auto name : {"vertexZweighted", "centrality", "centralityWeighted", "pthatWeighted"})
    OOMCCheckHistogram(*file, name, false, failures);

  OOMCCheckSparse(*file, "inclusiveJet", !response, failures);
  OOMCCheckSparse(*file, "energyEnergyCorrelator", false, failures);
  if (response) {
    OOMCCheckSparse(*file, "jetPtClosure", false, failures);
    OOMCCheckSparse(*file, "jetPtUnfoldingMeasured", false, failures);
    OOMCCheckSparse(*file, "jetPtUnfoldingTruth", false, failures);
    OOMCCheckSparse(*file, "jetPtUnfoldingResponse", false, failures);
    OOMCCheckSparse(*file, "particleDeltaRResponseMatrix", false, failures);
    OOMCCheckSparse(*file, "particlePtResponseMatrix", false, failures);
  }

  file->Close();
  std::cout << "OOMC validation " << (failures ? "FAIL" : "PASS")
            << " mode=" << selectedMode << " file=" << filename
            << " errors=" << failures << std::endl;
  gSystem->Exit(failures ? 1 : 0);
}

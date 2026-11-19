#include "../../src/Run3JetQA.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
const double pi = TMath::Pi();

void Require(bool condition, const std::string& message) {
  if(!condition) throw std::runtime_error(message);
}

void Equal(double actual, double expected, const std::string& message) {
  if(!std::isfinite(actual) ||
     std::abs(actual - expected) > 1e-10 * std::max(1., std::abs(expected)))
    throw std::runtime_error(message + ": got " + std::to_string(actual) +
                             ", expected " + std::to_string(expected));
}

struct TemporaryFiles {
  std::string directory;
  std::vector<std::string> paths;
  TemporaryFiles() {
    char pattern[] = "/tmp/eec-run3-jet-qa-XXXXXX";
    const char* result = mkdtemp(pattern);
    Require(result != nullptr, "Cannot create temporary directory");
    directory = result;
  }
  ~TemporaryFiles() {
    for(const auto& path : paths) std::remove(path.c_str());
    std::remove(directory.c_str());
  }
  std::string Path(const std::string& name) {
    paths.push_back(directory + "/" + name + ".root");
    return paths.back();
  }
};

struct Jet {
  double pt, eta, phi, wrappedPhi, centrality, weight;
  // The three bits encode trackMax pass, ID pass, and event-veto pass.
  unsigned flags;
};

// Explicit truth tables keep the expectations independent of the Fill predicates.
// Columns: both enabled, ID disabled, veto disabled, both disabled.
const unsigned selectionMasks[8][4] = {
  {  1,  1,  1, 1}, { 65, 65,  1, 1},
  {  5,  1,  5, 1}, { 85, 65,  5, 1},
  {  3,  3,  3, 3}, { 67, 99,  3, 3},
  { 15,  3, 15, 3}, {127, 99, 15, 3}
};
const bool legacySelected[8][4] = {
  {false, false, false, false}, {false, false, false, false},
  {false, false, false, false}, {false, false, false, false},
  {false, false, false, false}, {false, true,  false, true },
  {false, false, false, false}, {true,  true,  true,  true }
};

struct BinExpectation {
  double sum = 0;
  double sumw2 = 0;
};
struct ExpectedHistogram {
  std::map<int, BinExpectation> bins;
  int entries = 0;
  void Fill(TH3D& histogram, double x, double y, double z, double weight) {
    auto& bin = bins[histogram.FindBin(x, y, z)];
    bin.sum += weight;
    bin.sumw2 += weight * weight;
    ++entries;
  }
  void Check(const TH3D& histogram) const {
    const std::string name = histogram.GetName();
    Require(histogram.GetSumw2N() == histogram.GetNcells(), name + " missing Sumw2");
    Equal(histogram.GetEntries(), entries, name + " entries");
    for(int cell = 0; cell < histogram.GetNcells(); ++cell) {
      const auto found = bins.find(cell);
      const BinExpectation expected = found == bins.end() ? BinExpectation{} : found->second;
      Equal(histogram.GetBinContent(cell), expected.sum, name + " content bin " + std::to_string(cell));
      Equal(histogram.GetSumw2()->At(cell), expected.sumw2, name + " Sumw2 bin " + std::to_string(cell));
    }
  }
};

TH3D& Histogram(TFile& file, const std::string& name) {
  auto* histogram = dynamic_cast<TH3D*>(file.Get(("run3JetQA/" + name).c_str()));
  Require(histogram != nullptr, "Missing saved TH3D " + name);
  return *histogram;
}

void CheckGeometry(const TH3D& histogram, int coordinateBins, double low, double high) {
  const std::string name = histogram.GetName();
  Require(histogram.GetNbinsX() == coordinateBins, name + " coordinate bins");
  Equal(histogram.GetXaxis()->GetXmin(), low, name + " coordinate minimum");
  Equal(histogram.GetXaxis()->GetXmax(), high, name + " coordinate maximum");
  Require(histogram.GetNbinsY() == 7, name + " selection bins");
  Equal(histogram.GetYaxis()->GetXmin(), -.5, name + " selection minimum");
  Equal(histogram.GetYaxis()->GetXmax(), 6.5, name + " selection maximum");
  Require(histogram.GetNbinsZ() == 101, name + " centrality bins");
  Equal(histogram.GetZaxis()->GetXmin(), -1., name + " centrality minimum");
  Equal(histogram.GetZaxis()->GetXmax(), 100., name + " centrality maximum");
}

void CheckStageTotals(const TH3D& histogram, const TH3D& pt) {
  // Include both coordinate and centrality flow bins when comparing exposure.
  for(int stage = 1; stage <= 7; ++stage) {
    for(int activity = 0; activity <= pt.GetNbinsZ() + 1; ++activity) {
      double sum = 0, sumw2 = 0, ptSum = 0, ptSumw2 = 0;
      for(int coordinate = 0; coordinate <= histogram.GetNbinsX() + 1; ++coordinate) {
        const int bin = histogram.GetBin(coordinate, stage, activity);
        sum += histogram.GetBinContent(bin);
        sumw2 += histogram.GetSumw2()->At(bin);
      }
      for(int coordinate = 0; coordinate <= pt.GetNbinsX() + 1; ++coordinate) {
        const int bin = pt.GetBin(coordinate, stage, activity);
        ptSum += pt.GetBinContent(bin);
        ptSumw2 += pt.GetSumw2()->At(bin);
      }
      Equal(sum, ptSum, std::string(histogram.GetName()) + " stage exposure versus pT");
      Equal(sumw2, ptSumw2, std::string(histogram.GetName()) + " stage Sumw2 versus pT");
    }
  }
}

void TestConfiguration(TemporaryFiles& temporary, const std::string& mapPath, int configuration) {
  const bool idEnabled = configuration == 0 || configuration == 2;
  const bool vetoEnabled = configuration == 0 || configuration == 1;
  ConfigurationCard card;
  card.ReadInputLine("DataType 10");
  card.ReadInputLine("Run3JEC 1");
  card.ReadInputLine("JetPtBinEdgesEEC 0 100 200 1000");
  card.ReadInputLine(idEnabled ? "Run3JetID 1" : "Run3JetID 0");
  card.ReadInputLine(vetoEnabled ? "ApplyJetVeto 1" : "ApplyJetVeto 0");
  if(vetoEnabled) {
    card.ReadInputLine(("JetVetoMapFile " + mapPath).c_str());
    card.ReadInputLine("JetVetoMapHistogram mask");
  }
  Run3JetQA qa(&card);
  const std::vector<Jet> jets = {
    {-5., -2.5, pi, -pi, -2., 2., 0},
    {1205., 2.5, -pi, -pi, 100., -0.5, 1},
    {101., -1.2, 5*pi, -pi, 25., 3., 2},
    {121., 1.2, -5*pi, -pi, 25., 1.5, 3},
    {141., -.2, 2*pi+.4, .4, -1., 4., 4},
    {161., .2, -2*pi-.4, -.4, -1., -2., 5},
    {181., -.7, .4, .4, 50., 0.25, 6},
    {201., .7, -.4, -.4, 50., 2.5, 7},
    {-3., -2.5, pi, -pi, -2., 0.75, 7},
    {1001., 2.5, -pi, -pi, 100., 1.25, 7},
    {131., -2., 5*pi, -pi, 25., 1., 7},
    {131., 2., -5*pi, -pi, 25., -0.25, 7},
    {131., .2, 2*pi+.4, .4, 25., 0., 7},
    {131., .2, -2*pi-.4, -.4, 25., 2., 7},
    // Shared bins with unequal signed weights distinguish Sumw2 from counts.
    {131., .2, -2*pi-.4, -.4, 25., -0.5, 7}
  };
  for(const auto& jet : jets)
    qa.Fill(jet.pt, 100., 20., jet.centrality, jet.weight,
            jet.flags & 4, jet.flags & 2, jet.flags & 1,
            {100., jet.pt}, jet.eta, jet.phi);

  // Invalid axes must throw even for rejected jets, without partially filling
  // the original pT or track-fraction diagnostics. Saved expectations exclude
  // every invalid call, so the checks below also catch partial state changes.
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  for(bool selected : {false, true}) {
    for(double invalid : {nan, infinity, -infinity}) {
      for(bool invalidEta : {false, true}) {
        bool threw = false;
        try {
          qa.Fill(150., 100., 20., 25., 17., selected, selected, selected,
                  {100., 150.}, invalidEta ? invalid : .2, invalidEta ? .4 : invalid);
        } catch(const std::invalid_argument&) { threw = true; }
        Require(threw, "Nonfinite axis accepted");
      }
    }
  }

  const std::string path = temporary.Path("qa_" + std::to_string(configuration));
  {
    TFile output(path.c_str(), "RECREATE");
    Require(!output.IsZombie(), "Cannot create QA output");
    qa.Write();
  }
  TFile input(path.c_str(), "READ");
  Require(!input.IsZombie(), "Cannot reopen QA output");
  auto& ptWeighted = Histogram(input, "selectionPtWeighted");
  auto& ptCounts = Histogram(input, "selectionPtCounts");
  const char* names[] = {"selectionPtWeighted", "selectionPtCounts",
                        "selectionEtaWeighted", "selectionEtaCounts",
                        "selectionPhiWeighted", "selectionPhiCounts"};
  const char* labels[] = {"No jet-quality selection", "trackMax/rawPt only", "Jet ID only",
                         "trackMax/rawPt AND Jet ID", "Jet ID AND event veto",
                         "trackMax/rawPt AND Jet ID AND event veto", "event veto only"};
  for(int index = 0; index < 6; ++index) {
    auto& histogram = Histogram(input, names[index]);
    const int coordinate = index / 2;
    const bool weighted = index % 2 == 0;
    CheckGeometry(histogram, coordinate == 0 ? 500 : (coordinate == 1 ? 80 : 72),
                  coordinate == 0 ? 0. : (coordinate == 1 ? -2. : -pi),
                  coordinate == 0 ? 1000. : (coordinate == 1 ? 2. : pi));
    for(int stage = 0; stage < 7; ++stage) {
      const std::string expected = !idEnabled && stage >= 2 && stage <= 4 ? "INACTIVE: Jet ID disabled" :
        (!idEnabled && stage == 5 ? "trackMax/rawPt AND event veto (Jet ID disabled)" : labels[stage]);
      Require(histogram.GetYaxis()->GetBinLabel(stage + 1) == expected,
              std::string(names[index]) + " saved stage label " + std::to_string(stage));
    }
    ExpectedHistogram expected;
    for(const auto& jet : jets) {
      for(int stage = 0; stage < 7; ++stage) {
        if(selectionMasks[jet.flags][configuration] & (1u << stage))
          expected.Fill(histogram, coordinate == 0 ? jet.pt : (coordinate == 1 ? jet.eta : jet.wrappedPhi),
                        stage, jet.centrality, weighted ? jet.weight : 1.);
      }
    }
    expected.Check(histogram);
    CheckStageTotals(histogram, weighted ? ptWeighted : ptCounts);
  }

  // The retained kinematic histograms still use pT on axis Y and contain only
  // final-selected jets. Their complete bin contents and errors are preserved.
  for(bool eta : {false, true}) {
    auto& histogram = Histogram(input, eta ? "selectedJetEtaWeighted" : "selectedJetPhiWeighted");
    Require(histogram.GetNbinsX() == (eta ? 80 : 72), "Legacy coordinate bins changed");
    Require(histogram.GetNbinsY() == 3, "Legacy pT bins changed");
    const double ptEdges[] = {0., 100., 200., 1000.};
    for(int edge = 0; edge < 4; ++edge)
      Equal(histogram.GetYaxis()->GetBinLowEdge(edge + 1), ptEdges[edge], "Legacy pT edge");
    ExpectedHistogram expected;
    for(const auto& jet : jets)
      if(legacySelected[jet.flags][configuration])
        expected.Fill(histogram, eta ? jet.eta : jet.wrappedPhi, jet.pt, jet.centrality, jet.weight);
    expected.Check(histogram);
  }

  auto& trackFraction = Histogram(input, "trackMaxOverRawPt");
  auto& calibration = Histogram(input, "calibrationPtWeighted");
  ExpectedHistogram trackExpected, calibrationExpected;
  for(const auto& jet : jets) {
    trackExpected.Fill(trackFraction, .2, 100., jet.centrality, jet.weight);
    if(legacySelected[jet.flags][configuration]) {
      calibrationExpected.Fill(calibration, 100., 0., jet.centrality, jet.weight);
      calibrationExpected.Fill(calibration, jet.pt, 2., jet.centrality, jet.weight);
    }
  }
  trackExpected.Check(trackFraction);
  calibrationExpected.Check(calibration);
}
} // namespace

int main() {
  try {
    TemporaryFiles temporary;
    const std::string mapPath = temporary.Path("veto_map");
    {
      TFile file(mapPath.c_str(), "RECREATE");
      TH2D map("mask", "synthetic veto map", 2, -5.2, 5.2, 2, -pi, pi);
      map.SetBinContent(1, 1, 1.);
      map.Write();
    }
    for(int configuration = 0; configuration < 4; ++configuration)
      TestConfiguration(temporary, mapPath, configuration);
    std::cout << "Run3JetQA kinematics regression passed (4 configurations).\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr << "Run3JetQA kinematics regression failed: " << error.what() << '\n';
    return 1;
  }
}

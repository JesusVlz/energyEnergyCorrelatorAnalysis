#include "../../src/Run3JetSelection.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
const double pi = std::acos(-1.0);
const double nan = std::numeric_limits<double>::quiet_NaN();

EECJetSelection::JetInputs GoodJet() {
  EECJetSelection::JetInputs jet{};
  jet.eta = 1;
  jet.phi = 1;
  jet.rawPt = 20;
  jet.storedPt = 25;
  jet.neutralHadronFraction = 0.1;
  jet.neutralEmFraction = 0.1;
  jet.chargedHadronFraction = 0.5;
  jet.chargedEmFraction = 0.1;
  jet.muonFraction = 0.05;
  jet.chargedMultiplicity = 3;
  jet.neutralMultiplicity = 12;
  jet.numberOfConstituents = 15;
  return jet;
}

template<typename Exception = std::exception>
bool Throws(const std::function<void()>& action) {
  try { action(); } catch(const Exception&) { return true; }
  return false;
}

struct TemporaryMaps {
  std::string directory;
  std::vector<std::string> paths;
  TemporaryMaps() {
    char pattern[] = "/tmp/eec-run3-jet-selection-XXXXXX";
    char* path = mkdtemp(pattern);
    if(!path) throw std::runtime_error("Cannot create test directory");
    directory = path;
  }
  ~TemporaryMaps() {
    for(const auto& path : paths) std::remove(path.c_str());
    std::remove(directory.c_str());
  }
  std::string Path(const std::string& name) {
    const std::string path = directory + "/" + name + ".root";
    paths.push_back(path);
    return path;
  }
};

void MakeMap(const std::string& path, const std::string& kind = "valid") {
  TFile file(path.c_str(), "RECREATE");
  if(kind == "one-dimensional") {
    TH1D histogram("mask", "not a veto map", 4, -5.2, 5.2);
    histogram.Write();
    return;
  }
  const double xMax = kind == "short-eta" ? 4.7 : 5.2;
  const double yMax = kind == "wrong-phi" ? 5.2 : pi;
  TH2D histogram("mask", "synthetic map;#eta;#phi", 4, -5.2, xMax, 4, -yMax, yMax);
  histogram.SetBinContent(3, 3, 1);
  histogram.SetBinContent(2, 2, -2);
  histogram.SetBinContent(3, 1, 0.5);
  // The event decision is defined only by regular bins.
  histogram.SetBinContent(0, 0, 100);
  if(kind == "nan") histogram.SetBinContent(1, 1, nan);
  histogram.Write();
}

void TestJetID() {
  using EECJetSelection::PassJetID;
  auto jet = GoodJet();
  assert(PassJetID(jet));
  assert(PassJetID(jet, true));
  for(double eta : {0., 2.6, 2.7, 3., 5., -2.6, -2.7, -3., -5.}) {
    jet = GoodJet();
    jet.eta = eta;
    assert(PassJetID(jet));
  }
  jet.eta = std::nextafter(5., 6.);
  assert(!PassJetID(jet));

  // The eta-band boundaries matter: each restrictive test changes at the
  // published inclusive boundary and must remain symmetric in detector eta.
  for(double sign : {-1., 1.}) {
    jet = GoodJet();
    jet.neutralHadronFraction = 0.95;
    jet.eta = sign*2.6;
    assert(PassJetID(jet));
    jet.eta = sign*std::nextafter(2.6, 3.);
    assert(!PassJetID(jet));
    jet.eta = sign*2.7;
    assert(!PassJetID(jet));
    jet.eta = sign*std::nextafter(2.7, 3.);
    assert(PassJetID(jet));

    jet = GoodJet();
    jet.neutralEmFraction = 0.5;
    jet.eta = sign*3.;
    assert(PassJetID(jet));
    jet.eta = sign*std::nextafter(3., 4.);
    assert(!PassJetID(jet));
  }

  jet = GoodJet(); jet.neutralHadronFraction = 0.99; assert(!PassJetID(jet));
  jet = GoodJet(); jet.neutralEmFraction = 0.90; assert(!PassJetID(jet));
  jet = GoodJet(); jet.chargedHadronFraction = 0.01; assert(!PassJetID(jet));
  jet = GoodJet(); jet.chargedMultiplicity = 0; assert(!PassJetID(jet));
  jet = GoodJet(); jet.numberOfConstituents = 1; assert(!PassJetID(jet));
  jet = GoodJet(); jet.muonFraction = 0.8;
  assert(PassJetID(jet)); assert(!PassJetID(jet, true));
  jet = GoodJet(); jet.chargedEmFraction = 0.8;
  assert(PassJetID(jet)); assert(!PassJetID(jet, true));

  jet = GoodJet(); jet.eta = 2.65;
  jet.chargedHadronFraction = 0;
  jet.numberOfConstituents = 0;
  assert(PassJetID(jet));
  jet.neutralHadronFraction = 0.9; assert(!PassJetID(jet));
  jet = GoodJet(); jet.eta = 2.65; jet.neutralEmFraction = 0.99; assert(!PassJetID(jet));
  jet = GoodJet(); jet.eta = 2.65; jet.chargedMultiplicity = 0; assert(!PassJetID(jet));
  jet = GoodJet(); jet.eta = 2.65; jet.muonFraction = 0.8;
  assert(PassJetID(jet)); assert(!PassJetID(jet, true));
  jet = GoodJet(); jet.eta = 2.65; jet.chargedEmFraction = 0.8;
  assert(PassJetID(jet)); assert(!PassJetID(jet, true));

  jet = GoodJet(); jet.eta = 2.8; jet.neutralMultiplicity = 1; assert(!PassJetID(jet));
  jet = GoodJet(); jet.eta = 2.8; jet.neutralHadronFraction = 0.99; assert(!PassJetID(jet));
  jet = GoodJet(); jet.eta = 2.8; jet.neutralEmFraction = 0.99; assert(!PassJetID(jet));
  jet = GoodJet(); jet.eta = 2.8; jet.chargedMultiplicity = 0;
  jet.chargedHadronFraction = 0; jet.muonFraction = 1; jet.chargedEmFraction = 1;
  assert(PassJetID(jet, true));
  jet = GoodJet(); jet.eta = 4; jet.neutralMultiplicity = 10; assert(!PassJetID(jet));
  jet = GoodJet(); jet.eta = 4; jet.neutralEmFraction = 0.4; assert(!PassJetID(jet));
  jet = GoodJet(); jet.eta = 4; jet.neutralHadronFraction = 1; assert(PassJetID(jet));

  jet = GoodJet(); jet.eta = nan; assert(Throws<std::invalid_argument>([&] { PassJetID(jet); }));
  jet = GoodJet(); jet.neutralHadronFraction = nan;
  assert(Throws<std::invalid_argument>([&] { PassJetID(jet); }));
  jet = GoodJet(); jet.chargedEmFraction = -0.01;
  assert(Throws<std::invalid_argument>([&] { PassJetID(jet); }));
  jet = GoodJet(); jet.muonFraction = 1.1;
  assert(Throws<std::invalid_argument>([&] { PassJetID(jet); }));
  jet = GoodJet(); jet.chargedMultiplicity = -1;
  assert(Throws<std::invalid_argument>([&] { PassJetID(jet); }));

  // The five forest PF categories omit HF counts. A proven lower bound
  // passing a multiplicity cut is sufficient; an unresolved cut must throw.
  jet = GoodJet(); jet.multiplicitiesAreLowerBounds = true;
  assert(PassJetID(jet));
  jet.numberOfConstituents = 1;
  assert(Throws<std::runtime_error>([&] { PassJetID(jet); }));
  jet.neutralHadronFraction = 0.99;
  assert(!PassJetID(jet));  // Other definite failures resolve the decision.
  jet = GoodJet(); jet.multiplicitiesAreLowerBounds = true;
  jet.eta = 2.65; jet.numberOfConstituents = 0; jet.neutralMultiplicity = 0;
  assert(PassJetID(jet));  // Neither unknown total is used in this band.
  jet.eta = 2.8; jet.neutralMultiplicity = 1;
  assert(Throws<std::runtime_error>([&] { PassJetID(jet); }));
  jet.neutralMultiplicity = 2; assert(PassJetID(jet));
  jet.eta = 4; jet.neutralMultiplicity = 10;
  assert(Throws<std::runtime_error>([&] { PassJetID(jet); }));
  jet.neutralMultiplicity = 11; assert(PassJetID(jet));
}

void TestHINPrescriptions() {
  using EECJetSelection::JetIDPrescription;
  using EECJetSelection::PassJetID;
  const auto pbpb = JetIDPrescription::HINPbPb;
  const auto pp = JetIDPrescription::HINpp;
  auto baseline = GoodJet();
  baseline.neutralHadronFraction = 0.3;
  auto jet = baseline;
  for(double eta : {0., 2.6, 2.7, 3., 5., -2.6, -2.7, -3., -5.}) {
    jet = baseline; jet.eta = eta;
    assert(PassJetID(jet, pbpb));
    assert(PassJetID(jet, pp));
  }
  for(const auto prescription : {pbpb, pp}) {
    jet = baseline; jet.eta = std::nextafter(5., 6.);
    assert(!PassJetID(jet, prescription));
    jet.eta = -std::nextafter(5., 6.);
    assert(!PassJetID(jet, prescription));
  }

  // PbPb deliberately omits hadronic and multiplicity terms.
  jet = baseline;
  jet.neutralHadronFraction = 1; jet.chargedHadronFraction = 0;
  jet.chargedMultiplicity = 0; jet.neutralMultiplicity = 0;
  jet.numberOfConstituents = 0; jet.multiplicitiesAreLowerBounds = true;
  assert(PassJetID(jet, pbpb));
  assert(!PassJetID(jet, pp));
  for(double sign : {-1., 1.}) {
    jet = baseline; jet.neutralEmFraction = 0.95; jet.eta = sign*2.6;
    assert(!PassJetID(jet, pbpb));
    assert(!PassJetID(jet, pp));
    jet.eta = sign*std::nextafter(2.6, 3.);
    assert(PassJetID(jet, pbpb));
    assert(PassJetID(jet, pp));

    jet = baseline; jet.muonFraction = 0.8; jet.eta = sign*2.7;
    assert(!PassJetID(jet, pbpb));
    assert(!PassJetID(jet, pp));
    jet.eta = sign*std::nextafter(2.7, 3.);
    assert(PassJetID(jet, pbpb));
    assert(PassJetID(jet, pp));

    jet = baseline; jet.neutralEmFraction = 0.95; jet.eta = sign*3.;
    assert(PassJetID(jet, pbpb));
    assert(PassJetID(jet, pp));
    jet.eta = sign*std::nextafter(3., 4.);
    assert(!PassJetID(jet, pbpb));
    assert(!PassJetID(jet, pp));
  }
  jet = baseline; jet.neutralEmFraction = 0.9; assert(!PassJetID(jet, pbpb));
  jet = baseline; jet.muonFraction = 0.8; assert(!PassJetID(jet, pbpb));
  jet = baseline; jet.chargedEmFraction = 0.8; assert(!PassJetID(jet, pbpb));
  jet = baseline; jet.eta = 2.65; jet.neutralEmFraction = 0.99; assert(!PassJetID(jet, pbpb));
  jet = baseline; jet.eta = 2.65; jet.chargedEmFraction = 0.8; assert(!PassJetID(jet, pbpb));
  jet = baseline; jet.eta = 2.8; jet.neutralEmFraction = 0.99; assert(!PassJetID(jet, pbpb));
  jet = baseline; jet.eta = 4; jet.neutralEmFraction = 0.9; assert(!PassJetID(jet, pbpb));
  jet = baseline; jet.eta = 4; jet.muonFraction = 1; jet.chargedEmFraction = 1;
  assert(PassJetID(jet, pbpb));

  // HINpp differs from generic CHS at central NHF, transition NHF, and
  // forward NHF/NEF. Its multiplicities are the helper's category sums.
  jet = baseline; jet.neutralHadronFraction = 0.9; assert(!PassJetID(jet, pp));
  assert(PassJetID(jet, JetIDPrescription::Tight));
  jet = baseline; jet.neutralEmFraction = 0.9; assert(!PassJetID(jet, pp));
  jet = baseline; jet.chargedHadronFraction = 0.01; assert(!PassJetID(jet, pp));
  jet = baseline; jet.chargedMultiplicity = 0; assert(!PassJetID(jet, pp));
  jet = baseline; jet.numberOfConstituents = 1; jet.multiplicitiesAreLowerBounds = true;
  assert(!PassJetID(jet, pp));
  assert(Throws<std::runtime_error>([&] { PassJetID(jet, JetIDPrescription::Tight); }));
  jet = baseline; jet.chargedEmFraction = 0.8; assert(!PassJetID(jet, pp));
  jet = baseline; jet.muonFraction = 0.8; assert(!PassJetID(jet, pp));
  jet = baseline; jet.eta = 2.65; jet.neutralHadronFraction = 0.9; assert(!PassJetID(jet, pp));
  jet = baseline; jet.eta = 2.65; jet.neutralEmFraction = 0.99; assert(!PassJetID(jet, pp));
  jet = baseline; jet.eta = 2.65; jet.chargedMultiplicity = 0; assert(!PassJetID(jet, pp));
  jet = baseline; jet.eta = 2.65; jet.chargedEmFraction = 0.8; assert(!PassJetID(jet, pp));
  jet = baseline; jet.eta = 2.8; jet.neutralHadronFraction = 1;
  assert(PassJetID(jet, pp));
  assert(!PassJetID(jet, JetIDPrescription::Tight));
  jet.neutralMultiplicity = 1; jet.multiplicitiesAreLowerBounds = true;
  assert(!PassJetID(jet, pp));
  jet.neutralMultiplicity = 2; assert(PassJetID(jet, pp));
  jet.neutralEmFraction = 0.99; assert(!PassJetID(jet, pp));
  jet = baseline; jet.eta = 4; jet.neutralMultiplicity = 10; jet.multiplicitiesAreLowerBounds = true;
  assert(!PassJetID(jet, pp));
  jet.neutralMultiplicity = 11; assert(PassJetID(jet, pp));
  jet.neutralHadronFraction = 0.2; assert(!PassJetID(jet, pp));
  jet.neutralHadronFraction = std::nextafter(0.2, 1.); assert(PassJetID(jet, pp));
  jet.neutralEmFraction = 0.9; assert(!PassJetID(jet, pp));
  jet = baseline; jet.eta = 4; jet.neutralMultiplicity = -1;
  assert(Throws<std::invalid_argument>([&] { PassJetID(jet, pp); }));

  // Actual retained forests contain CHF=1.0000001192092896. Permit
  // Float_t endpoint rounding while comparing the unmodified stored fraction.
  const double fractionEpsilon = std::numeric_limits<float>::epsilon();
  jet = baseline; jet.chargedHadronFraction = 1+fractionEpsilon;
  assert(PassJetID(jet, pp));
  assert(PassJetID(jet, pbpb));
  jet.chargedHadronFraction = 1+1e-4;
  assert(Throws<std::invalid_argument>([&] { PassJetID(jet, pp); }));
  jet.chargedHadronFraction = nan;
  assert(Throws<std::invalid_argument>([&] { PassJetID(jet, pp); }));
  jet = baseline; jet.muonFraction = -fractionEpsilon;
  assert(PassJetID(jet, pp));
  jet.muonFraction = -1e-4;
  assert(Throws<std::invalid_argument>([&] { PassJetID(jet, pp); }));
  // The domain tolerance does not add a hadronic upper cut to forward
  // HINpp or relax its exact neutral EM threshold.
  jet = baseline; jet.eta = 4; jet.neutralHadronFraction = 1+fractionEpsilon;
  assert(PassJetID(jet, pp));
  jet.neutralEmFraction = 0.9;
  assert(!PassJetID(jet, pp));
  jet = baseline; jet.neutralEmFraction = 1+fractionEpsilon;
  assert(!PassJetID(jet, pbpb));

  jet = baseline;
  assert(Throws<std::invalid_argument>([&] { PassJetID(jet, static_cast<JetIDPrescription>(99)); }));
  assert(Throws<std::invalid_argument>([&] { EECJetSelection::IsEligibleVetoJet(jet, 0, static_cast<JetIDPrescription>(99)); }));
  assert(PassJetID(jet, false) == PassJetID(jet, JetIDPrescription::Tight));
  jet.muonFraction = 0.8;
  assert(PassJetID(jet, true) == PassJetID(jet, JetIDPrescription::TightLeptonVeto));
}

void TestVeto() {
  using EECJetSelection::JetVetoMap;
  using EECJetSelection::RejectEventForJet;
  TemporaryMaps temporary;
  const auto validPath = temporary.Path("valid");
  MakeMap(validPath);
  JetVetoMap map(validPath, "mask");
  assert(map.SourcePath() == validPath && map.HistogramKey() == "mask");
  assert(map.IsVetoed(1, 1));
  assert(map.IsVetoed(-1, -1));  // Negative mask content also vetoes.
  assert(!map.IsVetoed(1, -1));
  assert(map.IsVetoed(1, pi));
  assert(map.IsVetoed(1, -pi));
  const double floatPi = static_cast<double>(static_cast<float>(pi));
  assert(map.IsVetoed(1, floatPi) == map.IsVetoed(1, std::remainder(floatPi, 2*pi)));
  assert(map.IsVetoed(1, -floatPi) == map.IsVetoed(1, std::remainder(-floatPi, 2*pi)));
  assert(!map.IsVetoed(-5, 1));
  assert(!map.IsVetoed(5, 1));
  assert(Throws<std::invalid_argument>([&] { map.IsVetoed(5.01, 1); }));
  assert(Throws<std::invalid_argument>([&] { map.IsVetoed(1, pi+0.01); }));
  assert(Throws<std::invalid_argument>([&] { map.IsVetoed(1, nan); }));
  assert(Throws<std::invalid_argument>([&] { map.IsVetoed(nan, 1); }));

  auto jet = GoodJet();
  for(const auto prescription : {EECJetSelection::JetIDPrescription::HINPbPb,
                                  EECJetSelection::JetIDPrescription::HINpp}) {
    assert(RejectEventForJet(jet, std::nextafter(15., 16.), map, prescription));
    assert(!RejectEventForJet(jet, 15, map, prescription));
    jet.neutralEmFraction = 0.8; jet.chargedEmFraction = 0.1;
    assert(!RejectEventForJet(jet, 20, map, prescription));
    jet = GoodJet(); jet.chargedEmFraction = 0.8;
    assert(!RejectEventForJet(jet, 20, map, prescription));
    jet = GoodJet();
  }
  jet.chargedHadronFraction = 0; jet.chargedMultiplicity = 0;
  jet.numberOfConstituents = 0; jet.multiplicitiesAreLowerBounds = true;
  assert(RejectEventForJet(jet, 20, map, EECJetSelection::JetIDPrescription::HINPbPb));
  assert(!RejectEventForJet(jet, 20, map, EECJetSelection::JetIDPrescription::HINpp));
  jet = GoodJet();
  assert(RejectEventForJet(jet, std::nextafter(15., 16.), map));
  jet.phi = floatPi;
  assert(RejectEventForJet(jet, 20, map));
  jet.phi = -floatPi;
  assert(!RejectEventForJet(jet, 20, map));
  jet.phi = 1;
  assert(!RejectEventForJet(jet, 15., map));
  assert(!RejectEventForJet(jet, 14.9, map));
  jet.neutralEmFraction = 0.8;
  jet.chargedEmFraction = 0.1;
  assert(!RejectEventForJet(jet, 20, map));  // Sum == 0.9 is ineligible.
  // Use a clearly separated floating value for the below-threshold case.
  jet.chargedEmFraction = 0.099;
  assert(RejectEventForJet(jet, 20, map));
  jet = GoodJet(); jet.chargedMultiplicity = 0;
  assert(!RejectEventForJet(jet, 20, map));
  jet = GoodJet(); jet.phi = -1;
  assert(!RejectEventForJet(jet, 20, map));
  jet = GoodJet(); jet.multiplicitiesAreLowerBounds = true;
  jet.numberOfConstituents = 1;
  assert(!RejectEventForJet(jet, 15, map));
  jet.neutralEmFraction = 0.8; jet.chargedEmFraction = 0.1;
  assert(!RejectEventForJet(jet, 20, map));
  jet.neutralEmFraction = 0.1;
  assert(Throws<std::runtime_error>([&] { RejectEventForJet(jet, 20, map); }));

  // A sub-analysis-threshold jet can veto the complete event. The caller
  // must inspect every reco jet, rather than only its analysis-selected jet.
  std::vector<EECJetSelection::JetInputs> event(2, GoodJet());
  event[0].phi = -1;
  event[1].storedPt = 20;
  bool rejected = false;
  for(const auto& candidate : event) {
    rejected = rejected || RejectEventForJet(candidate, candidate.storedPt, map);
  }
  assert(rejected);

  jet = GoodJet();
  assert(Throws<std::invalid_argument>([&] { RejectEventForJet(jet, nan, map); }));
  assert(Throws<std::invalid_argument>([&] { RejectEventForJet(jet, -1, map); }));
  jet.neutralEmFraction = nan;
  assert(Throws<std::invalid_argument>([&] { RejectEventForJet(jet, 0, map); }));
  jet = GoodJet(); jet.phi = 4;
  assert(Throws<std::invalid_argument>([&] { RejectEventForJet(jet, 0, map); }));

  assert(Throws<std::invalid_argument>([&] { JetVetoMap missing("", "mask"); }));
  assert(Throws<std::invalid_argument>([&] { JetVetoMap missing(validPath, ""); }));
  assert(Throws<std::runtime_error>([&] { JetVetoMap missing(validPath, "absent"); }));
  assert(Throws<std::runtime_error>([&] { JetVetoMap missing(temporary.directory + "/absent.root", "mask"); }));
  for(const std::string kind : {"one-dimensional", "short-eta", "wrong-phi", "nan"}) {
    const auto path = temporary.Path(kind);
    MakeMap(path, kind);
    assert(Throws<std::runtime_error>([&] { JetVetoMap invalid(path, "mask"); }));
  }
  JetVetoMap moved(std::move(map));
  assert(moved.IsVetoed(1, 1));
  assert(Throws<std::logic_error>([&] { map.IsVetoed(1, 1); }));
  // Removal of the source must not change the owned, detached map.
  assert(std::remove(validPath.c_str()) == 0);
  assert(moved.IsVetoed(1, 1));
}
}  // namespace

int main() {
  TestJetID();
  TestHINPrescriptions();
  TestVeto();
  std::cout << "PASS: supplied tight/lepton-veto and pinned HIN PbPb/pp Jet ID, eta and threshold boundaries, "
               "event-veto selection, explicit detached ROOT map, negative masks, "
               "invalid payload and coordinate rejection\n";
}

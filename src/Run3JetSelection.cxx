#include "Run3JetSelection.h"

#include <TAxis.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace EECJetSelection {
namespace {
const double kPi = std::acos(-1.0);
const double kAxisTolerance = 1e-6;
// Native PF energy fractions are Float_t ratios. Values one rounding unit
// beyond physical endpoints occur in retained forests. Accept only this small
// representational neighborhood; never clamp values used by the ID cuts.
const double kFractionEndpointTolerance = 2*std::numeric_limits<float>::epsilon();
// Native forest coordinates are Float_t. An angle rounding across +/-pi by
// one Float_t rounding unit remains a valid periodic detector coordinate.
const double kPhiEndpointTolerance = std::numeric_limits<float>::epsilon()*kPi;

double CanonicalDetectorPhi(double phi) {
  if(!std::isfinite(phi) || std::abs(phi) > kPi+kPhiEndpointTolerance) {
    throw std::invalid_argument("Run-3 veto map: expected finite detector phi in [-pi,pi] radians (allowing Float_t endpoint rounding)");
  }
  if(phi == kPi) return -kPi;
  return std::remainder(phi, 2*kPi);
}

void ValidatePrescription(JetIDPrescription prescription) {
  switch(prescription) {
    case JetIDPrescription::Tight:
    case JetIDPrescription::TightLeptonVeto:
    case JetIDPrescription::HINPbPb:
    case JetIDPrescription::HINpp: return;
  }
  throw std::invalid_argument("Run-3 Jet ID: unknown explicit prescription");
}

void RequireFraction(double value, const char* name) {
  if(!std::isfinite(value) || value < -kFractionEndpointTolerance || value > 1+kFractionEndpointTolerance) {
    throw std::invalid_argument(std::string("Run-3 Jet ID: invalid ") + name +
      "; expected a finite energy fraction in [0,1], allowing 2 Float_t epsilon endpoint rounding");
  }
}

void ValidateJetIDInputs(const JetInputs& jet) {
  if(!std::isfinite(jet.eta)) {
    throw std::invalid_argument("Run-3 Jet ID: detector eta is not finite");
  }
  RequireFraction(jet.neutralHadronFraction, "neutral hadron fraction");
  RequireFraction(jet.neutralEmFraction, "neutral EM fraction");
  RequireFraction(jet.chargedHadronFraction, "charged hadron fraction");
  RequireFraction(jet.chargedEmFraction, "charged EM fraction");
  RequireFraction(jet.muonFraction, "muon fraction");
}

bool PassMultiplicity(int count, int threshold, bool lowerBound, const char* name) {
  if(count < 0) {
    throw std::invalid_argument(std::string("Run-3 Jet ID: negative ") + name);
  }
  if(count > threshold) return true;
  if(lowerBound) {
    throw std::runtime_error(std::string("Run-3 Jet ID: the forest's ") + name +
      " is only a lower bound and cannot resolve the ID threshold; provide the total multiplicity including HF constituents");
  }
  return false;
}

void ValidateAxis(const TAxis& axis, const char* name) {
  if(axis.GetNbins() < 1) {
    throw std::runtime_error(std::string("Run-3 veto map: empty ") + name + " axis");
  }
  double previous = axis.GetBinLowEdge(1);
  if(!std::isfinite(previous)) {
    throw std::runtime_error(std::string("Run-3 veto map: nonfinite ") + name + " axis");
  }
  for(int bin = 1; bin <= axis.GetNbins(); ++bin) {
    const double edge = axis.GetBinUpEdge(bin);
    if(!std::isfinite(edge) || !(edge > previous)) {
      throw std::runtime_error(std::string("Run-3 veto map: invalid ") + name + " bin edges");
    }
    previous = edge;
  }
}

void ValidateMap(const TH2& map) {
  ValidateAxis(*map.GetXaxis(), "eta (X)");
  ValidateAxis(*map.GetYaxis(), "phi (Y)");
  if(map.GetXaxis()->GetXmin() > -5 || map.GetXaxis()->GetXmax() <= 5) {
    throw std::runtime_error("Run-3 veto map: X axis must cover detector eta in [-5,5], including eta=5 in a regular bin");
  }
  if(std::abs(map.GetYaxis()->GetXmin() + kPi) > kAxisTolerance ||
     std::abs(map.GetYaxis()->GetXmax() - kPi) > kAxisTolerance ||
     map.GetYaxis()->GetXmin() > -kPi || map.GetYaxis()->GetXmax() < kPi) {
    throw std::runtime_error("Run-3 veto map: Y axis must span one detector phi period [-pi,pi] in radians");
  }
  for(int etaBin = 0; etaBin <= map.GetNbinsX()+1; ++etaBin) {
    for(int phiBin = 0; phiBin <= map.GetNbinsY()+1; ++phiBin) {
      if(!std::isfinite(map.GetBinContent(etaBin, phiBin))) {
        throw std::runtime_error("Run-3 veto map: nonfinite histogram content");
      }
    }
  }
}
}  // namespace

double DetectorPhi(double phi) { return CanonicalDetectorPhi(phi); }

bool PassJetID(const JetInputs& jet, bool leptonVeto) {
  return PassJetID(jet, leptonVeto ? JetIDPrescription::TightLeptonVeto : JetIDPrescription::Tight);
}

bool PassJetID(const JetInputs& jet, JetIDPrescription prescription) {
  ValidatePrescription(prescription);
  ValidateJetIDInputs(jet);
  const double eta = std::abs(jet.eta);
  if(eta > 5) return false;
  if(prescription == JetIDPrescription::HINPbPb) {
    if(eta <= 2.6) return jet.neutralEmFraction < 0.9 && jet.muonFraction < 0.8 && jet.chargedEmFraction < 0.8;
    if(eta <= 2.7) return jet.neutralEmFraction < 0.99 && jet.muonFraction < 0.8 && jet.chargedEmFraction < 0.8;
    if(eta <= 3.0) return jet.neutralEmFraction < 0.99;
    return jet.neutralEmFraction < 0.9;
  }
  if(prescription == JetIDPrescription::HINpp) {
    // The pinned helper explicitly defines m=CHM+NHM+CEM+NEM+MUM,
    // cm=CHM+CEM+MUM and nm=NHM+NEM. These category sums are the
    // prescription itself; do not substitute unprovided full-HF counts.
    if(eta <= 2.6) {
      return jet.neutralHadronFraction < 0.9 && jet.neutralEmFraction < 0.9 &&
        jet.muonFraction < 0.8 && jet.chargedHadronFraction > 0.01 && jet.chargedEmFraction < 0.8 &&
        PassMultiplicity(jet.chargedMultiplicity, 0, false, "charged category count") &&
        PassMultiplicity(jet.numberOfConstituents, 1, false, "total category count");
    }
    if(eta <= 2.7) {
      return jet.neutralHadronFraction < 0.9 && jet.neutralEmFraction < 0.99 &&
        jet.muonFraction < 0.8 && jet.chargedEmFraction < 0.8 &&
        PassMultiplicity(jet.chargedMultiplicity, 0, false, "charged category count");
    }
    if(eta <= 3.0) {
      return jet.neutralEmFraction < 0.99 &&
        PassMultiplicity(jet.neutralMultiplicity, 1, false, "neutral category count");
    }
    return jet.neutralHadronFraction > 0.2 && jet.neutralEmFraction < 0.9 &&
      PassMultiplicity(jet.neutralMultiplicity, 10, false, "neutral category count");
  }
  const bool leptonVeto = prescription == JetIDPrescription::TightLeptonVeto;
  if(eta <= 2.6) {
    return jet.neutralHadronFraction < 0.99 && jet.neutralEmFraction < 0.90 &&
      jet.chargedHadronFraction > 0.01 &&
      (!leptonVeto || (jet.muonFraction < 0.80 && jet.chargedEmFraction < 0.80)) &&
      PassMultiplicity(jet.chargedMultiplicity, 0, false, "charged multiplicity") &&
      PassMultiplicity(jet.numberOfConstituents, 1, jet.multiplicitiesAreLowerBounds, "constituent count");
  }
  if(eta <= 2.7) {
    return jet.neutralHadronFraction < 0.90 && jet.neutralEmFraction < 0.99 &&
      (!leptonVeto || (jet.muonFraction < 0.80 && jet.chargedEmFraction < 0.80)) &&
      PassMultiplicity(jet.chargedMultiplicity, 0, false, "charged multiplicity");
  }
  if(eta <= 3.0) {
    return jet.neutralHadronFraction < 0.99 && jet.neutralEmFraction < 0.99 &&
      PassMultiplicity(jet.neutralMultiplicity, 1, jet.multiplicitiesAreLowerBounds, "neutral multiplicity");
  }
  return jet.neutralEmFraction < 0.4 &&
    PassMultiplicity(jet.neutralMultiplicity, 10, jet.multiplicitiesAreLowerBounds, "neutral multiplicity");
}

JetVetoMap::JetVetoMap(const std::string& path, const std::string& histogramKey) :
  fSourcePath(path), fHistogramKey(histogramKey) {
  if(path.empty() || histogramKey.empty()) {
    throw std::invalid_argument("Run-3 veto map: explicit nonempty file path and histogram key are required");
  }
  // Keep file and clone registration local; the detached clone survives file close.
  TDirectory::TContext directoryContext(nullptr);
  std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
  if(!file || file->IsZombie()) {
    throw std::runtime_error("Run-3 veto map: cannot open " + path);
  }
  TH2* source = dynamic_cast<TH2*>(file->Get(histogramKey.c_str()));
  if(!source) {
    throw std::runtime_error("Run-3 veto map: missing two-dimensional histogram " + histogramKey + " in " + path);
  }
  ValidateMap(*source);
  TH2* detached = dynamic_cast<TH2*>(source->Clone());
  if(!detached) throw std::runtime_error("Run-3 veto map: histogram clone failed");
  detached->SetDirectory(nullptr);
  fHistogram.reset(detached);
}

JetVetoMap::~JetVetoMap() = default;
JetVetoMap::JetVetoMap(JetVetoMap&&) noexcept = default;
JetVetoMap& JetVetoMap::operator=(JetVetoMap&&) noexcept = default;

bool JetVetoMap::IsVetoed(double eta, double phi) const {
  if(!fHistogram) throw std::logic_error("Run-3 veto map: moved-from map cannot be queried");
  if(!std::isfinite(eta) || std::abs(eta) > 5) {
    throw std::invalid_argument("Run-3 veto map: expected finite detector eta in [-5,5]");
  }
  phi = DetectorPhi(phi);
  const int etaBin = fHistogram->GetXaxis()->FindFixBin(eta);
  const int phiBin = fHistogram->GetYaxis()->FindFixBin(phi);
  if(etaBin < 1 || etaBin > fHistogram->GetNbinsX() ||
     phiBin < 1 || phiBin > fHistogram->GetNbinsY()) {
    throw std::out_of_range("Run-3 veto map: detector coordinates fall outside regular map bins; no overflow fallback");
  }
  const double value = fHistogram->GetBinContent(etaBin, phiBin);
  if(!std::isfinite(value)) throw std::runtime_error("Run-3 veto map: queried nonfinite content");
  return value != 0;
}

bool IsEligibleVetoJet(const JetInputs& jet, double selectedPt, bool leptonVeto) {
  return IsEligibleVetoJet(jet, selectedPt,
      leptonVeto ? JetIDPrescription::TightLeptonVeto : JetIDPrescription::Tight);
}

bool IsEligibleVetoJet(const JetInputs& jet, double selectedPt, JetIDPrescription prescription) {
  ValidatePrescription(prescription);
  if(!std::isfinite(selectedPt) || selectedPt < 0) {
    throw std::invalid_argument("Run-3 event veto: selected pT must be finite and nonnegative");
  }
  DetectorPhi(jet.phi);
  ValidateJetIDInputs(jet);
  if(selectedPt <= 15 || jet.chargedEmFraction + jet.neutralEmFraction >= 0.9) return false;
  return PassJetID(jet, prescription);
}

bool RejectEventForJet(const JetInputs& jet, double selectedPt,
                       const JetVetoMap& map, bool leptonVeto) {
  return RejectEventForJet(jet, selectedPt, map,
      leptonVeto ? JetIDPrescription::TightLeptonVeto : JetIDPrescription::Tight);
}

bool RejectEventForJet(const JetInputs& jet, double selectedPt,
                       const JetVetoMap& map, JetIDPrescription prescription) {
  if(!IsEligibleVetoJet(jet, selectedPt, prescription)) return false;
  return map.IsVetoed(jet.eta, jet.phi);
}

}  // namespace EECJetSelection

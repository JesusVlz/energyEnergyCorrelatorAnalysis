#ifndef EEC_RUN3_JET_SELECTION_H
#define EEC_RUN3_JET_SELECTION_H

#include "Run3JetSelectionInputs.h"

#include <memory>
#include <string>

class TH2;

// Run-3 PF Jet ID and the event-veto predicate. All coordinates in this API are
// the detector jet axis (jteta/jtphi), independently of the EEC analysis axis.
namespace EECJetSelection {

// Canonical detector phi used both for map queries and QA. Exact +pi maps to
// -pi; Float_t rounding across endpoints wraps periodically. Invalid or clearly
// out-of-domain input throws, rather than mapping arbitrary angles silently.
double DetectorPhi(double phi);

// Explicit alternatives remain separate: the HIN helpers have different
// collision-specific criteria from the originally supplied generic CHS table.
enum class JetIDPrescription { Tight = 1, TightLeptonVeto = 2, HINPbPb = 3, HINpp = 4 };

// HIN Run-3 tight AK4/AK8 CHS criteria. leptonVeto additionally applies the
// charged EM and muon fraction requirements in |eta| <= 2.7. Nonfinite or
// unphysical inputs throw; fraction-domain validation allows 2 Float_t epsilon
// endpoint rounding without changing the raw values used by any ID threshold.
// Finite jets outside |eta| <= 5 fail the ID. Native
// five-category counts may omit HF: lower bounds can establish a passing
// multiplicity cut, but an unresolved threshold throws instead of guessing.
bool PassJetID(const JetInputs& jet, bool leptonVeto = false);
// HIN prescriptions reproduce the ID sections of HIN-Jets include/JetSelection_
// {PbPb,pp}.h at c64293a132af388af1827044779af36e57f884d2. HINpp defines its
// neutral and total counts as the five-category sums, independently of whether
// those sums include all physical HF constituents. No map default is adopted.
bool PassJetID(const JetInputs& jet, JetIDPrescription prescription);

// An explicit map and histogram key are always required. The histogram must
// have detector eta on X and phi in radians on Y, and cover |eta| <= 5 and a
// full [-pi,pi] period. This object owns a detached histogram after loading.
class JetVetoMap {
public:
  JetVetoMap(const std::string& path, const std::string& histogramKey);
  ~JetVetoMap();
  JetVetoMap(JetVetoMap&&) noexcept;
  JetVetoMap& operator=(JetVetoMap&&) noexcept;
  JetVetoMap(const JetVetoMap&) = delete;
  JetVetoMap& operator=(const JetVetoMap&) = delete;

  // Nonzero values, including negative ones, veto. Phi +pi is the same
  // periodic coordinate as -pi. Float_t rounding across either endpoint is
  // wrapped periodically; other out-of-domain coordinates throw.
  // No underflow, overflow, edge clamping or interpolation is used.
  bool IsVetoed(double eta, double phi) const;
  const std::string& SourcePath() const { return fSourcePath; }
  const std::string& HistogramKey() const { return fHistogramKey; }

private:
  std::string fSourcePath;
  std::string fHistogramKey;
  std::unique_ptr<TH2> fHistogram;
};

// The minimum selection for the veto, with the same validation as the full
// predicate. Jets failing pT or the EM-sum cut do not require multiplicity
// decisions. Useful for recording QA for every eligible detector jet.
bool IsEligibleVetoJet(const JetInputs& jet, double selectedPt, bool leptonVeto = false);
bool IsEligibleVetoJet(const JetInputs& jet, double selectedPt, JetIDPrescription prescription);

// Call for EVERY reconstructed jet in an event, before EEC analysis jet cuts.
// selectedPt is an explicitly chosen caller-side convention; this helper does
// not infer a correction chain, collection, era, or data/MC applicability.
// Reject the whole event when any call returns true. Invalid inputs throw.
bool RejectEventForJet(const JetInputs& jet, double selectedPt,
                       const JetVetoMap& map, bool leptonVeto = false);
bool RejectEventForJet(const JetInputs& jet, double selectedPt,
                       const JetVetoMap& map, JetIDPrescription prescription);

}  // namespace EECJetSelection

#endif

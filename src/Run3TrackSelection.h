#ifndef EEC_RUN3_TRACK_SELECTION_H
#define EEC_RUN3_TRACK_SELECTION_H

#include <cmath>
#include <stdexcept>
#include <string>

// Run-3 OO2025/ppRef2024 prescription from the supplied tracking selections.
// Deliberately separate from Run-2: do not change PbPb2018/pp2017 cuts.
namespace EECRun3Tracking {
inline const char* Version() { return "run3-tracking-2026-11-v1"; }
inline std::string CanonicalSelection(const std::string& value) {
  if(value.empty() || value == "Nominal") return "Nominal";
  if(value == "Loose" || value == "Tight") return value;
  throw std::invalid_argument("Run3TrackSelection must be Nominal, Loose or Tight");
}
inline double VertexLimit(const std::string& value) {
  const std::string selection = CanonicalSelection(value);
  return selection == "Loose" ? 5.0 : selection == "Tight" ? 2.0 : 3.0;
}
inline std::string InferSelection(double relativeError, double vertexLimit) {
  if(!std::isfinite(relativeError) || std::abs(relativeError - 0.1) > 1e-6)
    throw std::invalid_argument("Run-3 tracks require MaxTrackPtRelativeError 0.1 for all selections; Run-2 0.05/0.15 variations are not applicable");
  for(const char* name : {"Nominal", "Loose", "Tight"})
    if(std::isfinite(vertexLimit) && std::abs(vertexLimit - VertexLimit(name)) < 1e-6) return name;
  throw std::invalid_argument("Run-3 VertexMaxDistance must be 3 (Nominal), 5 (Loose) or 2 (Tight)");
}
inline std::string ResolveSelection(const std::string& configured, double relativeError, double vertexLimit) {
  const std::string inferred = InferSelection(relativeError, vertexLimit);
  if(!configured.empty() && CanonicalSelection(configured) != inferred)
    throw std::invalid_argument("Run3TrackSelection and VertexMaxDistance disagree; change the selection and its matching vertex limit together");
  return inferred;
}
inline std::string Payload(bool oxygen, const std::string& value) {
  const std::string selection = CanonicalSelection(value);
  if(oxygen)
    return "trackCorrectionTables/OO2025/Eff_OO_2025_Hijing_MB_Centrality_fromHihfpf_NoPU_3D_" + selection + "_Official_18Nov2025.root";
  return "trackCorrectionTables/pp2024/Eff_ppref_2024_Pythia_QCD_pThat15_NopU_2D_" + selection + "_Official.root";
}
inline bool PassMomentum(double pt, double error) {
  // At exactly 10 GeV the upper cut is NOT applied. Equality at 0.1 passes.
  if(!std::isfinite(pt) || pt <= 0 || !std::isfinite(error) || error < 0) return false;
  return pt <= 10.0 || error / pt <= 0.1;
}
inline bool PassVertex(double distance, double error, double limit) {
  // A missing/undefined significance is not an accepted measurement.
  return std::isfinite(distance) && std::isfinite(error) && error > 0 &&
         std::isfinite(limit) && limit > 0 && std::abs(distance / error) <= limit;
}
} // namespace EECRun3Tracking
#endif

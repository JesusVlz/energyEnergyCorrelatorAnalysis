#ifndef EEC_RUN3_JET_GEOMETRY_H
#define EEC_RUN3_JET_GEOMETRY_H

#include <cmath>
#include <stdexcept>

// Collection/radius and circular-coordinate checks shared by preflight tests.
// The analysis cone radius must identify the native jet collection.
namespace EECJetGeometry {
inline double WrapPhi(double phi) {
  return std::remainder(phi, 2.0 * std::acos(-1.0));
}

inline void ValidateRun3Collection(bool isOO, bool isPpRef, int jetType, double radius) {
  if(!isOO && !isPpRef) return;
  if(jetType < 0 || jetType > 3 || !std::isfinite(radius)) {
    throw std::invalid_argument("Run-3: invalid jet collection or radius");
  }
  const double nativeRadius = isOO ? 0.2 * (jetType + 1) : (jetType == 3 ? 0.8 : 0.4);
  if(std::abs(radius - nativeRadius) > 1e-6) {
    throw std::invalid_argument("Run-3: JetRadius must match the native JetType collection; widening an AK4 cone does not create AK8 jets");
  }
}

inline void RequireRun3Calibration(bool isOO, bool isPpRef, int jetType) {
  if((isOO && jetType != 1) || (isPpRef && jetType == 3)) {
    throw std::invalid_argument("Run-3: radius-matched jet calibration is unavailable for this collection. AK4 payloads must not be applied to AK8 (or other OO radii); use the reader-only diagnostic until applicable payloads are supplied");
  }
}
}

#endif

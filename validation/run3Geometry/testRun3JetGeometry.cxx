#include "../../src/Run3JetGeometry.h"
#include <cassert>
#include <iostream>
#include <functional>

bool rejects(const std::function<void()>& operation) {
  try { operation(); } catch(const std::invalid_argument&) { return true; }
  return false;
}
int main() {
  const double pi = std::acos(-1.0);
  for(double phi : {-pi, -3.0, -1.0, 0.0, 1.0, 3.0, pi}) {
    for(int sign : {-1, 1}) {
      const double control = EECJetGeometry::WrapPhi(phi + sign*pi/2);
      assert(std::abs(control) <= pi);
      assert(std::abs(EECJetGeometry::WrapPhi(control-phi) - sign*pi/2) < 1e-12);
    }
  }
  EECJetGeometry::ValidateRun3Collection(true, false, 1, 0.4);
  EECJetGeometry::ValidateRun3Collection(true, false, 3, 0.8);
  EECJetGeometry::ValidateRun3Collection(false, true, 3, 0.8);
  assert(rejects([] { EECJetGeometry::ValidateRun3Collection(true, false, 1, 0.8); }));
  assert(rejects([] { EECJetGeometry::ValidateRun3Collection(false, true, 1, 0.8); }));
  assert(rejects([] { EECJetGeometry::ValidateRun3Collection(true, false, 3, 0.4); }));
  assert(rejects([] { EECJetGeometry::RequireRun3Calibration(true, false, 3); }));
  assert(rejects([] { EECJetGeometry::RequireRun3Calibration(false, true, 3); }));
  EECJetGeometry::RequireRun3Calibration(true, false, 1);
  EECJetGeometry::RequireRun3Calibration(false, true, 1);
  std::cout << "PASS: canonical phi, native collection/radius checks, AK8 calibration rejection\n";
}

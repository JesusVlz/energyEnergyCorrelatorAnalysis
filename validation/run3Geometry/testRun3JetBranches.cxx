#include "../../src/Run3JetBranchValidation.h"
#include <cassert>
#include <functional>
#include <iostream>

bool rejects(const std::function<void()>& operation, const std::string& branch) {
  try { operation(); } catch(const std::runtime_error& error) {
    return std::string(error.what()).find(branch) != std::string::npos;
  }
  return false;
}
int main() {
  TTree gen("gen", "synthetic branch-presence fixture");
  float value = 0;
  for(const char* name : {"ngen", "genpt", "geneta", "genphi"}) gen.Branch(name, &value);
  EECJetGeometry::RequireRun3JetBranches(&gen, true, true, 0, false);
  assert(rejects([&] { EECJetGeometry::RequireRun3JetBranches(&gen, true, true, 1, false); }, "WTAgeneta"));
  gen.Branch("WTAgeneta", &value);
  assert(rejects([&] { EECJetGeometry::RequireRun3JetBranches(&gen, true, true, 1, false); }, "WTAgenphi"));
  gen.Branch("WTAgenphi", &value);
  EECJetGeometry::RequireRun3JetBranches(&gen, true, true, 1, false);
  assert(rejects([&] { EECJetGeometry::RequireRun3JetBranches(&gen, true, true, 1, true); }, "nref"));
  TTree reco("reco", "synthetic branch-presence fixture");
  for(const char* name : {"nref", "jtpt", "rawpt", "jteta", "jtphi", "trackMax"}) reco.Branch(name, &value);
  assert(rejects([&] { EECJetGeometry::RequireRun3JetBranches(&reco, false, true, 1, false); }, "WTAeta"));
  reco.Branch("WTAeta", &value);
  reco.Branch("WTAphi", &value);
  EECJetGeometry::RequireRun3JetBranches(&reco, false, true, 1, false);
  assert(rejects([&] { EECJetGeometry::RequireRun3JetBranches(&reco, false, true, 1, true); }, "refpt"));
  std::cout << "PASS: missing Run-3 WTA/gen axes and matching fields fail before binding\n";
}

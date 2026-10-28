#ifndef EEC_RUN3_JET_BRANCH_VALIDATION_H
#define EEC_RUN3_JET_BRANCH_VALIDATION_H

#include <TTree.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace EECJetGeometry {
inline void RequireRun3JetBranches(TTree* tree, bool generatorJets, bool isMC, int axis, bool matchJets) {
  if(tree == nullptr) throw std::runtime_error("Run-3: missing native jet tree");
  std::vector<const char*> required;
  if(generatorJets) {
    required = {"ngen", "genpt", "geneta", "genphi"};
    if(axis == 1) required.insert(required.end(), {"WTAgeneta", "WTAgenphi"});
  } else {
    required = {"nref", "jtpt", "rawpt", "jteta", "jtphi", "trackMax"};
    if(axis == 1) required.insert(required.end(), {"WTAeta", "WTAphi"});
  }
  if(isMC && matchJets) {
    required.insert(required.end(), {"nref", "refpt", "refeta", "refphi", "refparton_flavorForB"});
    if(generatorJets) {
      required.insert(required.end(), {"jtpt", "rawpt"});
      if(axis == 1) required.insert(required.end(), {"WTAeta", "WTAphi"});
      else required.insert(required.end(), {"jteta", "jtphi"});
    } else {
      required.insert(required.end(), {"ngen", "genpt", "geneta", "genphi"});
      if(axis == 1) required.insert(required.end(), {"WTAgeneta", "WTAgenphi"});
    }
  }
  for(const char* branch : required) {
    if(tree->GetBranch(branch) == nullptr) {
      throw std::runtime_error(std::string("Run-3 native jet tree ") + tree->GetName() + " is missing required branch " + branch);
    }
  }
}
}

#endif

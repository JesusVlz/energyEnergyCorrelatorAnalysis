#ifndef EEC_RUN3_JET_SELECTION_INPUTS_H
#define EEC_RUN3_JET_SELECTION_INPUTS_H

#include <array>
#include <limits>
#include <stdexcept>
#include <string>

#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TTree.h>

namespace EECJetSelection {

// Detector eta/phi are the standard jet coordinates, never the analysis WTA axis.
struct JetInputs {
  double eta, phi, rawPt, storedPt;
  double neutralHadronFraction, neutralEmFraction;
  double chargedHadronFraction, chargedEmFraction, muonFraction;
  int chargedMultiplicity, neutralMultiplicity, numberOfConstituents;
  // The native HiForest five-category counts omit HF hadrons and HF EM objects.
  // A neutral/constituent lower bound above an ID threshold proves a pass; a
  // lower bound below it cannot establish a failing Jet ID decision.
  bool multiplicitiesAreLowerBounds = false;
};

namespace Detail {

inline void RequireJetLeaf(TTree* tree, const char* name, const char* type,
                           const char* count = nullptr) {
  TBranch* branch = tree ? tree->GetBranch(name) : nullptr;
  TLeaf* leaf = branch ? branch->GetLeaf(name) : nullptr;
  if(!leaf || !branch->GetListOfLeaves() || branch->GetListOfLeaves()->GetEntries() != 1 ||
     std::string(leaf->GetTypeName()) != type || leaf->GetLenStatic() != 1 ||
     (count ? (!leaf->GetLeafCount() || std::string(leaf->GetLeafCount()->GetName()) != count)
            : leaf->GetLeafCount() != nullptr)) {
    throw std::runtime_error(std::string("Run-3 jet selection: incompatible or missing ") +
                             (tree ? tree->GetName() : "jet tree") + "/" + name +
                             "; expected " + type + (count ? std::string("[") + count + "]" : " scalar"));
  }
}

// Native branch types are pinned by HiInclusiveJetAnalyzer: fractions /F,
// multiplicities /I, all indexed by nref. Validate before any payload is read.
inline void ValidateNativeJetInputs(TTree* tree) {
  RequireJetLeaf(tree, "nref", "Int_t");
  for(const char* name : {"jtpt", "rawpt", "jteta", "jtphi", "jtPfCHF", "jtPfNHF",
                          "jtPfCEF", "jtPfNEF", "jtPfMUF"})
    RequireJetLeaf(tree, name, "Float_t", "nref");
  for(const char* name : {"jtPfCHM", "jtPfNHM", "jtPfCEM", "jtPfNEM", "jtPfMUM"})
    RequireJetLeaf(tree, name, "Int_t", "nref");
}

template <size_t Capacity>
class NativePFJetInputs {
public:
  void Bind(TTree* tree) {
    ValidateNativeJetInputs(tree);
    const char* floatNames[] = {"jtPfCHF", "jtPfNHF", "jtPfCEF", "jtPfNEF", "jtPfMUF"};
    const char* intNames[] = {"jtPfCHM", "jtPfNHM", "jtPfCEM", "jtPfNEM", "jtPfMUM"};
    for(size_t i = 0; i < 5; ++i) {
      tree->SetBranchStatus(floatNames[i], 1);
      tree->SetBranchStatus(intNames[i], 1);
      if(tree->SetBranchAddress(floatNames[i], fFractions[i].data()) < 0 ||
         tree->SetBranchAddress(intNames[i], fMultiplicities[i].data()) < 0)
        throw std::runtime_error("Run-3 jet selection: cannot bind native PF composition");
    }
  }

  JetInputs Get(size_t i, double eta, double phi, double rawPt, double storedPt) const {
    if(i >= Capacity) throw std::out_of_range("Run-3 jet selection: jet index exceeds native buffer");
    long long charged = 0, neutral = 0;
    for(size_t j = 0; j < 5; ++j) {
      if(fMultiplicities[j][i] < 0) throw std::runtime_error("Run-3 jet selection: negative native PF multiplicity");
      if(j == 0 || j == 2 || j == 4) charged += fMultiplicities[j][i];
      else neutral += fMultiplicities[j][i];
    }
    if(charged + neutral > std::numeric_limits<int>::max())
      throw std::runtime_error("Run-3 jet selection: native PF multiplicity sum overflows int");
    return JetInputs{eta, phi, rawPt, storedPt, fFractions[1][i], fFractions[3][i],
                     fFractions[0][i], fFractions[2][i], fFractions[4][i],
                     static_cast<int>(charged), static_cast<int>(neutral),
                     static_cast<int>(charged + neutral), true};
  }

private:
  std::array<std::array<Float_t, Capacity>, 5> fFractions{};
  std::array<std::array<Int_t, Capacity>, 5> fMultiplicities{};
};

} // namespace Detail
} // namespace EECJetSelection

#endif

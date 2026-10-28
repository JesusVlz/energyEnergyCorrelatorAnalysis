#include "../../src/GeneratorLevelForestReader.h"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
  if(argc != 3) throw std::invalid_argument("usage: testNativeGeneratorReader OO.root pp.root");
  for(int jetType : {1, 3}) {
    TFile input(argv[1], "READ");
    GeneratorLevelForestReader reader(ForestReader::kOOMC, 0, jetType, 1, true);
    reader.ReadForestFromFile(&input);
    auto tree = dynamic_cast<TTree*>(input.Get(jetType == 1 ? "akCs4PFJetAnalyzer/t" : "akCs8PFJetAnalyzer/t"));
    assert(tree != nullptr);
    int compared = 0;
    int matched = 0;
    for(int event = 0; event < 100 && event < tree->GetEntries(); ++event) {
      tree->GetEntry(event);
      reader.GetEvent(event);
      assert(reader.GetNJets() == tree->GetLeaf("ngen")->GetValue());
      for(int j = 0; j < reader.GetNJets(); ++j) {
        ++compared;
        assert(std::abs(reader.GetJetPt(j)-tree->GetLeaf("genpt")->GetValue(j)) < 1e-5);
        assert(std::abs(reader.GetJetPhi(j)-tree->GetLeaf("WTAgenphi")->GetValue(j)) < 1e-5);
        assert(std::abs(reader.GetJetEta(j)-tree->GetLeaf("WTAgeneta")->GetValue(j)) < 1e-5);
        const int match = static_cast<int>(tree->GetLeaf("genmatchindex")->GetValue(j));
        assert(reader.HasMatchingJet(j) == (match >= 0));
        if(match < 0) continue;
        ++matched;
        assert(std::abs(reader.GetMatchedPt(j)-tree->GetLeaf("rawpt")->GetValue(match)) < 1e-5);
        assert(std::abs(reader.GetMatchedPhi(j)-tree->GetLeaf("WTAphi")->GetValue(match)) < 1e-5);
        assert(std::abs(reader.GetMatchedEta(j)-tree->GetLeaf("WTAeta")->GetValue(match)) < 1e-5);
      }
    }
    assert(compared > 0);
    assert(matched > 0);
    std::cout << "PASS: OO native generator/WTA JetType " << jetType
              << ", jets checked " << compared << ", matched " << matched << "\n";
  }
  bool missingRejected = false;
  try {
    TFile input(argv[2], "READ");
    GeneratorLevelForestReader reader(ForestReader::kPpRef5p36TeVMC, 0, 3, 1, true);
    reader.ReadForestFromFile(&input);
  } catch(const std::runtime_error& error) {
    missingRejected = std::string(error.what()).find("ak8PFJetAnalyzer/t") != std::string::npos;
  }
  assert(missingRejected);
  std::cout << "PASS: unavailable pp AK8 is rejected without AK4 fallback\n";
}

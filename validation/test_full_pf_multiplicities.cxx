// From the repository root:
// Compile this file and src/Run3JetSelection.cxx with root-config --cflags --libs.
// /tmp/test_full_pf_multiplicities
#include "../src/Run3JetSelection.h"
#include <TTree.h>
#include <cassert>
#include <climits>
#include <map>
#include <stdexcept>
#include <string>

namespace {
struct Fixture {
  int nref=1;
  std::map<std::string,float> fractions;
  std::map<std::string,int> counts;
  TTree tree;
  explicit Fixture(int mode, int charged=1, int neutral=0, int fullCharged=1, int fullNeutral=1)
    : tree("t","PF multiplicity fixture") {
    tree.Branch("nref",&nref,"nref/I");
    for(const char* n:{"jtpt","rawpt","jteta","jtphi","jtPfCHF","jtPfNHF","jtPfCEF","jtPfNEF","jtPfMUF"}) {
      fractions[n]=0.1f;
      tree.Branch(n,&fractions[n],(std::string(n)+"[nref]/F").c_str());
    }
    fractions["jtPfCHF"]=0.5f;
    for(const char* n:{"jtPfCHM","jtPfNHM","jtPfCEM","jtPfNEM","jtPfMUM"}) {
      counts[n]=0;
      tree.Branch(n,&counts[n],(std::string(n)+"[nref]/I").c_str());
    }
    counts["jtPfCHM"]=charged;
    counts["jtPfNHM"]=neutral;
    if(mode) {
      counts["jtPfCM"]=fullCharged;
      tree.Branch("jtPfCM",&counts["jtPfCM"],"jtPfCM[nref]/I");
    }
    if(mode==1) {
      counts["jtPfNM"]=fullNeutral;
      tree.Branch("jtPfNM",&counts["jtPfNM"],"jtPfNM[nref]/I");
    } else if(mode==3) {
      fractions["jtPfNM"]=float(fullNeutral);
      tree.Branch("jtPfNM",&fractions["jtPfNM"],"jtPfNM[nref]/F");
    }
    tree.Fill();
  }
};

template<class F> void Throws(F operation) {
  bool threw=false;
  try {operation();} catch(const std::exception&) {threw=true;}
  assert(threw);
}

EECJetSelection::JetInputs Read(Fixture& f) {
  EECJetSelection::Detail::NativePFJetInputs<2> input;
  input.Bind(&f.tree);
  assert(f.tree.GetEntry(0)>0);
  return input.Get(0,0.2,0.1,50,60);
}
}

int main() {
  using namespace EECJetSelection;
  Fixture oldUnresolved(0);
  auto legacy=Read(oldUnresolved);
  Throws([&]{PassJetID(legacy,JetIDPrescription::MiniPOGpp);});
  Fixture oldSufficient(0,1,1);
  assert(PassJetID(Read(oldSufficient),JetIDPrescription::MiniPOGpp));

  // Same category lower bound, now with one measured extra neutral (e.g. HF).
  Fixture withHF(1,1,0,1,1);
  auto exact=Read(withHF);
  assert(PassJetID(exact,JetIDPrescription::MiniPOGpp));
  // HINpp's separate pinned category-sum prescription must not change.
  assert(!PassJetID(exact,JetIDPrescription::HINpp));
  Fixture oneConstituent(1,1,0,1,0);
  assert(!PassJetID(Read(oneConstituent),JetIDPrescription::MiniPOGpp));

  for(double eta:{2.8,3.5}) {
    auto forward=exact;forward.eta=eta;forward.neutralHadronFraction=0.3;
    forward.fullNeutralMultiplicity=eta<3?2:11;
    assert(PassJetID(forward,JetIDPrescription::TightLeptonVeto));
    forward.fullNeutralMultiplicity=eta<3?1:10;
    assert(!PassJetID(forward,JetIDPrescription::TightLeptonVeto));
  }
  Fixture partial(2),wrongType(3),belowCategories(1,1,0,0,1),overflow(1,1,0,INT_MAX,1);
  Throws([&]{Read(partial);});Throws([&]{Read(wrongType);});
  Throws([&]{Read(belowCategories);});Throws([&]{Read(overflow);});
  Detail::NativePFJetInputs<2> chain;
  chain.Bind(&withHF.tree);
  chain.ValidateSchema(&withHF.tree);
  Throws([&]{chain.ValidateSchema(&oldUnresolved.tree);});
  assert(PassJetID(legacy,JetIDPrescription::MiniPOGOO));
}

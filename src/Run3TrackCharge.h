#ifndef EEC_RUN3_TRACK_CHARGE_H
#define EEC_RUN3_TRACK_CHARGE_H
#include <algorithm>
#include <cmath>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <TBranch.h>
#include <TChain.h>
#include <TLeaf.h>
#include <TTree.h>

namespace EECRun3Tracking {
// ROOT writes into buffers whose addresses stay stable, including when a reader
// is shallow-copied in the same way as its existing TChain pointers.
class TrackCharge {
  struct Buffers {
    std::vector<char> chars; std::vector<char>* pc = &chars;
    std::vector<signed char> schars; std::vector<signed char>* ps = &schars;
    std::vector<unsigned char> uchars; std::vector<unsigned char>* pu = &uchars;
    std::vector<short> shorts; std::vector<short>* ph = &shorts;
    std::vector<int> ints; std::vector<int>* pi = &ints;
  };
  TTree* tree_ = nullptr;
  std::shared_ptr<Buffers> buffers_;
  std::string kind_;
  static std::string Kind(TBranch* branch) {
    if(!branch) throw std::runtime_error("Run-3 tracking requires physical trkCharge; re-forest input that omitted it");
    std::string value = branch->GetClassName();
    for(std::size_t p; (p=value.find("std::"))!=std::string::npos;) value.erase(p,5);
    value.erase(std::remove_if(value.begin(),value.end(),[](unsigned char c){return std::isspace(c);}),value.end());
    return value;
  }
public:
  void Bind(TTree* tree) {
    if(!tree) throw std::runtime_error("Run-3 track charge: missing track tree");
    if(tree->GetEntries() > 0) tree->LoadTree(0);
    tree_ = tree; buffers_ = std::make_shared<Buffers>();
    TBranch* branch = tree_->GetBranch("trkCharge");
    kind_ = Kind(branch);
    tree_->SetBranchStatus("trkCharge",1);
    int code = 0;
    if(kind_ == "vector<char>") code = tree_->SetBranchAddress("trkCharge",&buffers_->pc);
    else if(kind_ == "vector<signedchar>") code = tree_->SetBranchAddress("trkCharge",&buffers_->ps);
    else if(kind_ == "vector<unsignedchar>") code = tree_->SetBranchAddress("trkCharge",&buffers_->pu);
    else if(kind_ == "vector<short>") code = tree_->SetBranchAddress("trkCharge",&buffers_->ph);
    else if(kind_ == "vector<int>") code = tree_->SetBranchAddress("trkCharge",&buffers_->pi);
    else if(kind_.empty()) {
      // AOD numeric leaves use ROOT's internal leaf buffer, not a guessed C++ type.
      TLeaf* leaf = branch->GetLeaf("trkCharge");
      const std::string type = leaf ? leaf->GetTypeName() : "";
      if(type != "Int_t" && type != "Short_t" && type != "Char_t")
        throw std::runtime_error("Unsupported Run-3 AOD trkCharge leaf type: " + type);
      if(!leaf->GetLeafCount() || std::string(leaf->GetLeafCount()->GetName()) != "nTrk")
        throw std::runtime_error("Run-3 AOD trkCharge must be indexed by nTrk");
    } else throw std::runtime_error("Unsupported Run-3 trkCharge class: " + kind_);
    if(code < 0) throw std::runtime_error("Cannot bind physical Run-3 trkCharge");
    branch->SetAutoDelete(false);
  }
  int Get(std::size_t index) const {
    if(!tree_ || !buffers_) throw std::runtime_error("Run-3 track charge requested before binding");
    TTree* current = tree_;
    if(auto* chain = dynamic_cast<TChain*>(tree_)) current = chain->GetTree();
    TBranch* branch = current ? current->GetBranch("trkCharge") : nullptr;
    if(Kind(branch) != kind_) throw std::runtime_error("Run-3 trkCharge schema changed inside a chain");
    int value = 0;
    if(kind_ == "vector<char>") value = static_cast<signed char>(buffers_->pc->at(index));
    else if(kind_ == "vector<signedchar>") value = buffers_->ps->at(index);
    else if(kind_ == "vector<unsignedchar>") value = static_cast<signed char>(buffers_->pu->at(index));
    else if(kind_ == "vector<short>") value = buffers_->ph->at(index);
    else if(kind_ == "vector<int>") value = buffers_->pi->at(index);
    else {
      TLeaf* leaf = branch->GetLeaf("trkCharge");
      if(!leaf || index >= static_cast<std::size_t>(leaf->GetLen()))
        throw std::out_of_range("Run-3 track charge index outside AOD leaf");
      const double charge = leaf->GetValue(index);
      if(!std::isfinite(charge) || charge != std::nearbyint(charge))
        throw std::runtime_error("Invalid Run-3 physical track charge");
      value = static_cast<int>(charge);
    }
    // MiniAOD track trees can retain neutral packed candidates. Return their
    // stored zero so PassTrackCuts can reject them via its existing charge cut.
    // Missing charge remains an error; never replace it with a fabricated sign.
    if(value < -1 || value > 1) throw std::runtime_error("Invalid Run-3 stored track charge; expected -1, 0 or +1");
    return value;
  }
};
} // namespace EECRun3Tracking
#endif

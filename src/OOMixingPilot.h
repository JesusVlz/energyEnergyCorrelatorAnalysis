#ifndef OO_MIXING_PILOT_H
#define OO_MIXING_PILOT_H

// Opt-in MC GenGen pilot primitives. No reconstruction corrections or ROOT types.
#include "OOEventIdentity.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace oo_eec {
struct PhysicalId {
  std::uint64_t run=0, lumi=0, event=0;
  bool operator<(const PhysicalId& other) const { return std::tie(run,lumi,event)<std::tie(other.run,other.lumi,other.event); }
  bool operator==(const PhysicalId& other) const { return run==other.run && lumi==other.lumi && event==other.event; }
};
inline std::uint64_t PilotKey(const PhysicalId& id, bool donor){
  return EventKey(donor ? "OO2025_HIJING_existing_owner_pool_pilot" : "OO2025_embedded_existing_owner_pilot", id.run,id.lumi,id.event);
}
inline unsigned int PilotPool(const PhysicalId& id, std::uint64_t seed){ return (PilotKey(id,true)^seed)%2; }
inline bool IsPilotMatchingParticle(double pt,double eta,int charge,int status,int subevent,bool pureHijing){
  return std::isfinite(pt) && std::isfinite(eta) && pt>0.7 && std::abs(eta)<2.4 && charge!=0 && status==1 && (pureHijing || subevent>0);
}
inline double PilotHfOff(double raw){ return raw; }
inline double PilotRelativePhi(double phi,double axis){ return std::atan2(std::sin(phi-axis),std::cos(phi-axis)); }
struct BackgroundFingerprint {
  std::uint64_t count=0,sum=0,xorValue=0;
  void add(float pt,float eta,float phi,int charge){
    std::uint32_t a,b,c;
    std::memcpy(&a,&pt,4); std::memcpy(&b,&eta,4); std::memcpy(&c,&phi,4);
    const auto h=SplitMix64(a)^SplitMix64(std::uint64_t(b)<<1)^SplitMix64(std::uint64_t(c)<<2)^SplitMix64(std::uint64_t(std::uint32_t(charge))<<3);
    ++count;sum+=h;xorValue^=SplitMix64(h);
  }
  std::uint64_t key() const { return EventKey("OO2025_HIJING_stored_kinematics_v1",count,sum,xorValue); }
  std::string text() const { std::ostringstream out;out<<count<<':'<<std::hex<<std::setfill('0')<<std::setw(16)<<sum<<':'<<std::setw(16)<<xorValue;return out.str(); }
};
struct PilotEvent {
  PhysicalId id;
  int entry=-1, nh=0, hiBin=-1;
  double vz=0, weight=1, hf=0;
  bool eligible=false, backgroundKnown=false, backgroundExcluded=false, physicalExcluded=false;
  std::uint64_t background=0;
  BackgroundFingerprint fingerprint;
};
inline unsigned int PilotPool(const PilotEvent& event,std::uint64_t seed){
  return event.backgroundKnown ? (event.background^seed)%2 : PilotPool(event.id,seed);
}
inline bool PilotSharesBackground(const PilotEvent& a,const PilotEvent& b){
  return a.id==b.id || (a.backgroundKnown && b.backgroundKnown && a.background==b.background);
}
inline void ValidatePilotPool(const std::vector<PilotEvent>& pool){
  std::set<PhysicalId> seen;
  for(const auto& event:pool){
    if(!seen.insert(event.id).second) throw std::runtime_error("Duplicate physical run/lumi/event identity in OO pilot pool");
    if(!std::isfinite(event.vz) || !std::isfinite(event.weight) || event.nh<0)
      throw std::runtime_error("Invalid OO pilot pool matching variable or generator weight");
  }
}
inline void ExcludePilotHoldout(PilotEvent& donor,const std::set<PhysicalId>& signalIds,const std::set<std::uint64_t>& signalBackgrounds){
  donor.physicalExcluded=signalIds.count(donor.id)!=0;
  donor.backgroundExcluded=donor.backgroundKnown && signalBackgrounds.count(donor.background)!=0;
  if(donor.physicalExcluded || donor.backgroundExcluded) donor.eligible=false;
}
struct PilotAssignment {
  std::array<int,2> partners{{-1,-1}}; // indices into eligible metadata, never entries/sentinels
  std::size_t inspected=0;
  std::string reason="no_pool_A";
  bool matched() const { return partners[0]>=0 && partners[1]>=0; }
};
inline PilotAssignment SelectPilotPartners(const PilotEvent& signal,const std::vector<PilotEvent>& pool,std::uint64_t seed){
  if(!std::isfinite(signal.vz) || signal.nh<0) throw std::invalid_argument("Invalid OO signal matching coordinate");
  const double tolerance=std::max(1.,.01*signal.nh);
  PilotAssignment result;
  // Exactly one finite scan of the explicitly capped pool; tolerances never widen.
  std::array<std::tuple<double,std::uint64_t,PhysicalId>,2> best;
  std::array<bool,2> found{{false,false}};
  for(std::size_t index=0; index<pool.size(); ++index){
    ++result.inspected;
    const auto& candidate=pool[index];
    if(!candidate.eligible || PilotSharesBackground(signal,candidate)) continue;
    const double dz=std::abs(candidate.vz-signal.vz), dn=std::abs(candidate.nh-signal.nh);
    if(dz>.5 || dn>tolerance) continue;
    const auto fold=PilotPool(candidate,seed);
    const auto rank=std::make_tuple(dz/.5+dn/tolerance,SplitMix64(PilotKey(signal.id,false)^PilotKey(candidate.id,true)^seed),candidate.id);
    if(!found[fold] || rank<best[fold]){ found[fold]=true; best[fold]=rank; result.partners[fold]=index; }
  }
  if(!found[0] || !found[1]){
    result.reason=!found[0] ? "no_pool_A" : "no_pool_B";
    result.partners={{-1,-1}}; // no partial estimator
  } else if(PilotSharesBackground(pool[result.partners[0]],pool[result.partners[1]])){
    result.reason="shared_partner_background"; result.partners={{-1,-1}};
  } else result.reason="matched";
  return result;
}
// Campaign policy: reserve one immutable donor prefix across all jobs. Signals
// sharing any reserved raw identity or observed nonempty background are excluded
// before the common exposure, even if the reserved donor itself fails event cuts.
inline std::string PilotReservedSignalReason(const PilotEvent& signal,
    const std::set<PhysicalId>& reservedIds,const std::set<std::uint64_t>& reservedBackgrounds){
  if(reservedIds.count(signal.id)) return "signal_in_reserved_donor_physical_id";
  if(signal.backgroundKnown && reservedBackgrounds.count(signal.background))
    return "signal_in_reserved_donor_background";
  return "";
}
struct ClosurePilotState {
  bool enabled=false;
  std::uint64_t seed=0;
  int maxSignal=0,maxPool=0,signalRead=0,jetIndex=-1,partitionPolicy=0;
  PilotEvent signal;
  PilotAssignment assignment;
  std::vector<PilotEvent> pool;
  std::set<PhysicalId> signalIds;
  std::set<PhysicalId> signalHoldoutIds;
  std::set<PhysicalId> reservedPoolIds;
  std::set<std::uint64_t> reservedPoolBackgrounds;
  std::set<std::uint64_t> signalBackgrounds;
  std::map<int,unsigned int> uses;
  std::ofstream partners,poolLedger,jets,contributions;
};
} // namespace oo_eec
#endif

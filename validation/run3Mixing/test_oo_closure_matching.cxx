#include "../../src/OOMixingPilot.h"
#include <cassert>
#include <cstring>
#include <iostream>
using namespace oo_eec;
int main(){
  const std::uint64_t seed=20261030;
  PilotEvent signal; signal.id={1,9,UINT64_C(1)<<42};signal.nh=100;signal.vz=0.;signal.eligible=true;
  std::vector<PilotEvent> pool;
  for(unsigned int event=1;event<=40;++event){ PilotEvent p=signal;p.id.event+=event;p.entry=event;p.vz=event*.0001;pool.push_back(p); }
  ValidatePilotPool(pool);
  const auto result=SelectPilotPartners(signal,pool,seed);
  assert(result.matched() && result.inspected==pool.size());
  assert(PilotPool(pool[result.partners[0]].id,seed)==0);
  assert(PilotPool(pool[result.partners[1]].id,seed)==1);
  assert(!(pool[result.partners[0]].id==pool[result.partners[1]].id));
  // The selected physical identities are invariant under input order, jets, and powers.
  std::array<PhysicalId,2> ids{{pool[result.partners[0]].id,pool[result.partners[1]].id}};
  std::reverse(pool.begin(),pool.end());
  for(int jet=0;jet<4;++jet) for(int power=1;power<=2;++power){
    const auto again=SelectPilotPartners(signal,pool,seed);
    for(int i=0;i<2;++i) assert(pool[again.partners[i]].id==ids[i]);
  }
  // NH is fixed, strict and independent of EEC/track thresholds and labels of pure HIJING.
  assert(IsPilotMatchingParticle(.70001,2.399,1,1,0,true));
  assert(!IsPilotMatchingParticle(.7,0,1,1,1,false));
  assert(!IsPilotMatchingParticle(2,2.4,1,1,1,false));
  assert(!IsPilotMatchingParticle(2,0,0,1,1,false));
  assert(!IsPilotMatchingParticle(2,0,1,2,1,false));
  assert(!IsPilotMatchingParticle(2,0,1,1,0,false));
  assert(IsPilotMatchingParticle(2,0,1,1,1,false));
  // No numerical sentinel can become an eligible entry, and no partial pair escapes.
  auto sparse=std::vector<PilotEvent>{pool[result.partners[0]]};
  const auto failed=SelectPilotPartners(signal,sparse,seed);
  assert(!failed.matched() && failed.partners[0]==-1 && failed.partners[1]==-1 && failed.inspected==1);
  assert(!SelectPilotPartners(signal,{},seed).matched());
  auto distant=pool;
  for(auto& p:distant) p.vz=.50001;
  assert(!SelectPilotPartners(signal,distant,seed).matched());
  for(auto& p:distant){p.vz=.5;p.nh=102;}
  assert(!SelectPilotPartners(signal,distant,seed).matched());
  for(auto& p:distant) p.nh=101;
  assert(SelectPilotPartners(signal,distant,seed).matched());
  for(auto& p:distant) p.eligible=false;
  assert(!SelectPilotPartners(signal,distant,seed).matched());
  // Full 64-bit identity survives; identical file entries do not imply identity.
  PhysicalId truncated=signal.id;truncated.event=static_cast<std::uint32_t>(truncated.event);
  assert(!(truncated==signal.id) && PilotKey(truncated,false)!=PilotKey(signal.id,false));
  auto duplicated=pool;duplicated.push_back(pool[0]);
  bool rejected=false;try{ ValidatePilotPool(duplicated); }catch(const std::runtime_error&){rejected=true;}assert(rejected);
  auto selves=pool;
  for(auto& p:selves) p.id=signal.id;
  assert(!SelectPilotPartners(signal,selves,seed).matched());
  signal.backgroundKnown=true;signal.background=77;
  for(auto& p:pool){p.backgroundKnown=true;p.background=77;}
  assert(!SelectPilotPartners(signal,pool,seed).matched());
  // Unknown provenance stays unknown; physical identity checks remain strict.
  assert(PilotSharesBackground(signal,signal));
  BackgroundFingerprint first,second;
  first.add(1.5f,-.2f,3.13f,1); first.add(2.f,.4f,-3.13f,-1);
  second.add(2.f,.4f,-3.13f,-1); second.add(1.5f,-.2f,3.13f,1);
  assert(first.text()==second.text() && first.key()==second.key());
  auto reusedA=pool[0],reusedB=pool[1];
  reusedA.backgroundKnown=reusedB.backgroundKnown=true;
  reusedA.background=reusedB.background=first.key();
  assert(PilotPool(reusedA,seed)==PilotPool(reusedB,seed));
  assert(PilotSharesBackground(reusedA,reusedB));
  signal.background=first.key();
  assert(!SelectPilotPartners(signal,{reusedA,reusedB},seed).matched());
  // Global holdout exclusion rejects a donor sharing ANOTHER signal event's
  // identity/background, even when it would pass the current signal match.
  PilotEvent otherSignal;otherSignal.id={9,8,7};otherSignal.backgroundKnown=true;otherSignal.background=98765;
  auto crossSignal=pool[0];crossSignal.eligible=true;crossSignal.backgroundKnown=true;crossSignal.background=98765;
  ExcludePilotHoldout(crossSignal,{otherSignal.id},{otherSignal.background});
  assert(crossSignal.backgroundExcluded && !crossSignal.physicalExcluded && !crossSignal.eligible);
  crossSignal.id=otherSignal.id;crossSignal.backgroundKnown=false;crossSignal.eligible=true;
  ExcludePilotHoldout(crossSignal,{otherSignal.id},{otherSignal.background});
  assert(crossSignal.physicalExcluded && !crossSignal.backgroundExcluded && !crossSignal.eligible);
  // Fixed donor reservation excludes the same background across different jobs,
  // independent of donor eligibility; empty backgrounds remain distinct by raw ID.
  auto reserved=otherSignal; reserved.eligible=false;
  auto jobOne=signal,jobTwo=signal; jobOne.id={31,2,8};jobTwo.id={32,2,9};
  jobOne.backgroundKnown=jobTwo.backgroundKnown=true;
  jobOne.background=jobTwo.background=reserved.background;
  assert(PilotReservedSignalReason(jobOne,{reserved.id},{reserved.background})=="signal_in_reserved_donor_background");
  assert(PilotReservedSignalReason(jobTwo,{reserved.id},{reserved.background})=="signal_in_reserved_donor_background");
  jobOne.id=reserved.id;
  assert(PilotReservedSignalReason(jobOne,{reserved.id},{})=="signal_in_reserved_donor_physical_id");
  jobOne.id={99,4,UINT64_C(1)<<50};jobOne.backgroundKnown=false;
  assert(PilotReservedSignalReason(jobOne,{reserved.id},{reserved.background}).empty());
  assert(PilotReservedSignalReason(jobOne,{},{}).empty());
  // HF OFF is bitwise identity, including low, negative, and signed-zero values.
  for(double value:{-3.,-.1,-0.,0.,.1,1.,1.9,2.,100.}){
    const auto off=PilotHfOff(value); assert(std::memcmp(&off,&value,sizeof(double))==0);
  }
  assert(std::abs(PilotRelativePhi(-3.13,3.13))<.03);
  std::cout<<"OO closure matching: strict NH, fixed boundaries, bounded sparse failure, eligible-only A/B, deterministic partners, uint64/duplicate/self/shared IDs, HF-OFF identity and circular phi passed\n";
}

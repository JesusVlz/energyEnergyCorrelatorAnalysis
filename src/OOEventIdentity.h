#ifndef OO_EVENT_IDENTITY_H
#define OO_EVENT_IDENTITY_H

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

// Shared physical event identity for folds/replicas across radii, powers and methods.
// sampleNamespace identifies the collision sample (e.g. OO2025_MC_campaignX),
// never radius, power, file path or entry. All physical IDs must be available;
// callers must reject missing/duplicate identities rather than fall back to entries.
namespace oo_eec {
inline std::uint64_t EventKey(const std::string& sampleNamespace, std::uint64_t run,
                              std::uint64_t lumi, std::uint64_t event){
  if(sampleNamespace.empty()) throw std::invalid_argument("A physical sample namespace is required");
  std::uint64_t key = UINT64_C(14695981039346656037);
  const auto append = [&key](unsigned char byte){ key ^= byte; key *= UINT64_C(1099511628211); };
  for(unsigned char byte : sampleNamespace) append(byte);
  append(0);
  for(std::uint64_t value : {run, lumi, event}){
    for(int byte=0; byte<8; ++byte){ append(static_cast<unsigned char>(value & 255)); value >>= 8; }
  }
  return key;
}
inline unsigned int EventFold(std::uint64_t key, unsigned int folds=2){
  if(folds < 2) throw std::invalid_argument("At least two independent event folds are required");
  return key % folds;
}
inline std::uint64_t SplitMix64(std::uint64_t value){
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}
inline unsigned int ReplicaWeight(std::uint64_t key, std::uint64_t replica){
  const std::uint64_t random = SplitMix64(key ^ SplitMix64(replica + UINT64_C(0x9e3779b97f4a7c15)));
  const double candidate = (static_cast<double>(random >> 11) + 0.5) / 9007199254740992.0;
  const double uniform = candidate < 1.0 ? candidate : std::nextafter(1.0, 0.0);
  // Deterministic Poisson(1) multiplier shared by every contribution of this event.
  double probability = std::exp(-1.0), cumulative = probability;
  unsigned int count=0;
  while(uniform > cumulative){ ++count; probability /= count; cumulative += probability; }
  return count;
}
} // namespace oo_eec
#endif

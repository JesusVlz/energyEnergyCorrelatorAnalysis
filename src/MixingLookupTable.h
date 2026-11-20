#ifndef MIXING_LOOKUP_TABLE_H
#define MIXING_LOOKUP_TABLE_H

#include <Rtypes.h>
#include <TRandom3.h>
#include <TString.h>

#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

struct MixingLookupConfig {
  Int_t dataType = -1;
  Bool_t megaSkimMode = false;
  Double_t vzCut = 15.0;
  Int_t triggerSelection = 0;
  Bool_t applyJetVeto = false;
  std::string jetVetoMapFile;
  std::string jetVetoMapHistogram;
  std::string jetVetoPtConvention;
  Int_t run3JetID = 0;
  std::string run3JetIDTable;
  Bool_t useRun3JEC = false;
  Bool_t usePpL1FastJet = false;
  Bool_t ooClosurePilot = false;
  Bool_t requireGenMultiplicity = false;
};

struct MixingDonorRecord {
  Int_t globalIndex = -1;
  Int_t fileIndex = -1;
  Int_t localEntry = -1;
  UInt_t run = 0;
  UInt_t lumi = 0;
  ULong64_t eventNumber = 0;
  Bool_t hasPhysicalId = false;
  Double_t vz = 0.0;
  Int_t storedHiBin = -1;
  Double_t hfSum = 0.0;
  Double_t rawHfSum = 0.0;
  Double_t eventWeight = 1.0;
  Int_t genMultiplicity = -1;
  Bool_t passStandardEventCuts = false;
  Bool_t passEventCuts = false;
  Bool_t passMBTrigger = true;
  Bool_t pvBit = false;
  Bool_t hfCoincBit = false;
  Bool_t pileupBit = false;
  Bool_t vzBit = false;

  // Opt-in OOClosurePilot fields
  Int_t ooNh = -1;
  Bool_t backgroundKnown = false;
  std::uint64_t backgroundKey = 0;
  std::string backgroundSignature;

  Int_t NormalizedHiBin() const { return storedHiBin < 0 ? 1 : storedHiBin; }

  Bool_t IsEligible(Int_t dataType, Bool_t megaSkimMode) const {
    if (megaSkimMode) return true;
    const Bool_t mbOk = (dataType != 9 /* ForestReader::kOO */) || passMBTrigger;
    return passEventCuts && mbOk;
  }
};

class MixingLookupTable {
public:
  enum MatchMode {
    kMatchHiBin = 0,
    kMatchHFEnergy = 1,
    kMatchMultiplicity = 2
  };

  MixingLookupTable() = default;

  static std::string CanonicalizeMixingPath(const std::string& rawPath);
  static std::vector<std::string> CandidateLookupPaths(const std::string& mixingListPath);
  static std::uint64_t ComputeShardHash(Int_t jobShardId, const std::vector<TString>& signalFiles);

  void Clear();
  void SetFileMetadata(const std::vector<TString>& mixingFiles, const std::vector<Int_t>& fileEntryCounts);
  void AddRecord(const MixingDonorRecord& record);
  void ReserveRecords(std::size_t count) { fRecords.reserve(count); }

  bool WriteToFile(const std::string& outputPath, const MixingLookupConfig& config) const;
  bool TryLoadFromFile(const std::string& inputPath,
                       const MixingLookupConfig& expected,
                       const std::vector<TString>& mixingFiles,
                       std::string& rejectReason);

  // Select a deterministic per-job file window across mixingFiles.
  // When all records are already loaded, automatically expands if the initial window has < 2 eligible events.
  std::vector<Int_t> SelectActiveFiles(std::size_t nTotalFiles,
                                       Int_t maxOpenFiles,
                                       std::uint64_t shardHash,
                                       Int_t dataType,
                                       Bool_t megaSkimMode) const;

  // Build the 2D (vz, activity) lookup index over the selected active files (empty activeFiles = all loaded records).
  void BuildIndex(Int_t matchMode,
                  Bool_t binnedMixing,
                  Int_t dataType,
                  Bool_t megaSkimMode,
                  const std::set<Int_t>& activeFiles,
                  const std::vector<Double_t>& vzBins,
                  const std::vector<Double_t>& hiBinBins,
                  const std::vector<Double_t>& hfBins,
                  const std::vector<Double_t>& multBins);

  // O(1) / O(log N) bucketed partner lookup preserving exact tolerances, widening rules, and non-self checks.
  void FindMatchedDonors(std::vector<Int_t>& matchedGlobalIndices,
                         Int_t nEventsToMatch,
                         Double_t signalVz,
                         Double_t signalActivity,
                         Double_t hfEnergyShift,
                         ULong64_t signalEventNumber,
                         Int_t signalLoopIndex,
                         Int_t debugLevel,
                         TRandom3* rng,
                         Bool_t hasSignalPhysicalId = false,
                         UInt_t signalRun = 0,
                         UInt_t signalLumi = 0) const;

  const std::vector<MixingDonorRecord>& GetRecords() const { return fRecords; }
  const MixingDonorRecord& GetRecord(Int_t globalIndex) const { return fRecords.at(globalIndex); }
  const std::vector<Int_t>& GetActiveRecordIndices() const { return fActiveRecordIndices; }
  std::size_t GetEligibleActiveCount() const { return fEligibleActiveCount; }
  std::size_t GetTotalRecordCount() const { return fRecords.size(); }

  static Int_t ComputeMixingBin(Double_t value, const std::vector<Double_t>& borders);

private:
  static Int_t VzCell(Double_t vz) {
    return static_cast<Int_t>(std::floor(vz / 0.5));
  }

  std::vector<std::string> fCanonicalFiles;
  std::vector<Int_t> fFileEntryCounts;
  std::vector<MixingDonorRecord> fRecords;

  // Active-shard state & 2D spatial/binned lookup structures
  Int_t fMatchMode = kMatchHiBin;
  Bool_t fBinnedMixing = false;
  Int_t fDataType = -1;
  Bool_t fMegaSkimMode = false;
  std::vector<Double_t> fVzBins;
  std::vector<Double_t> fHiBinBins;
  std::vector<Double_t> fHfBins;
  std::vector<Double_t> fMultBins;

  std::vector<Int_t> fActiveRecordIndices;
  std::size_t fEligibleActiveCount = 0;
  Double_t fMinVz = 0.0;
  Double_t fMaxVz = 0.0;
  Double_t fMinActivity = 0.0;
  Double_t fMaxActivity = 0.0;

  // Unbinned (vzCell, hiBin) -> global record indices
  std::map<std::pair<Int_t, Int_t>, std::vector<Int_t>> fUnbinnedHiBinBuckets;
  // Unbinned vzCell -> global record indices sorted ascending by hfSum or genMultiplicity
  std::map<Int_t, std::vector<Int_t>> fUnbinnedSortedActivityBuckets;
  // Binned (vzBin, activityBin) -> global record indices
  std::map<std::pair<Int_t, Int_t>, std::vector<Int_t>> fBinnedBuckets;
  Int_t fMinVzBin = 0;
  Int_t fMaxVzBin = 0;
  Int_t fMinActBin = 0;
  Int_t fMaxActBin = 0;
};

#endif

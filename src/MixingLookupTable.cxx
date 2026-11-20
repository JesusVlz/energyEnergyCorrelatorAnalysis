#include "MixingLookupTable.h"
#include "OOEventIdentity.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

namespace {

std::string TrimWhitespace(const std::string& text) {
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

bool NearlyEqual(double a, double b, double eps = 1e-6) {
  return std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= eps;
}

std::uint64_t Fnv1a64(const std::string& text) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (unsigned char byte : text) {
    hash ^= byte;
    hash *= UINT64_C(1099511628211);
  }
  return hash;
}

void PickRandomDistinct(std::vector<Int_t>& candidates,
                        std::vector<Int_t>& chosen,
                        Int_t targetCount,
                        TRandom3* rng) {
  while (!candidates.empty() && static_cast<Int_t>(chosen.size()) < targetCount) {
    const std::size_t pick = rng ? static_cast<std::size_t>(rng->Integer(candidates.size())) : 0;
    chosen.push_back(candidates[pick]);
    candidates[pick] = candidates.back();
    candidates.pop_back();
  }
}

} // namespace

std::string MixingLookupTable::CanonicalizeMixingPath(const std::string& rawPath) {
  std::string p = TrimWhitespace(rawPath);
  const auto schemePos = p.find("://");
  if (schemePos != std::string::npos) {
    const auto slashAfterHost = p.find('/', schemePos + 3);
    if (slashAfterHost != std::string::npos) {
      p = p.substr(slashAfterHost);
    }
  }
  while (p.size() > 1 && p[0] == '/' && p[1] == '/') {
    p.erase(0, 1);
  }
  const auto storePos = p.find("/store/");
  if (storePos != std::string::npos) {
    return p.substr(storePos);
  }
  const auto eosPos = p.find("/eos/");
  if (eosPos != std::string::npos) {
    return p.substr(eosPos);
  }
  return std::filesystem::path(p).lexically_normal().string();
}

std::vector<std::string> MixingLookupTable::CandidateLookupPaths(const std::string& mixingListPath) {
  std::vector<std::string> candidates;
  if (mixingListPath.empty() || mixingListPath == "none") return candidates;
  candidates.push_back(mixingListPath + ".lookup.tsv");
  if (mixingListPath.size() > 4 && mixingListPath.substr(mixingListPath.size() - 4) == ".txt") {
    candidates.push_back(mixingListPath.substr(0, mixingListPath.size() - 4) + ".lookup.tsv");
  }
  return candidates;
}

std::uint64_t MixingLookupTable::ComputeShardHash(Int_t jobShardId, const std::vector<TString>& signalFiles) {
  std::uint64_t h = Fnv1a64("EEC_MIXING_SHARD_V1");
  if (jobShardId >= 0) {
    h ^= oo_eec::SplitMix64(static_cast<std::uint64_t>(jobShardId) + 0x9e3779b97f4a7c15ULL);
  }
  for (const auto& file : signalFiles) {
    const std::string canonical = CanonicalizeMixingPath(file.Data());
    h ^= oo_eec::SplitMix64(Fnv1a64(canonical));
    h = oo_eec::SplitMix64(h);
  }
  return h;
}

void MixingLookupTable::Clear() {
  fCanonicalFiles.clear();
  fFileEntryCounts.clear();
  fRecords.clear();
  fActiveRecordIndices.clear();
  fEligibleActiveCount = 0;
  fUnbinnedHiBinBuckets.clear();
  fUnbinnedSortedActivityBuckets.clear();
  fBinnedBuckets.clear();
}

void MixingLookupTable::SetFileMetadata(const std::vector<TString>& mixingFiles,
                                       const std::vector<Int_t>& fileEntryCounts) {
  if (mixingFiles.size() != fileEntryCounts.size()) {
    throw std::invalid_argument("MixingLookupTable: file count and entry-count vector sizes differ");
  }
  fCanonicalFiles.clear();
  fCanonicalFiles.reserve(mixingFiles.size());
  for (const auto& path : mixingFiles) {
    fCanonicalFiles.push_back(CanonicalizeMixingPath(path.Data()));
  }
  fFileEntryCounts = fileEntryCounts;
}

void MixingLookupTable::AddRecord(const MixingDonorRecord& record) {
  MixingDonorRecord copy = record;
  copy.globalIndex = static_cast<Int_t>(fRecords.size());
  fRecords.push_back(std::move(copy));
}

bool MixingLookupTable::WriteToFile(const std::string& outputPath, const MixingLookupConfig& config) const {
  if (outputPath.empty()) throw std::invalid_argument("MixingLookupTable: empty output path");
  const std::filesystem::path target(outputPath);
  if (target.has_parent_path()) {
    std::filesystem::create_directories(target.parent_path());
  }
  const std::string tempPath = outputPath + ".tmp." + std::to_string(::getpid());
  {
    std::ofstream out(tempPath, std::ios::out | std::ios::trunc);
    if (!out) {
      throw std::runtime_error("MixingLookupTable: cannot open temporary output file " + tempPath);
    }
    out << "# EEC_MIXING_LOOKUP_V1\n";
    out << "# META DataType=" << config.dataType << "\n";
    out << "# META MegaSkimMixing=" << (config.megaSkimMode ? 1 : 0) << "\n";
    out << std::setprecision(17);
    out << "# META ZVertexCut=" << config.vzCut << "\n";
    out << "# META TriggerSelection=" << config.triggerSelection << "\n";
    out << "# META ApplyJetVeto=" << (config.applyJetVeto ? 1 : 0) << "\n";
    out << "# META JetVetoMapFile=" << config.jetVetoMapFile << "\n";
    out << "# META JetVetoMapHistogram=" << config.jetVetoMapHistogram << "\n";
    out << "# META JetVetoPtConvention=" << config.jetVetoPtConvention << "\n";
    out << "# META Run3JetID=" << config.run3JetID << "\n";
    out << "# META Run3JetIDTable=" << config.run3JetIDTable << "\n";
    out << "# META Run3JEC=" << (config.useRun3JEC ? 1 : 0) << "\n";
    out << "# META Run3UsePpL1FastJet=" << (config.usePpL1FastJet ? 1 : 0) << "\n";
    out << "# META OOClosurePilot=" << (config.ooClosurePilot ? 1 : 0) << "\n";
    out << "# META HasGenMultiplicity=" << (config.requireGenMultiplicity ? 1 : 0) << "\n";
    out << "# META NumFiles=" << fCanonicalFiles.size() << "\n";
    out << "# META NumEvents=" << fRecords.size() << "\n";
    for (std::size_t i = 0; i < fCanonicalFiles.size(); ++i) {
      out << "# FILE " << i << "\t" << fFileEntryCounts.at(i) << "\t" << fCanonicalFiles.at(i) << "\n";
    }
    out << "# COLUMNS fileIdx\tlocalEntry\trun\tlumi\tevt\thasPhysId\tvz\tstoredHiBin\thfSum\trawHfSum\tweight\tgenMult\tpassStdCuts\tpassCuts\tpassMB\tpv\thfCoinc\tpileup\tvzPass\tooNh\tbgKnown\tbgKey\tbgSig\n";
    out << std::setprecision(17);
    for (const auto& r : fRecords) {
      out << r.fileIndex << '\t'
          << r.localEntry << '\t'
          << r.run << '\t'
          << r.lumi << '\t'
          << r.eventNumber << '\t'
          << (r.hasPhysicalId ? 1 : 0) << '\t'
          << r.vz << '\t'
          << r.storedHiBin << '\t'
          << r.hfSum << '\t'
          << r.rawHfSum << '\t'
          << r.eventWeight << '\t'
          << r.genMultiplicity << '\t'
          << (r.passStandardEventCuts ? 1 : 0) << '\t'
          << (r.passEventCuts ? 1 : 0) << '\t'
          << (r.passMBTrigger ? 1 : 0) << '\t'
          << (r.pvBit ? 1 : 0) << '\t'
          << (r.hfCoincBit ? 1 : 0) << '\t'
          << (r.pileupBit ? 1 : 0) << '\t'
          << (r.vzBit ? 1 : 0) << '\t'
          << r.ooNh << '\t'
          << (r.backgroundKnown ? 1 : 0) << '\t'
          << r.backgroundKey << '\t'
          << (r.backgroundSignature.empty() ? "-" : r.backgroundSignature) << '\n';
    }
    out.flush();
    if (!out) {
      std::error_code ec;
      std::filesystem::remove(tempPath, ec);
      throw std::runtime_error("MixingLookupTable: write failed for " + tempPath);
    }
  }
  std::filesystem::rename(tempPath, outputPath);
  return true;
}

bool MixingLookupTable::TryLoadFromFile(const std::string& inputPath,
                                        const MixingLookupConfig& expected,
                                        const std::vector<TString>& mixingFiles,
                                        std::string& rejectReason) {
  Clear();
  rejectReason.clear();
  std::ifstream in(inputPath);
  if (!in) {
    rejectReason = "cannot open file " + inputPath;
    return false;
  }

  std::string line;
  if (!std::getline(in, line) || TrimWhitespace(line) != "# EEC_MIXING_LOOKUP_V1") {
    rejectReason = "missing # EEC_MIXING_LOOKUP_V1 header";
    return false;
  }

  std::map<std::string, std::string> meta;
  std::vector<std::string> tableFiles;
  std::vector<Int_t> tableCounts;
  std::size_t expectedNumEvents = 0;
  bool headerValidated = false;

  auto validateHeader = [&]() -> bool {
    try {
      if (std::stoi(meta.at("DataType")) != expected.dataType) {
        rejectReason = "DataType mismatch";
        return false;
      }
      if ((std::stoi(meta.at("MegaSkimMixing")) != 0) != expected.megaSkimMode) {
        rejectReason = "MegaSkimMixing mismatch";
        return false;
      }
      if (!NearlyEqual(std::stod(meta.at("ZVertexCut")), expected.vzCut)) {
        rejectReason = "ZVertexCut mismatch";
        return false;
      }
      if ((expected.dataType == 9 || expected.ooClosurePilot) &&
          std::stoi(meta.at("TriggerSelection")) != expected.triggerSelection) {
        rejectReason = "TriggerSelection mismatch";
        return false;
      }
      const bool tableVeto = std::stoi(meta.at("ApplyJetVeto")) != 0;
      if (tableVeto != expected.applyJetVeto) {
        rejectReason = "ApplyJetVeto mismatch";
        return false;
      }
      if (expected.applyJetVeto) {
        if (meta["JetVetoMapFile"] != expected.jetVetoMapFile ||
            meta["JetVetoMapHistogram"] != expected.jetVetoMapHistogram ||
            meta["JetVetoPtConvention"] != expected.jetVetoPtConvention ||
            std::stoi(meta.at("Run3JetID")) != expected.run3JetID ||
            (expected.run3JetID != 0 && meta["Run3JetIDTable"] != expected.run3JetIDTable)) {
          rejectReason = "Run-3 jet veto configuration mismatch";
          return false;
        }
        if (expected.jetVetoPtConvention == "recomputed" &&
            ((std::stoi(meta.at("Run3JEC")) != 0) != expected.useRun3JEC ||
             (std::stoi(meta.at("Run3UsePpL1FastJet")) != 0) != expected.usePpL1FastJet)) {
          rejectReason = "Recomputed jet veto JEC/L1 configuration mismatch";
          return false;
        }
      }
      if (expected.ooClosurePilot && std::stoi(meta.at("OOClosurePilot")) == 0) {
        rejectReason = "Lookup table lacks OOClosurePilot metadata";
        return false;
      }
      if (expected.requireGenMultiplicity && std::stoi(meta.at("HasGenMultiplicity")) == 0) {
        rejectReason = "Lookup table lacks generator multiplicity metadata";
        return false;
      }
      const std::size_t numFiles = static_cast<std::size_t>(std::stoull(meta.at("NumFiles")));
      expectedNumEvents = static_cast<std::size_t>(std::stoull(meta.at("NumEvents")));
      if (numFiles != mixingFiles.size() || tableFiles.size() != numFiles) {
        rejectReason = "MixingFileList file count differs from lookup table";
        return false;
      }
      for (std::size_t i = 0; i < numFiles; ++i) {
        const std::string expectedCanonical = CanonicalizeMixingPath(mixingFiles[i].Data());
        if (tableFiles[i] != expectedCanonical) {
          rejectReason = "MixingFileList entry " + std::to_string(i) + " differs (" +
                         tableFiles[i] + " vs " + expectedCanonical + ")";
          return false;
        }
      }
    } catch (const std::exception& err) {
      rejectReason = std::string("malformed lookup header: ") + err.what();
      return false;
    }
    fCanonicalFiles = tableFiles;
    fFileEntryCounts = tableCounts;
    fRecords.reserve(expectedNumEvents);
    headerValidated = true;
    return true;
  };

  while (std::getline(in, line)) {
    const std::string trimmed = TrimWhitespace(line);
    if (trimmed.empty()) continue;
    if (trimmed[0] == '#') {
      if (trimmed.rfind("# META ", 0) == 0) {
        const std::string kv = trimmed.substr(7);
        const auto eq = kv.find('=');
        if (eq != std::string::npos) {
          meta[TrimWhitespace(kv.substr(0, eq))] = TrimWhitespace(kv.substr(eq + 1));
        }
      } else if (trimmed.rfind("# FILE ", 0) == 0) {
        std::istringstream ss(trimmed.substr(7));
        std::size_t idx = 0;
        Int_t count = 0;
        std::string fpath;
        if (!(ss >> idx >> count) || !std::getline(ss, fpath)) {
          rejectReason = "malformed # FILE line in lookup table";
          Clear();
          return false;
        }
        fpath = TrimWhitespace(fpath);
        if (idx != tableFiles.size()) {
          rejectReason = "out-of-order # FILE index in lookup table";
          Clear();
          return false;
        }
        tableFiles.push_back(fpath);
        tableCounts.push_back(count);
      }
      continue;
    }

    if (!headerValidated) {
      if (!validateHeader()) {
        Clear();
        return false;
      }
    }

    std::istringstream row(trimmed);
    MixingDonorRecord rec;
    int hasPhys = 0, passStd = 0, passCuts = 0, passMB = 0, pv = 0, hfCoinc = 0, pileup = 0, vzPass = 0, bgKnown = 0;
    std::string bgSig;
    if (!(row >> rec.fileIndex >> rec.localEntry >> rec.run >> rec.lumi >> rec.eventNumber
              >> hasPhys >> rec.vz >> rec.storedHiBin >> rec.hfSum >> rec.rawHfSum
              >> rec.eventWeight >> rec.genMultiplicity >> passStd >> passCuts >> passMB
              >> pv >> hfCoinc >> pileup >> vzPass
              >> rec.ooNh >> bgKnown >> rec.backgroundKey >> bgSig)) {
      rejectReason = "malformed event row in lookup table";
      Clear();
      return false;
    }
    rec.hasPhysicalId = (hasPhys != 0);
    rec.passStandardEventCuts = (passStd != 0);
    rec.passEventCuts = (passCuts != 0);
    rec.passMBTrigger = (passMB != 0);
    rec.pvBit = (pv != 0);
    rec.hfCoincBit = (hfCoinc != 0);
    rec.pileupBit = (pileup != 0);
    rec.vzBit = (vzPass != 0);
    rec.backgroundKnown = (bgKnown != 0);
    rec.backgroundSignature = (bgSig == "-") ? "" : bgSig;
    AddRecord(rec);
  }

  if (!headerValidated && !validateHeader()) {
    Clear();
    return false;
  }

  if (fRecords.size() != expectedNumEvents) {
    rejectReason = "event count mismatch in lookup table body";
    Clear();
    return false;
  }
  std::size_t sumCounts = 0;
  for (Int_t c : fFileEntryCounts) sumCounts += static_cast<std::size_t>(std::max(0, c));
  if (sumCounts != fRecords.size()) {
    rejectReason = "per-file entry counts do not sum to total events";
    Clear();
    return false;
  }
  return true;
}

std::vector<Int_t> MixingLookupTable::SelectActiveFiles(std::size_t nTotalFiles,
                                                        Int_t maxOpenFiles,
                                                        std::uint64_t shardHash,
                                                        Int_t dataType,
                                                        Bool_t megaSkimMode) const {
  std::vector<Int_t> selected;
  if (nTotalFiles == 0) return selected;
  const std::size_t initialCount =
    (maxOpenFiles <= 0 || static_cast<std::size_t>(maxOpenFiles) >= nTotalFiles)
      ? nTotalFiles
      : static_cast<std::size_t>(maxOpenFiles);
  const std::size_t startFile = static_cast<std::size_t>(shardHash % nTotalFiles);
  selected.reserve(nTotalFiles);
  for (std::size_t offset = 0; offset < initialCount; ++offset) {
    selected.push_back(static_cast<Int_t>((startFile + offset) % nTotalFiles));
  }

  // If records are already loaded (from lookup table), ensure at least 2 eligible donors
  // exist in the selected shard window by expanding across remaining files if needed.
  if (!fRecords.empty() && selected.size() < nTotalFiles && dataType >= 0) {
    std::vector<std::size_t> eligiblePerFile(nTotalFiles, 0);
    for (const auto& rec : fRecords) {
      if (rec.fileIndex >= 0 && static_cast<std::size_t>(rec.fileIndex) < nTotalFiles &&
          rec.IsEligible(dataType, megaSkimMode)) {
        ++eligiblePerFile[static_cast<std::size_t>(rec.fileIndex)];
      }
    }
    std::size_t totalEligible = 0;
    for (Int_t f : selected) totalEligible += eligiblePerFile[static_cast<std::size_t>(f)];
    for (std::size_t offset = initialCount; offset < nTotalFiles && totalEligible < 2; ++offset) {
      const Int_t nextFile = static_cast<Int_t>((startFile + offset) % nTotalFiles);
      selected.push_back(nextFile);
      totalEligible += eligiblePerFile[static_cast<std::size_t>(nextFile)];
    }
  }
  return selected;
}

Int_t MixingLookupTable::ComputeMixingBin(Double_t value, const std::vector<Double_t>& borders) {
  if (borders.size() > 2) {
    Int_t bin = -1;
    for (Double_t border : borders) {
      if (value < border) return bin;
      ++bin;
    }
    return bin;
  }
  const Double_t anchor = borders.at(0);
  const Double_t width = borders.at(1);
  Double_t anchoredValue = value - anchor;
  if (anchoredValue < 0) anchoredValue -= width;
  return static_cast<Int_t>(anchoredValue / width);
}

void MixingLookupTable::BuildIndex(Int_t matchMode,
                                   Bool_t binnedMixing,
                                   Int_t dataType,
                                   Bool_t megaSkimMode,
                                   const std::set<Int_t>& activeFiles,
                                   const std::vector<Double_t>& vzBins,
                                   const std::vector<Double_t>& hiBinBins,
                                   const std::vector<Double_t>& hfBins,
                                   const std::vector<Double_t>& multBins) {
  fMatchMode = matchMode;
  fBinnedMixing = binnedMixing;
  fDataType = dataType;
  fMegaSkimMode = megaSkimMode;
  fVzBins = vzBins;
  fHiBinBins = hiBinBins;
  fHfBins = hfBins;
  fMultBins = multBins;

  fActiveRecordIndices.clear();
  fEligibleActiveCount = 0;
  fUnbinnedHiBinBuckets.clear();
  fUnbinnedSortedActivityBuckets.clear();
  fBinnedBuckets.clear();

  fMinVz = std::numeric_limits<Double_t>::infinity();
  fMaxVz = -std::numeric_limits<Double_t>::infinity();
  fMinActivity = std::numeric_limits<Double_t>::infinity();
  fMaxActivity = -std::numeric_limits<Double_t>::infinity();
  fMinVzBin = std::numeric_limits<Int_t>::max();
  fMaxVzBin = std::numeric_limits<Int_t>::min();
  fMinActBin = std::numeric_limits<Int_t>::max();
  fMaxActBin = std::numeric_limits<Int_t>::min();

  for (const auto& rec : fRecords) {
    if (!activeFiles.empty() && activeFiles.count(rec.fileIndex) == 0) continue;
    fActiveRecordIndices.push_back(rec.globalIndex);
    if (!rec.IsEligible(fDataType, fMegaSkimMode)) continue;

    ++fEligibleActiveCount;
    fMinVz = std::min(fMinVz, rec.vz);
    fMaxVz = std::max(fMaxVz, rec.vz);

    Double_t actVal = 0.0;
    if (fMatchMode == kMatchHiBin) actVal = static_cast<Double_t>(rec.NormalizedHiBin());
    else if (fMatchMode == kMatchHFEnergy) actVal = rec.hfSum;
    else if (fMatchMode == kMatchMultiplicity) actVal = static_cast<Double_t>(rec.genMultiplicity);
    else throw std::invalid_argument("MixingLookupTable: unsupported EventMatchForMixing " + std::to_string(fMatchMode));

    fMinActivity = std::min(fMinActivity, actVal);
    fMaxActivity = std::max(fMaxActivity, actVal);

    if (fBinnedMixing) {
      const Int_t vzBin = ComputeMixingBin(rec.vz, fVzBins);
      Int_t actBin = 0;
      if (fMatchMode == kMatchHiBin) actBin = ComputeMixingBin(actVal, fHiBinBins);
      else if (fMatchMode == kMatchHFEnergy) actBin = ComputeMixingBin(actVal, fHfBins);
      else actBin = ComputeMixingBin(actVal, fMultBins);

      fMinVzBin = std::min(fMinVzBin, vzBin);
      fMaxVzBin = std::max(fMaxVzBin, vzBin);
      fMinActBin = std::min(fMinActBin, actBin);
      fMaxActBin = std::max(fMaxActBin, actBin);
      fBinnedBuckets[{vzBin, actBin}].push_back(rec.globalIndex);
    } else {
      const Int_t vc = VzCell(rec.vz);
      if (fMatchMode == kMatchHiBin) {
        fUnbinnedHiBinBuckets[{vc, rec.NormalizedHiBin()}].push_back(rec.globalIndex);
      } else {
        fUnbinnedSortedActivityBuckets[vc].push_back(rec.globalIndex);
      }
    }
  }

  const auto byCanonicalOrder = [&](Int_t a, Int_t b) {
    if (fRecords[a].fileIndex != fRecords[b].fileIndex) {
      return fRecords[a].fileIndex < fRecords[b].fileIndex;
    }
    return fRecords[a].localEntry < fRecords[b].localEntry;
  };
  std::sort(fActiveRecordIndices.begin(), fActiveRecordIndices.end(), byCanonicalOrder);

  if (fBinnedMixing) {
    for (auto& kv : fBinnedBuckets) {
      std::sort(kv.second.begin(), kv.second.end(), byCanonicalOrder);
    }
  } else if (fMatchMode == kMatchHiBin) {
    for (auto& kv : fUnbinnedHiBinBuckets) {
      std::sort(kv.second.begin(), kv.second.end(), byCanonicalOrder);
    }
  } else {
    for (auto& kv : fUnbinnedSortedActivityBuckets) {
      auto& vec = kv.second;
      std::sort(vec.begin(), vec.end(), [&](Int_t a, Int_t b) {
        const Double_t va = (fMatchMode == kMatchHFEnergy)
          ? fRecords[a].hfSum
          : static_cast<Double_t>(fRecords[a].genMultiplicity);
        const Double_t vb = (fMatchMode == kMatchHFEnergy)
          ? fRecords[b].hfSum
          : static_cast<Double_t>(fRecords[b].genMultiplicity);
        if (va != vb) return va < vb;
        return byCanonicalOrder(a, b);
      });
    }
  }
}

void MixingLookupTable::FindMatchedDonors(std::vector<Int_t>& matchedGlobalIndices,
                                          Int_t nEventsToMatch,
                                          Double_t signalVz,
                                          Double_t signalActivity,
                                          Double_t hfEnergyShift,
                                          ULong64_t signalEventNumber,
                                          Int_t signalLoopIndex,
                                          Int_t debugLevel,
                                          TRandom3* rng,
                                          Bool_t hasSignalPhysicalId,
                                          UInt_t signalRun,
                                          UInt_t signalLumi) const {
  matchedGlobalIndices.clear();
  if (nEventsToMatch <= 0) return;
  if (static_cast<Int_t>(fEligibleActiveCount) < nEventsToMatch) {
    throw std::runtime_error("MixingLookupTable: fewer eligible donor events than required matches");
  }

  const auto alreadyPicked = [&](Int_t idx) {
    return std::find(matchedGlobalIndices.begin(), matchedGlobalIndices.end(), idx) != matchedGlobalIndices.end();
  };
  const auto isSamePhysicalEvent = [&](const MixingDonorRecord& rec) {
    if (hasSignalPhysicalId && rec.hasPhysicalId) {
      return rec.run == signalRun && rec.lumi == signalLumi && rec.eventNumber == signalEventNumber;
    }
    return rec.eventNumber == signalEventNumber;
  };

  if (fBinnedMixing) {
    const Int_t vzBin = ComputeMixingBin(signalVz, fVzBins);
    Int_t actBin = 0;
    if (fMatchMode == kMatchHiBin) actBin = ComputeMixingBin(signalActivity, fHiBinBins);
    else if (fMatchMode == kMatchHFEnergy) actBin = ComputeMixingBin(signalActivity, fHfBins);
    else actBin = ComputeMixingBin(signalActivity, fMultBins);

    const Int_t maxVzDist = std::max(std::abs(fMaxVzBin - vzBin), std::abs(vzBin - fMinVzBin));
    const Int_t maxActDist = std::max(std::abs(fMaxActBin - actBin), std::abs(actBin - fMinActBin));
    const Int_t maxStep = std::max(maxVzDist, maxActDist);

    std::vector<Int_t> stepCandidates;
    for (Int_t step = 0; step <= maxStep && static_cast<Int_t>(matchedGlobalIndices.size()) < nEventsToMatch; ++step) {
      if (step > 0 && debugLevel > 0) {
        std::cout << "Could not find matching mixed events for event " << signalLoopIndex
                  << " at bin step " << (step - 1) << "; widening bin tolerance to " << step << std::endl;
      }
      stepCandidates.clear();
      auto it = fBinnedBuckets.lower_bound({vzBin - step, std::numeric_limits<Int_t>::min()});
      for (; it != fBinnedBuckets.end() && it->first.first <= vzBin + step; ++it) {
        const Int_t dv = std::abs(it->first.first - vzBin);
        const Int_t da = std::abs(it->first.second - actBin);
        if (da > step || std::max(dv, da) != step) continue;
        for (Int_t idx : it->second) {
          if (isSamePhysicalEvent(fRecords[idx])) continue;
          if (alreadyPicked(idx)) continue;
          stepCandidates.push_back(idx);
        }
      }
      PickRandomDistinct(stepCandidates, matchedGlobalIndices, nEventsToMatch, rng);
    }

    if (static_cast<Int_t>(matchedGlobalIndices.size()) < nEventsToMatch) {
      throw std::runtime_error("MixingLookupTable: fewer than required non-self eligible donor events in binned pool");
    }
    return;
  }

  // Unbinned 2D lookup
  Double_t targetActivity = signalActivity;
  if (fMatchMode == kMatchHFEnergy) {
    if (targetActivity < 2.0 * hfEnergyShift) {
      targetActivity /= 2.0;
    } else {
      targetActivity -= hfEnergyShift;
    }
    if (targetActivity < 2.0) targetActivity = 2.0;
  }

  const Double_t maxVzSpan = std::max(std::abs(fMaxVz - signalVz), std::abs(signalVz - fMinVz));
  const Double_t maxActSpan = std::max(std::abs(fMaxActivity - targetActivity), std::abs(targetActivity - fMinActivity));

  std::vector<Int_t> stepCandidates;
  for (Int_t step = 0; static_cast<Int_t>(matchedGlobalIndices.size()) < nEventsToMatch; ++step) {
    if (step > 0 && debugLevel > 0) {
      std::cout << "Could not find matching mixed events for event " << signalLoopIndex
                << "; increasing unbinned tolerance to step " << step << std::endl;
    }
    const Double_t vzTolerance = 0.5 * (step + 1);
    const Double_t prevVzTolerance = (step == 0) ? -1.0 : (0.5 * step);

    const Int_t minVzCell = VzCell(signalVz - vzTolerance - 1e-4);
    const Int_t maxVzCell = VzCell(signalVz + vzTolerance + 1e-4);

    stepCandidates.clear();

    if (fMatchMode == kMatchHiBin) {
      const Int_t signalHiBin = static_cast<Int_t>(std::llround(targetActivity));
      const Int_t hiBinTolerance = step;
      const Int_t prevHiBinTolerance = step - 1;
      for (Int_t vc = minVzCell; vc <= maxVzCell; ++vc) {
        for (Int_t hb = signalHiBin - hiBinTolerance; hb <= signalHiBin + hiBinTolerance; ++hb) {
          auto it = fUnbinnedHiBinBuckets.find({vc, hb});
          if (it == fUnbinnedHiBinBuckets.end()) continue;
          for (Int_t idx : it->second) {
            const auto& rec = fRecords[idx];
            const Double_t dz = std::abs(rec.vz - signalVz);
            const Int_t dhb = std::abs(rec.NormalizedHiBin() - signalHiBin);
            if (dz > vzTolerance + 1e-4 || dhb > hiBinTolerance) continue;
            if (step > 0 && dz <= prevVzTolerance + 1e-4 && dhb <= prevHiBinTolerance) continue;
            if (isSamePhysicalEvent(rec) || alreadyPicked(idx)) continue;
            stepCandidates.push_back(idx);
          }
        }
      }
      PickRandomDistinct(stepCandidates, matchedGlobalIndices, nEventsToMatch, rng);
      if (vzTolerance >= maxVzSpan && static_cast<Double_t>(hiBinTolerance) >= maxActSpan) break;
    } else {
      const bool highRegime = (fMatchMode == kMatchHFEnergy) ? (targetActivity > 100.0)
                                                             : (std::llround(targetActivity) > 100);
      const Double_t actTolerance = highRegime ? (targetActivity * 0.01 * (step + 1)) : (1.0 * (step + 1));
      const Double_t prevActTolerance = (step == 0) ? -1.0
        : (highRegime ? (targetActivity * 0.01 * step) : (1.0 * step));

      for (Int_t vc = minVzCell; vc <= maxVzCell; ++vc) {
        auto it = fUnbinnedSortedActivityBuckets.find(vc);
        if (it == fUnbinnedSortedActivityBuckets.end()) continue;
        const auto& vec = it->second;
        const Double_t lowBound = targetActivity - actTolerance - 1e-9;
        const Double_t highBound = targetActivity + actTolerance + 1e-9;

        auto pos = std::lower_bound(vec.begin(), vec.end(), lowBound, [&](Int_t idx, Double_t bound) {
          const Double_t val = (fMatchMode == kMatchHFEnergy)
            ? fRecords[idx].hfSum
            : static_cast<Double_t>(fRecords[idx].genMultiplicity);
          return val < bound;
        });
        for (; pos != vec.end(); ++pos) {
          const auto& rec = fRecords[*pos];
          const Double_t val = (fMatchMode == kMatchHFEnergy)
            ? rec.hfSum
            : static_cast<Double_t>(rec.genMultiplicity);
          if (val > highBound) break;
          const Double_t dz = std::abs(rec.vz - signalVz);
          const Double_t da = std::abs(val - targetActivity);
          if (dz > vzTolerance + 1e-4 || da > actTolerance) continue;
          if (step > 0 && dz <= prevVzTolerance + 1e-4 && da <= prevActTolerance) continue;
          if (isSamePhysicalEvent(rec) || alreadyPicked(*pos)) continue;
          stepCandidates.push_back(*pos);
        }
      }
      PickRandomDistinct(stepCandidates, matchedGlobalIndices, nEventsToMatch, rng);
      if (vzTolerance >= maxVzSpan && actTolerance >= maxActSpan) break;
    }
  }

  if (static_cast<Int_t>(matchedGlobalIndices.size()) < nEventsToMatch) {
    throw std::runtime_error("MixingLookupTable: fewer than required non-self eligible donor events in pool");
  }
}

#ifndef TRACKPAIREFFICIENCYCORRECTOR_H
#define TRACKPAIREFFICIENCYCORRECTOR_H

// Root includes
#include <TH1.h>
#include <TFile.h>

// Own Includes
#include "BinningCard.h"

/*
 * TrackPairEfficiencyCorrector class
 *
 * Class whose purpose is to provide weight for reflected cone histograms to take into account jet reconstructions effects to background particles
 */
class TrackPairEfficiencyCorrector {
  
public:
 
  // Maximum dimensions for the arrays
  static const int kMaxCentralityBins = 5;      // Maximum allowed number of centrality bins
  static const int kMaxTrackPtBins = 12;        // Maximum allowed number of track pT bins
  static const int kMaxJetPtBins = 20;          // Maximum allowed number of jet pT bins
  
  TrackPairEfficiencyCorrector();                   // Contructor
  TrackPairEfficiencyCorrector(TString inputFileName, bool useSmoothedCorrection, bool isPbPbData);   // Custom constructor
  ~TrackPairEfficiencyCorrector() = default;        // Destructor
  
  std::pair<double,double> GetTrackPairEfficiencyCorrection(const double deltaR, const double centrality, const double triggerPt, const double associatedPt, const double jetPt = -1) const;

  void ReadCorrectionTables(); // Initialize the correction tables from file
  
  bool GetDisableCorrection(); // Getter for flag to disable the track pair efficiency correction
  void SetDisableCorrection(const bool disable); // Setter for disabling the correction
  
private:
  
  // Card with binning information
  BinningCard* fCard;
  
  // File from which the histograms are read
  TFile* fInputFile;
  
  // Binning information read to be read from the input file
  double fCentralityBinBorders[kMaxCentralityBins+1] = {0};
  double fTrackPtBinBorders[kMaxTrackPtBins+1] = {0};
  double fJetPtBinBorders[kMaxJetPtBins+1] = {0};
  int fnCentralityBins;
  int fnTrackPtBins;
  int fnJetPtBins;
  
  // Flags for selecting correct corrections
  bool fDisableCorrection;      // Disable the correction and just return 1
  bool fUseSmoothedCorrection;  // True: Use smoothed correction. False: Use raw correction
  bool fIsPbPbData;             // PbPb and pp have different fluctuation reduction settings. Thus we need a flag for data type
  
  // Histograms from which the corrections are read
  TH1D* fCorrectionTable[kMaxCentralityBins][kMaxTrackPtBins][kMaxTrackPtBins];
  TH1D* fCorrectionTableCloseToJet[kMaxCentralityBins][kMaxTrackPtBins][kMaxTrackPtBins][kMaxJetPtBins];
  
  // Methods
  int FindBinIndex(const double* array, const int nBins, const double value) const; // Find a bin index from a given array
  int FindCentralityBin(const double centrality) const; // Find centrality bin index
  int FindTrackPtBin(const double trackPt) const; // Find track pT bin index
  int FindJetPtBin(const double jetPt) const; // Find jet pT bin index
  
};

#endif
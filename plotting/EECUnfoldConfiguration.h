#ifndef EECUNFOLDCONFIGURATION_H
#define EECUNFOLDCONFIGURATION_H

// Own includes
#include "EECCard.h"

/*
 * EECUnfoldConfiguration class
 *
 * Class whose purpose is to provide unfolding configuration for the data
 */
class EECUnfoldConfiguration {
  
public:
 
  // Helper enumeration for different defined Monte Carlo splits and systematic parameter sets
  enum enumUnfoldingParameterSet{kNominalSmear, kJetPtResolutionUncertaintyDown, kJetPtResolutionUncertaintyUp, kJetEnergyScaleUncertaintyDown, kJetEnergyScaleUncertaintyUp, kJetPtPriorUncertainty, kCentralityShiftDown, kCentralityShiftUp, kNumberOfIterationsDown, kNumberOfIterationsUp, kMCStatistics, kNParameterSets};
  enum enumMonteCarloSplit{kWholeDataset, kMonteCarloSplit1, kMonteCarloSplit2, kNDatasetSplits};

  // Dimensions for the arrays are defined by the files used the obtain the iteration numbers
  static const int kNCentralityBins = 4;      // Number of centrality bins for which background scale is determined
  static const int kNTrackPtBins = 8;         // Number of track pT bins for which the background scale is determined
  
  EECUnfoldConfiguration();               // Constructor
  EECUnfoldConfiguration(EECCard* card, const int iSplit = 0, const int iSystematic = 0, const int iWeightExponent = 1);  // Custom constructor with a card
  ~EECUnfoldConfiguration();              // Destructor

  // Getter for unfolding configuration
  double GetNumberOfIterations(const std::pair<double,double> centralityBinBorders, double trackPtBorderLow) const; // Getter for the number of iterations
  TString GetResponseFileName() const; // Getter for the name of the file containing the response matrices
  
private:
  
  // Private variables
  bool fIsPbPbData;                // Flag for PbPb data, read from card
  int fSystematicIndex;            // Index for systematic uncertainty study
  int fSplitIndex;                 // Index for the MC split. 0 = Whole dataset. 1 = Split 1, 2 = Split 2
  int fWeightExponent;             // Exponent used for the energy weights for energy-energy correlators
  TString fResponseMatrixFileName; // Response matrix file name for which the configuration is determined

  // Binning information for the scaling tables
  const double fCentralityBinBorderLow[kNCentralityBins] = {0, 10, 30, 50};
  const double fCentralityBinBorderHigh[kNCentralityBins] = {10, 30, 50, 90};
  const double fCentralityBinBorderLowShifted[kNCentralityBins] = {4, 14, 34, 54};
  const double fCentralityBinBorderHighShifted[kNCentralityBins] = {14, 34, 54, 94};
  const double fTrackPtBinBorderLow[kNTrackPtBins] = {0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4};

  // Array that holds the best number of iterations for each response matrix
  double fBestNumberOfIterations[kNCentralityBins][kNTrackPtBins];
  
  // Initialize the number of iterations array based on the predefined configuration index
  void InitializeArrays();
  
};

#endif

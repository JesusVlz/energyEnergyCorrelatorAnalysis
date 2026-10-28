#ifndef EECCARD_H
#define EECCARD_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TObjString.h>
#include <TVectorT.h>

/*
 * EECCard class
 *
 * This class reads the ConfigurationCard from the input root file and decodes the
 * necessary information for the analysis.
 */
class EECCard {
  
public:
 
  // Indices for card entries
  enum enumCardEntries{
    kDataType,                    // Data type in the data file (pp, PbPb, pp MC, PbPb MC)
    kMcCorrelationType,           // Monte Carlo correlation type (RecoReco, RecoGen, GenReco, GenGen)
    kMatchJets,                   // 0 = Reco and Gen jets are not matched, 1 = They are matched
    kTriggerSelection,            // 0 = Do not use any triggers, 1 = Require CaloJet80, 2 = Require CaloJet100, 3 = Require CaloJet80 or CaloJet100
    kJetType,                     // 0 = Calorimeter jets, 1 = PF jets
    kJetAxis,                     // 0 = E-scheme axis, 1 = WTA axis
    kJetEtaCut,                   // Eta cut for jets
    kMinPtCut,                    // Minimum allowed pT for the inclusive jets
    kMaxPtCut,                    // Maximum allowed pT for the inclusive jets
    kCutBadPhiRegion,             // Cut the phi region with bad tracked performance from the analysis
    kMinMaxTrackPtFraction,       // Minimum fraction of jet pT taken by the highest pT track in jet
    kMaxMaxTrackPtFraction,       // Maximum fraction of jet pT taken by the highest pT track in jet
    kJetUncertaintyMode,          // 0 = Nominal jet pT, 1 = pT minus uncertainty, 2 = pT plus uncertainty
    kTrackEtaCut,                 // Eta cut for tracks
    kMinTrackPtCut,               // Minimum accepted track pT
    kMaxTrackPtCut,               // Maximum accepted track pT
    kMaxTrackPtRelativeError,     // Maximum relative error allowed for track pT
    kVertexMaxDistance,           // Maximum allowed distance of tracks from reconstructed vertex
    kCalorimeterSignalLimitPt,    // Limit for track pT above which a signal in calorimeters is required
    kHighPtEtFraction,            // Minimum fraction between pT and Et for high pT tracks
    kChi2QualityCut,              // Maximum accepted chi2 for reconstructed tracks
    kMinimumTrackHits,            // Minimum number of hits in tracking for a track
    kSubeventCut,                 // 0 = Subevent 0 (Pythia), 1 = Subevent > 0 (Hydjet), 2 = No subevent selection
    kTrackEfficiencyVariation,    // Amount of variation applied to track efficiency corrections to evaluate systematic uncertainties
    kJetPtWeightConfiguration,    // In MC: 0 = Do not weight jet pT. 1 = Weight jet pT. 2 = Weight only true jet pT for response matrices
    kDisableTrackPairEfficiencyCorrection, // 0 = Keep track pair efficiency correction enabled. 1 = Disable track pair efficiency correction
    kZVertexCut,                  // Maximum accepted vz in the event
    kLowPtHatCut,                 // Minimum accepted pT hat
    kHighPtHatCut,                // Maximum accepted pT hat
    kMultiplicityMode,            // Select between multiplicity weighting and centrality weighting
    kJetRadius,                   // Radius around the jet axis which is used to select the tracks for energy-energy correlation analysis
    kWeightExponent,              // Exponent given to the energy-energy correlator weight
    kJetPtBinEdgesEEC,            // Jet pT binning for energy-energy correlation analysis
    kTrackPtBinEdgesEEC,          // Track pT binning for energy-energy correlation analysis
    kSkipCovarianceMatrix,        // 0 = Fill covariance matrix, 1 = Do not fill covariance matrix
    kMinJetPtUnfoldingReco,       // Minimum jet pT for reconstructed jets in unfolding response matrix
    kMinJetPtUnfoldingTruth,      // Minimum jet pT for generator level jets in unfolding response matrix
    kCentralityBinEdges,          // Centrality bin edges
    kTrackPtBinEdges,             // Track pT bin edges
    kPtHatBinEdges,               // pT hat bin edges
    kBackgroundMethods,           // Bit 0 = Estimate reflected cone, Bit 1 = Estimate perpendicular cone, Bit 2 = Estimate mixed cone
    kMegaSkimMixing,              // 0 = Regular mixing files, 1 = Megaskimmed micing files 
    kAllowJetsInReflectedCone,    // 0 = Do not allow jets in reflected cone, 1 = Allow jets in reflected cone
    kFirstUnfoldedCentralityBin,  // Index of the first centrality bin that has been unfolded
    kLastUnfoldedCentralityBin,   // Index of the last centrality bin that has been unfolded
    kFirstUnfoldedTrackPtBin,     // Index of the first track pT bin that has been unfolded
    kLastUnfoldedTrackPtBin,      // Index of the last track pT bin that has been unfolded
    kFirstUnfoldedJetPtBin,       // Index of the first jet pT bin that has been unfolded
    kLastUnfoldedJetPtBin,        // Index of the last jet pT bin that has been unfolded
    kBackgroundMethod,            // Index for the method used to perform background subtraction
    kBackgroundSystematic,        // Index for the systematic uncertainty estimation for background subtraction
    knEntries};                   // Number of entries in the card
  
  // Enumeration for input files used in postprocessing
  enum enumFileNames{kInputFileName, kResponseMatrixFile, knFileNames};
  
private:
  
  // Names for each entry read from the configuration card
  const char* fCardEntryNames[knEntries] = {"DataType","McCorrelationType","MatchJets","TriggerSelection","JetType","JetAxis","JetEtaCut","MinJetPtCut","MaxJetPtCut","CutBadPhi","MinMaxTrackPtFraction","MaxMaxTrackPtFraction","JetUncertainty","TrackEtaCut","MinTrackPtCut","MaxTrackPtCut","MaxTrackPtRelativeError","VertexMaxDistance","CalorimeterSignalLimitPt","HighPtEtFraction","Chi2QualityCut","MinimumTrackHits","SubeventCut","TrackEfficiencyVariation","JetPtWeight","DisableTrackPairEfficiencyCorrection","ZVertexCut","LowPtHatCut","HighPtHatCut","MultiplicityMode","JetRadius","WeightExponent","JetPtBinEdgesEEC","TrackPtBinEdgesEEC","SkipCovarianceMatrix","MinJetPtUnfoldingReco","MinJetPtUnfoldingTruth","CentralityBinEdges","TrackPtBinEdges","PtHatBinEdges","BackgroundMethods","MegaSkimMixing","AllowJetsInReflectedCone","FirstUnfoldedCentralityBin","LastUnfoldedCentralityBin","FirstUnfoldedTrackPtBin","LastUnfoldedTrackPtBin","FirstUnfoldedJetPtBin","LastUnfoldedJetPtBin","BackgroundSubtractionMethod","BackgroundSubtractionSystematic"};
  const char* fFileNameType[knFileNames] = {"input", "response matrix"};
  const char* fFileNameSaveName[knFileNames] = {"InputFile", "ResponseMatrixFile"};
  
  TFile* fInputFile;          // Input file from which all the data is read
  TString fCardDirectory;     // Path to the ConfigurationCard directory
  int fDataType;              // Total number of centrality bins in the analysis
  int fMonteCarloType;        // Type of Monte Carlo used for jet-track correlations
  TString fDataTypeString;    // Total number of eta gaps in the analysis
  TString fAlternativeDataTypeString; // Alternative data type string
  TString fDataTypeStringWithoutMCType; // Data type string without the MC type appended to it
  bool fBackgroundLegacyMode; // Legacy naming for background subtraction methods
  
  void FindDataTypeString(); // Construct a data type string based on information on the card
  void ReadVectors();        // Read the vectors from the file
  
  // Strings for git hash
  TObjString* fGitHash;
  TObjString* fProjectionGitHash;
  TObjString* fProcessGitHash;
  TObjString* fUnfoldingGitHash;
  
  // Vectors for all the lines inside the card
  TVectorT<float>* fCardEntries[knEntries];   // Array of all the vectors in the card
  TObjString* fFileNames[knFileNames];        // Array for filenames used in postprocessing
  
  // Private methods
  int GetNBins(const int index) const;                                // Get the number of bins for internal index
  double GetLowBinBorder(const int index, const int iBin) const;      // Get the low border of i:th bin from internal index
  double GetHighBinBorder(const int index, const int iBin) const;     // Get the high border of i:th bin from internal index
  int GetBinIndex(const int index, const double value) const;         // Get the bin index in the i:th bin from internal index based on given value
  int FindBinIndex(const int index, const double lowBinBorder, const double highBinBorder) const; // Get the bin index in the i:th bin from internal index based on the provided bin borders
   
public:
  
  EECCard(TFile* inFile); // Contructor with input file
  ~EECCard();             // Destructor
  
  TString GetDataType() const;            // Getter for data type string
  TString GetAlternativeDataType(const bool includeMCtype = true) const; // Getter for alternative data type string
  void Write(TDirectory* file);            // Write the contents of the card to a file
  void WriteProcessHash(TDirectory* file); // Write the git hash used for processing histograms to the file
  void WriteUnfoldInfo(TDirectory* file);  // Write the information about parameters used in unfolding energy-energy correlators
  void WriteBackground(TDirectory* file);  // Write the information related to background subtraction
  void Print() const;                      // Print the contents of the card to the console
  
  int GetNCentralityBins() const; // Get the number of centrality bins
  int GetNTrackPtBins() const;    // Get the number of track pT bins
  int GetNJetPtBinsEEC() const;   // Get the number of jet pT bins in energy-energy correlator analysis
  int GetNTrackPtBinsEEC() const; // Get the number of track pT bins in energy-energy correlator analysis
  double GetLowBinBorderCentrality(const int iBin) const;  // Get the low border of i:th centrality bin
  double GetLowBinBorderTrackPt(const int iBin) const;     // Get the low border of i:th track pT bin
  double GetLowBinBorderJetPtEEC(const int iBin) const;    // Get the low border of i:th jet pT bin in energy-energy correlator analysis
  double GetLowBinBorderTrackPtEEC(const int iBin) const;  // Get the low border of i:th track pT bin in energy-energy correlator analysis
  double GetHighBinBorderCentrality(const int iBin) const; // Get the high border of i:th centrality bin
  double GetHighBinBorderTrackPt(const int iBin) const;    // Get the high border of i:th track pT bin
  double GetHighBinBorderJetPtEEC(const int iBin) const;   // Get the high border of i:th jet pT bin in energy-energy correlator analysis
  double GetHighBinBorderTrackPtEEC(const int iBin) const; // Get the high border of i:th track pT bin in energy-energy correlator analysis
  std::pair<double,double> GetBinBordersCentrality(const int iBin) const; // Get the bin borders of the i:th centrality bin
  std::pair<double,double> GetBinBordersTrackPt(const int iBin) const; // Get the bin borders of the i:th track pT bin
  std::pair<double,double> GetBinBordersJetPtEEC(const int iBin) const; // Get the bin borders of the i:th jet pT bin in energy-energy correlator analysis
  std::pair<double,double> GetBinBordersTrackPtEEC(const int iBin) const; // Get the bin borders of the i:th track pT bin in energy-energy correlator analysis
  int GetBinIndexCentrality(const double value) const;     // Get the bin index for a given centrality value
  int GetBinIndexTrackPt(const double value) const;        // Get the bin index for a given track pT value
  int GetBinIndexJetPtEEC(const double value) const;       // Get the bin index for a given jet pT value in energy-energy correlator analysis
  int GetBinIndexTrackPtEEC(const double value) const;     // Get the bin index for a given track pT value in energy-energy correlator analysis
  int FindBinIndexCentrality(const double lowBorder, const double highBorder) const; // Find if a centrality bin with given borders exists and return its index
  int FindBinIndexCentrality(const std::pair<double,double> binBorders) const; // Find if a centrality bin with given borders exists and return its index
  int FindBinIndexTrackPt(const double lowBorder, const double highBorder) const;    // Find if a track pT bin with given borders exists and return its index
  int FindBinIndexTrackPt(const std::pair<double,double> binBorders) const;    // Find if a track pT bin with given borders exists and return its index
  int FindBinIndexJetPtEEC(const double lowBorder, const double highBorder) const;   // Find if a jet pT bin in energy-energy correlator analysis with given borders exists and return its index
  int FindBinIndexJetPtEEC(const std::pair<double,double> binBorders) const;   // Find if a jet pT bin in energy-energy correlator analysis with given borders exists and return its index
  int FindBinIndexTrackPtEEC(const double lowBorder, const double highBorder) const; // Find if a track pT bin in energy-energy correlator analysis with given borders exists and return its index
  int FindBinIndexTrackPtEEC(const std::pair<double,double> binBorders) const; // Find if a track pT bin in energy-energy correlator analysis with given borders exists and return its index
  int GetFirstUnfoldedCentralityBin() const; // Getter for the first unfolded centrality bin index
  int GetLastUnfoldedCentralityBin() const;  // Getter for the last unfolded centrality bin index
  int GetFirstUnfoldedTrackPtBin() const;    // Getter for the first unfolded track pT bin index
  int GetLastUnfoldedTrackPtBin() const;     // Getter for the last unfolded track pT bin index
  int GetFirstUnfoldedJetPtBin() const;      // Getter for the first unfolded jet pT bin index
  int GetLastUnfoldedJetPtBin() const;       // Getter for the last unfolded jet pT bin index
  int GetSubeventCut() const;       // Get the index for used subevent cut
  int GetJetType() const;           // Get the jet type index
  double GetJetPtCut() const;       // Get the minimum jet pT cut
  double GetJetRadius() const;   // Getter for the jet radius R
  int GetBackgroundMethods() const; // Get infprmation which backgorund methods are included in the file
  int GetWeightExponent(int index = 1) const; // Get the weight exponent used in energy-energy correlators
  int FindWeightExponentIndex(double weightExponent) const; // Find the index in card for the input weight exponent
  int GetNWeightExponents() const;  // Get the number of weight exponents that are defined in the file
  
  void AddOneDimensionalVector(int entryIndex, float entryContent); // Add one dimensional vector to the card
  void AddVector(int entryIndex, int dimension, double* contents); // Add a vector to the card
  void AddFileName(int entryIndex, TString fileName); // Add a file name to the card
  void AddProjectionGitHash(const char* gitHash); // Add a git hash used to project the histograms to the file
  void AddProcessGitHash(const char* gitHash); // Add a git hash used to process the histograms in the file
  void AddUnfoldingGitHash(const char* gitHash); // Add a git hash used to do jet pT unfolding for the energy-energy correlators
  
};

#endif

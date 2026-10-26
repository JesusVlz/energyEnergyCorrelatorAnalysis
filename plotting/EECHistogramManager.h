#ifndef EECHISTOGRAMMANAGER_H
#define EECHISTOGRAMMANAGER_H

// Deprecation message
#include <iostream>
#include <bitset>

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <THnBase.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>

// Own includes
#include "EECCard.h"
#include "../src/EECHistograms.h"
#include "EECBackgroundScale.h"
#include "EECSignalToBackgroundUnfoldingScale.h"
#include "AlgorithmLibrary.h"

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class EECHistogramManager {
 
public:
  
  // Indices for different track histogram categories
  enum enumTrackHistograms{kTrack, kUncorrectedTrack, knTrackCategories};
  
  // Indices for different multiplicity within the jet cone histogram categories
  enum enumMultiplicityInJetCone{kMultiplicityInJetCone, kMultiplicityInReflectedCone, kMultiplicityInJetConeUncorrected, kMultiplicityInReflectedConeUncorrected, knMultiplicityInJetConeTypes};
  
  // Indices for different particle density types measured around the jet axis
  enum enumParticleDensityAroundJets{kParticleDensityAroundJetAxis, kParticlePtDensityAroundJetAxis, kParticleDensityAroundJetAxisPtBinned, kParticlePtDensityAroundJetAxisPtBinned, knParticleDensityAroundJetAxisTypes};
  
  // Indices for maximum particle pT within the jet cone types
  enum enumMaxParticlePtWithinJetConeType{kMaxSignalParticlePt, kMaxBackgroundParticlePt, knMaxParticlePtWithinJetConeTypes};
  
  // Indices for different energy-energy correlator categories
  enum enumEnergyEnergyCorrelators{kEnergyEnergyCorrelator, kEnergyEnergyCorrelatorEfficiencyVariationPlus, kEnergyEnergyCorrelatorEfficiencyVariationMinus, kEnergyEnergyCorrelatorPairEfficiencyVariationPlus, kEnergyEnergyCorrelatorPairEfficiencyVariationMinus, knEnergyEnergyCorrelatorTypes};
  
  // Indices for different energy-energy correlator processing levels
  enum enumEnergyEnergyCorrelatorProcessing{kEnergyEnergyCorrelatorNormalized, kEnergyEnergyCorrelatorBackground, kEnergyEnergyCorrelatorSignal, kEnergyEnergyCorrelatorUnfolded, kEnergyEnergyCorrelatorBackgroundAfterUnfolding, kEnergyEnergyCorrelatorUnfoldedSignal, knEnergyEnergyCorrelatorProcessingLevels};

  // Indices for unfolding distributions
  enum enumUnfoldingDistribution{kUnfoldingMeasured, kUnfoldingTruth, knUnfoldingDistributionTypes};

  // Indices for covariance matrices
  enum enumCovarianceMatrix{kCovarianceMatrixMeasured, kCovarianceMatrixUnfolded, knCovarianceMatrixTypes};

  // Indices for track/particle matching histograms
  enum enumTrackParticleMatchingQA{kNumberOfParticlesCloseToTrack, kHasMatchingParticle, knTrackParticleMatchingQAHistograms};
  enum enumTrackParticleMatchingResponse{kTrackParticleMatchingDeltaRRresponse, kTrackParticleMatchingPtResponse, knTrackParticleMatchingResponseTypes};
  
  // Dimensions for histogram arrays
  static const int kMaxCentralityBins = 5;       // Maximum allowed number of centrality bins
  static const int kMaxTrackPtBins = 10;         // Maximum allowed number of track pT bins
  static const int knGenJetPtBins = 45;          // Number of generator level jet pT bins for jet pT closures
  static const int knJetEtaBins = 50;            // Number of jet eta bins for jet pT closures
  static const int knJetPhiBins = 64;            // Number of jet phi bins for jet pT closures
  static const int kMaxJetPtBinsEEC = 60;       // Maximum allowed number of jet pT bins for energy-energy correlators
  static const int kMaxTrackPtBinsEEC = 20;     // Maximum allowed number of track pT bins for energy-energy correlators
  static const int knProjectedMaxParticlePtBins = 6; // Number of pT bins projected from the max particle pT within the jets histograms
  static const int knMaxJetsBins = 10;           // Maximum allowed number of high-pT jets in an event
  
private:
  
  // Naming for track histograms
  const char* fTrackHistogramNames[knTrackCategories] = {"track","trackUncorrected"}; // Names that different track histograms have in the input file
  const char* fTrackAxisNames[knTrackCategories] = {"Track","Uncorrected track"}; // Names attached to the figure axes
  
  // Naming for jet histograms
  const char* fJetHistogramName = "inclusiveJet";
  const char* fJetAxisName = "Jet";
  
  // Naming for multiplicity within the jet cone histograms
  const char* fMultiplicityInJetsHistogramNames[knMultiplicityInJetConeTypes] = {"multiplicityInJetCone", "multiplicityInReflectedCone", "multiplicityInJetConeUncorrected", "multiplicityInReflectedConeUncorrected"};
  const char* fMultiplicityInJetsAxisNames[knMultiplicityInJetConeTypes] = {"Multiplicity in jet cone", "Multiplicity in reflected cone", "UC multiplicity in jet cone", "UC multiplicity in reflected cone"};
  
  // Naming for particle density around the jet axis histograms
  const char* fParticleDensityAroundJetsHistogramNames[knParticleDensityAroundJetAxisTypes] = {"particleDensity", "particlePtDensity", "particleDensity", "particlePtDensity"};
  const char* fParticleDensityAroundJetsSaveNames[knParticleDensityAroundJetAxisTypes] = {"particleDensity", "particlePtDensity", "particleDensityPtBinned", "particlePtDensityPtBinned"};
  const char* fParticleDensityAroundJetsAxisNames[knParticleDensityAroundJetAxisTypes] = {"#rho(N_{ch})", "#rho(p_{T}^{ch})", "#rho(N_{ch})", "#rho(p_{T}^{ch})"};
  
  // Naming for energy-energy correlator histograms
  const char* fEnergyEnergyCorrelatorHistogramNames[knEnergyEnergyCorrelatorTypes] = {"energyEnergyCorrelator", "energyEnergyCorrelatorEfficiencyVariationPlus", "energyEnergyCorrelatorEfficiencyVariationMinus", "energyEnergyCorrelatorPairEfficiencyVariationPlus", "energyEnergyCorrelatorPairEfficiencyVariationMinus"};
  const char* fEnergyEnergyCorrelatorAxisNames[knEnergyEnergyCorrelatorTypes] = {"EEC", "EEC efficiency variation +", "EEC efficiency variation -", "EEC pair eff variation +", "EEC pair edd variation -"};
  const char* fEnergyEnergyCorrelatorProcessedSaveString[knEnergyEnergyCorrelatorProcessingLevels] = {"Normalized", "Background", "Signal", "Unfolded", "BackgroundAfterUnfolding", "UnfoldedSignal"};

  // Naming for reflected cone QA histograms
  const char* fReflectedConeQAFolderName = "reflectedConeQA";
  const char* fNumberOfJetsWithinReflectedConeName = "jetNumberInReflectedCone"; 
  const char* fJetPtWithinReflectedConeName = "jetPtInReflectedCone"; 
  
  // Naming for closure particle
  const char* fClosureParticleName[EECHistograms::knClosureParticleTypes+1] = {"_quark","_gluon",""};
  
  // Naming for subevent types
  const char* fSubeventTypeName[EECHistograms::knSubeventTypes] = {"Pythia", "Hydjet"};
  const char* fSubeventCombinationName[EECHistograms::knSubeventCombinations] = {"Pythia-Pythia", "Pythia-Hydjet", "Hydjet-Pythia", "Hydjet-Hydjet"};
  
  // Naming for pairing types
  const char* fPairingTypeSaveName[EECHistograms::knPairingTypes] = {"SignalConePair", "SignalReflectedConePair", "ReflectedConePair", "SignalMixedConePair", "MixedConePair", "SignalSecondMixedConePair", "MixedMixedConePair", "SecondMixedConePair", "SignalPerpendicularConePair", "PerpendicularConePair", "SignalSecondPerpendicularConePair", "PerpendicularPerpendicularConePair", "SecondPerpendicularConePair"};
  
  // Naming for jet cone types
  const char* fJetConeTypeSaveName[EECHistograms::knJetConeTypes] = {"", "ReflectedCone", "MixedCone"};
  
  // Maximum particle pT within the jet cone study
  const char* fMaxParticlePtInJetConeHistogramName = "maxParticlePtInJet";
  const char* fMaxParticlePtInJetConeSaveName[knMaxParticlePtWithinJetConeTypes] = {"maxParticlePtInJet", "maxBackgroundParticlePtInJet"};
  const char* fMaxParticlePtInJetConeAxisName[knMaxParticlePtWithinJetConeTypes] = {"Max particle p_{T}  (GeV)", "Max background p_{T}  (GeV)"};
  const double fProjectedMaxParticlePtBinBorders[knProjectedMaxParticlePtBins+1] = {10,15,20,30,40,50,500};

  // Naming for jet pT unfolding distributions
  const char* fJetPtUnfoldingDistributionName[knUnfoldingDistributionTypes] = {"jetPtUnfoldingMeasured", "jetPtUnfoldingTruth"};
  const char* fJetPtResponseMatrixName = "jetPtUnfoldingResponse";
  const char* fJetPtOneDimensionalUnfoldingDistributionName[knUnfoldingDistributionTypes] = {"oneDimensionalJetPtUnfoldingMeasured", "oneDimensionalJetPtUnfoldingTruth"};
  const char* fJetPtOneDimensionalResponseMatrixName = "oneDimensionalJetPtUnfoldingResponse";
  const char* fJetPtCovarianceMatrixName[knCovarianceMatrixTypes] = {"jetPtUnfoldingCovariance", "jetPtCovarianceAfterUnfolding"};

  // Naming for track/particle matching histograms
  const char* fTrackParticleMatchingQAName[knTrackParticleMatchingQAHistograms] = {"particlesCloseToTracks","tracksWithMatchedParticle"};
  const char* fTrackParticleMatchingResponseName[knTrackParticleMatchingResponseTypes] = {"particleDeltaRResponseMatrix","particlePtResponseMatrix"};
  const char* fTrackParticlePtClosureSaveName = "particlePairPtClosure";
  
public:
  
  EECHistogramManager();                                    // Default constructor
  EECHistogramManager(TFile* inputFile);                    // Constructor with input file
  EECHistogramManager(TFile* inputFile, EECCard* card);     // Constructor with input file and card
  EECHistogramManager(EECCard* card);                       // Constructor with card
  EECHistogramManager(const EECHistogramManager& in);       // Copy constructor
  ~EECHistogramManager();                                   // Destructor
  
  void LoadHistograms();          // Load the histograms from the inputfile
  void Write(const char* fileName, const char* fileOption);          // Write all the loaded histograms into a file
  void WriteProcessed(const char* fileName, const char* fileOption); // Write the processed histograms into a file
  void WriteCombinedMixedConeHistograms(const char* fileName, const char* fileOption); // Write the mixed cone histograms after combining them
  void WriteUnfoldedEnergyEnergyCorrelators(const char* fileName, const char* fileOption); // Write the unfolded energy-energy correlators to a file
  void WriteProcessedAfterUnfolding(const char* fileName, const char* fileOption); // Write the processed energy-energy correlators after unfolding into a file
  void WriteCovarianceMatrixAfterUnfolding(const char* fileName, const char* fileOption); // Write the unfolded covariance matrices into a file
  void LoadProcessedHistograms(); // Load processed histograms from the inputfile
  
  // Setters for binning information
  void SetCentralityBins(const bool readBinsFromFile, const int nBins, const double* binBorders, bool setIndices = true); // Set up centrality bin indices according to provided bin borders
  void SetTrackPtBins(const bool readBinsFromFile, const int nBins, const double* binBorders, bool setIndices = true);    // Set up track pT bin indices according to provided bin borders
  void SetJetPtBinsEEC(const bool readBinsFromFile, const int nBins, const double* binBorders, bool setIndices = true);   // Set up jet pT bin indices for energy-energy correlator according to provided bin borders
  void SetTrackPtBinsEEC(const bool readBinsFromFile, const int nBins, const double* binBorders, bool setIndices = true); // Set up track pT bin indices for energy-energy correlator according to provided bin borders
  
  // Setters for event information and dijets
  void SetLoadEventInformation(const bool loadOrNot); // Setter for loading event information
  
  // Setters for jets
  void SetLoadJetHistograms(const bool loadOrNot);        // Setter for loading all jet histograms
  
  // Setters for tracks
  void SetLoadTracks(const bool loadOrNot);            // Setter for loading tracks
  void SetLoadTracksUncorrected(const bool loadOrNot); // Setter for loading uncorrected tracks
  void SetLoadAllTracks(const bool loadTracks, const bool loadUncorrected); // Setter for loading all track histograms
  
  // Setters for multiplicity and particle density around jets
  void SetLoadMultiplicityInJets(const bool loadOrNot);          // Setter for loading track multiplicity within the jet cone
  void SetLoadParticleDensityAroundJets(const bool loadOrNot);   // Setter for loading particle density histograms around the jet axis
  void SetLoadParticlePtDensityAroundJets(const bool loadOrNot); // Setter for loading particle pT density histograms around the jet axis
  void SetLoadParticleDensityAroundJetsPtBinned(const bool loadOrNot);   // Setter for loading pT binned particle density histograms around the jet axis
  void SetLoadParticlePtDensityAroundJetsPtBinned(const bool loadOrNot); // Setter for loading pT binned particle pT density histograms around the jet axis
  void SetLoadAllParticleDensitiesAroundJets(const bool loadRegular, const bool loadPtWeighted, const bool loadPtBinned, const bool loadPtBinnedPtWeighted); // Setter for loading all particle density histograms around the jet axis
  
  // Setter for loading maximum particle pT within the jet cone histograms
  void SetLoadMaxParticlePtWithinJetCone(const bool loadOrNot);
  
  // Setters for energy-energy correlator histograms
  void SetLoadEnergyEnergyCorrelators(const bool loadOrNot);                         // Setter for loading energy-energy correlators
  void SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(const bool loadOrNot);  // Setter for loading energy-energy correlators with positive track efficiency variation 
  void SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(const bool loadOrNot); // Setter for loading energy-energy correlators with negative track efficiency variation 
  void SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(const bool loadOrNot);  // Setter for loading energy-energy correlators with positive track pair efficiency variation 
  void SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(const bool loadOrNot); // Setter for loading energy-energy correlators with negative track pair efficiency variation 
  void SetLoadAllEnergyEnergyCorrelators(const bool loadRegular, const bool loadEfficiencyVariationPlus, const bool loadEfficiencyVariationMinus, const bool loadPairEfficiencyVariationPlus, const bool loadPairEfficiencyVariationMinus); // Setter for loading all energy-energy correlators
  void SetLoadOnlyUnfoldedSignal(const bool loadOrNot);  // Setter for loading only unfolded signals for energy-energy correlators
  
  // Setters for reflected cone QA
  void SetLoadReflectedConeQAHistograms(const bool loadOrNot);

  // Setters for jet pT unfolding study
  void SetLoadJetPtUnfoldingHistograms(const bool loadOrNot);               // Setter for loading histograms needed in the jet pT unfolding study
  void SetLoadJetPtUnfoldingCovariance(const bool loadOrNot);               // Setter for loading the covariance histograms for input to jet pT unfolding
  void SetLoadJetPtOneDimensionalUnfoldingHistograms(const bool loadOrNot); // Setter for loading histograms needed in the one dimensional jet pT unfolding study

  // Setters for track/particle matching study
  void SetLoadTrackParticleMatchingHistograms(const bool loadOrNot); // Setter for loading histograms needed in track/particle matching study

  // Setters for jet flavor and loaded weight exponent
  void SetJetFlavor(const int iFlavor);  // For Monte Carlo, can select if we are looking for quark or gluon initiated jets
  void SetLoadedWeightExponent(const double weightExponent); // Define the weight exponent value that is searched from the file
  
  // Setter for loading additional histograms
  void SetLoad2DHistograms(const bool loadOrNot);           // Setter for loading two-dimensional histograms
  void SetLoadJetPtClosureHistograms(const bool loadOrNot); // Setter for loading jet pT closure histograms
  void SetLoadJetPtResponseMatrix(const bool loadOrNot);    // Setter for loading jet pT response matrix
  
  // Setters for ranges for different bins
  void SetCentralityBinRange(const int first, const int last);          // Setter for centrality bin range
  void SetTrackPtBinRange(const int first, const int last);             // Setter for track pT bin range
  void SetJetPtBinRangeEEC(const int first, const int last);            // Setter for jet pT bin range in energy-energy correlator histograms
  void SetTrackPtBinRangeEEC(const int first, const int last);          // Setter for track pT bin range in energy-energy correlator histograms
  
  // Unfolding is done in a separate macro. Thus provide setter for unfolded energy-energy correlators so they can be stored in the histogram manager
  void SetUnfoldedEnergyEnergyCorrelator(const TH1D* unfoldedEnergyEnergyCorrelator, const int iEnergyEnergyCorrelatorType, const int iCentrality, const int iJetPt, const int iTrackPt);
  void SetUnfoldedCoverianceMatrix(const TH2D* unfoldedCovarianceMatrix, const int iCentrality, const int iTrackPt);

  // Combine mixed cone histograms from other histogram manager to those in this histogram manager
  void CombineMixedConeBackgrounds(EECHistogramManager* anotherManager);

  // Getters for number of bins in histograms
  int GetNCentralityBins() const;          // Getter for the number of centrality bins
  int GetNTrackPtBins() const;             // Getter for the number of track pT bins
  int GetNJetPtBinsEEC() const;            // Getter for the number of jet pT bins in energy-energy correlator histograms
  int GetNTrackPtBinsEEC() const;          // Getter for the number of track pT bins in energy-energy correlator histograms
  double GetCentralityBinBorder(const int iCentrality) const;  // Getter for i:th centrality bin border
  double GetTrackPtBinBorder(const int iTrackPt) const;        // Getter for i:th track pT bin border
  double GetJetPtBinBorderEEC(const int iJetPt) const;         // Getter for i:th jet pT bin border in energy-energy correlator histograms
  double GetTrackPtBinBorderEEC(const int iTrackPt) const;     // Getter for i:th track pT bin border in energy-energy correlator histograms
  double GetMaxTrackPtWithinJetConeBinBorder(const int iTrackPt) const; // Getter for i:th track pT bin border in projections for maximum particle pT within the jet cone
  
  // Getters for histogram and axis naming
  const char* GetTrackHistogramName(int iTrackType) const; // Getter for track histogram name
  const char* GetTrackAxisName(int iTrackType) const;      // Getter for name suitable for x-axis in a given track histogram
  
  const char* GetJetHistogramName() const; // Getter for the jet histogram name
  const char* GetJetAxisName() const;      // Getter for name suitable for x-axis in a jet histogram
  
  const char* GetMultiplicityInJetConeHistogramName(int iMultiplicityType) const; // Getter for multiplicity in jet cone histogram name
  const char* GetMultiplicityInJetConeAxisName(int iMultiplicityType) const;      // Getter for multiplicity in jet cone axis name
  
  const char* GetParticleDensityAroundJetAxisHistogramName(int iParticleDensityType) const; // Getter for the particle density around jet axis histogram name
  const char* GetParticleDensityAroundJetAxisSaveName(int iParticleDensityType) const; // Getter for the particle density around jet axis save name
  const char* GetParticleDensityAroundJetAxisAxisName(int iParticleDensityType) const; // Getter for the particle density around jet axis axis name
  
  const char* GetEnergyEnergyCorrelatorHistogramName(int iEnergyEnergyCorrelatorType) const; // Getter for energy-energy correlator histogram name
  const char* GetEnergyEnergyCorrelatorAxisName(int iEnergyEnergyCorrelatorType) const;      // Getter for energy-energy correlator axis name
  const char* GetEnergyEnergyCorrelatorProcessSaveName(int iProcessingLevel) const;          // Getter for energy-energy correlator processing save name
  
  const char* GetSubeventType(const int iSubeventType) const;            // Getter for subevent type
  const char* GetSubeventCombination(const int iSubeventType) const;     // Getter for subevent combination
  TString GetSubeventCombinationSaveName(const int iSubeventType) const; // Getter for a well thought save name for subevent combination
  
  const char* GetPairingTypeSaveName(const int iPairingType) const; // Getter for pairing type save names
  
  const char* GetJetConeTypeSaveName(const int iJetConeType) const; // Getter for jet cone type save name
  
  const char* GetMaxParticlePtWithinJetConeSaveName(const int iMaxParticlePtWithinJetConeType) const; // Getter for maximum particle pT within jet cone save name
  const char* GetMaxParticlePtWithinJetConeAxisName(const int iMaxParticlePtWithinJetConeType) const; // Getter for maximum particle pT within jet cone axis name

  const char* GetJetPtUnfoldingCovarianceSaveName(const int iCovarianceType) const; // Getter for the save name for covariance matrices
  
  TString GetSystem() const;  // Getter for collision system
  
  // Getters for event information histograms
  TH1D* GetHistogramVertexZ();            // Getter for z-vertex histogram
  TH1D* GetHistogramVertexZWeighted();    // Getter for weighted z-vertex histogram
  TH1D* GetHistogramEvents();             // Getter for histogram for number of events surviving different event cuts
  TH1D* GetHistogramTriggers();           // Getter for histogram for trigger selection
  TH1D* GetHistogramTrackCuts();          // Getter for histogram for number of tracks surviving different track cuts
  TH1D* GetHistogramCentrality();         // Getter for centrality histogram in all events
  TH1D* GetHistogramCentralityWeighted(); // Getter for weighted centrality histogram in all events
  TH1D* GetHistogramPtHat();              // Getter for pThat histogram
  TH1D* GetHistogramPtHatWeighted();      // Getter for weighted pThat histogram
  
  TH1D* GetHistogramMultiplicity(int iCentrality);               // Getter for multiplicity from all events
  TH1D* GetHistogramMultiplicityWeighted(int iCentrality);       // Getter for efficiency weighted multiplicity from all events
  TH2D* GetHistogramMultiplicityMap();                           // Getter for multiplicity vs. centrality map
  TH2D* GetHistogramWeightedMultiplicityMap();                   // Getter for efficiency weighted multiplicity vs. centrality map

  // Getters for HF energy histograms
  TH1D* GetHistogramHFPlus(int iJetPt = -1, int nJets = -1);             // Getter for energy in HF plus calorimeters
  TH1D* GetHistogramHFMinus(int iJetPt = -1, int nJets = -1);            // Getter for energy in HF minus calorimeters
  TH1D* GetHistogramHFSum(int iJetPt = -1, int nJets = -1);              // Getter for energy in HF calorimeters
  TH2D* GetHistogramHFPlusVsHFMinus(int iJetPt = -1, int nJets = -1);    // Getter for 2D map between HF plus and HF minus

  TH1D* GetHistogramHFPlus(std::pair<double,double> jetPtBin, int nJets = -1);           // Getter for energy in HF plus calorimeters
  TH1D* GetHistogramHFMinus(std::pair<double,double> jetPtBin, int nJets = -1);          // Getter for energy in HF minus calorimeters
  TH1D* GetHistogramHFSum(std::pair<double,double> jetPtBin, int nJets = -1);            // Getter for energy in HF calorimeters
  TH2D* GetHistogramHFPlusVsHFMinus(std::pair<double,double> jetPtBin, int nJets = -1);  // Getter for 2D map between HF plus and HF minus
  
  // Getters for jet histograms
  TH1D* GetHistogramJetPt(int iCentrality);     // Jet pT histograms
  TH1D* GetHistogramJetPhi(int iCentrality);    // Jet phi histograms
  TH1D* GetHistogramJetEta(int iCentrality);    // Jet eta histograms
  TH1D* GetHistogramJetEtaCM(int iCentrality);  // Jet eta_CM histograms
  TH2D* GetHistogramJetEtaPhi(int iCentrality); // 2D eta-phi histogram for jets
  
  // Getters for histograms for tracks
  TH1D* GetHistogramTrackPt(const int iTrackType, const int iCentrality);                         // Track pT histograms
  TH1D* GetHistogramTrackPhi(const int iTrackType, const int iCentrality, const int iTrackPt);    // Track phi histograms
  TH1D* GetHistogramTrackEta(const int iTrackType, const int iCentrality, const int iTrackPt);    // Track eta histograms
  TH2D* GetHistogramTrackEtaPhi(const int iTrackType, const int iCentrality, const int iTrackPt); // 2D eta-phi histogram for track
  
  // Getters for multiplicity and particle density histograms within the jet cones
  TH1D* GetHistogramMultiplicityInJetCone(const int iCentrality, const int iJetPt, const int iTrackPt, const int MultiplicityType = kMultiplicityInJetCone, const int iSubevent = EECHistograms::knSubeventTypes); // Multiplicity within the jet cone
  TH1D* GetHistogramParticleDensityAroundJetAxis(const int iCentrality, const int iJetPt, const int iTrackPt, const int iJetConeType = EECHistograms::kSignalCone, const int iParticleDensityType = kParticleDensityAroundJetAxis, const int iSubevent = EECHistograms::knSubeventTypes); // Particle density around the jet axis
  
  // Getters for the maximum particle pT histograms within the jet cone
  TH1D* GetHistogramMaxParticlePtInJetCone(const int iMaxParticlePtWithinJetConeType, const int iCentrality, const int iJetPt, const int iTrackPt = knProjectedMaxParticlePtBins); // Maximum particle pT in jet cone
  TH1D* GetHistogramMaxParticlePtInJetConePtCut(const int iMaxParticlePtWithinJetConeType, const int iCentrality, const int iJetPt, const int iTrackPt); // Maximum particle pT in jet cone with pT cut for background particles
  TH1D* GetHistogramMaxSignalParticlePtInJetCone(const int iCentrality, const int iJetPt, const int iTrackPt = knProjectedMaxParticlePtBins); // Maximum signal particle pT in jet cone
  TH1D* GetHistogramMaxSignalParticlePtInJetConePtCut(const int iCentrality, const int iJetPt, const int iTrackPt); // Maximum signal particle pT in jet cone with pT cut for background particles
  TH1D* GetHistogramMaxBackgroundParticlePtInJetCone(const int iCentrality, const int iJetPt, const int iTrackPt = knProjectedMaxParticlePtBins); // Maximum background particle pT in jet cone
  TH1D* GetHistogramMaxBackgroundParticlePtInJetConePtCut(const int iCentrality, const int iJetPt, const int iTrackPt); // Maximum background particle pT in jet cone with pT cut for signal particles
  
  // Getters for energy-energy correlator histograms
  TH1D* GetHistogramEnergyEnergyCorrelator(const int iEnergyEnergyCorrelatorType, const int iCentrality, const int iJetPt, const int iTrackPt, const int iPairingType = EECHistograms::kSameJetPair, const int iSubevent = EECHistograms::knSubeventCombinations);  // Energy-energy correlator histograms
  TH1D* GetHistogramEnergyEnergyCorrelator(const int iEnergyEnergyCorrelatorType, const std::pair<int,int> centralityBin, const std::pair<int,int> jetPtBin, const double trackPtCut, const int iPairingType = EECHistograms::kSameJetPair, const int iSubevent = EECHistograms::knSubeventCombinations);  // Energy-energy correlator histograms
  TH1D* GetHistogramEnergyEnergyCorrelatorProcessed(const int iEnergyEnergyCorrelatorType, const int iCentrality, const int iJetPt, const int iTrackPt, const int iProcessingLevel);  // Processed energy-energy correlator histograms

  // Getters for reflected cone QA histograms
  TH1D* GetHistogramNumberOfJetsWithinReflectedCone(const int iCentrality);
  TH1D* GetHistogramJetPtWithinReflectedCone(const int iCentrality);
  
  // Getter for jet pT response matrix
  TH2D* GetHistogramJetPtResponseMatrix(const int iCentrality);

  // Getter for jet pT closure histograms
  TH1D* GetHistogramJetPtClosure(const int iGenPtBin, const int iEtaBin, const int iPhiBin, const int iCentrality, const int iClosureParticle); // Jet pT closure
  
  // Getters for jet pT unfolding histograms
  TH1D* GetHistogramJetPtUnfoldingMeasured(const int iCentrality, const int iTrackPt); // Getter for measured jet pT unfolding distribution
  TH1D* GetHistogramJetPtUnfoldingTruth(const int iCentrality, const int iTrackPt);    // Getter for truth jet pT unfolding distribution
  TH2D* GetHistogramJetPtUnfoldingResponse(const int iCentrality, const int iTrackPt); // Getter for jet pT unfolding response
  TH2D* GetHistogramJetPtUnfoldingCovariance(const int iCovarianceMatrixType, const int iCentrality, const int iTrackPt); // Getter for jet pT unfolding covariance

  TH1D* GetHistogramJetPtOneDimensionalUnfoldingMeasured(const int iCentrality); // Getter for measured jet pT one dimensional unfolding distribution
  TH1D* GetHistogramJetPtOneDimensionalUnfoldingTruth(const int iCentrality);    // Getter for truth jet pT one dimensional unfolding distribution
  TH2D* GetHistogramJetPtOneDimensionalUnfoldingResponse(const int iCentrality); // Getter for one dimensional jet pT unfolding response

  // Getters for track/particle matching histograms
  TH1D* GetHistogramParticlesCloseToTrack(const int iCentrality, const int iJetPt, const int iTrackPt); // Getter for number of particles close to tracks
  TH1D* GetHistogramHasMatchingParticle(const int iCentrality, const int iJetPt, const int iTrackPt); // Getter for flag if a matching particle is found
  TH2D* GetHistogramTrackParticleDeltaRResponse(const int iCentrality, const int iJetPt, const int iTrackPt); // Getter for deltaR response matrix between track pairs and matched particle pairs
  TH2D* GetHistogramTrackParticlePtResponse(const int iCentrality, const int iJetPt, const int iTrackPt); // Getter for pT response matrix between pT1*pT2 from track pairs and particle pairs
  TH1D* GetHistogramTrackParticlePtClosure(const int iCentrality, const int iJetPt, const int iTrackPt); // Getter for track pT1*pT2 / particle pT1*pT2 histograms

  // Generic getters for one and two dimensional histograms
  TH1D* GetOneDimensionalHistogram(TString name, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0, int bin6 = 0, int bin7 = 0); // Getter for any one-dimensional histogram based on input string
  TH2D* GetTwoDimensionalHistogram(TString name, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0); // Getter for any two-dimensional histogram based on input string
  
  
  // Getters for the loaded centrality and track pT bins
  int GetFirstCentralityBin() const;          // Get the first loaded centrality bin
  int GetLastCentralityBin() const;           // Get the last loaded centrality bin
  int GetFirstTrackPtBin() const;             // Get the first loaded track pT bin
  int GetLastTrackPtBin() const;              // Get the last loaded track pT bin
  int GetFirstJetPtBinEEC() const;            // Get the first loaded energy-energy correlator jet pT bin
  int GetLastJetPtBinEEC() const;             // Get the last loaded energy-energy correlator jet pT bin
  int GetFirstTrackPtBinEEC() const;          // Get the first loaded energy-energy correlator track pT bin
  int GetLastTrackPtBinEEC() const;           // Get the last loaded energy-energy correlator track pT bin
  
  // Getters for normalization information
  int GetNEvents();                      // Getter for the number of events passing the cuts
  double GetJetPtIntegral(const int iCentrality); // Getter for integral over inclusive jet pT in a given centrality
  double GetJetPtIntegral(int iCentrality, const double minPt, const double maxPt); // Getter for integral over inclusive jet pT in a given pT range within a given centrality bin
  
  // Getter for the card
  EECCard* GetCard() const;  // Getter for the JCard
  
  // Post-processing the energy-energy correlator histograms
  void SubtractBackground(int iMethod, const int iSystematic);             // Subtract the background from the energy-energy correlator histograms
  void SubtractBackgroundFromUnfolded(int iMethod, const int iSystematic); // Subtract the background from unfolded energy-energy correlator histograms
  
private:
  
  // Data members
  TFile* fInputFile;                  // File from which the histograms are read
  EECCard* fCard;                     // Card inside the data file for binning, cut collision system etc. information
  TString fSystemAndEnergy;           // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy;    // Same a before but without white spaces and dots
  
  // ==============================================
  // ======== Flags for histograms to load ========
  // ==============================================
  
  bool fLoadEventInformation;                              // Load the event information histograms
  bool fLoadJets;                                          // Load the jet histograms
  bool fLoadTracks[knTrackCategories];                     // Load the track histograms
  bool fLoad2DHistograms;                                  // Load also two-dimensional (eta,phi) histograms
  bool fLoadJetPtClosureHistograms;                        // Load the jet pT closure histograms
  bool fLoadJetPtResponseMatrix;                           // Load the jet pT response matrix
  bool fLoadMultiplicityInJetHistograms;                   // Load the multiplicity histograms within jet cone
  bool fLoadParticleDensityAroundJetsHistograms[knParticleDensityAroundJetAxisTypes];  // Load the particle density histograms around the jet axis
  bool fLoadMaxParticlePtWithinJetConeHistograms;          // Load the maximum particle pT within the jet cone histograms
  bool fLoadEnergyEnergyCorrelatorHistograms[knEnergyEnergyCorrelatorTypes];   // Load the energy-energy correlator histograms
  bool fLoadReflectedConeQAHistograms;                     // Load the reflected cone QA histograms
  bool fLoadJetPtUnfoldingHistograms;                      // Load the histograms needed in jet pT unfolding study
  bool fLoadJetPtUnfoldingCovariance;                      // Load the covariance distribution for unfolding
  bool fLoadJetPtOneDimensionalUnfoldingHistograms;        // Load the histograms needed in the one-dimensional jet pT unfolding study
  bool fLoadTrackParticleMatchingHistograms;               // Load the histograms for track/particle matching study
  int  fJetFlavor;                                         // Select the flavor for loaded jets (1 = Quark, 2 = Gluon)
  double fLoadedWeightExponent;                            // Value for weight exponent in energy-energy correlators that is searched from the files
  bool fLoadOnlyUnfoldedSignal;                            // In many cases we only need unfolded signal. To save memory, we can only load these histograms from the file
  
  // ==============================================
  // ======== Ranges of histograms to load ========
  // ==============================================
  
  int fFirstLoadedCentralityBin;  // First centrality bin that is loaded
  int fLastLoadedCentralityBin;   // Last centrality bin that is loaded
  int fFirstLoadedTrackPtBin;     // First track pT bin that is loaded
  int fLastLoadedTrackPtBin;      // Last track pT bin that is loaded
  int fFirstLoadedJetPtBinEEC;    // First loaded jet pT bin for the energy-energy correlator histograms
  int fLastLoadedJetPtBinEEC;     // Last loaded jet pT bin for the energy-energy correlator histograms
  int fFirstLoadedTrackPtBinEEC;  // First loaded track pT bin for the energy-energy correlator histograms
  int fLastLoadedTrackPtBinEEC;   // Last loaded track pT bin for the energy-energy correlator histograms
  
  // =============================================
  // ============ Binning information ============
  // =============================================
  int fCentralityBinIndices[kMaxCentralityBins+1];           // Indices for centrality bins in centrality binned histograms
  double fCentralityBinBorders[kMaxCentralityBins+1];        // Centrality bin borders, from which bin indices are obtained
  int fTrackPtBinIndices[kMaxTrackPtBins+1];                 // Indices for track pT bins in track eta and phi histograms
  double fTrackPtBinBorders[kMaxTrackPtBins+1];              // Track pT bin borders, from which bin indices are obtained
  int fJetPtIndicesEEC[kMaxJetPtBinsEEC+1];                  // Indices for jet pT bins in energy-energy correlator histograms
  double fJetPtBinBordersEEC[kMaxJetPtBinsEEC+1];            // Jet pT bin borders in energy-energy correlator histograms
  int fTrackPtIndicesEEC[kMaxTrackPtBinsEEC+1];              // Indices for track pT bins in energy-energy correlator histograms
  double fTrackPtBinBordersEEC[kMaxTrackPtBinsEEC+1];        // Track pT bin borders in energy-energy correlator histograms
  int fnCentralityBins;                                      // Number of centrality bins in the JCard of the data file
  int fnTrackPtBins;                                         // Number of track pT bins in the JCard of the data file
  int fnJetPtBinsEEC;                                        // Number of jet pT bins for the energy-energy correlator histograms
  int fnTrackPtBinsEEC;                                      // Number of track pT bins for the energy-energy correlator histograms
  std::vector<int>* fIncludedPairingTypes;                   // Vector for pairing types that are included in the data file

  // =============================================
  // ===== Histograms for the dijet analysis =====
  // =============================================

  // Event information histograms
  TH1D* fhVertexZ;            // Vertex z position
  TH1D* fhVertexZWeighted;    // Weighted vertex z-position (only meaningfull for MC)
  TH1D* fhEvents;             // Number of events surviving different event cuts
  TH1D* fhTriggers;           // Trigger selection for the events
  TH1D* fhTrackCuts;          // Number of tracks surviving different track cuts
  TH1D* fhCentrality;         // Centrality of all events
  TH1D* fhCentralityWeighted; // Weighted centrality distribution in all events (only meaningful for MC)
  TH1D* fhPtHat;              // pT hat for MC events (only meaningful for MC)
  TH1D* fhPtHatWeighted;      // Weighted pT hat distribution (only meaningful for MC)

  TH1D* fhMultiplicity[kMaxCentralityBins+1];          // Multiplicity form all events
  TH1D* fhMultiplicityWeighted[kMaxCentralityBins+1];  // Efficiency weighted multiplicity form all events
  TH2D* fhMultiplicityMap;                             // Multiplicity vs. centrality map
  TH2D* fhMultiplicityMapWeighted;                     // Efficiency weighted multiplicity vs. centrality map

  TH1D* fhHFPlus[kMaxJetPtBinsEEC][knMaxJetsBins+1];           // Energy in HF plus calorimeters
  TH1D* fhHFMinus[kMaxJetPtBinsEEC][knMaxJetsBins+1];          // Energy in HF minus calorimeters
  TH1D* fhHFSum[kMaxJetPtBinsEEC][knMaxJetsBins+1];            // Energy in HF calorimeters
  TH2D* fhHFPlusVsHFMinus[kMaxJetPtBinsEEC][knMaxJetsBins+1];  // 2D-map between energy in plus and minus sides of HF calorimeters

  // Histograms for jets
  TH1D* fhJetPt[kMaxCentralityBins];      // Jet pT histograms
  TH1D* fhJetPhi[kMaxCentralityBins];     // Jet phi histograms
  TH1D* fhJetEta[kMaxCentralityBins];     // Jet eta histograms
  TH1D* fhJetEtaCM[kMaxCentralityBins];   // Jet eta_CM histograms
  TH2D* fhJetEtaPhi[kMaxCentralityBins];  // 2D eta-phi histogram for jets

  // Histograms for tracks
  TH1D* fhTrackPt[knTrackCategories][kMaxCentralityBins];                         // Track pT histograms
  TH1D* fhTrackPhi[knTrackCategories][kMaxCentralityBins][kMaxTrackPtBins+1];     // Track phi histograms
  TH1D* fhTrackEta[knTrackCategories][kMaxCentralityBins][kMaxTrackPtBins+1];     // Track eta histograms
  TH2D* fhTrackEtaPhi[knTrackCategories][kMaxCentralityBins][kMaxTrackPtBins+1];  // 2D eta-phi histogram for track

  // Histograms for multiplicity and particle density within the jet cones
  TH1D* fhMultiplicityInJetCone[kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC][knMultiplicityInJetConeTypes][EECHistograms::knSubeventTypes+1];
  TH1D* fhParticleDensityAroundJetAxis[kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC][EECHistograms::knJetConeTypes][knParticleDensityAroundJetAxisTypes][EECHistograms::knSubeventTypes + 1];

  // Histograms for maximum particle pT within the jet cone
  TH1D* fhMaxParticlePtInJetConePtCut[knMaxParticlePtWithinJetConeTypes][kMaxCentralityBins][kMaxJetPtBinsEEC][knProjectedMaxParticlePtBins];
  TH1D* fhMaxParticlePtInJetConePtBin[knMaxParticlePtWithinJetConeTypes][kMaxCentralityBins][kMaxJetPtBinsEEC][knProjectedMaxParticlePtBins+1];

  // Histograms for energy-energy correlators
  TH1D* fhEnergyEnergyCorrelator[knEnergyEnergyCorrelatorTypes][kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC][EECHistograms::knPairingTypes][EECHistograms::knSubeventCombinations+1];  // Raw correlators read from data file
  TH1D* fhEnergyEnergyCorrelatorProcessed[knEnergyEnergyCorrelatorTypes][kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC][knEnergyEnergyCorrelatorProcessingLevels];   // Postprocessed energy-energy correlators

  // Quality assurance histograms for reflected cone
  TH1D* fhNumberOfJetsWithinReflectedCone[kMaxCentralityBins];
  TH1D* fhJetPtWithinReflectedCone[kMaxCentralityBins];

  // Jet pT response matrix
  TH2D* fhJetPtResponseMatrix[kMaxCentralityBins];

  // Histograms for jet pT closure
  TH1D* fhJetPtClosure[knGenJetPtBins+1][knJetEtaBins+1][knJetPhiBins+1][kMaxCentralityBins][EECHistograms::knClosureParticleTypes+1];  // Jet pT closure

  // Histograms for jet pT unfolding study
  TH1D* fhJetPtUnfoldingDistribution[knUnfoldingDistributionTypes][kMaxCentralityBins][kMaxTrackPtBinsEEC];
  TH2D* fhJetPtUnfoldingResponse[kMaxCentralityBins][kMaxTrackPtBinsEEC];
  TH2D* fhJetPtUnfoldingCovariance[knCovarianceMatrixTypes][kMaxCentralityBins][kMaxTrackPtBinsEEC];

  // Histograms for one-dimensional jet pT unfolding study
  TH1D* fhOneDimensionalJetPtUnfoldingDistribution[knUnfoldingDistributionTypes][kMaxCentralityBins];
  TH2D* fhOneDimensionalJetPtUnfoldingResponse[kMaxCentralityBins];

  // Histograms to study track/particle matching
  TH1D* fhTrackParticleMatchQA[knTrackParticleMatchingQAHistograms][kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC];
  TH2D* fhTrackParticleResponse[knTrackParticleMatchingResponseTypes][kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC];
  TH1D* fhTrackParticlePtClosure[kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC];

  // Private methods
  void InitializeFromCard(); // Initialize several member variables from EECCard
  
  // Methods related to processing the histograms
  void NormalizeToDeltaRBinArea(TH1D* histogramInNeedOfNormalization); // Normalize each bin in a histogram histogram to DeltaR bin area
  
  // Binning related methods
  void SetBinIndices(const char* histogramName, const int nBins, int* binIndices, const double* binBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetBinBordersAndIndices(const char* histogramName, const int nBins, double *copyBinBorders, int* binIndices, const double* binBorders, const int iAxis, const bool setIndices); // Read the bin indices for given bin borders
  
  // Finders for histograms with different amount of restrictions
  TH2D* FindHistogram2D(THnBase* histogramArray, int xAxis, int yAxis, const bool normalizeToBinWidth = true); // Extract a 2D histogram without axis restrictions from THnSparseD
  TH2D* FindHistogram2D(THnBase* histogramArray, int xAxis, int yAxis, int nAxes, int* axisNumber, int* lowBinIndex, int* highBinIndex, const bool normalizeToBinWidth = true); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH2D* FindHistogram2D(THnBase* histogramArray, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0, const bool normalizeToBinWidth = true); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(THnBase* histogramArray, int xAxis, const bool normalizeToBinWidth = true); // Extract a histogram without axis restrictions from THnSparseD
  TH1D* FindHistogram(THnBase* histogramArray, int xAxis, int nAxes, int* axisNumber, int* lowBinIndex, int* highBinIndex, const bool normalizeToBinWidth = true); // Extract a histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(THnBase* histogramArray, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0, const bool normalizeToBinWidth = true); // Extract a histogram using given axis restrictions from THnSparseD
  
  // Loaders for different groups of histograms
  void LoadHFEnergyHistograms(); // Loader for HF energy histograms
  void LoadMultiplicityHistograms(); // Loader for multiplicity histograms
  void LoadJetHistograms(); // Loader for jet histograms
  void LoadTrackHistograms(); // Loader for track histograms
  void LoadMultiplicityInJetConeHistograms(); // Loader for multiplicity histograms within the jet cone
  void LoadParticleDensityHistograms(); // Loader for particle density histograms around the jet cone
  void LoadMaxParticlePtInJetConeHistograms(); // Loader for maximum particle pT in jet cone histograms
  void LoadEnergyEnergyCorrelatorHistograms(); // Loader for energy-energy correlator histograms
  void LoadReflectedConeQAHistograms();        // Loader for reflected cone QA histograms
  void LoadJetPtResponseMatrix();    // Loader for the jet pT response matrices
  void LoadJetPtClosureHistograms(); // Loader for jet pT closure histograms
  void LoadJetPtUnfoldingHistograms(); // Loader for jet pT unfolding histograms
  void LoadJetPtUnfoldingCovariance(); // Loader for covariance histograms used in jet pT unfolding algorithm
  void LoadJetPtOneDimensionalUnfoldingHistograms(); // Loader for one dimensional jet pT unfolding histograms
  void LoadTrackParticleMatchingHistograms(); // Loader for track/particle matching histograms
  void StabilizeBackground(); // Stabilize the background histograms by combining all jet pT bins for histograms that do not depend on jet pT

  void LoadEventsHistogram();           // Loader for the nEvents histogram
  void LoadJetPtHistogram(const int iCentrality); // Loader for jet pT histogram

  // Handling of weight exponents
  void CheckWeightExponent(); // Check that the weight exponent requested is present in the input data
  void UpdateWeightExponent(); // If one weight exponent from many is projected, update the information in card
  
  // Generic setter for bin indice and borders
  void SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double* binBorders, const char* errorMessage, const int maxBins, const bool setIndices); // Generic bin setter
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int& first, int& last); // Sanity check for given binning
  int BinIndexCheck(const int nBins, const int binIndex) const; // Check that given index is in defined range
  
  // Methods for histogram writing
  void WriteJetHistograms();                          // Write the jet histograms to the file that is currently open
  void WriteTrackHistograms();                        // Write the track histograms to the file that is currently open
  void WriteMultiplicityInJetConeHistograms();        // Write the multiplicity histograms within the jet cone
  void WriteParticleDensityAroundJetsHistograms();    // Write the particle density histograms around the jet axes
  void WriteMaxParticlePtWithinJetConeHistograms();   // Write the maximum particle pT within the jet cone histograms
  void WriteEnergyEnergyCorrelatorHistograms();       // Write the energy-energy correlator histograms to the file that is currently open
  void WriteReflectedConeQAHistograms();              // Write the reflected cone QA histograms to the file that is currently open
  void WriteJetPtResponseMatrix();                    // Write the jet pT response matrices
  void WriteClosureHistograms();                      // Write the closure histograms to the file that is currently open
  void WriteJetPtUnfoldingHistograms();               // Write the jet pT unfolding histograms to the output file
  void WriteJetPtUnfoldingCovariance();               // Write the covariance histograms used in jet pT unfolding
  void WriteJetPtOneDimensionalUnfoldingHistograms(); // Write the jet pT one-dimensional unfolding histograms to the output file
  void WriteTrackParticleMatchingHistograms();        // Write the track/particle matching histograms to the output file
  
};

#endif

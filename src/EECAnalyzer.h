// Class for the main analysis algorithms for leading-subleading jet analysis

#ifndef EECANALYZER_H
#define EECANALYZER_H

// C++ includes
#include <vector>
#include <bitset>
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <tuple>      // For returning several arguments in a transparent manner
#include <fstream>
#include <string>
#include <map>        // To avoid using variable length arrays, that are clang extension

// Root includes
#include <TString.h>
#include <TRandom3.h>
#include <TMath.h>

// Own includes
#include "ConfigurationCard.h"
#include "EECHistograms.h"
#include "HighForestReader.h"
#include "GeneratorLevelForestReader.h"
#include "UnfoldingForestReader.h"
#include "JetCorrector.h"
#include "JetUncertainty.h"
#include "trackingEfficiency2018PbPb.h"
#include "trackingEfficiency2017pp.h"
#include "trackingEfficiency2016pPb.h"
#include "TrackingEfficiencyInterface.h"
#include "TrackPairEfficiencyCorrector.h"
#include "JetMetScalingFactorManager.h"
#include "SmearingProvider.h"

class EECAnalyzer{
  
private:
  
  enum enumFilledHistograms{kFillEventInformation, kFillJets, kFillTracks, kFillJetConeHistograms, kFillEnergyEnergyCorrelators, kFillEnergyEnergyCorrelatorsSystematics, kFillJetPtClosure, kFillJetPtUnfoldingResponse, kFillTrackParticleMatchingHistograms, knFillTypes}; // Histograms to fill
  enum enumSubeventCuts{kSubeventZero,kSubeventNonZero,kSubeventAny,knSubeventCuts}; // Cuts for subevent index
  enum enumMcCorrelationType{kRecoReco,kRecoGen,kGenReco,kGenGen,knMcCorrelationTypes}; // How to correlate jets and tracks in MC
  enum enumUnfoldingLevel{kUnfoldingReconstructed, kUnfoldingTruth, kNUnfoldingAxes}; // Select the axis on unfolding response matrix
  enum enumTupleDecoder{kTrackPt, kTrackEta, kTrackPhi, kTrackEfficiencyCorrection, kTrackCharge, kMatchIndex}; // Components of the n-tuple in track vector
  enum enumMatchForMixing{kMatchHiBin, kMatchHFEnergy, kMatchMultiplicity, knMatchMethods}; // Different methods to match the events for mixing
  enum enumSmallTupleDecoder{kPossibleMatchPtDifference, kPossibleMatchIndex}; // Components of the n-tuple in track vector

  // Definitions for nonsensical bin values
  static inline constexpr int kNonsenseValue = -50;           // Non-sensical mixing variable value
  static inline constexpr int kNonsenseBin = -500;            // Non-sensical mixing bin index
  
public:
  
  // Constructors and destructor
  EECAnalyzer(); // Default constructor
  EECAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard* newCard, Int_t runLocal, Int_t mixingListIndex); // Custom constructor
  EECAnalyzer(const EECAnalyzer& in); // Copy constructor
  virtual ~EECAnalyzer(); // Destructor
  EECAnalyzer& operator=(const EECAnalyzer& obj); // Equal sign operator
  
  // Methods
  void RunAnalysis();                     // Run the dijet analysis
  EECHistograms* GetHistograms() const;   // Getter for histograms

 private:
  
  // Private methods
  void CalculateEnergyEnergyCorrelator(const vector<double> selectedTrackPt[EECHistograms::knJetConeTypes], const vector<double> relativeTrackEta[EECHistograms::knJetConeTypes], const vector<double> relativeTrackPhi[EECHistograms::knJetConeTypes], const vector<int> selectedTrackSubevent[EECHistograms::knJetConeTypes], const double jetPt);  // Calculate energy-energy correlators
  void CalculateEnergyEnergyCorrelatorForUnfolding(const vector<double> selectedTrackPt, const vector<double> relativeTrackEta, const vector<double> relativeTrackPhi, const double jetPt, const double genPt); // Calculate energy-energy correlators for unfolding
  void FillOneDimensionalJetPtUnfoldingHistograms(const double jetPt, const double genPt); // Fill histograms for one dimensional jet pT unfolding
  void FillJetPtResponseMatrix(const Int_t jetIndex); // Fill jet pT response matrix
  void FillJetPtClosureHistograms(const Int_t jetIndex); // Fill jet pT closure histograms
  void FillUnfoldingResponse(); // Fill the histograms needed for unfolding study
  void ConstructParticleResponses(); // Construct DeltaR and pT1*pT2 response matrices
  void ReadConfigurationFromCard(); // Read all the configuration from the input card

  // Methods for event mixing
  void PrepareMixingVectors(); // Find vz and hiBin values from mixed event in preparation for event mixing
  void PrepareBinnedMixingVectors(); // Prepare variables needed for mixing in a binned manner
  void FindMultiplicityMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, Int_t const currentMultiplicity, const Int_t eventNumber, const Int_t iEvent, const Bool_t binnedMixing); // Find the mixed events that are matched with the signal event using multiplicity
  void FindHiBinMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, const Int_t hibin, const Int_t eventNumber, const Int_t iEvent, const Bool_t binnedMixing); // Find the mixed events that are matched with the signal event using hiBin
  void FindHFEnergyMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, Double_t hfEnergy, const Int_t eventNumber, const Int_t iEvent, const Bool_t binnedMixing); // Find the mixed events that are matched with the signal event using HF energy
  void FindBinnedMultiplicityMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, Int_t const currentMultiplicity, const Int_t eventNumber, const Int_t iEvent); // Find the mixed events that are matched with the signal event using multiplicity bins
  void FindBinnedHiBinMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, const Int_t hibin, const Int_t eventNumber, const Int_t iEvent); // Find the mixed events that are matched with the signal event using hiBin bins
  void FindBinnedHFEnergyMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, const Double_t hfEnergy, const Int_t eventNumber, const Int_t iEvent); // Find the mixed events that are matched with the signal event using HF energy bins
  Int_t FindMixingBin(Double_t value, Double_t anchor, Double_t width); // Find a generic mixing bin
  Int_t FindMixingBin(Double_t value, std::vector<Double_t> borders); // Find a generic mixing bin from bin borders
  Int_t FindMixingVzBin(Double_t vz); // Find mixing vz bin from a vz value
  Int_t FindMixingHFEnergyBin(Double_t hfEnergy); // Find mixing HF energy bin from HF energy value
  Int_t FindMixingHiBinBin(Int_t hiBin); // Find mixing HiBin bin from HiBin value
  Int_t FindMixingMultiplicityBin(Int_t multiplicity); // Find mixing multiplicity bin from multiplicity value
  std::vector<Double_t> GetMixingBinBorders(std::vector<Double_t> binBorderVector, std::vector<Double_t> valueVector); // Get vector with bin borders used in mixing
  std::vector<Double_t> GetBinnedMixingBinBorders(std::vector<Double_t> binBorderVector, std::vector<Int_t> binVector); // Get vector with bin borders used in mixing
  Double_t FindExampleValue(Int_t iBin, std::vector<Double_t> borders); // Find an example value from certain mixing bin
  
  // Other methods
  Bool_t PassSubeventCut(const Int_t subeventIndex) const;  // Check if the track passes the set subevent cut
  Bool_t PassTrackCuts(ForestReader* trackReader, const Int_t iTrack, TH1F* trackCutHistogram, const Bool_t bypassFill = false); // Check if a track passes all the track cuts
  Bool_t PassTrackCuts(UnfoldingForestReader* trackReader, const Int_t iTrack); // Check if a track passes all the track cuts
  Bool_t PassGenParticleSelection(ForestReader* trackReader, const Int_t iTrack); // Check if a generator particle passes the defined selections
  Bool_t PassGenParticleSelection(UnfoldingForestReader* trackReader, const Int_t iTrack); // Check if a generator particle passes the defined selections
  Bool_t PassEventCuts(ForestReader* eventReader, const Bool_t fillHistograms); // Check if the event passes the event cuts
  Double_t GetTrackEfficiencyCorrection(const Int_t iTrack); // Get the track efficiency correction for a given track
  Double_t  GetTrackEfficiencyCorrection(const Float_t trackPt, const Float_t trackEta, const Int_t hiBin); // Get the track efficiency correction for given track and event information
  Double_t GetVzWeight(const Double_t vz) const;  // Get the proper vz weighting depending on analyzed system
  Double_t GetCentralityWeight(const Int_t hiBin) const; // Get the proper centrality weighting depending on analyzed system
  Double_t GetMultiplicityWeight(const Double_t multiplicity) const; // Get the proper multiplicity weight for MC
  Double_t GetJetPtWeight(const Double_t jetPt) const; // Get the proper jet pT weighting for 2017 and 2018 MC
  Double_t GetSmearingFactor(Double_t jetPt, Double_t jetEta, const Double_t centrality); // Getter for jet pT smearing factor
  Int_t GetCentralityBin(const Double_t centrality) const; // Getter for centrality bin
  Double_t GetMultiplicity(); // Get the track multiplicity in the current event
  Int_t GetGenMultiplicity(ForestReader* trackReader, const Int_t subeventIndex, const Bool_t isMegaSkim); // Get the generator level multiplicity
  Double_t GetCentralityFromMultiplicity(const Double_t multiplicity) const; // Get the analysis centrality bin corresponding to the given multiplicity
  Double_t GetDeltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2) const; // Get deltaR between two objects
  Int_t GetSubeventCombination(const Int_t subevent1, const Int_t subevent2) const; // Get the subevent combination type from two track subevents
  Int_t GetSubeventIndex(const Int_t subevent) const; // Get the subevent index for a track
  Double_t GetReflectedEta(const Double_t eta) const; // Get jet eta reflected around zero, avoiding overlapping jet cones
  Double_t GetPerpendicularPhi(const Double_t phi, const Int_t direction) const; // Get perpendicular jet phi coordinate in a defined direction
  Double_t TransformToUnfoldingAxis(const Double_t deltaR, const Double_t jetPt) const; // Transform the deltaR value to the unfolding axis
  Double_t SimpleSmearDeltaR(const Double_t deltaR); // Simple smearing for deltaR to see if it affects the final distributions

  
  // Private data members
  ForestReader* fJetReader;                      // Reader for jets in the event
  ForestReader* fRecoJetReader;                  // Reader for reconstructed jets in the event
  ForestReader* fTrackReader;                    // Reader for tracks in the event
  ForestReader* fMixedEventReader;               // Reader for mixed events
  UnfoldingForestReader* fUnfoldingForestReader; // Reader for unfolding study
  std::vector<TString> fFileNames;               // Vector for all the files to loop over
  ConfigurationCard* fCard;                      // Configuration card for the analysis
  EECHistograms* fHistograms;                    // Filled histograms
  TF1* fVzWeightFunction;                        // Weighting function for vz. Needed for MC.
  TF1* fCentralityWeightFunctionCentral;         // Weighting function for central centrality classes. Needed for MC.
  TF1* fCentralityWeightFunctionPeripheral;      // Weighting function for peripheral centrality classes. Needed for MC.
  TF1* fMultiplicityWeightFunction;              // Track multiplicity based weighting function. Can be done instead of centrality weight.
  TF1* fPtWeightFunction;                        // Weighting function for jet pT. Needed for MC.
  TF1* fSmearingFunction;                        // Additional smearing for jets. Needed in systematic uncertainty study.
  TrackingEfficiencyInterface* fTrackEfficiencyCorrector;  // Tracking efficiency corrector
  JetCorrector* fJetCorrector;                   // Class for making jet energy correction
  JetUncertainty* fJetUncertainty;               // Class for finding uncertainty for jet pT
  TrackPairEfficiencyCorrector* fTrackPairEfficiencyCorrector; // Track pair efficiency corrector
  JetMetScalingFactorManager* fEnergyResolutionSmearingFinder; // Manager to find proper jet energy resolution scaling factors provided by the JetMet group
  SmearingProvider* fDeltaRSmearer;              // Realistic smearing for DeltaR done from response matrices
  SmearingProvider* fEnergyWeightSmearer;        // Realistic smearing for energy weights from reco/gen ratios
  TRandom3* fRng;                                // Random number generator

  // Histograms to determine the covariance matrix
  TH1D* fThisEventCorrelator[5][10];        // Energy-energy correlators calculated exclusively from this event
  
  // Analyzed data and forest types
  Int_t fDataType;                   // Analyzed data type
  Int_t fTriggerSelection;           // 0 = Do not use any triggers, 1 = Require CaloJet80, 2 = Require CaloJet100, 3 = Require CaloJet80 or CaloJet100
  Int_t fJetType;                    // Type of jets used for analysis. 0 = Calo jets, 1 = PF jets
  Int_t fMatchJets;                  // Jet matching flag. 0 = Do not match jets, 1 = Match jets, 2 = Anti-match jets
  Int_t fDebugLevel;                 // Amount of debug messages printed to console
  Bool_t fIsRealData;                // Simplification from fDataType: True = any measured data, False = any simulation
  Bool_t fIsPPbData;                 // Flag for any pPb dataset: True = any pPb dataset, False = anything else
  Bool_t fIsPpData;                  // Flag for any pp dataset: True = any pp dataset, False = anything else
  Bool_t fIsOOData;                  // Flag for any OO dataset: True = any OO dataset, False = anything else
  
  // Weights for filling the MC histograms
  Double_t fVzWeight;                // Weight for vz in MC
  Double_t fCentralityWeight;        // Weight for centrality in MC
  Double_t fPtHatWeight;             // Weight for pT hat in MC
  Double_t fTotalEventWeight;        // Combined weight factor for MC
  
  // Jet and track selection cuts
  Int_t fJetAxis;                      // Used jet axis type. 0 = Anti-kT jet axis, 1 = Axis from leading PF candidate
  Double_t fVzCut;                     // Cut for vertez z-position in an event
  Double_t fMinimumPtHat;              // Minimum accepted pT hat value
  Double_t fMaximumPtHat;              // Maximum accepted pT hat value
  Double_t fJetEtaCut;                 // Eta cut around midrapidity
  Double_t fJetMinimumPtCut;           // Minimum pT cut for jets
  Double_t fJetMaximumPtCut;           // Maximum pT accepted for jets (and tracks)
  Bool_t fCutBadPhiRegion;             // Cut the phi region with bad tracker performance from the analysis
  Bool_t fCutWithJetEtaCM;             // Cut symmetric region in center-of-mass jet eta that is within defined jet eta acceptence
  Double_t fMinimumMaxTrackPtFraction; // Cut for jets consisting only from soft particles
  Double_t fMaximumMaxTrackPtFraction; // Cut for jets consisting only from one high pT
  Int_t fJetUncertaintyMode;           // Use uncertainty for jet pT. 0 = Nominal, 1 = Minus uncertainty, 2 = Plus uncertainty
  Double_t fTrackEtaCut;               // Eta cut around midrapidity
  Double_t fTrackMinPtCut;             // Minimum pT cut for tracks
  Double_t fTrackMaxPtCut;             // Maximum pT cut for tracks
  Double_t fMaxTrackPtRelativeError;   // Maximum relative error for pT
  Double_t fMaxTrackDistanceToVertex;  // Maximum distance to primary vetrex
  Double_t fCalorimeterSignalLimitPt;  // Require signal in calorimeters for track above this pT
  Double_t fHighPtEtFraction;          // For high pT tracks, minimum required Et as a fraction of track pT
  Double_t fChi2QualityCut;            // Quality cut for track reconstruction
  Double_t fMinimumTrackHits;          // Quality cut for track hits
  Int_t fSubeventCut;                  // Cut for the subevent index

  // Systematic variations
  Double_t fTrackEfficiencyVariation;   // Relative amount with which the tracking efficiency corrections are varied to estimate systematic uncertainties
  std::string fTrackSelectionVariation; // Automatically determine from track cuts which variation is used, and choose proper efficiency table accordingly
  Int_t fJetPtWeightConfiguration;      // Jet pT weighting configuration in MC

  // Extra jet cuts for jet pT unfolding
  Double_t fReconstructedJetMinimumPtCut;  // Minimum jet pT cut for reconstructed jets for jet pT unfolding
  Double_t fGeneratorJetMinimumPtCut;      // Minimum jet pT cut for generator level jets for jet pT unfolding
  Bool_t fSkipCovarianceMatrix;            // Do not fill the covariance matrix. Speeds up the code significantly.
  
  // Correlation type for Monte Carlo
  Int_t fMcCorrelationType;            // Correlation type for Monte Carlo. See enumeration enumMcCorrelationType
  
  // Configuration for energy-energy correlators
  Double_t fJetRadius;       // Jet radius parameter
  std::vector<Double_t> fWeightExponent;  // Exponents for the energy weight in energy-energy correlators
  Bool_t fSmearDeltaR;       // Flag for smearing DeltaR in energy-energy correlators
  Bool_t fSmearEnergyWeight; // Flag for smearing energy weight in energy-energy correlators
  Bool_t fDoReflectedCone;   // Estimate background from eta-reflected cones
  Bool_t fDoMixedCone;       // Estimate background from cones dropped in mixed events
  Bool_t fMegaSkimMode;      // Use mega skimmed mixing files, that contain bare minimum information for mixing
  Bool_t fDoPerpendicularCone;         // Estimate background from perpendicular cones

  // Event mixing variables
  Int_t fLocalRun;                      // Flag for local vs. CRAB/slurm run
  Int_t fLocalMixing;                   // Mixing is done locally for slurm, via xrootd for CRAB
  Int_t fMixingListIndex;               // Index of the used mixing file list
  Int_t fMixingStartIndex;              // Event index from which event mixing is started
  Int_t fRunningMixingIndex;            // Mixed event index 
  Int_t fnEventsInMixingFile;           // Number of mixed events available
  std::vector<Double_t> fMixedEventVz;           // vz values in mixed events
  std::vector<Int_t> fMixedEventHiBin;           // HiBin values in mixed events
  std::vector<Int_t> fMixedEventMultiplicity;    // Multiplicity in the mixed event
  std::vector<Double_t> fMixedEventHFEnergy;     // HF energy in the mixed event
  std::vector<ULong64_t> fMixedEventEventNumber; // Event number of the mixed event
  std::vector<Int_t> fMixedEventVzBin;           // vz bin indices for event mixing
  std::vector<Int_t> fMixedEventHiBinBin;        // HiBin bin indices for event mixing
  std::vector<Int_t> fMixedEventMultiplicityBin; // Multiplicity bin indices for event mixing
  std::vector<Int_t> fMixedEventHFEnergyBin;     // HF energy bin indices for event mixing
  Int_t fEventMatchForMixing;           // Method used to match events for mixing analysis
  Bool_t fBinnedMixing;                 // Flag for doing event mixing in binned or unbinned way
  std::vector<Double_t> fVzMixingBinBorders;           // Bin borders for vz used for event matching in mixing
  std::vector<Double_t> fHiBinMixingBinBorders;        // Bin borders for HiBin used for event matching in mixing
  std::vector<Double_t> fHFEnergyMixingBinBorders;     // Bin borders for HF energy used for event matching in mixing
  std::vector<Double_t> fMultiplicityMixingBinBorders; // Bin borders for multiplicity used for event matching in mixing
  Double_t fHFEnergyMatchingShift;      // Shift the matched HF energy value to find a reference value without hard interaction
  
  // Which histograms are filled. Do not fill all in order to save memory and not to crash jobs.
  Bool_t fFillEventInformation;                   // Fill event information histograms
  Bool_t fFillJetHistograms;                      // Fill single and dijet histograms
  Bool_t fFillTrackHistograms;                    // Fill inclusive tracks and tracks in dijet events
  Bool_t fFillJetConeHistograms;                  // Fill particle multiplicity and density histograms around the jet cone
  Bool_t fFillEnergyEnergyCorrelators;            // Fill energy-energy correlator histograms
  Bool_t fFillEnergyEnergyCorrelatorsSystematics; // Fill energy-energy correlator histograms for systematic uncertainty analysis
  Bool_t fFillJetPtClosure;                       // Fill jet pT closure histograms
  Bool_t fFillJetPtUnfoldingResponse;             // Fill the jet pT unfolding response
  Bool_t fFillTrackParticleMatchingHistograms;    // Fill the matching histograms between particles and tracks
  
  // Weighting mode flag for MC
  Bool_t fMultiplicityMode;       // True: Weight multiplicity to match the data. False: Weight centrality to match the data

};

#endif

// Class for histograms needed in the dijet analysis

#ifndef EECHISTOGRAMS_H
#define EECHISTOGRAMS_H

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

// Own includes
#include "ConfigurationCard.h"

class EECHistograms{
  
public:
  
  // Enumeration for event types to event histogram and track cuts for track cut histogram
  enum enumEventTypes {kAll, kTriggered, kPrimaryVertex, kHfCoincidence, kClusterCompatibility, kHBHENoise, kBeamScraping, kPileup, kVzCut, kJetOver80, kJetOver120, knEventTypes};
  enum enumTriggerBits {kNoTrigger, kOnlyCaloJet60, kOnlyCaloJet80, kOnlyCaloJet100, kCaloJet60And80, kCaloJet60And100, kCaloJet80And100, kCaloJet60And80And100, knTriggerCombinations};
  enum enumTrackCuts {kAllTracks, kMcCharge, kMcSube, kMcStatus, kPtCuts, kEtaCut, kTrackAlgorithm, kHighPurity, kPtError, kVertexDistance, kCaloSignal, kReconstructionQuality, knTrackCuts};
  enum enumClosureParticleType {kQuark,kGluon,knClosureParticleTypes};
  enum enumJetConeTypes{kSignalCone, kReflectedCone, kMixedCone, kSecondMixedCone, kPerpendicularCone, kSecondPerpendicularCone, knJetConeTypes};
  enum enumPairingType{kSameJetPair, kSignalReflectedConePair, kReflectedConePair, kSignalMixedConePair, kMixedConePair, kSignalSecondMixedConePair, kMixedMixedConePair, kSecondMixedConePair, kSignalPerpendicularConePair, kPerpendicularConePair, kSignalSecondPerpendicularConePair, kPerpendicularPerpendicularConePair, kSecondPerpendicularConePair, knPairingTypes}; // Pair tracks between different cones
  enum enumSubeventTypes{kPythia, kHydjet, knSubeventTypes};
  enum enumSubeventCombinations{kPythiaPythia, kPythiaHydjet, kHydjetPythia, kHydjetHydjet, knSubeventCombinations};
    
  // Constructors and destructor
  EECHistograms(); // Default constructor
  EECHistograms(ConfigurationCard* newCard); // Custom constructor
  EECHistograms(const EECHistograms& in); // Copy constructor
  virtual ~EECHistograms(); // Destructor
  EECHistograms& operator=(const EECHistograms& obj); // Equal sign operator
  
  // Methods
  void CreateHistograms();                   // Create all histograms
  void Write() const;                        // Write the histograms to a file that is opened somewhere else
  void Write(TString outputFileName) const;  // Write the histograms to a file
  void SetCard(ConfigurationCard* newCard);  // Set a new configuration card for the histogram class

  void FillEventCoverage(Int_t hiBin, Int_t stage, Bool_t isOO);
  void FillSelectedEvent(Int_t hiBin, Double_t hfSum, Double_t vz, Double_t eventWeight);
  void MarkTruncatedInput();

  // Getters
  Int_t GetNDeltaRBinsEEC(); // Get the number of deltaR bins for energy-energy correlators
  Double_t GetDeltaRBinBorderLowEEC(Int_t iBin);  // Getter for the low bin border of the i:th deltaR bin
  Double_t GetDeltaRBinBorderHighEEC(Int_t iBin); // Getter for the high bin border of the i:th deltaR bin
  
  // Histograms defined public to allow easier access to them. Should not be abused
  // Notation in comments: l = leading jet, s = subleading jet, inc - inclusive jet, uc = uncorrected, ptw = pT weighted
  TH1F* fhVertexZ;                 // Vertex z-position
  TH1F* fhVertexZWeighted;         // Weighted vertex z-position (only meaningfull for MC)
  TH1F* fhEvents;                  // Number of events. For binning see enumEventTypes.
  TH1D* fhJetVetoEventCounts = nullptr;
  TH2D* fhJetVetoEtaPhiBefore[2] = {nullptr, nullptr}; // signal and mixing; all eligible reco jets
  TH2D* fhJetVetoEtaPhiAfter[2] = {nullptr, nullptr}; // only events passing the full event veto
  TH1F* fhTriggers;                // Number of events selected for each trigger. For binning, see enumTriggerBits
  TH1F* fhTriggersAfterSelection;  // Same as before, but after trigger selection. For binning, see enumTriggerBits
  TH1F* fhTrackCuts;               // Number of tracks. For binning see enumTrackCuts.
  TH1F* fhCentrality;              // Centrality information. -0.5 for pp or PYTHIA.
  TH1F* fhCentralityWeighted;      // Weighted centrality distribution (only meaningful for MC)
  TH1F* fhPtHat;                   // pT hat for MC events (only meaningful for MC)
  TH1F* fhPtHatWeighted;           // Weighted pT hat distribution
  TH2I* fhMixedEventCounts;        // Count of events available for mixing for each matched variable
  THnSparseF* fhHFEnergy;          // HF energy as a function of leading jet pT
  THnSparseF* fhMultiplicity;      // Track multiplicity from all events [multiplicity][centrality]
  THnSparseF* fhInclusiveJet;      // Inclusive jet information. Axes: [jet pT][jet phi][jet eta][cent]
  THnSparseF* fhTrack;             // Track histogram. Axes: [pT][phi][eta][cent]
  THnSparseF* fhTrackUncorrected;  // Track histogram for uncorrected tracks. Axes: [uc pT][uc phi][uc eta][cent]
  THnSparseF* fhParticleDensityAroundJet;   // Particle density around the studied jets
  THnSparseF* fhParticlePtDensityAroundJet; // pT weighted particle density around the studied jets
  THnSparseF* fhParticleMultiplicityInJet;                      // Multiplicity of particles within the studied jet cones
  THnSparseF* fhParticleMultiplicityInReflectedCone;            // Multiplicity of particles within the reflected jet cones
  THnSparseF* fhParticleMultiplicityInJetUncorrected;           // Uncorrected multiplicity of particles within the studied jet cones
  THnSparseF* fhParticleMultiplicityInReflectedConeUncorrected; // Uncorrected multiplicity of particles within the reflected jet cones
  THnSparseF* fhMaxPtParticleInJet;                // Histogram for maximum pT particle found in the jet
  THnSparseF* fhEnergyEnergyCorrelator;                  // Histogram for energy-energy correlator. Axes: [dR][jetPt][trackPt][centrality]
  THnSparseF* fhEnergyEnergyCorrelatorEfficiencyVariationPlus;  // Histogram for energy-energy correlator with track efficiency correction increased by systematic variation amount. Axes: [dR][jetPt][trackPt][centrality]
  THnSparseF* fhEnergyEnergyCorrelatorEfficiencyVariationMinus; // Histogram for energy-energy correlator with track efficiency correction decreased by systematic variation amount. Axes: [dR][jetPt][trackPt][centrality]
  THnSparseF* fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus;  // Histogram for energy-energy correlator with track pair efficiency correction increased by systematic variation amount. Axes: [dR][jetPt][trackPt][centrality]
  THnSparseF* fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus; // Histogram for energy-energy correlator with track pair efficiency correction decreased by systematic variation amount. Axes: [dR][jetPt][trackPt][centrality]
  THnSparseF* fhJetPtClosure; // Jet pT closure histograms. Also information for response matrix. [gen pT][reco pT][centrality][q/g][reco/gen]

  // Quality assurance histograms for the reflected cone
  THnSparseF* fhJetNumberInReflectedCone;  // Number of jets above 25 GeV found in reflected cone
  THnSparseF* fhJetPtInReflectedCone;      // pT of the jets found in reflected cone
  
  // Extra histogram for unfolding study
  THnSparseF* fhUnfoldingMeasured;  // Measured disribution for unfolding
  THnSparseF* fhUnfoldingTruth;     // Truth distribution for unfolding
  THnSparseF* fhUnfoldingResponse;  // Unfolding response matrix

  // Histograms for 1-dimensional jet pT unfolding
  THnSparseF* fhJetPtUnfoldingMeasured;   // Measured jet pT distribution for unfolding test
  THnSparseF* fhJetPtUnfoldingTruth;      // Truth jet pT distribution for unfolding test
  THnSparseF* fhJetPtUnfoldingResponse;   // Jet pT unfolding response matrix
  THnSparseF* fhJetPtUnfoldingCovariance; // Covariance matrix for jet pT unfolding

  // Reconstructed to generator level particle matching
  THnSparseF* fhParticlesCloseToTracks;     // Histogram counting particles close to tracks
  THnSparseF* fhTracksWithMatchedParticle;  // Histogram telling the fraction of particles for which a match is found
  THnSparseF* fhParticleDeltaRResponse;     // Response matrix for particle pair deltaR values
  THnSparseF* fhParticlePtResponse;         // Response matrix for particle pair pT1*pT2 values

  // Generic event-coverage QA, independent of the ordinary centrality output binning.
  TH1D* fhJetExposure; // Sum standard event weights once per selected jet
  TH2D* fhEventCoverage; // [stored activity category][selection stage], unweighted
  TH1D* fhSelectedActivityWeighted; // Accepted standard event weights, Sumw2 retained
  THnSparseD* fhEventActivity; // Actual stored [hiBin][hfSum][vz], including flow bins
  TH1D* fhInputStatus; // Read files; files truncated by explicit pilot limit

private:
  
  ConfigurationCard* fCard;    // Card for binning info
  const TString kEventTypeStrings[knEventTypes] = {"All", "Triggered", "PrimVertex", "HfCoin2Th4", "ClustCompt", "HBHENoise", "BeamScrape", "Pile-up", "v_{z} cut", "Jet>80", "Jet>120"}; // Strings corresponding to event types
  const TString kTriggerStrings[knTriggerCombinations] = {"No trigger", "60 && !80 && !100", "!60 && 80 && !100", "!60 && !80 && 100", "60 && 80 && !100", "60 && !80 && 100", "!60 && 80 && 100", "60 && 80 && 100"};
  const TString kTrackCutStrings[knTrackCuts] = {"All", "MC Charge", "MC sube", "MC status", "p_{T} cut", "#eta cut", "Track algo", "HighPurity", "p_{T} error", "vertexDist", "caloSignal", "RecoQuality"}; // String corresponding to track cuts
  
};

#endif

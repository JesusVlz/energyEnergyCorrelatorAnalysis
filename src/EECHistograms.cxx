// Histograms needed in the dijet analysis

// C++ includes
#include <assert.h>

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "EECHistograms.h"

/*
 * Default constructor
 */
EECHistograms::EECHistograms() :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhTriggers(0),
  fhTriggersAfterSelection(0),
  fhTrackCuts(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhMixedEventCounts(0),
  fhHFEnergy(0),
  fhMultiplicity(0),
  fhInclusiveJet(0),
  fhTrack(0),
  fhTrackUncorrected(0),
  fhParticleDensityAroundJet(0),
  fhParticlePtDensityAroundJet(0),
  fhParticleMultiplicityInJet(0),
  fhParticleMultiplicityInReflectedCone(0),
  fhParticleMultiplicityInJetUncorrected(0),
  fhParticleMultiplicityInReflectedConeUncorrected(0),
  fhMaxPtParticleInJet(0),
  fhEnergyEnergyCorrelator(0),
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus(0),
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus(0), 
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus(0), 
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus(0),
  fhJetPtClosure(0),
  fhJetNumberInReflectedCone(0),
  fhJetPtInReflectedCone(0),
  fhUnfoldingMeasured(0),
  fhUnfoldingTruth(0),
  fhUnfoldingResponse(0),
  fhJetPtUnfoldingMeasured(0),
  fhJetPtUnfoldingTruth(0),
  fhJetPtUnfoldingResponse(0),
  fhJetPtUnfoldingCovariance(0),
  fhParticlesCloseToTracks(0),
  fhTracksWithMatchedParticle(0),
  fhParticleDeltaRResponse(0),
  fhParticlePtResponse(0),
  fCard(0)
{
  // Default constructor
  
}

/*
 * Custom constructor
 */
EECHistograms::EECHistograms(ConfigurationCard* newCard) :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhTriggers(0),
  fhTriggersAfterSelection(0),
  fhTrackCuts(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhMixedEventCounts(0),
  fhHFEnergy(0),
  fhMultiplicity(0),
  fhInclusiveJet(0),
  fhTrack(0),
  fhTrackUncorrected(0),
  fhParticleDensityAroundJet(0),
  fhParticlePtDensityAroundJet(0),
  fhParticleMultiplicityInJet(0),
  fhParticleMultiplicityInReflectedCone(0),
  fhParticleMultiplicityInJetUncorrected(0),
  fhParticleMultiplicityInReflectedConeUncorrected(0),
  fhMaxPtParticleInJet(0),
  fhEnergyEnergyCorrelator(0),
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus(0),
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus(0),
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus(0), 
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus(0),  
  fhJetPtClosure(0),
  fhJetNumberInReflectedCone(0),
  fhJetPtInReflectedCone(0),
  fhUnfoldingMeasured(0),
  fhUnfoldingTruth(0),
  fhUnfoldingResponse(0),
  fhJetPtUnfoldingMeasured(0),
  fhJetPtUnfoldingTruth(0),
  fhJetPtUnfoldingResponse(0),
  fhJetPtUnfoldingCovariance(0),
  fhParticlesCloseToTracks(0),
  fhTracksWithMatchedParticle(0),
  fhParticleDeltaRResponse(0),
  fhParticlePtResponse(0),
  fCard(newCard)
{
  // Custom constructor

}

/*
 * Copy constructor
 */
EECHistograms::EECHistograms(const EECHistograms& in) :
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhTriggers(in.fhTriggers),
  fhTriggersAfterSelection(in.fhTriggersAfterSelection),
  fhTrackCuts(in.fhTrackCuts),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted),
  fhPtHat(in.fhPtHat),
  fhPtHatWeighted(in.fhPtHatWeighted),
  fhMixedEventCounts(in.fhMixedEventCounts),
  fhHFEnergy(in.fhHFEnergy),
  fhMultiplicity(in.fhMultiplicity),
  fhInclusiveJet(in.fhInclusiveJet),
  fhTrack(in.fhTrack),
  fhTrackUncorrected(in.fhTrackUncorrected),
  fhParticleDensityAroundJet(in.fhParticleDensityAroundJet),
  fhParticlePtDensityAroundJet(in.fhParticlePtDensityAroundJet),
  fhParticleMultiplicityInJet(in.fhParticleMultiplicityInJet),
  fhParticleMultiplicityInReflectedCone(in.fhParticleMultiplicityInReflectedCone),
  fhParticleMultiplicityInJetUncorrected(in.fhParticleMultiplicityInJetUncorrected),
  fhParticleMultiplicityInReflectedConeUncorrected(in.fhParticleMultiplicityInReflectedConeUncorrected),
  fhMaxPtParticleInJet(in.fhMaxPtParticleInJet),
  fhEnergyEnergyCorrelator(in.fhEnergyEnergyCorrelator),
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus(in.fhEnergyEnergyCorrelatorEfficiencyVariationPlus),
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus(in.fhEnergyEnergyCorrelatorEfficiencyVariationMinus),
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus(in.fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus), 
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus(in.fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus),
  fhJetPtClosure(in.fhJetPtClosure),
  fhJetNumberInReflectedCone(in.fhJetNumberInReflectedCone),
  fhJetPtInReflectedCone(in.fhJetPtInReflectedCone),
  fhUnfoldingMeasured(in.fhUnfoldingMeasured),
  fhUnfoldingTruth(in.fhUnfoldingTruth),
  fhUnfoldingResponse(in.fhUnfoldingResponse),
  fhJetPtUnfoldingMeasured(in.fhJetPtUnfoldingMeasured),
  fhJetPtUnfoldingTruth(in.fhJetPtUnfoldingTruth),
  fhJetPtUnfoldingResponse(in.fhJetPtUnfoldingResponse),
  fhJetPtUnfoldingCovariance(in.fhJetPtUnfoldingCovariance),
  fhParticlesCloseToTracks(in.fhParticlesCloseToTracks),
  fhTracksWithMatchedParticle(in.fhTracksWithMatchedParticle),
  fhParticleDeltaRResponse(in.fhParticleDeltaRResponse),
  fhParticlePtResponse(in.fhParticlePtResponse),
  fCard(in.fCard)
{
  // Copy constructor

}

/*
 * Assingment operator
 */
EECHistograms& EECHistograms::operator=(const EECHistograms& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fhVertexZ = in.fhVertexZ;
  fhVertexZWeighted = in.fhVertexZWeighted;
  fhEvents = in.fhEvents;
  fhTriggers = in.fhTriggers;
  fhTriggersAfterSelection = in.fhTriggersAfterSelection;
  fhTrackCuts = in.fhTrackCuts;
  fhCentrality = in.fhCentrality;
  fhCentralityWeighted = in.fhCentralityWeighted;
  fhPtHat = in.fhPtHat;
  fhPtHatWeighted = in.fhPtHatWeighted;
  fhMixedEventCounts = in.fhMixedEventCounts;
  fhHFEnergy = in.fhHFEnergy;
  fhMultiplicity = in.fhMultiplicity;
  fhInclusiveJet = in.fhInclusiveJet;
  fhTrack = in.fhTrack;
  fhTrackUncorrected = in.fhTrackUncorrected;
  fhParticleDensityAroundJet = in.fhParticleDensityAroundJet;
  fhParticlePtDensityAroundJet = in.fhParticlePtDensityAroundJet;
  fhParticleMultiplicityInJet = in.fhParticleMultiplicityInJet;
  fhParticleMultiplicityInReflectedCone = in.fhParticleMultiplicityInReflectedCone;
  fhParticleMultiplicityInJetUncorrected = in.fhParticleMultiplicityInJetUncorrected;
  fhParticleMultiplicityInReflectedConeUncorrected = in.fhParticleMultiplicityInReflectedConeUncorrected;
  fhMaxPtParticleInJet = in.fhMaxPtParticleInJet;
  fhEnergyEnergyCorrelator = in.fhEnergyEnergyCorrelator;
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus = in.fhEnergyEnergyCorrelatorEfficiencyVariationPlus;
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus = in.fhEnergyEnergyCorrelatorEfficiencyVariationMinus;
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus = in.fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus;
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus = in.fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus;
  fhJetPtClosure = in.fhJetPtClosure;
  fhJetNumberInReflectedCone = in.fhJetNumberInReflectedCone;
  fhJetPtInReflectedCone = in.fhJetPtInReflectedCone;
  fhUnfoldingMeasured = in.fhUnfoldingMeasured;
  fhUnfoldingTruth = in.fhUnfoldingTruth;
  fhUnfoldingResponse = in.fhUnfoldingResponse;
  fhJetPtUnfoldingMeasured = in.fhJetPtUnfoldingMeasured;
  fhJetPtUnfoldingTruth = in.fhJetPtUnfoldingTruth;
  fhJetPtUnfoldingResponse = in.fhJetPtUnfoldingResponse;
  fhJetPtUnfoldingCovariance = in.fhJetPtUnfoldingCovariance;
  fhParticlesCloseToTracks = in.fhParticlesCloseToTracks;
  fhTracksWithMatchedParticle = in.fhTracksWithMatchedParticle;
  fhParticleDeltaRResponse = in.fhParticleDeltaRResponse;
  fhParticlePtResponse = in.fhParticlePtResponse;
  fCard = in.fCard;
  
  return *this;
}

/*
 * Destructor
 */
EECHistograms::~EECHistograms(){
  // destructor
  delete fhVertexZ;
  delete fhVertexZWeighted;
  delete fhEvents;
  delete fhTriggers;
  delete fhTriggersAfterSelection;
  delete fhTrackCuts;
  delete fhCentrality;
  delete fhCentralityWeighted;
  delete fhPtHat;
  delete fhPtHatWeighted;
  delete fhMixedEventCounts;
  delete fhHFEnergy;
  delete fhMultiplicity;
  delete fhInclusiveJet;
  delete fhTrack;
  delete fhTrackUncorrected;
  delete fhParticleDensityAroundJet;
  delete fhParticlePtDensityAroundJet;
  delete fhParticleMultiplicityInJet;
  delete fhParticleMultiplicityInReflectedCone;
  delete fhParticleMultiplicityInJetUncorrected;
  delete fhParticleMultiplicityInReflectedConeUncorrected;
  delete fhMaxPtParticleInJet;
  delete fhEnergyEnergyCorrelator;
  delete fhEnergyEnergyCorrelatorEfficiencyVariationPlus;
  delete fhEnergyEnergyCorrelatorEfficiencyVariationMinus;
  delete fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus;
  delete fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus;
  delete fhJetPtClosure;
  delete fhJetNumberInReflectedCone;
  delete fhJetPtInReflectedCone;
  delete fhUnfoldingMeasured;
  delete fhUnfoldingTruth;
  delete fhUnfoldingResponse;
  delete fhJetPtUnfoldingMeasured;
  delete fhJetPtUnfoldingTruth;
  delete fhJetPtUnfoldingResponse;
  delete fhJetPtUnfoldingCovariance;
  delete fhParticlesCloseToTracks;
  delete fhTracksWithMatchedParticle;
  delete fhParticleDeltaRResponse;
  delete fhParticlePtResponse;
}

/*
 * Set the configuration card used for the histogram class
 */
void EECHistograms::SetCard(ConfigurationCard* newCard){
  fCard = newCard;
}

/*
 * Create the necessary histograms
 */
void EECHistograms::CreateHistograms(){
  
  // ======== Common binning information for histograms =========
  
  // Centrality
  const Double_t minCentrality = -1;   // Minimum centrality bin is negative since hiBin is -0.5 for pp
  const Double_t maxCentrality = 100;  // Maximum centrality bin
  const Int_t nCentralityBins = 202;      // Number of centrality bins
  
  // Jet pT
  const Double_t minPtJet = 0;     // Minimum jet pT
  const Double_t maxPtJet = 500;   // Maximum jet pT
  const Int_t nPtBinsJet = 100;    // Number of jet pT bins
  
  //Track pT
  const Double_t minPtTrack = 0;   // Minimum track pT for track histograms
  const Double_t maxPtTrack = 50;  // Maximum track pT for track histograms (Hallie's analysis = 20)
  const Int_t nPtBinsTrack = 500;  // Number of track pT bins for track histograms (Hallie's analysis = 500)
  
  // Track pT for maximum particle pT inside a jet
  const Double_t minParticlePt = 0;
  const Double_t maxParticlePt = 100;
  const Int_t nPtBinsParticle = 100;
  
  // Phi
  const Double_t minPhi = -TMath::Pi();  // Minimum phi
  const Double_t maxPhi = TMath::Pi();   // Maximum phi
  const Int_t nPhiBins = 64;             // Number of phi bins
  
  // Eta
  const Double_t minEta = -2.5;    // Minimum eta (current eta cut for tracks = 2.4)
  const Double_t maxEta = 2.5;     // Maximum eta (current eta cut for tracks = 2.4)
  const Int_t nEtaBins = 50;       // Number of eta bins
  
  // Vertex z-position
  const Double_t minVz = -20;   // Minimum vz
  const Double_t maxVz = 20;    // Maximum vz
  const Int_t nVzBins = 80;     // Number of vz bins
  
  // pT hat
  const Double_t minPtHat = 0;     // Minimum pT hat
  const Double_t maxPtHat = 460;   // Maximum pT hat
  const Int_t nFinePtHatBins = 230; // Number of fine pT hat bins
  
  // Generator level pT binning for closure histograms
  const Double_t minClosurePt = 50;                             // Minimum gen jet pT for closure plots
  const Double_t maxClosurePt = 500;                            // Maximum gen jet pT for closure plots
  const Int_t nClosurePtBins = (maxClosurePt-minClosurePt)/10;  // Bin width of 10 for the Gen pT in closure plots
  
  // Particle type for closure plots (0 = quark, 1 = gluon)
  const Double_t minClosureParticleType = -0.5;                        // Closure particle type indexing starts from zero
  const Double_t maxClosureParticleType = knClosureParticleTypes-0.5;  // Maximum closure particle type index
  const Int_t nClosureParticleTypeBins = knClosureParticleTypes;       // Bin width for particle type is 1
  
  // Track pairing type for energy-energy correlators (0 = same jet, 1 = signal-reflected cone 2 = only eta-reflected jet)
  const Double_t minTrackPairingType = -0.5;                // Track pairing type indexing starts from zero
  const Double_t maxTrackPairingType = knPairingTypes-0.5;  // Maximum track pairing type index
  const Int_t nTrackPairingTypeBins = knPairingTypes;       // Bin width for track pairing type is 1
  
  // Subevent index for particles (0 = Particle is from embedded Pythia, 1 = Particle is from Hydjet)
  const Double_t minSubeventType = -0.5;
  const Double_t maxSubeventType = knSubeventTypes-0.5;
  const Int_t nSubeventTypeBins = knSubeventTypes;
  
  // Subevent pairing index for particles (0 = both paritcles from Pythia, 1 = one particle from Pythia, one from Hydjet, 2 = Both particles from Hydjet)
  const Double_t minSubeventCombination = -0.5;
  const Double_t maxSubeventCombination = knSubeventCombinations-0.5;
  const Int_t nSubeventCombinationBins = knSubeventCombinations;
  
  // Jet cone type (0 = signal cone, 1 = reflected cone)
  const Double_t minJetConeType = -0.5;
  const Double_t maxJetConeType = knJetConeTypes-0.5;
  const Int_t nJetConeTypes = knJetConeTypes;
  
  // Binning for reco/gen ratio for closure histograms
  const Double_t minClosureRatio = 0.025;    // Minimum ratio for the closure plots
  const Double_t maxClosureRatio = 2.025;    // Maximum ratio for the closure plots
  const Int_t nClosureRatioBins = 40;    // Number of closure ratio bins
  
  // Binning for multiplicity
  const Double_t minMultiplicity = 0;
  const Double_t maxMultiplicity = 5000;
  const Double_t maxMultiplicityWeighted = 4000;
  const Int_t nMultiplicityBins = 500;
  const Int_t nMultiplicityBinsWeighted = 400;
  
  // Binning for multiplicity in within jet cone
  const Double_t minMultiplicityInJets = -0.5;
  const Double_t maxMultiplicityInJets = 149.5;
  const Int_t nMultiplicityInJetsBins = 150;

  // Number of particles close to a track
  const Double_t minNumberOfTracks = -0.5;
  const Double_t maxNumberOfTracks = 9.5;
  const Int_t nNumberOftracksBins = 10;

  // Bins for pT1*pT2
  const Double_t minPt1TimesPt2 = 4;
  const Double_t maxPt1TimesPt2 = 90000;
  const Int_t nPt1TimesPt2Bins = 24;
  const Double_t pt1TimesPt2Bins[nPt1TimesPt2Bins+1] = {4,6,9,12,16,20,25,30,35,40,45,50,60,70,80,90,100,120,140,180,240,400,600,2000,90000};

  // Bins for jet pT within the reflected cone
  const Double_t minJetPtReflectedCone = 25;
  const Double_t maxJetPtReflectedCone = 5020;
  const Int_t nJetPtInReflectedConeBins = 5;
  const Double_t jetPtReflectedConeBins[nJetPtInReflectedConeBins+1] = {25, 40, 60, 80, 100, 5020};

  const Double_t minNumberOfJetsReflectedCone = -0.5;
  const Double_t maxNumberOfJetsReflectedCone = 4.5;
  const Int_t nNumberOfJetsReflectedConeBins = 5;

  // Bins for HF energy
  const Double_t minHFenergy = 0;
  const Double_t maxHFenergy = 300;
  const Int_t nHFenergyBins = 60;

  // Bins for number of jets in the collision
  const Double_t minJets = -0.5;
  const Double_t maxJets = 9.5;
  const Int_t nJetNumberBins = 10;
  
  // Centrality bins for THnSparses (We run into memory issues, if have all the bins)
  const Int_t nWideCentralityBins = fCard->GetNBin("CentralityBinEdges");
  Double_t* wideCentralityBins = new Double_t[nWideCentralityBins+1];
  for(Int_t iCentrality = 0; iCentrality < nWideCentralityBins+1; iCentrality++){
    wideCentralityBins[iCentrality] = fCard->Get("CentralityBinEdges",iCentrality);
  }
  
  // Bins for the pT hat histogram
  const Int_t nPtHatBins = fCard->GetNBin("PtHatBinEdges");
  Double_t* ptHatBins = new Double_t[nPtHatBins+1];
  for(Int_t iPtHat = 0; iPtHat < nPtHatBins+1; iPtHat++){
    ptHatBins[iPtHat] = fCard->Get("PtHatBinEdges",iPtHat);
  }
  
  // Jet pT binning for energy-energy correlator histograms
  const Int_t nJetPtBinsEEC = fCard->GetNBin("JetPtBinEdgesEEC");
  Double_t* jetPtBinsEEC = new Double_t[nJetPtBinsEEC+1];
  for(Int_t iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
    jetPtBinsEEC[iJetPt] = fCard->Get("JetPtBinEdgesEEC",iJetPt);
  }
  const Double_t minJetPtEEC = jetPtBinsEEC[0];
  const Double_t maxJetPtEEC = jetPtBinsEEC[nJetPtBinsEEC];
  
  // Track pT binning for energy-energy correlator histograms
  const Int_t nTrackPtBinsEEC = fCard->GetNBin("TrackPtBinEdgesEEC");
  Double_t* trackPtBinsEEC = new Double_t[nTrackPtBinsEEC+1];
  for(Int_t iTrackPt = 0; iTrackPt < nTrackPtBinsEEC+1; iTrackPt++){
    trackPtBinsEEC[iTrackPt] = fCard->Get("TrackPtBinEdgesEEC",iTrackPt);
  }
  const Double_t minTrackPtEEC = trackPtBinsEEC[0];
  const Double_t maxTrackPtEEC = trackPtBinsEEC[nTrackPtBinsEEC];

  // Create a set of energy-energy correlator histograms for each defined energy weight
  const Int_t nEnergyWeightsEEC = fCard->GetN("WeightExponent");
  const Double_t minEnergyWeightEEC = -0.5;
  const Double_t maxEnergyWeightEEC = nEnergyWeightsEEC - 0.5;

  // DeltaR binning for track density. Range must cover the full cone diameter (2*R),
  // since two tracks inside a cone of radius R can be separated by up to 2R.
  const Int_t nDeltaRBinsTrackDensity = 80;
  const Double_t minDeltaRTrackDensity = 0;
  const Double_t maxDeltaRTrackDensity = 2*fCard->Get("JetRadius");
  
  // Logarithmic deltaR binning for energy-energy correlator histograms. Range must cover
  // the full cone diameter (2*R) for the same reason as above.
  const Int_t nDeltaRBinsEEC = 32;
  const Double_t minDeltaREEC = 0;
  const Double_t maxDeltaREEC = 2*fCard->Get("JetRadius");
  const Double_t binnerShift = 0.01;
  const Double_t deltaRlogBinWidth = (TMath::Log(maxDeltaREEC+binnerShift) - TMath::Log(minDeltaREEC+binnerShift)) / nDeltaRBinsEEC;
  Double_t deltaRBinsEEC[nDeltaRBinsEEC+1];
  for(int iDeltaR = 0; iDeltaR <= nDeltaRBinsEEC; iDeltaR++){
    deltaRBinsEEC[iDeltaR] = (minDeltaREEC+binnerShift)*TMath::Exp(iDeltaR*deltaRlogBinWidth)-binnerShift;
  }

  // Binning for the two-dimensional unfolding response
  const Int_t nUnfoldingBins = nDeltaRBinsEEC*nJetPtBinsEEC;
  const Double_t minUnfoldingBin = 0;
  const Double_t maxUnfoldingBin = nJetPtBinsEEC*maxDeltaREEC;
  Double_t* fullUnfoldingBinning = new Double_t[nUnfoldingBins+1];
  for(int iDeltaR = 0; iDeltaR < nDeltaRBinsEEC; iDeltaR++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      fullUnfoldingBinning[iDeltaR+iJetPt*nDeltaRBinsEEC] = deltaRBinsEEC[iDeltaR]+maxDeltaREEC*iJetPt;
    }
  }
  fullUnfoldingBinning[nUnfoldingBins] = maxUnfoldingBin;

  // Arrays for creating THnSparses
  const Int_t nAxesMultiplicity = 3;
  Int_t nBinsMultiplicity[nAxesMultiplicity];
  Double_t lowBinBorderMultiplicity[nAxesMultiplicity];
  Double_t highBinBorderMultiplicity[nAxesMultiplicity];
  
  const Int_t nAxesMultiplicityInJets = 5;
  Int_t nBinsMultiplicityInJets[nAxesMultiplicityInJets];
  Double_t lowBinBorderMultiplicityInJets[nAxesMultiplicityInJets];
  Double_t highBinBorderMultiplicityInJets[nAxesMultiplicityInJets];
  
  const Int_t nAxesJet = 6;
  Int_t nBinsJet[nAxesJet];
  Double_t lowBinBorderJet[nAxesJet];
  Double_t highBinBorderJet[nAxesJet];
  
  const Int_t nAxesTrack = 4;
  Int_t nBinsTrack[nAxesTrack];
  Double_t lowBinBorderTrack[nAxesTrack];
  Double_t highBinBorderTrack[nAxesTrack];
  
  const Int_t nAxesTrackDensity = 6;
  Int_t nBinsTrackDensity[nAxesTrackDensity];
  Double_t lowBinBorderTrackDensity[nAxesTrackDensity];
  Double_t highBinBorderTrackDensity[nAxesTrackDensity];
  
  const Int_t nAxesMaxParticlePtInJet = 4;
  Int_t nBinsMaxParticlePtInJet[nAxesMaxParticlePtInJet];
  Double_t lowBinBorderMaxParticlePtInJet[nAxesMaxParticlePtInJet];
  Double_t highBinBorderMaxParticlePtInJet[nAxesMaxParticlePtInJet];
  
  const Int_t nAxesEnergyEnergyCorrelator = 7;
  Int_t nBinsEnergyEnergyCorrelator[nAxesEnergyEnergyCorrelator];
  Double_t lowBinBorderEnergyEnergyCorrelator[nAxesEnergyEnergyCorrelator];
  Double_t highBinBorderEnergyEnergyCorrelator[nAxesEnergyEnergyCorrelator];
  
  const Int_t nAxesJetClosure = 7;
  Int_t nBinsJetClosure[nAxesJetClosure];
  Double_t lowBinBorderJetClosure[nAxesJetClosure];
  Double_t highBinBorderJetClosure[nAxesJetClosure];

  const Int_t nAxesJetPtUnfoldDistribution = 4;
  Int_t nBinsJetPtUnfoldDistribution[nAxesJetPtUnfoldDistribution];
  Double_t lowBinBorderJetPtUnfoldDistribution[nAxesJetPtUnfoldDistribution];
  Double_t highBinBorderJetPtUnfoldDistribution[nAxesJetPtUnfoldDistribution];

  const Int_t nAxesJetPtUnfoldResponse = 5;
  Int_t nBinsJetPtUnfoldResponse[nAxesJetPtUnfoldResponse];
  Double_t lowBinBorderJetPtUnfoldResponse[nAxesJetPtUnfoldResponse];
  Double_t highBinBorderJetPtUnfoldResponse[nAxesJetPtUnfoldResponse];

  const Int_t nAxesParticleMatchingQA = 4;
  Int_t nBinsParticleMatchingQA[nAxesParticleMatchingQA];
  Double_t lowBinBorderParticleMatchingQA[nAxesParticleMatchingQA];
  Double_t highBinBorderParticleMatchingQA[nAxesParticleMatchingQA];

  const Int_t nAxesParticleDeltaRResponseMatrix = 5;
  Int_t nBinsParticleDeltaRResponseMatrix[nAxesParticleDeltaRResponseMatrix];
  Double_t lowBinBorderParticleDeltaRResponseMatrix[nAxesParticleDeltaRResponseMatrix];
  Double_t highBinBorderParticleDeltaRResponseMatrix[nAxesParticleDeltaRResponseMatrix];

  const Int_t nAxesParticlePtResponseMatrix = 7;
  Int_t nBinsParticlePtResponseMatrix[nAxesParticlePtResponseMatrix];
  Double_t lowBinBorderParticlePtResponseMatrix[nAxesParticlePtResponseMatrix];
  Double_t highBinBorderParticlePtResponseMatrix[nAxesParticlePtResponseMatrix];

  const Int_t nAxesHFEnergy = 5;
  Int_t nBinsHFEnergy[nAxesHFEnergy];
  Double_t lowBinBorderHFEnergy[nAxesHFEnergy];
  Double_t highBinBorderHFEnergy[nAxesHFEnergy];

  const Int_t nAxesReflectedConeQA = 2;
  Int_t nBinsReflectedConeQA[nAxesReflectedConeQA];
  Double_t lowBinBorderReflectedConeQA[nAxesReflectedConeQA];
  Double_t highBinBorderReflectedConeQA[nAxesReflectedConeQA];

  const Int_t nAxesOneDimensionalJetPtUnfoldDistribution = 2;
  Int_t nBinsOneDimensionalJetPtUnfoldDistribution[nAxesOneDimensionalJetPtUnfoldDistribution];
  Double_t lowBinBorderOneDimensionalJetPtUnfoldDistribution[nAxesOneDimensionalJetPtUnfoldDistribution];
  Double_t highBinBorderOneDimensionalJetPtUnfoldDistribution[nAxesOneDimensionalJetPtUnfoldDistribution];

  const Int_t nAxesOneDimensionalJetPtUnfoldResponse = 3;
  Int_t nBinsOneDimensionalJetPtUnfoldResponse[nAxesOneDimensionalJetPtUnfoldResponse];
  Double_t lowBinBorderOneDimensionalJetPtUnfoldResponse[nAxesOneDimensionalJetPtUnfoldResponse];
  Double_t highBinBorderOneDimensionalJetPtUnfoldResponse[nAxesOneDimensionalJetPtUnfoldResponse];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1F("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhVertexZWeighted = new TH1F("vertexZweighted","vertexZweighted",nVzBins,minVz,maxVz); fhVertexZWeighted->Sumw2();
  fhEvents = new TH1F("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhTriggers = new TH1F("triggers","triggers",knTriggerCombinations,-0.5,knTriggerCombinations-0.5); fhTriggers->Sumw2();
  fhTriggersAfterSelection = new TH1F("triggersAfterSelection","triggersAfterSelection",knTriggerCombinations,-0.5,knTriggerCombinations-0.5); fhTriggersAfterSelection->Sumw2();
  fhTrackCuts = new TH1F("trackCuts","trackCuts",knTrackCuts,-0.5,knTrackCuts-0.5); fhTrackCuts->Sumw2();
  fhCentrality = new TH1F("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  fhCentralityWeighted = new TH1F("centralityWeighted","centralityWeighted",nCentralityBins,minCentrality,maxCentrality); fhCentralityWeighted->Sumw2();
  fhPtHat = new TH1F("pthat","pthat",nPtHatBins,ptHatBins); fhPtHat->Sumw2();
  fhPtHatWeighted = new TH1F("pthatWeighted","pthatWeighted",nFinePtHatBins,minPtHat,maxPtHat); fhPtHatWeighted->Sumw2();
  fhMixedEventCounts = new TH2I("mixedEventCounts", "mixedEventCounts", 1, 0, 1, 1, 0, 1);
  
  // For the event histogram, label each bin corresponding to an event cut
  for(Int_t i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // For the trigger histogram, label each bin corresponding to a trigger selection
  for(Int_t i = 0; i < knTriggerCombinations; i++){
    fhTriggers->GetXaxis()->SetBinLabel(i+1,kTriggerStrings[i]);
    fhTriggersAfterSelection->GetXaxis()->SetBinLabel(i+1,kTriggerStrings[i]);
  }
  
  // For the track cut histogram, label each bin corresponding to a track cut
  for(Int_t i = 0; i < knTrackCuts; i++){
    fhTrackCuts->GetXaxis()->SetBinLabel(i+1,kTrackCutStrings[i]);
  }

  // ======== THnSparses for HF energy ========
  
  // Axis 0 for the HF energy histogram: HF plus
  nBinsHFEnergy[0] = nHFenergyBins;          // nBins for energy in HF plus
  lowBinBorderHFEnergy[0] = minHFenergy;     // low bin border for energy in HF plus
  highBinBorderHFEnergy[0] = maxHFenergy;    // high bin border for energy in HF plus
  
  // Axis 1 for the HF energy histogram: HF minus
  nBinsHFEnergy[1] = nHFenergyBins;          // nBins for energy in HF minus
  lowBinBorderHFEnergy[1] = minHFenergy;     // low bin border for energy in HF minus
  highBinBorderHFEnergy[1] = maxHFenergy;    // high bin border for energy in HF minus
  
  // Axis 2 for the HF energy histogram: HF sum
  nBinsHFEnergy[2] = nHFenergyBins*2;        // nBins for energy in HF calorimeters
  lowBinBorderHFEnergy[2] = minHFenergy;     // low bin border for energy in HF calorimeters
  highBinBorderHFEnergy[2] = maxHFenergy*2;  // high bin border for energy in HF calorimeters

  // Axis 3 for the HF energy histogram: leading jet pT
  nBinsHFEnergy[3] = nJetPtBinsEEC;          // nBins for leading jet pT
  lowBinBorderHFEnergy[3] = minJetPtEEC;     // low bin border for leading jet pT
  highBinBorderHFEnergy[3] = maxJetPtEEC;    // high bin border for leading jet pT

  // Axis 4 for the HF energy histogram: number of jets in the event
  nBinsHFEnergy[4] = nJetNumberBins;         // nBins for leading jet pT
  lowBinBorderHFEnergy[4] = minJets;         // low bin border for leading jet pT
  highBinBorderHFEnergy[4] = maxJets;        // high bin border for leading jet pT
  
  // Create the HF energy histograms using the above binning information
  fhHFEnergy = new THnSparseF("hfEnergy", "hfEnergy", nAxesHFEnergy, nBinsHFEnergy, lowBinBorderHFEnergy, highBinBorderHFEnergy); fhHFEnergy->Sumw2();

  // Custom bin borders for leading jet pT
  fhHFEnergy->SetBinEdges(3,jetPtBinsEEC);   // Jet pT bins
  
  // ======== THnSparses for multiplicity ========
  
  // Axis 0 for the multiplicity histogram: multiplicity
  nBinsMultiplicity[0] = nMultiplicityBins;       // nBins for multiplicity
  lowBinBorderMultiplicity[0] = minMultiplicity;  // low bin border for multiplicity
  highBinBorderMultiplicity[0] = maxMultiplicity; // high bin border for multiplicity
  
  // Axis 1 for the multiplicity histogram: weighted multiplicity
  nBinsMultiplicity[1] = nMultiplicityBinsWeighted;       // nBins for weighted multiplicity
  lowBinBorderMultiplicity[1] = minMultiplicity;          // low bin border for weighted multiplicity
  highBinBorderMultiplicity[1] = maxMultiplicityWeighted; // high bin border for weighted multiplicity
  
  // Axis 2 for the multiplicity histogram: centrality
  nBinsMultiplicity[2] = nCentralityBins;     // nBins for wide centrality bins
  lowBinBorderMultiplicity[2] = minCentrality;    // low bin border for centrality
  highBinBorderMultiplicity[2] = maxCentrality;   // high bin border for centrality
  
  // Create the multiplicity histograms using the above binning information
  fhMultiplicity = new THnSparseF("multiplicity", "multiplicity", nAxesMultiplicity, nBinsMultiplicity, lowBinBorderMultiplicity, highBinBorderMultiplicity); fhMultiplicity->Sumw2();
  
  // ======== THnSparse for multiplicity within the jet cone ========
  
  // Axis 0 for the multiplicity in jets histogram: multiplicity
  nBinsMultiplicityInJets[0] = nMultiplicityInJetsBins;       // nBins for non-corrected multiplicity
  lowBinBorderMultiplicityInJets[0] = minMultiplicityInJets;  // low bin border for non-corrected multiplicity
  highBinBorderMultiplicityInJets[0] = maxMultiplicityInJets; // high bin border for non-corrected multiplicity
  
  // Axis 1 for the multiplicity in jets histogram: jet pT
  nBinsMultiplicityInJets[1] = nJetPtBinsEEC;         // nBins for jet pT in energy-energy correlator histograms
  lowBinBorderMultiplicityInJets[1] = minJetPtEEC;    // low bin border for jet pT in energy-energy correlator histograms
  highBinBorderMultiplicityInJets[1] = maxJetPtEEC;   // high bin border for jet pT in energy-energy correlator histograms
  
  // Axis 2 for the multiplicity in jets histogram: track pT
  nBinsMultiplicityInJets[2] = nTrackPtBinsEEC;        // nBins for track pT in energy-energy correlator histograms
  lowBinBorderMultiplicityInJets[2] = minTrackPtEEC;   // low bin border for track pT in energy-energy correlator histograms
  highBinBorderMultiplicityInJets[2] = maxTrackPtEEC;  // high bin border for track pT in energy-energy correlator histograms
  
  // Axis 3 for the multiplicity in jets histogram: centrality
  nBinsMultiplicityInJets[3] = nWideCentralityBins;    // nBins for wide centrality bins
  lowBinBorderMultiplicityInJets[3] = minCentrality;   // low bin border for centrality
  highBinBorderMultiplicityInJets[3] = maxCentrality;  // high bin border for centrality
  
  // Axis 4 for the multiplicity in jets histogram: track subevent (only relevant for Monte Carlo)
  nBinsMultiplicityInJets[4] = nSubeventTypeBins+1;    // nBins for subevent types (Pythia/Hydjet/Combined)
  lowBinBorderMultiplicityInJets[4] = minSubeventType;       // low bin border for subevent types
  highBinBorderMultiplicityInJets[4] = maxSubeventType+1;    // high bin border for subevent types
  
  // Create the multiplicity in jets histrograms using the above binning information
  fhParticleMultiplicityInJet = new THnSparseF("multiplicityInJetCone", "multiplicityInJetCone", nAxesMultiplicityInJets, nBinsMultiplicityInJets, lowBinBorderMultiplicityInJets, highBinBorderMultiplicityInJets); fhParticleMultiplicityInJet->Sumw2();
  fhParticleMultiplicityInReflectedCone = new THnSparseF("multiplicityInReflectedCone", "multiplicityInReflectedCone", nAxesMultiplicityInJets, nBinsMultiplicityInJets, lowBinBorderMultiplicityInJets, highBinBorderMultiplicityInJets); fhParticleMultiplicityInReflectedCone->Sumw2();
  fhParticleMultiplicityInJetUncorrected = new THnSparseF("multiplicityInJetConeUncorrected", "multiplicityInJetConeUncorrected", nAxesMultiplicityInJets, nBinsMultiplicityInJets, lowBinBorderMultiplicityInJets, highBinBorderMultiplicityInJets); fhParticleMultiplicityInJetUncorrected->Sumw2();
  fhParticleMultiplicityInReflectedConeUncorrected = new THnSparseF("multiplicityInReflectedConeUncorrected", "multiplicityInReflectedConeUncorrected", nAxesMultiplicityInJets, nBinsMultiplicityInJets, lowBinBorderMultiplicityInJets, highBinBorderMultiplicityInJets); fhParticleMultiplicityInReflectedConeUncorrected->Sumw2();
  
  fhParticleMultiplicityInJet->SetBinEdges(1,jetPtBinsEEC);                              // Jet pT bins
  fhParticleMultiplicityInReflectedCone->SetBinEdges(1,jetPtBinsEEC);                    // Jet pT bins
  fhParticleMultiplicityInJetUncorrected->SetBinEdges(1,jetPtBinsEEC);                   // Jet pT bins
  fhParticleMultiplicityInReflectedConeUncorrected->SetBinEdges(1,jetPtBinsEEC);         // Jet pT bins
  
  fhParticleMultiplicityInJet->SetBinEdges(2,trackPtBinsEEC);                            // Track pT bins
  fhParticleMultiplicityInReflectedCone->SetBinEdges(2,trackPtBinsEEC);                  // Track pT bins
  fhParticleMultiplicityInJetUncorrected->SetBinEdges(2,trackPtBinsEEC);                 // Track pT bins
  fhParticleMultiplicityInReflectedConeUncorrected->SetBinEdges(2,trackPtBinsEEC);       // Track pT bins
  
  fhParticleMultiplicityInJet->SetBinEdges(3,wideCentralityBins);                        // Centrality bins
  fhParticleMultiplicityInReflectedCone->SetBinEdges(3,wideCentralityBins);              // Centrality bins
  fhParticleMultiplicityInJetUncorrected->SetBinEdges(3,wideCentralityBins);             // Centrality bins
  fhParticleMultiplicityInReflectedConeUncorrected->SetBinEdges(3,wideCentralityBins);   // Centrality bins
  
  // ======== THnSparse for all jets ========
  
  // Axis 0 for the inclusive jet histogram: jet pT
  nBinsJet[0] = nPtBinsJet;         // nBins for inclusive jet pT
  lowBinBorderJet[0] = minPtJet;    // low bin border for inclusive jet pT
  highBinBorderJet[0] = maxPtJet;   // high bin border for inclusive jet pT
  
  // Axis 1 for the inclusive jet histogram: jet phi
  nBinsJet[1] = nPhiBins;        // nBins for inclusive jet phi
  lowBinBorderJet[1] = minPhi;   // low bin border for inclusive jet phi
  highBinBorderJet[1] = maxPhi;  // high bin border for inclusive jet phi
  
  // Axis 2 for the inclusive jet histogram: jet eta
  nBinsJet[2] = nEtaBins;        // nBins for inclusive jet eta
  lowBinBorderJet[2] = minEta;   // low bin border for inclusive jet eta
  highBinBorderJet[2] = maxEta;  // high bin border for inclusive jet eta

  // Axis 3 for the inclusive jet histogram: jet eta_CM
  nBinsJet[3] = nEtaBins;        // nBins for inclusive jet eta_CM
  lowBinBorderJet[3] = minEta;   // low bin border for inclusive jet eta_CM
  highBinBorderJet[3] = maxEta;  // high bin border for inclusive jet eta_CM
  
  // Axis 3 for the inclusive jet histogram: centrality
  nBinsJet[4] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderJet[4] = minCentrality;  // low bin border for centrality
  highBinBorderJet[4] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the jet histogram: jet flavor (quark/gluon)
  nBinsJet[5] = nClosureParticleTypeBins;        // nBins for jet flavor
  lowBinBorderJet[5] = minClosureParticleType;   // low bin border for jet flavor
  highBinBorderJet[5] = maxClosureParticleType;  // high bin border for jet flavor
  
  // Create the histogram for all jets using the above binning information
  fhInclusiveJet = new THnSparseF("inclusiveJet","inclusiveJet",nAxesJet,nBinsJet,lowBinBorderJet,highBinBorderJet); fhInclusiveJet->Sumw2();

  // Set custom centrality bins for histograms
  fhInclusiveJet->SetBinEdges(4,wideCentralityBins);
  
  // ======== THnSparses for tracks and uncorrected tracks ========
  
  // Axis 0 for the track histogram: track pT
  nBinsTrack[0] = nPtBinsTrack;         // nBins for track pT
  lowBinBorderTrack[0] = minPtTrack;    // low bin border for track pT
  highBinBorderTrack[0] = maxPtTrack;   // high bin border for track pT
  
  // Axis 1 for the track histogram: track phi
  nBinsTrack[1] = nPhiBins;         // nBins for track phi
  lowBinBorderTrack[1] = minPhi;    // low bin border for track phi
  highBinBorderTrack[1] = maxPhi;   // high bin border for track phi
  
  // Axis 2 for the track histogram: track eta
  nBinsTrack[2] = nEtaBins;         // nBins for track eta
  lowBinBorderTrack[2] = minEta;    // low bin border for track eta
  highBinBorderTrack[2] = maxEta;   // high bin border for track eta
  
  // Axis 3 for the track histogram: centrality
  nBinsTrack[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderTrack[3] = minCentrality;  // low bin border for centrality
  highBinBorderTrack[3] = maxCentrality; // high bin border for centrality
  
  // Create the histograms for tracks and uncorrected tracks using the above binning information
  fhTrack = new THnSparseF("track","track",nAxesTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrack->Sumw2();
  fhTrackUncorrected = new THnSparseF("trackUncorrected","trackUncorrected",nAxesTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrackUncorrected->Sumw2();

  // Set custom centrality bins for histograms
  fhTrack->SetBinEdges(3,wideCentralityBins);
  fhTrackUncorrected->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparse for maximum particle pT within a jet ========
  
  // Axis 0 for the track density histogram: deltaR
  nBinsMaxParticlePtInJet[0] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderMaxParticlePtInJet[0] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderMaxParticlePtInJet[0] = maxJetPtEEC;   // high bin border for jet pT
  
  // Axis 1 for the track density histogram: particle pT (signal)
  nBinsMaxParticlePtInJet[1] = nPtBinsParticle;        // nBins for max particle pT
  lowBinBorderMaxParticlePtInJet[1] = minParticlePt;   // low bin border for max particle pT
  highBinBorderMaxParticlePtInJet[1] = maxParticlePt;  // high bin border for max particle pT
  
  // Axis 2 for the track density histogram: particle pT (background)
  nBinsMaxParticlePtInJet[2] = nPtBinsParticle;        // nBins for max particle pT
  lowBinBorderMaxParticlePtInJet[2] = minParticlePt;   // low bin border for max particle pT
  highBinBorderMaxParticlePtInJet[2] = maxParticlePt;  // high bin border for max particle pT
  
  // Axis 3 for track density histogram: centrality
  nBinsMaxParticlePtInJet[3] = nWideCentralityBins;    // nBins for wide centrality bins
  lowBinBorderMaxParticlePtInJet[3] = minCentrality;   // low bin border for centrality
  highBinBorderMaxParticlePtInJet[3] = maxCentrality;  // high bin border for centrality
  
  // Create the histograms for track density and pT weighted track density using the above binning information
  fhMaxPtParticleInJet = new THnSparseF("maxParticlePtInJet", "maxParticlePtInJet", nAxesMaxParticlePtInJet, nBinsMaxParticlePtInJet, lowBinBorderMaxParticlePtInJet, highBinBorderMaxParticlePtInJet); fhMaxPtParticleInJet->Sumw2();
  
  // Set custom bin borders for histograms
  fhMaxPtParticleInJet->SetBinEdges(0,jetPtBinsEEC);              // Jet pT bins
  fhMaxPtParticleInJet->SetBinEdges(3,wideCentralityBins);        // Centrality bins
  
  // ======== THnSparse for track density around the jets ========
  
  // Axis 0 for the track density histogram: deltaR
  nBinsTrackDensity[0] = nDeltaRBinsTrackDensity;        // nBins for deltaR between the track and the jet axis
  lowBinBorderTrackDensity[0] = minDeltaRTrackDensity;   // low bin border for deltaR
  highBinBorderTrackDensity[0] = maxDeltaRTrackDensity;  // high bin border for deltaR
  
  // Axis 1 for the track density histogram: jet pT
  nBinsTrackDensity[1] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderTrackDensity[1] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderTrackDensity[1] = maxJetPtEEC;   // high bin border for jet pT
  
  // Axis 2 for the track density histogram: track pT
  nBinsTrackDensity[2] = nTrackPtBinsEEC;        // nBins for track pT
  lowBinBorderTrackDensity[2] = minTrackPtEEC;   // low bin border for track pT
  highBinBorderTrackDensity[2] = maxTrackPtEEC;  // high bin border for track pT
  
  // Axis 3 for track density histogram: centrality
  nBinsTrackDensity[3] = nWideCentralityBins;    // nBins for wide centrality bins
  lowBinBorderTrackDensity[3] = minCentrality;   // low bin border for centrality
  highBinBorderTrackDensity[3] = maxCentrality;  // high bin border for centrality
  
  // Axis 4 for track density histogram: jet cone type (signal cone/reflected cone)
  nBinsTrackDensity[4] = nJetConeTypes;           // nBins for wide centrality bins
  lowBinBorderTrackDensity[4] = minJetConeType;   // low bin border for centrality
  highBinBorderTrackDensity[4] = maxJetConeType;  // high bin border for centrality
  
  // Axis 5 for track density histogram: track subevent (only relevant for Monte Carlo)
  nBinsTrackDensity[5] = nSubeventTypeBins;       // nBins for subevent types (Pythia/Hydjet)
  lowBinBorderTrackDensity[5] = minSubeventType;  // low bin border for subevent combinations
  highBinBorderTrackDensity[5] = maxSubeventType; // high bin border for subevent combinations
  
  // Create the histograms for track density and pT weighted track density using the above binning information
  fhParticleDensityAroundJet = new THnSparseF("particleDensity", "particleDensity", nAxesTrackDensity, nBinsTrackDensity ,lowBinBorderTrackDensity, highBinBorderTrackDensity); fhParticleDensityAroundJet->Sumw2();
  fhParticlePtDensityAroundJet = new THnSparseF("particlePtDensity", "particlePtDensity", nAxesTrackDensity, nBinsTrackDensity, lowBinBorderTrackDensity, highBinBorderTrackDensity); fhParticlePtDensityAroundJet->Sumw2();
  
  // Set custom bin borders for histograms
  fhParticleDensityAroundJet->SetBinEdges(1,jetPtBinsEEC);              // Jet pT bins
  fhParticlePtDensityAroundJet->SetBinEdges(1,jetPtBinsEEC);            // Jet pT bins
  
  fhParticleDensityAroundJet->SetBinEdges(2,trackPtBinsEEC);            // Track pT bins
  fhParticlePtDensityAroundJet->SetBinEdges(2,trackPtBinsEEC);          // Track pT bins
  
  fhParticleDensityAroundJet->SetBinEdges(3,wideCentralityBins);        // Centrality bins
  fhParticlePtDensityAroundJet->SetBinEdges(3,wideCentralityBins);      // Centrality bins
  
  // ======== THnSparses for energy-energy correlators ========
  
  // Axis 0 for the energy-energy correlator histogram: deltaR
  nBinsEnergyEnergyCorrelator[0] = nDeltaRBinsEEC;        // nBins for deltaR between the two tracks
  lowBinBorderEnergyEnergyCorrelator[0] = minDeltaREEC;   // low bin border for deltaR
  highBinBorderEnergyEnergyCorrelator[0] = maxDeltaREEC;  // high bin border for deltaR
  
  // Axis 1 for the energy-energy correlator histogram: jet pT
  nBinsEnergyEnergyCorrelator[1] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderEnergyEnergyCorrelator[1] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderEnergyEnergyCorrelator[1] = maxJetPtEEC;   // high bin border for jet pT
  
  // Axis 2 for the energy-energy correlator histogram: track pT
  nBinsEnergyEnergyCorrelator[2] = nTrackPtBinsEEC;        // nBins for track pT
  lowBinBorderEnergyEnergyCorrelator[2] = minTrackPtEEC;   // low bin border for track pT
  highBinBorderEnergyEnergyCorrelator[2] = maxTrackPtEEC;  // high bin border for track pT
  
  // Axis 3 for the energy-energy correlator histogram: centrality
  nBinsEnergyEnergyCorrelator[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderEnergyEnergyCorrelator[3] = minCentrality;  // low bin border for centrality
  highBinBorderEnergyEnergyCorrelator[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the energy-energy correlator histogram: track pairing type (same jet/reflected cone jet)
  nBinsEnergyEnergyCorrelator[4] = nTrackPairingTypeBins;       // nBins for track pairing types
  lowBinBorderEnergyEnergyCorrelator[4] = minTrackPairingType;  // low bin border for track pairing type
  highBinBorderEnergyEnergyCorrelator[4] = maxTrackPairingType; // high bin border for track pairing type
  
  // Axis 5 for the energy-energy correlator histogram: track subevent (only relevant for Monte Carlo)
  nBinsEnergyEnergyCorrelator[5] = nSubeventCombinationBins;       // nBins for subevent combinations (pythia-pythia, pythia-hydjet, hydjet-hydjet)
  lowBinBorderEnergyEnergyCorrelator[5] = minSubeventCombination;  // low bin border for subevent combinations
  highBinBorderEnergyEnergyCorrelator[5] = maxSubeventCombination; // high bin border for subevent combinations

  // Axis 6 for the energy-energy correlator histogram: energy weight index
  nBinsEnergyEnergyCorrelator[6] = nEnergyWeightsEEC;          // Number of energy weights defined in the card
  lowBinBorderEnergyEnergyCorrelator[6] = minEnergyWeightEEC;  // Index of the first defined energy weight
  highBinBorderEnergyEnergyCorrelator[6] = maxEnergyWeightEEC; // Index of the last defined energy weight

  
  // Create the histograms for energy-energy correlators with and without track efficiency corrections
  fhEnergyEnergyCorrelator = new THnSparseF("energyEnergyCorrelator", "energyEnergyCorrelator", nAxesEnergyEnergyCorrelator, nBinsEnergyEnergyCorrelator, lowBinBorderEnergyEnergyCorrelator, highBinBorderEnergyEnergyCorrelator); fhEnergyEnergyCorrelator->Sumw2();
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus = new THnSparseF("energyEnergyCorrelatorEfficiencyVariationPlus", "energyEnergyCorrelatorEfficiencyVariationPlus", nAxesEnergyEnergyCorrelator, nBinsEnergyEnergyCorrelator, lowBinBorderEnergyEnergyCorrelator, highBinBorderEnergyEnergyCorrelator); fhEnergyEnergyCorrelatorEfficiencyVariationPlus->Sumw2();
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus = new THnSparseF("energyEnergyCorrelatorEfficiencyVariationMinus", "energyEnergyCorrelatorEfficiencyVariationMinus", nAxesEnergyEnergyCorrelator, nBinsEnergyEnergyCorrelator, lowBinBorderEnergyEnergyCorrelator, highBinBorderEnergyEnergyCorrelator); fhEnergyEnergyCorrelatorEfficiencyVariationMinus->Sumw2();
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus = new THnSparseF("energyEnergyCorrelatorPairEfficiencyVariationPlus", "energyEnergyCorrelatorPairEfficiencyVariationPlus", nAxesEnergyEnergyCorrelator, nBinsEnergyEnergyCorrelator, lowBinBorderEnergyEnergyCorrelator, highBinBorderEnergyEnergyCorrelator); fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus->Sumw2();
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus = new THnSparseF("energyEnergyCorrelatorPairEfficiencyVariationMinus", "energyEnergyCorrelatorPairEfficiencyVariationMinus", nAxesEnergyEnergyCorrelator, nBinsEnergyEnergyCorrelator, lowBinBorderEnergyEnergyCorrelator, highBinBorderEnergyEnergyCorrelator); fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus->Sumw2();
  
  
  // Set custom bin borders for histograms
  fhEnergyEnergyCorrelator->SetBinEdges(0,deltaRBinsEEC);                                  // DeltaR bins
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus->SetBinEdges(0,deltaRBinsEEC);           // DeltaR bins
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus->SetBinEdges(0,deltaRBinsEEC);          // DeltaR bins
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus->SetBinEdges(0,deltaRBinsEEC);       // DeltaR bins
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus->SetBinEdges(0,deltaRBinsEEC);      // DeltaR bins
  
  fhEnergyEnergyCorrelator->SetBinEdges(1,jetPtBinsEEC);                                   // Jet pT bins
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus->SetBinEdges(1,jetPtBinsEEC);            // Jet pT bins
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus->SetBinEdges(1,jetPtBinsEEC);           // Jet pT bins
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus->SetBinEdges(1,jetPtBinsEEC);        // Jet pT bins
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus->SetBinEdges(1,jetPtBinsEEC);       // Jet pT bins
  
  fhEnergyEnergyCorrelator->SetBinEdges(2,trackPtBinsEEC);                                 // Track pT bins
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus->SetBinEdges(2,trackPtBinsEEC);          // Track pT bins
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus->SetBinEdges(2,trackPtBinsEEC);         // Track pT bins
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus->SetBinEdges(2,trackPtBinsEEC);      // Track pT bins
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus->SetBinEdges(2,trackPtBinsEEC);     // Track pT bins
  
  fhEnergyEnergyCorrelator->SetBinEdges(3,wideCentralityBins);                             // Centrality bins
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus->SetBinEdges(3,wideCentralityBins);      // Centrality bins
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus->SetBinEdges(3,wideCentralityBins);     // Centrality bins
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus->SetBinEdges(3,wideCentralityBins);  // Centrality bins
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus->SetBinEdges(3,wideCentralityBins); // Centrality bins
  
  // ======== THnSparses for jet pT closures ========
  
  // Axis 0 for the jet pT closure histogram: generator level jet pT
  nBinsJetClosure[0] = nClosurePtBins;       // nBins for generator level pT bins in closure plots
  lowBinBorderJetClosure[0] = minClosurePt;  // low bin border generator level pT in closure plots
  highBinBorderJetClosure[0] = maxClosurePt; // high bin border generator level pT in closure plots
  
  // Axis 1 for the jet pT closure histogram: reconstructed jet pT
  nBinsJetClosure[1] = nClosurePtBins;       // nBins for reconstructed jet pT bins in closure plots
  lowBinBorderJetClosure[1] = minClosurePt;  // low bin border for reconstructed jet pT in closure plots
  highBinBorderJetClosure[1] = maxClosurePt; // high bin border for reconstructed jet pT in closure plots
  
  // Axis 2 for the jet pT closure histogram: generator level jet eta
  nBinsJetClosure[2] = nEtaBins;             // nBins for jet eta
  lowBinBorderJetClosure[2] = minEta;        // low bin border for jet eta
  highBinBorderJetClosure[2] = maxEta;       // high bin border for jet eta
  
  // Axis 3 for the jet pT closure histogram: centrality
  nBinsJetClosure[3] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetClosure[3] = minCentrality;    // low bin border for centrality
  highBinBorderJetClosure[3] = maxCentrality;   // high bin border for centrality
  
  // Axis 4 for the jet pT closure histogram: ref parton = quark/gluon
  nBinsJetClosure[4] = nClosureParticleTypeBins;         // nBins for reference parton
  lowBinBorderJetClosure[4] = minClosureParticleType;    // low bin border for reference parton
  highBinBorderJetClosure[4] = maxClosureParticleType;   // high bin border for reference parton
  
  // Axis 5 for the jet pT closure histogram: reco/gen ratio for closure
  nBinsJetClosure[5] = nClosureRatioBins;        // nBins for closure ratio
  lowBinBorderJetClosure[5] = minClosureRatio;   // low bin border for closure ratio
  highBinBorderJetClosure[5] = maxClosureRatio;  // high bin border for closure ratio

  // Axis 6 for the jet pT closure histogram: generator level jet eta
  nBinsJetClosure[6] = nPhiBins;             // nBins for jet eta
  lowBinBorderJetClosure[6] = minPhi;        // low bin border for jet eta
  highBinBorderJetClosure[6] = maxPhi;       // high bin border for jet eta
  
  // Create histograms for jet pT closure
  fhJetPtClosure = new THnSparseF("jetPtClosure", "jetPtClosure", nAxesJetClosure, nBinsJetClosure, lowBinBorderJetClosure, highBinBorderJetClosure); fhJetPtClosure->Sumw2();
  
  // Set custom centrality bins for histograms
  fhJetPtClosure->SetBinEdges(3,wideCentralityBins);
  
  // ======== RooUnfold compatible energy-energy correlator distributions. Assume only jet pT is unfolded ========

  // Axis 0 for the jet pT unfolding distribution histogram: deltaR in jet pT bins
  nBinsJetPtUnfoldDistribution[0] = nUnfoldingBins;          // nBins for the unfolded quantities
  lowBinBorderJetPtUnfoldDistribution[0] = minUnfoldingBin;  // low bin border for the unfolded quantities
  highBinBorderJetPtUnfoldDistribution[0] = maxUnfoldingBin; // high bin border for the unfolded quantities

  // Axis 1 for the jet pT unfolding distribution histogram: track pT
  nBinsJetPtUnfoldDistribution[1] = nTrackPtBinsEEC;        // nBins for track pT in energy-energy correlator histograms
  lowBinBorderJetPtUnfoldDistribution[1] = minTrackPtEEC;   // low bin border for track pT in energy-energy correlator histograms
  highBinBorderJetPtUnfoldDistribution[1] = maxTrackPtEEC;  // high bin border for track pT in energy-energy correlator histograms

  // Axis 2 for the jet pT unfolding distribution histogram: centrality
  nBinsJetPtUnfoldDistribution[2] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetPtUnfoldDistribution[2] = minCentrality;    // low bin border for centrality
  highBinBorderJetPtUnfoldDistribution[2] = maxCentrality;   // high bin border for centrality

  // Axis 3 for the jet pT unfolding distribution histogram: energy weight index
  nBinsJetPtUnfoldDistribution[3] = nEnergyWeightsEEC;          // Number of energy weights defined in the card
  lowBinBorderJetPtUnfoldDistribution[3] = minEnergyWeightEEC;  // Index of the first defined energy weight
  highBinBorderJetPtUnfoldDistribution[3] = maxEnergyWeightEEC; // Index of the last defined energy weight

  // Create the histogram for reconstructed jet pT unfolding distributions
  fhUnfoldingMeasured = new THnSparseF("jetPtUnfoldingMeasured", "jetPtUnfoldingMeasured", nAxesJetPtUnfoldDistribution, nBinsJetPtUnfoldDistribution, lowBinBorderJetPtUnfoldDistribution, highBinBorderJetPtUnfoldDistribution); fhUnfoldingMeasured->Sumw2();

  // Set custom bin axes for the histograms
  fhUnfoldingMeasured->SetBinEdges(0, fullUnfoldingBinning);
  fhUnfoldingMeasured->SetBinEdges(1, trackPtBinsEEC);
  fhUnfoldingMeasured->SetBinEdges(2, wideCentralityBins);

  // The unfolding axis can have different binning for truth level

  // Axis 0 for the jet pT unfolding distribution histogram: deltaR in jet pT bins
  nBinsJetPtUnfoldDistribution[0] = nUnfoldingBins;          // nBins for the unfolded quantities
  lowBinBorderJetPtUnfoldDistribution[0] = minUnfoldingBin;  // low bin border for the unfolded quantities
  highBinBorderJetPtUnfoldDistribution[0] = maxUnfoldingBin; // high bin border for the unfolded quantities

  // Create the histogram for generator level jet pT unfolding distributions
  fhUnfoldingTruth = new THnSparseF("jetPtUnfoldingTruth", "jetPtUnfoldingTruth", nAxesJetPtUnfoldDistribution, nBinsJetPtUnfoldDistribution, lowBinBorderJetPtUnfoldDistribution, highBinBorderJetPtUnfoldDistribution); fhUnfoldingTruth->Sumw2();

  // Set custom bin axes for the histograms
  fhUnfoldingTruth->SetBinEdges(0, fullUnfoldingBinning);
  fhUnfoldingTruth->SetBinEdges(1, trackPtBinsEEC);
  fhUnfoldingTruth->SetBinEdges(2, wideCentralityBins);

  // ======== RooUnfold compatible response matrix. Assume only jet pT is unfolded ========

  // Axis 0 for the jet pT unfolding response histogram: deltaR in reconstructed jet pT bins
  nBinsJetPtUnfoldResponse[0] = nUnfoldingBins;          // nBins for the unfolded quantities
  lowBinBorderJetPtUnfoldResponse[0] = minUnfoldingBin;  // low bin border for the unfolded quantities
  highBinBorderJetPtUnfoldResponse[0] = maxUnfoldingBin; // high bin border for the unfolded quantities
  
  // Axis 1 for the jet pT unfolding response histogram: deltaR in generator level jet pT bins
  nBinsJetPtUnfoldResponse[1] = nUnfoldingBins;          // nBins for the unfolded quantities
  lowBinBorderJetPtUnfoldResponse[1] = minUnfoldingBin;  // low bin border for the unfolded quantities
  highBinBorderJetPtUnfoldResponse[1] = maxUnfoldingBin; // high bin border for the unfolded quantities
  
  // Axis 2 for the jet pT unfolding response histogram: track pT
  nBinsJetPtUnfoldResponse[2] = nTrackPtBinsEEC;        // nBins for track pT in energy-energy correlator histograms
  lowBinBorderJetPtUnfoldResponse[2] = minTrackPtEEC;   // low bin border for track pT in energy-energy correlator histograms
  highBinBorderJetPtUnfoldResponse[2] = maxTrackPtEEC;  // high bin border for track pT in energy-energy correlator histograms

  // Axis 3 for the jet pT unfolding response histogram: centrality
  nBinsJetPtUnfoldResponse[3] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetPtUnfoldResponse[3] = minCentrality;    // low bin border for centrality
  highBinBorderJetPtUnfoldResponse[3] = maxCentrality;   // high bin border for centrality

  // Axis 4 for the jet pT unfolding response histogram: energy weight index
  nBinsJetPtUnfoldResponse[4] = nEnergyWeightsEEC;          // Number of energy weights defined in the card
  lowBinBorderJetPtUnfoldResponse[4] = minEnergyWeightEEC;  // Index of the first defined energy weight
  highBinBorderJetPtUnfoldResponse[4] = maxEnergyWeightEEC; // Index of the last defined energy weight

  // Create the histogram for jet pT unfolding response
  fhUnfoldingResponse = new THnSparseF("jetPtUnfoldingResponse", "jetPtUnfoldingResponse", nAxesJetPtUnfoldResponse, nBinsJetPtUnfoldResponse, lowBinBorderJetPtUnfoldResponse, highBinBorderJetPtUnfoldResponse); fhUnfoldingResponse->Sumw2();
  fhJetPtUnfoldingCovariance = new THnSparseF("jetPtUnfoldingCovariance", "jetPtUnfoldingCovariance", nAxesJetPtUnfoldResponse, nBinsJetPtUnfoldResponse, lowBinBorderJetPtUnfoldResponse, highBinBorderJetPtUnfoldResponse); fhJetPtUnfoldingCovariance->Sumw2();

  // Set custom bin axes for the histograms
  fhUnfoldingResponse->SetBinEdges(0, fullUnfoldingBinning);
  fhJetPtUnfoldingCovariance->SetBinEdges(0, fullUnfoldingBinning);
  fhUnfoldingResponse->SetBinEdges(1, fullUnfoldingBinning);
  fhJetPtUnfoldingCovariance->SetBinEdges(1, fullUnfoldingBinning);
  fhUnfoldingResponse->SetBinEdges(2, trackPtBinsEEC);
  fhJetPtUnfoldingCovariance->SetBinEdges(2, trackPtBinsEEC);
  fhUnfoldingResponse->SetBinEdges(3, wideCentralityBins);
  fhJetPtUnfoldingCovariance->SetBinEdges(3, wideCentralityBins);

  // ======== RooUnfold compatible jet pT distributions for 1-dimensional unfolding ========

  // Axis 0 for the 1-dimensional jet pT unfolding histogram: jet pT
  nBinsOneDimensionalJetPtUnfoldDistribution[0] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderOneDimensionalJetPtUnfoldDistribution[0] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderOneDimensionalJetPtUnfoldDistribution[0] = maxJetPtEEC;   // high bin border for jet pT
  
  // Axis 1 for the 1-dimensional jet pT unfolding histogram: centrality
  nBinsOneDimensionalJetPtUnfoldDistribution[1] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderOneDimensionalJetPtUnfoldDistribution[1] = minCentrality;  // low bin border for centrality
  highBinBorderOneDimensionalJetPtUnfoldDistribution[1] = maxCentrality; // high bin border for centrality
  
  // Create the histogram for all jets using the above binning information
  fhJetPtUnfoldingMeasured = new THnSparseF("oneDimensionalJetPtUnfoldingMeasured", "oneDimensionalJetPtUnfoldingMeasured", nAxesOneDimensionalJetPtUnfoldDistribution, nBinsOneDimensionalJetPtUnfoldDistribution, lowBinBorderOneDimensionalJetPtUnfoldDistribution, highBinBorderOneDimensionalJetPtUnfoldDistribution); fhJetPtUnfoldingMeasured->Sumw2();
  fhJetPtUnfoldingTruth = new THnSparseF("oneDimensionalJetPtUnfoldingTruth", "oneDimensionalJetPtUnfoldingTruth", nAxesOneDimensionalJetPtUnfoldDistribution, nBinsOneDimensionalJetPtUnfoldDistribution, lowBinBorderOneDimensionalJetPtUnfoldDistribution, highBinBorderOneDimensionalJetPtUnfoldDistribution); fhJetPtUnfoldingTruth->Sumw2();

  // Set custom bin boundaries for the histograms
  fhJetPtUnfoldingMeasured->SetBinEdges(0,jetPtBinsEEC);        // Custom jet pT bins
  fhJetPtUnfoldingTruth->SetBinEdges(0,jetPtBinsEEC);           // Custom jet pT bins
  fhJetPtUnfoldingMeasured->SetBinEdges(1,wideCentralityBins);  // Custom centrality bins
  fhJetPtUnfoldingTruth->SetBinEdges(1,wideCentralityBins);     // Custom centrality bins

  // ======== RooUnfold compatible jet pT response matrix for 1-dimensional unfolding ========

  // Axis 0 for the 1-dimensional jet pT unfolding response matrix: measured jet pT
  nBinsOneDimensionalJetPtUnfoldResponse[0] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderOneDimensionalJetPtUnfoldResponse[0] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderOneDimensionalJetPtUnfoldResponse[0] = maxJetPtEEC;   // high bin border for jet pT

  // Axis 1 for the 1-dimensional jet pT unfolding response matrix: truth jet pT
  nBinsOneDimensionalJetPtUnfoldResponse[1] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderOneDimensionalJetPtUnfoldResponse[1] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderOneDimensionalJetPtUnfoldResponse[1] = maxJetPtEEC;   // high bin border for jet pT
  
  // Axis 2 for the 1-dimensional jet pT unfolding response matrix: centrality
  nBinsOneDimensionalJetPtUnfoldResponse[2] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderOneDimensionalJetPtUnfoldResponse[2] = minCentrality;  // low bin border for centrality
  highBinBorderOneDimensionalJetPtUnfoldResponse[2] = maxCentrality; // high bin border for centrality
  
  // Create the histogram for all jets using the above binning information
  fhJetPtUnfoldingResponse = new THnSparseF("oneDimensionalJetPtUnfoldingResponse", "oneDimensionalJetPtUnfoldingResponse", nAxesOneDimensionalJetPtUnfoldResponse, nBinsOneDimensionalJetPtUnfoldResponse, lowBinBorderOneDimensionalJetPtUnfoldResponse, highBinBorderOneDimensionalJetPtUnfoldResponse); fhJetPtUnfoldingResponse->Sumw2();

  // Set custom bin boundaries for the histograms
  fhJetPtUnfoldingResponse->SetBinEdges(0,jetPtBinsEEC);        // Custom measured jet pT bins
  fhJetPtUnfoldingResponse->SetBinEdges(1,jetPtBinsEEC);        // Custom truth jet pT bins
  fhJetPtUnfoldingResponse->SetBinEdges(2,wideCentralityBins);  // Custom centrality bins

  // ======== Testing on how to do track to particle matching ========

  // Axis 0 of the particle matching QA histograms: number of generator level particles close to track
  nBinsParticleMatchingQA[0] = nNumberOftracksBins;
  lowBinBorderParticleMatchingQA[0] = minNumberOfTracks;
  highBinBorderParticleMatchingQA[0] = maxNumberOfTracks;

  // Axis 1 for the particle matching QA histograms: jet pT
  nBinsParticleMatchingQA[1] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderParticleMatchingQA[1] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderParticleMatchingQA[1] = maxJetPtEEC;   // high bin border for jet pT

  // Axis 2 for the particle matching QA histograms: track pT
  nBinsParticleMatchingQA[2] = nTrackPtBinsEEC;        // nBins for track pT
  lowBinBorderParticleMatchingQA[2] = minTrackPtEEC;   // low bin border for track pT
  highBinBorderParticleMatchingQA[2] = maxTrackPtEEC;  // high bin border for track pT

  // Axis 2 for the particle matching QA histograms: centrality
  nBinsParticleMatchingQA[3] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderParticleMatchingQA[3] = minCentrality;    // low bin border for centrality
  highBinBorderParticleMatchingQA[3] = maxCentrality;   // high bin border for centrality

  // Create the histogram for jet pT unfolding response
  fhParticlesCloseToTracks = new THnSparseF("particlesCloseToTracks", "particlesCloseToTracks", nAxesParticleMatchingQA, nBinsParticleMatchingQA, lowBinBorderParticleMatchingQA, highBinBorderParticleMatchingQA); fhParticlesCloseToTracks->Sumw2();
  fhTracksWithMatchedParticle = new THnSparseF("tracksWithMatchedParticle", "tracksWithMatchedParticle", nAxesParticleMatchingQA, nBinsParticleMatchingQA, lowBinBorderParticleMatchingQA, highBinBorderParticleMatchingQA); fhTracksWithMatchedParticle->Sumw2();

  // Set custom bin axes for the histograms
  fhParticlesCloseToTracks->SetBinEdges(1,jetPtBinsEEC);              // Jet pT bins
  fhTracksWithMatchedParticle->SetBinEdges(1,jetPtBinsEEC); 
  fhParticlesCloseToTracks->SetBinEdges(2,trackPtBinsEEC);            // Track pT bins
  fhTracksWithMatchedParticle->SetBinEdges(2,trackPtBinsEEC); 
  fhParticlesCloseToTracks->SetBinEdges(3, wideCentralityBins);       // Centrality bins
  fhTracksWithMatchedParticle->SetBinEdges(3, wideCentralityBins);

  // ======== Response matrix for deltaR between two particles ========

  // Axis 0 of the particle deltaR response matrix histogram: deltaR between two tracks matched to particles
  nBinsParticleDeltaRResponseMatrix[0] = nDeltaRBinsEEC;        // nBins for deltaR
  lowBinBorderParticleDeltaRResponseMatrix[0] = minDeltaREEC;   // low bin border for deltaR
  highBinBorderParticleDeltaRResponseMatrix[0] = maxDeltaREEC;  // high bin border for deltaR

  // Axis 1 of the particle deltaR response matrix histogram: deltaR between two particles matched to tracks
  nBinsParticleDeltaRResponseMatrix[1] = nDeltaRBinsEEC;        // nBins for deltaR
  lowBinBorderParticleDeltaRResponseMatrix[1] = minDeltaREEC;   // low bin border for deltaR
  highBinBorderParticleDeltaRResponseMatrix[1] = maxDeltaREEC;  // high bin border for deltaR

  // Axis 2 for the particle deltaR response matrix histogram: jet pT
  nBinsParticleDeltaRResponseMatrix[2] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderParticleDeltaRResponseMatrix[2] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderParticleDeltaRResponseMatrix[2] = maxJetPtEEC;   // high bin border for jet pT

  // Axis 3 for the particle deltaR response matrix histogram: lower reconstructed track pT of the pair
  nBinsParticleDeltaRResponseMatrix[3] = nTrackPtBinsEEC;         // nBins for track pT
  lowBinBorderParticleDeltaRResponseMatrix[3] = minTrackPtEEC;    // low bin border for track pT
  highBinBorderParticleDeltaRResponseMatrix[3] = maxTrackPtEEC;   // high bin border for track pT

  // Axis 4 for the particle deltaR response matrix histogram: centrality
  nBinsParticleDeltaRResponseMatrix[4] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderParticleDeltaRResponseMatrix[4] = minCentrality;    // low bin border for centrality
  highBinBorderParticleDeltaRResponseMatrix[4] = maxCentrality;   // high bin border for centrality

  // Create the histogram for the response matrix between two particles
  fhParticleDeltaRResponse = new THnSparseF("particleDeltaRResponseMatrix", "particleDeltaRResponseMatrix", nAxesParticleDeltaRResponseMatrix, nBinsParticleDeltaRResponseMatrix, lowBinBorderParticleDeltaRResponseMatrix, highBinBorderParticleDeltaRResponseMatrix); fhParticleDeltaRResponse->Sumw2();

  // Set custom bin axes for the histograms
  fhParticleDeltaRResponse->SetBinEdges(0,deltaRBinsEEC);         // DeltaR bins for reconstructed tracks
  fhParticleDeltaRResponse->SetBinEdges(1,deltaRBinsEEC);         // DeltaR bins for generator level particles
  fhParticleDeltaRResponse->SetBinEdges(2,jetPtBinsEEC);          // Jet pT bins
  fhParticleDeltaRResponse->SetBinEdges(3,trackPtBinsEEC);        // Track pT bins
  fhParticleDeltaRResponse->SetBinEdges(4, wideCentralityBins);   // Centrality bins

  // ======== Response matrix for pT1*pT2 between two particles ========

  // Axis 0 of the particle pT1*pT2 response matrix histogram: pT1*pT2 between two tracks matched to particles
  nBinsParticlePtResponseMatrix[0] = nPt1TimesPt2Bins;        // nBins for pT1*pT2
  lowBinBorderParticlePtResponseMatrix[0] = minPt1TimesPt2;   // low bin border for pT1*pT2
  highBinBorderParticlePtResponseMatrix[0] = maxPt1TimesPt2;  // high bin border for pT1*pT2

  // Axis 1 of the particle pT1*pT2 response matrix histogram: pT1*pT2 between two particles matched to tracks
  nBinsParticlePtResponseMatrix[1] = nPt1TimesPt2Bins;        // nBins for pT1*pT2
  lowBinBorderParticlePtResponseMatrix[1] = minPt1TimesPt2;   // low bin border for pT1*pT2
  highBinBorderParticlePtResponseMatrix[1] = maxPt1TimesPt2;  // high bin border for pT1*pT2

  // Axis 2 for the particle pT1*pT2 response matrix histogram: jet pT
  nBinsParticlePtResponseMatrix[2] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderParticlePtResponseMatrix[2] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderParticlePtResponseMatrix[2] = maxJetPtEEC;   // high bin border for jet pT

  // Axis 3 for the particle pT1*pT2 response matrix histogram: lower reconstructed track pT of the pair
  nBinsParticlePtResponseMatrix[3] = nTrackPtBinsEEC;         // nBins for track pT
  lowBinBorderParticlePtResponseMatrix[3] = minTrackPtEEC;    // low bin border for track pT
  highBinBorderParticlePtResponseMatrix[3] = maxTrackPtEEC;   // high bin border for track pT

  // Axis 4 for the particle pT1*pT2 response matrix histogram: centrality
  nBinsParticlePtResponseMatrix[4] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderParticlePtResponseMatrix[4] = minCentrality;    // low bin border for centrality
  highBinBorderParticlePtResponseMatrix[4] = maxCentrality;   // high bin border for centrality

  // Axis 5 for the particle pT1*pT2 response matrix histogram: reco/gen ratio for closure
  nBinsParticlePtResponseMatrix[5] = nClosureRatioBins;        // nBins for closure ratio
  lowBinBorderParticlePtResponseMatrix[5] = minClosureRatio;   // low bin border for closure ratio
  highBinBorderParticlePtResponseMatrix[5] = maxClosureRatio;  // high bin border for closure ratio

  // Axis 6 for the particle pT1*pT2 response matrix histogram: energy weight index
  nBinsParticlePtResponseMatrix[6] = nEnergyWeightsEEC;          // Number of energy weights defined in the card
  lowBinBorderParticlePtResponseMatrix[6] = minEnergyWeightEEC;  // Index of the first defined energy weight
  highBinBorderParticlePtResponseMatrix[6] = maxEnergyWeightEEC; // Index of the last defined energy weight

  // Create the histogram for the response matrix between two particles
  fhParticlePtResponse = new THnSparseF("particlePtResponseMatrix", "particlePtResponseMatrix", nAxesParticlePtResponseMatrix, nBinsParticlePtResponseMatrix, lowBinBorderParticlePtResponseMatrix, highBinBorderParticlePtResponseMatrix); fhParticlePtResponse->Sumw2();

  // Set custom bin axes for the histograms
  fhParticlePtResponse->SetBinEdges(0,pt1TimesPt2Bins);       // pT1*pT2 bins for reconstructed tracks
  fhParticlePtResponse->SetBinEdges(1,pt1TimesPt2Bins);       // pT1*pT2 bins for generator level particles
  fhParticlePtResponse->SetBinEdges(2,jetPtBinsEEC);          // Jet pT bins
  fhParticlePtResponse->SetBinEdges(3,trackPtBinsEEC);        // Track pT bins
  fhParticlePtResponse->SetBinEdges(4,wideCentralityBins);   // Centrality bins

  // ======== Jet pT histograms within the reflected cone ========

  // Axis 0 for the number of jets in reflected cone histograms: number of jets
  nBinsReflectedConeQA[0] = nNumberOfJetsReflectedConeBins;         // nBins for number of jets within the reflected cone
  lowBinBorderReflectedConeQA[0] = minNumberOfJetsReflectedCone;    // low bin border for centrality
  highBinBorderReflectedConeQA[0] = maxNumberOfJetsReflectedCone;   // high bin border for centrality

  // Axis 1 for the number of jets in reflected cone histograms: centrality
  nBinsReflectedConeQA[1] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderReflectedConeQA[1] = minCentrality;    // low bin border for centrality
  highBinBorderReflectedConeQA[1] = maxCentrality;   // high bin border for centrality

  // Create the histogram for the number of jets within the reflected cone
  fhJetNumberInReflectedCone = new THnSparseF("jetNumberInReflectedCone", "jetNumberInReflectedCone", nAxesReflectedConeQA, nBinsReflectedConeQA, lowBinBorderReflectedConeQA, highBinBorderReflectedConeQA); fhJetNumberInReflectedCone->Sumw2();

  // Set custom bin axes for the histograms
  fhJetNumberInReflectedCone->SetBinEdges(1,wideCentralityBins);  // Centrality bins 

  // Axis 0 for the jet pT in reflected cone histograms: jet pT
  nBinsReflectedConeQA[0] = nJetPtInReflectedConeBins;       // nBins for jet pT within the reflected cone
  lowBinBorderReflectedConeQA[0] = minJetPtReflectedCone;    // low bin border for jet pT within the reflected cone
  highBinBorderReflectedConeQA[0] = maxJetPtReflectedCone;   // high bin border for jet pT within the reflected cone

  // Axis 1 for the number of jets in reflected cone histograms: centrality
  nBinsReflectedConeQA[1] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderReflectedConeQA[1] = minCentrality;    // low bin border for centrality
  highBinBorderReflectedConeQA[1] = maxCentrality;   // high bin border for centrality

  // Create the histogram for the jet pT within the reflected cone
  fhJetPtInReflectedCone = new THnSparseF("jetPtInReflectedCone", "jetPtInReflectedCone", nAxesReflectedConeQA, nBinsReflectedConeQA, lowBinBorderReflectedConeQA, highBinBorderReflectedConeQA); fhJetPtInReflectedCone->Sumw2();

  // Set custom bin axes for the histograms
  fhJetPtInReflectedCone->SetBinEdges(0,jetPtReflectedConeBins);  // Jet pT bins 
  fhJetPtInReflectedCone->SetBinEdges(1,wideCentralityBins);      // Centrality bins 

  // Release memory allocated to new arrays
  delete[] wideCentralityBins;
  delete[] ptHatBins;
  delete[] jetPtBinsEEC;
  delete[] trackPtBinsEEC;
  delete[] fullUnfoldingBinning;

}

/*
 * Write the histograms to file
 */
void EECHistograms::Write() const{
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhVertexZWeighted->Write();
  fhEvents->Write();
  fhTriggers->Write();
  fhTriggersAfterSelection->Write();
  fhTrackCuts->Write();
  fhCentrality->Write();
  fhCentralityWeighted->Write();
  fhHFEnergy->Write();
  fhPtHat->Write();
  fhPtHatWeighted->Write();
  fhMixedEventCounts->Write();
  fhMultiplicity->Write();
  fhInclusiveJet->Write();
  fhTrack->Write();
  fhTrackUncorrected->Write();
  fhParticleDensityAroundJet->Write();
  fhParticlePtDensityAroundJet->Write();
  fhParticleMultiplicityInJet->Write();
  fhParticleMultiplicityInReflectedCone->Write();
  fhParticleMultiplicityInJetUncorrected->Write();
  fhParticleMultiplicityInReflectedConeUncorrected->Write();
  fhMaxPtParticleInJet->Write();
  fhEnergyEnergyCorrelator->Write();
  fhEnergyEnergyCorrelatorEfficiencyVariationPlus->Write();
  fhEnergyEnergyCorrelatorEfficiencyVariationMinus->Write();
  fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus->Write();
  fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus->Write();
  fhJetPtClosure->Write();
  fhJetNumberInReflectedCone->Write();
  fhJetPtInReflectedCone->Write();
  fhUnfoldingMeasured->Write();
  fhUnfoldingTruth->Write();
  fhUnfoldingResponse->Write();
  fhJetPtUnfoldingMeasured->Write();
  fhJetPtUnfoldingTruth->Write();
  fhJetPtUnfoldingResponse->Write();
  fhJetPtUnfoldingCovariance->Write();
  fhParticlesCloseToTracks->Write();
  fhTracksWithMatchedParticle->Write();
  fhParticleDeltaRResponse->Write();
  fhParticlePtResponse->Write();
}

/*
 * Write the histograms to a given file
 */
void EECHistograms::Write(TString outputFileName) const{
  
  // Define the output file
  TFile* outputFile = new TFile(outputFileName, "RECREATE");
  
  // Write the histograms to file
  Write();
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
}

// Get the number of deltaR bins for energy-energy correlators
Int_t EECHistograms::GetNDeltaRBinsEEC(){
  return fhEnergyEnergyCorrelator->GetAxis(0)->GetNbins();
}

// Getter for the low bin border of the i:th deltaR bin
Double_t EECHistograms::GetDeltaRBinBorderLowEEC(Int_t iBin){
  return fhEnergyEnergyCorrelator->GetAxis(0)->GetBinLowEdge(iBin);
} 

// Getter for the high bin border of the i:th deltaR bin
Double_t EECHistograms::GetDeltaRBinBorderHighEEC(Int_t iBin){
  return fhEnergyEnergyCorrelator->GetAxis(0)->GetBinUpEdge(iBin);
}

/*
 * Implementation of EECHistogramManager
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>

// Own includes
#include "EECHistogramManager.h"

/*
 * Default constructor
 */
EECHistogramManager::EECHistogramManager() :
  fInputFile(NULL),
  fCard(NULL),
  fSystemAndEnergy(""),
  fCompactSystemAndEnergy(""),
  fLoadEventInformation(false),
  fLoadJets(false),
  fLoad2DHistograms(false),
  fLoadJetPtClosureHistograms(false),
  fLoadJetPtResponseMatrix(false),
  fLoadMultiplicityInJetHistograms(false),
  fLoadMaxParticlePtWithinJetConeHistograms(false),
  fLoadReflectedConeQAHistograms(false),
  fLoadJetPtUnfoldingHistograms(false),
  fLoadJetPtUnfoldingCovariance(false),
  fLoadJetPtOneDimensionalUnfoldingHistograms(false),
  fLoadTrackParticleMatchingHistograms(false),
  fJetFlavor(0),
  fLoadedWeightExponent(1),
  fLoadOnlyUnfoldedSignal(false),
  fFirstLoadedCentralityBin(0),
  fLastLoadedCentralityBin(1),
  fFirstLoadedTrackPtBin(0),
  fLastLoadedTrackPtBin(1),
  fFirstLoadedJetPtBinEEC(0),
  fLastLoadedJetPtBinEEC(1),
  fFirstLoadedTrackPtBinEEC(0),
  fLastLoadedTrackPtBinEEC(1),
  fnCentralityBins(kMaxCentralityBins),
  fnTrackPtBins(kMaxTrackPtBins),
  fnJetPtBinsEEC(kMaxJetPtBinsEEC),
  fnTrackPtBinsEEC(kMaxTrackPtBinsEEC),
  fIncludedPairingTypes(0)
{
  
  // Do not draw anything by default
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    fLoadTracks[iTrackType] = false;
  }
  
  // Do not draw the particle density histograms around the jet axis by by default
  for(int iParticleDensityType = 0; iParticleDensityType < knParticleDensityAroundJetAxisTypes; iParticleDensityType++){
    fLoadParticleDensityAroundJetsHistograms[iParticleDensityType] = false;
  }
  
  // Do not draw energy-energy correlator histograms for default
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType] = false;
  }
  
  // Default binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = iCentrality+1;
    fCentralityBinBorders[iCentrality] = 0;
  }
  
  // Default binning for track pT
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    fTrackPtBinBorders[iTrackPt] = 0;
  }
  
  // Default binning for jet pT in energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC + 1; iJetPt++){
    fJetPtIndicesEEC[iJetPt] = iJetPt+1;
    fJetPtBinBordersEEC[iJetPt] = 0;
  }
  
  // Default binning for track pT in energy-energy correlator histograms
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC + 1; iTrackPt++){
    fTrackPtIndicesEEC[iTrackPt] = iTrackPt+1;
    fTrackPtBinBordersEEC[iTrackPt] = 0;
  }
  
  // Initialize all the other histograms to null
  fhVertexZ = NULL;            // Vertex z position
  fhVertexZWeighted = NULL;    // Weighted vertex z-position (only meaningfull for MC)
  fhEvents = NULL;             // Number of events surviving different event cuts
  fhTriggers = NULL;           // Trigger selection information
  fhTrackCuts = NULL;          // Number of tracks surviving different track cuts
  fhCentrality = NULL;         // Centrality of all events
  fhCentralityWeighted = NULL; // Weighted centrality distribution in all events (only meaningful for MC)
  fhPtHat = NULL;              // pT hat for MC events (only meaningful for MC)
  fhPtHatWeighted = NULL;      // Weighted pT hat distribution (only meaningful for MC)
  
  fhMultiplicityMap = NULL;              // Multiplicity vs. centrality map
  fhMultiplicityMapWeighted = NULL;      // Efficiency weighted multiplicity vs. centrality map

  for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
    for(int iNJet = 0 ; iNJet < knMaxJetsBins+1; iNJet++){
      fhHFPlus[iJetPt][iNJet] = NULL;             // Energy in HF plus calorimeters
      fhHFMinus[iJetPt][iNJet] = NULL;            // Energy in HF minus calorimeters
      fhHFSum[iJetPt][iNJet] = NULL;              // Energy in HF calorimeters
      fhHFPlusVsHFMinus[iJetPt][iNJet] = NULL;    // 2D-map of HF plus vs. HF minus
    }
  }
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    
    fhMultiplicity[iCentrality] = NULL;               // Track multiplicity from all events
    fhMultiplicityWeighted[iCentrality] = NULL;       // Efficiency weighted track multiplicity from all events
    
    // Jet histograms
    fhJetPt[iCentrality] = NULL;      // Jet pT histograms
    fhJetPhi[iCentrality] = NULL;     // Jet phi histograms
    fhJetEta[iCentrality] = NULL;     // Jet eta histograms
    fhJetEtaCM[iCentrality] = NULL;   // Jet eta_CM histograms
    fhJetEtaPhi[iCentrality] = NULL;  // 2D eta-phi histogram for jets
    
    // Track histograms
    for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
      fhTrackPt[iTrackType][iCentrality] = NULL;   // Track pT histograms
      
      // Loop over track pT bins
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
        fhTrackPhi[iTrackType][iCentrality][iTrackPt] = NULL;    // Track phi histograms
        fhTrackEta[iTrackType][iCentrality][iTrackPt] = NULL;    // Track eta histograms
        fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] = NULL; // 2D eta-phi histogram for track
      } // Track pT loop
      
    } // Track category loop
    
    // Multiplicity and particle density histograms within the jet cone
    for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
        for(int iMultiplicityType = 0; iMultiplicityType < knMultiplicityInJetConeTypes; iMultiplicityType++){
          for(int iSubeventType = 0; iSubeventType < EECHistograms::knSubeventTypes+1; iSubeventType++){
            fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubeventType] = NULL;
          } // Subevent loop
        } // Multiplicity type loop
        for(int iJetConeType = 0; iJetConeType < EECHistograms::knJetConeTypes; iJetConeType++){
          for(int iParticleDensityType = 0; iParticleDensityType < knParticleDensityAroundJetAxisTypes; iParticleDensityType++){
            for(int iSubeventType = 0; iSubeventType < EECHistograms::knSubeventTypes+1; iSubeventType++){
              fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubeventType] = NULL;
            } // Subevent loop
          } // Particle density type loop
        } // Jet cone type loop
      } // Track pT bins for multiplicity and particle densities
    } // Jet pT bins for energy-energy correlators
    
    // Maximum particle pT within the jet cone
    for(int iMaxParticlePtType = 0; iMaxParticlePtType < knMaxParticlePtWithinJetConeTypes; iMaxParticlePtType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < knProjectedMaxParticlePtBins; iTrackPt++){
          fhMaxParticlePtInJetConePtCut[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt] = NULL;
          fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt] = NULL;
        } // Particle pT selection for the maximum particle pT within the jet
        fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][knProjectedMaxParticlePtBins] = NULL;
      } // Jet pT bins for energy-energy correlators
    } // Maximum particle pT within jet cone type loop
    
    // Energy-energy correlator histograms
    for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
          for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations+1; iSubevent++){
              fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = NULL;
            } // Subevent loop
          } // Pairing type loop (same jet/reflected cone)
        } // Track pT bins for energy-energy correlators
      } // Jet pT bins for energy-energy correlators
    } // Energy-energy correlator type loop
    
    // Post-processed energy-energy correlator histograms
    for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
          for(int iProcessingLevel = 0; iProcessingLevel < knEnergyEnergyCorrelatorProcessingLevels; iProcessingLevel++){
              fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel] = NULL;
          } // Pairing type loop (same jet/reflected cone)
        } // Track pT bins for energy-energy correlators
      } // Jet pT bins for energy-energy correlators
    } // Energy-energy correlator type loop

    // QA histograms for energy-energy correlators
    fhNumberOfJetsWithinReflectedCone[iCentrality] = NULL;
    fhJetPtWithinReflectedCone[iCentrality] = NULL;
    
    // Jet pT closure histograms
    for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
      for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){
        for(int iJetPhi = 0; iJetPhi <= knJetPhiBins; iJetPhi++){
          for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
            fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iClosureParticle] = NULL;
          } // Closure particle loop
        } // Jet phi bin loop 
      } // Jet eta bin loop
    } // Gen jet pT loop

    // Jet pT response matrix
    fhJetPtResponseMatrix[iCentrality] = NULL;

    // Histograms required in the jet pT unfolding study
    for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
      fhJetPtUnfoldingResponse[iCentrality][iTrackPt] = NULL;
      for(int iCoverianceMatrixType = 0; iCoverianceMatrixType < knCovarianceMatrixTypes; iCoverianceMatrixType++){
        fhJetPtUnfoldingCovariance[iCoverianceMatrixType][iCentrality][iTrackPt] = NULL;
      }
      for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
        fhJetPtUnfoldingDistribution[iUnfoldType][iCentrality][iTrackPt] = NULL;
      }
    }

    // Histograms required in the one dimensional jet pT unfolding study
    fhOneDimensionalJetPtUnfoldingResponse[iCentrality] = NULL;
    for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
      fhOneDimensionalJetPtUnfoldingDistribution[iUnfoldType][iCentrality] = NULL;
    }

    // Histograms for the track/particle matching study
    for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
        for(int iTrackParticleQAHistogram = 0; iTrackParticleQAHistogram < knTrackParticleMatchingQAHistograms; iTrackParticleQAHistogram++){
          fhTrackParticleMatchQA[iTrackParticleQAHistogram][iCentrality][iJetPt][iTrackPt] = NULL;
        }
        for(int iTrackParticleResponseHistogram = 0; iTrackParticleResponseHistogram < knTrackParticleMatchingResponseTypes; iTrackParticleResponseHistogram++){
          fhTrackParticleResponse[iTrackParticleResponseHistogram][iCentrality][iJetPt][iTrackPt] = NULL;
        }
        fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop

  } // Centrality loop
}

/*
 * Constructor with input file
 */
EECHistogramManager::EECHistogramManager(TFile* inputFile) :
  EECHistogramManager()
{
  fInputFile = inputFile;
  
  // Read card from inputfile
  fCard = new EECCard(inputFile);
  
  // Initialize values using the information in card
  InitializeFromCard();
  
}

/*
 * Constructor with input file and input card
 */
EECHistogramManager::EECHistogramManager(TFile* inputFile, EECCard* card) :
  EECHistogramManager()
{
  fInputFile = inputFile;
  
  // Initialize values using the information in card
  fCard = card;
  InitializeFromCard();
  
}

/*
 * Constructor with input card
 */
EECHistogramManager::EECHistogramManager(EECCard* card) :
  EECHistogramManager()
{

  // Initialize values using the information in card
  fCard = card;
  InitializeFromCard();
  
}

/*
 * Initialize several member variables from EECCard
 */
void EECHistogramManager::InitializeFromCard(){
  
  // Read the collision system from the card
  TString collisionSystem = fCard->GetDataType();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Read bins for centrality, track pT, jet pT for energy-energy correlators, and track pT for energy-energy correlators from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  fnTrackPtBins = fCard->GetNTrackPtBins();
  fnJetPtBinsEEC = fCard->GetNJetPtBinsEEC();
  fnTrackPtBinsEEC = fCard->GetNTrackPtBinsEEC();
  
  // Centrality binning
  for(int iCentrality = 0; iCentrality <= fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fLastLoadedCentralityBin = fnCentralityBins-1;
  
  // Track pT binning
  for(int iTrackPt = 0; iTrackPt <= fnTrackPtBins; iTrackPt++){
    fTrackPtBinBorders[iTrackPt] = fCard->GetLowBinBorderTrackPt(iTrackPt);
  }
  fLastLoadedTrackPtBin = fnTrackPtBins-1;
  
  // Jet pT binning for energy-energy correlators
  for(int iJetPt = 0; iJetPt <= fnJetPtBinsEEC; iJetPt++){
    fJetPtBinBordersEEC[iJetPt] = fCard->GetLowBinBorderJetPtEEC(iJetPt);
  }
  fLastLoadedJetPtBinEEC = fnJetPtBinsEEC-1;
  
  // Track pT binnins for energy-energy correlators
  for(int iTrackPt = 0; iTrackPt <= fnTrackPtBinsEEC; iTrackPt++){
    fTrackPtBinBordersEEC[iTrackPt] = fCard->GetLowBinBorderTrackPtEEC(iTrackPt);
  }
  fLastLoadedTrackPtBinEEC = fnTrackPtBinsEEC-1;
  
  // Remove centrality selection from pp data
  if(collisionSystem.Contains("pp")){
    fLastLoadedCentralityBin = 0;
    fCentralityBinBorders[0] = -0.5;
  }

  // Find which pairings for background estimations that are defined in the card
  int backgroundMethods = fCard->GetBackgroundMethods();
  std::bitset<3> backgroundBitChecker(backgroundMethods);
  bool doReflectedCone = backgroundBitChecker.test(0);
  bool doPerpendicularCone = backgroundBitChecker.test(1);
  bool doMixedCone = backgroundBitChecker.test(2);

  // Jet cone+jet cone pairs are always included
  fIncludedPairingTypes = new std::vector<int>;
  fIncludedPairingTypes->push_back(EECHistograms::kSameJetPair);

  // Only load pairings for background estimations that are defined in the card
  if(doReflectedCone){
    fIncludedPairingTypes->push_back(EECHistograms::kSignalReflectedConePair);
    fIncludedPairingTypes->push_back(EECHistograms::kReflectedConePair);
  }

  if(doPerpendicularCone){
    fIncludedPairingTypes->push_back(EECHistograms::kSignalPerpendicularConePair);
    fIncludedPairingTypes->push_back(EECHistograms::kPerpendicularConePair);
    fIncludedPairingTypes->push_back(EECHistograms::kSignalSecondPerpendicularConePair);
    fIncludedPairingTypes->push_back(EECHistograms::kPerpendicularPerpendicularConePair);
    fIncludedPairingTypes->push_back(EECHistograms::kSecondPerpendicularConePair);
  }

  if(doMixedCone){
    fIncludedPairingTypes->push_back(EECHistograms::kSignalMixedConePair);
    fIncludedPairingTypes->push_back(EECHistograms::kMixedConePair);
    fIncludedPairingTypes->push_back(EECHistograms::kSignalSecondMixedConePair);
    fIncludedPairingTypes->push_back(EECHistograms::kMixedMixedConePair);
    fIncludedPairingTypes->push_back(EECHistograms::kSecondMixedConePair);
  }
  
}

/*
 * Copy constructor
 */
EECHistogramManager::EECHistogramManager(const EECHistogramManager& in) :
  fInputFile(in.fInputFile),
  fCard(in.fCard),
  fSystemAndEnergy(in.fSystemAndEnergy),
  fCompactSystemAndEnergy(in.fCompactSystemAndEnergy),
  fLoadEventInformation(in.fLoadEventInformation),
  fLoadJets(in.fLoadJets),
  fLoad2DHistograms(in.fLoad2DHistograms),
  fLoadJetPtClosureHistograms(in.fLoadJetPtClosureHistograms),
  fLoadJetPtResponseMatrix(in.fLoadJetPtResponseMatrix),
  fLoadMultiplicityInJetHistograms(in.fLoadMultiplicityInJetHistograms),
  fLoadMaxParticlePtWithinJetConeHistograms(in.fLoadMaxParticlePtWithinJetConeHistograms),
  fLoadReflectedConeQAHistograms(in.fLoadReflectedConeQAHistograms),
  fLoadJetPtUnfoldingHistograms(in.fLoadJetPtUnfoldingHistograms),
  fLoadJetPtUnfoldingCovariance(in.fLoadJetPtUnfoldingCovariance),
  fLoadJetPtOneDimensionalUnfoldingHistograms(in.fLoadJetPtOneDimensionalUnfoldingHistograms),
  fLoadTrackParticleMatchingHistograms(in.fLoadTrackParticleMatchingHistograms),
  fJetFlavor(in.fJetFlavor),
  fLoadedWeightExponent(in.fLoadedWeightExponent),
  fLoadOnlyUnfoldedSignal(in.fLoadOnlyUnfoldedSignal),
  fFirstLoadedCentralityBin(in.fFirstLoadedCentralityBin),
  fLastLoadedCentralityBin(in.fLastLoadedCentralityBin),
  fFirstLoadedTrackPtBin(in.fFirstLoadedTrackPtBin),
  fLastLoadedTrackPtBin(in.fLastLoadedTrackPtBin),
  fFirstLoadedJetPtBinEEC(in.fFirstLoadedJetPtBinEEC),
  fLastLoadedJetPtBinEEC(in.fLastLoadedJetPtBinEEC),
  fFirstLoadedTrackPtBinEEC(in.fFirstLoadedTrackPtBinEEC),
  fLastLoadedTrackPtBinEEC(in.fLastLoadedTrackPtBinEEC),
  fIncludedPairingTypes(in.fIncludedPairingTypes),
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhTriggers(in.fhTriggers),
  fhTrackCuts(in.fhTrackCuts),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted)
{
  // Copy constructor
  
  // Copy all values
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    fLoadTracks[iTrackType] = in.fLoadTracks[iTrackType];
  }
  
  // Copy loading particle density histograms
  for(int iParticleDensityType = 0; iParticleDensityType < knParticleDensityAroundJetAxisTypes; iParticleDensityType++){
    fLoadParticleDensityAroundJetsHistograms[iParticleDensityType] = in.fLoadParticleDensityAroundJetsHistograms[iParticleDensityType];
  }
  
  // Copy loading energy-energy correlator histograms
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType] = in.fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType];
  }
  
  // Copy binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = in.fCentralityBinIndices[iCentrality];
    fCentralityBinBorders[iCentrality] = in.fCentralityBinBorders[iCentrality];
  }
  
  // Copy binning for track pT
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = in.fTrackPtBinIndices[iTrackPt];
    fTrackPtBinBorders[iTrackPt] = in.fTrackPtBinBorders[iTrackPt];
  }
  
  // Copy binning for jet pT in energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC+1; iJetPt++){
    fJetPtIndicesEEC[iJetPt] = in.fJetPtIndicesEEC[iJetPt];
    fJetPtBinBordersEEC[iJetPt] = in.fJetPtBinBordersEEC[iJetPt];
  }
  
  // Copy binning for track pT in energy-energy correlator histograms
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC+1; iTrackPt++){
    fTrackPtIndicesEEC[iTrackPt] = in.fTrackPtIndicesEEC[iTrackPt];
    fTrackPtBinBordersEEC[iTrackPt] = in.fTrackPtBinBordersEEC[iTrackPt];
  }
  
  // Multiplicity maps
  fhMultiplicityMap = in.fhMultiplicityMap;
  fhMultiplicityMapWeighted = in.fhMultiplicityMapWeighted;

  // HF energy histograms
  for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
    for(int iNJet = 0 ; iNJet < knMaxJetsBins+1; iNJet++){
      fhHFPlus[iJetPt][iNJet] = in.fhHFPlus[iJetPt][iNJet];                    // Energy in HF plus calorimeters
      fhHFMinus[iJetPt][iNJet] = in.fhHFMinus[iJetPt][iNJet];                  // Energy in HF minus calorimeters
      fhHFSum[iJetPt][iNJet] = in.fhHFSum[iJetPt][iNJet];                      // Energy in HF calorimeters
      fhHFPlusVsHFMinus[iJetPt][iNJet] = in.fhHFPlusVsHFMinus[iJetPt][iNJet];  // 2D-map of HF plus vs. HF minus
    }
  }
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    
    fhMultiplicity[iCentrality] = in.fhMultiplicity[iCentrality];                           // Track multiplicity from all events
    fhMultiplicityWeighted[iCentrality] = in.fhMultiplicityWeighted[iCentrality];           // Efficiency weighted track multiplicity from all events
    
    // Jet histograms
    fhJetPt[iCentrality] = in.fhJetPt[iCentrality];         // Jet pT histograms
    fhJetPhi[iCentrality] = in.fhJetPhi[iCentrality];       // Jet phi histograms
    fhJetEta[iCentrality] = in.fhJetEta[iCentrality];       // Jet eta histograms
    fhJetEtaCM[iCentrality] = in.fhJetEtaCM[iCentrality];   // Jet eta_CM histograms
    fhJetEtaPhi[iCentrality] = in.fhJetEtaPhi[iCentrality]; // 2D eta-phi histogram for jets

    
    // Track histograms
    for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
      fhTrackPt[iTrackType][iCentrality] = in.fhTrackPt[iTrackType][iCentrality];   // Track pT histograms
      
      // Loop over track pT bins
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
        fhTrackPhi[iTrackType][iCentrality][iTrackPt] = in.fhTrackPhi[iTrackType][iCentrality][iTrackPt];    // Track phi histograms
        fhTrackEta[iTrackType][iCentrality][iTrackPt] = in.fhTrackEta[iTrackType][iCentrality][iTrackPt];    // Track eta histograms
        fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] = in.fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt]; // 2D eta-phi histogram for track
      } // Track pT loop
      
    } // Track category loop
    
    // Multiplicity and particle density histograms within the jet cone
    for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
        for(int iMultiplicityType = 0; iMultiplicityType < knMultiplicityInJetConeTypes; iMultiplicityType++){
          for(int iSubeventType = 0; iSubeventType < EECHistograms::knSubeventTypes+1; iSubeventType++){
            fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubeventType] = in.fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubeventType];
          } // Subevent loop
        } // Multiplicity type loop
        for(int iJetConeType = 0; iJetConeType < EECHistograms::knJetConeTypes; iJetConeType++){
          for(int iParticleDensityType = 0; iParticleDensityType < knParticleDensityAroundJetAxisTypes; iParticleDensityType++){
            for(int iSubeventType = 0; iSubeventType < EECHistograms::knSubeventTypes+1; iSubeventType++){
              fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubeventType] = in.fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubeventType];
            } // Subevent loop
          } // Particle density type loop
        } // Jet cone type loop
      } // Track pT bins for multiplicity and particle densities
    } // Jet pT bins for energy-energy correlators
    
    // Maximum particle pT within the jet cone
    for(int iMaxParticlePtType = 0; iMaxParticlePtType < knMaxParticlePtWithinJetConeTypes; iMaxParticlePtType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < knProjectedMaxParticlePtBins; iTrackPt++){
          fhMaxParticlePtInJetConePtCut[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt] = in.fhMaxParticlePtInJetConePtCut[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt];
          fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt] = in.fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt];
        } // Particle pT selection for the maximum particle pT within the jet
        fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][knProjectedMaxParticlePtBins] = in.fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][knProjectedMaxParticlePtBins];
      } // Jet pT bins for energy-energy correlators
    } // Maximum particle pT within jet cone type loop
    
    // Energy-energy correlator histograms
    for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
          for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations+1; iSubevent++){
              fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = in.fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent];
            } // Subevent loop
          } // Pairing type loop
        } // Track pT bins for energy-energy correlators
      } // Jet pT bins for energy-energy correlators
    } // Energy-energy correlator type loop
    
    // Post-processed energy-energy correlator histograms
    for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
          for(int iProcessingLevel = 0; iProcessingLevel < knEnergyEnergyCorrelatorProcessingLevels; iProcessingLevel++){
              fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel] = in.fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel];
          } // Pairing type loop (same jet/reflected cone)
        } // Track pT bins for energy-energy correlators
      } // Jet pT bins for energy-energy correlators
    } // Energy-energy correlator type loop

    // QA histograms for energy-energy correlators
    fhNumberOfJetsWithinReflectedCone[iCentrality] = in.fhNumberOfJetsWithinReflectedCone[iCentrality];
    fhJetPtWithinReflectedCone[iCentrality] = in.fhJetPtWithinReflectedCone[iCentrality];
    
    // Jet pT closure histograms
    for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
      for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){
        for(int iJetPhi = 0; iJetPhi <= knJetPhiBins; iJetPhi++){
          for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
            fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iClosureParticle] = in.fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iClosureParticle];
          } // Closure particle loop
        } // Jet phi bin loop
      } // Jet eta bin loop
    } // Gen jet pT loop
    
    // Jet pT response matrix
    fhJetPtResponseMatrix[iCentrality] = in.fhJetPtResponseMatrix[iCentrality];

    // Histograms required in the jet pT unfolding study
    for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
      fhJetPtUnfoldingResponse[iCentrality][iTrackPt] = in.fhJetPtUnfoldingResponse[iCentrality][iTrackPt];
      for(int iCoverianceMatrixType = 0; iCoverianceMatrixType < knCovarianceMatrixTypes; iCoverianceMatrixType++){
        fhJetPtUnfoldingCovariance[iCoverianceMatrixType][iCentrality][iTrackPt] = in.fhJetPtUnfoldingCovariance[iCoverianceMatrixType][iCentrality][iTrackPt];;
      }
      for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
        fhJetPtUnfoldingDistribution[iUnfoldType][iCentrality][iTrackPt] = in.fhJetPtUnfoldingDistribution[iUnfoldType][iCentrality][iTrackPt];
      }
    }

    // Histograms required in the one dimensional jet pT unfolding study
    fhOneDimensionalJetPtUnfoldingResponse[iCentrality] = in.fhOneDimensionalJetPtUnfoldingResponse[iCentrality];
    for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
      fhOneDimensionalJetPtUnfoldingDistribution[iUnfoldType][iCentrality] = in.fhOneDimensionalJetPtUnfoldingDistribution[iUnfoldType][iCentrality];
    }

    // Histograms for the track/particle matching study
    for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
        for(int iTrackParticleQAHistogram = 0; iTrackParticleQAHistogram < knTrackParticleMatchingQAHistograms; iTrackParticleQAHistogram++){
          fhTrackParticleMatchQA[iTrackParticleQAHistogram][iCentrality][iJetPt][iTrackPt] = in.fhTrackParticleMatchQA[iTrackParticleQAHistogram][iCentrality][iJetPt][iTrackPt];
        }
        for(int iTrackParticleResponseHistogram = 0; iTrackParticleResponseHistogram < knTrackParticleMatchingResponseTypes; iTrackParticleResponseHistogram++){
          fhTrackParticleResponse[iTrackParticleResponseHistogram][iCentrality][iJetPt][iTrackPt] = in.fhTrackParticleResponse[iTrackParticleResponseHistogram][iCentrality][iJetPt][iTrackPt];
        }
        fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt] = in.fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt];
      } // Track pT loop
    } // Jet pT loop

  } // Centrality loop
}

/*
 * Destructor
 */
EECHistogramManager::~EECHistogramManager(){
  delete fCard;
  delete fIncludedPairingTypes;
}

/*
 * Normalize the energy-energy correlator distributions and subtract the background from them
 *
 *  int iMethod = Background subtraction method
 *                0: Use mixed event background subtraction with two mixed events
 *                1: Use ALICE style perpendicular cone subtraction
 *                2: Use reflected cone background subtraction with MC based scaling
 *                3: Do not subtract background, just copy the raw distribution
 *                4: Use generator level truth signal for the signal histogram
 *                
 *  const int iSystematic = Index for systematic uncertainty estimation for background subtraction.
 *                          0: Nominal results, no systematic uncertainty estimation
 *                          1: Systematic uncertainty derived from 2% centrality shifted simulation
 *                          2: Systematic uncertainty derived from 6% centrality shifted simulation
 */
void EECHistogramManager::SubtractBackground(int iMethod, const int iSystematic){

  // Sanity check for selected background subtraction method
  if(iMethod < 0) iMethod = 0;
  if(iMethod > 4) iMethod = 4;
  
  double normalizationFactor;
  EECBackgroundScale* scaleProvider = new EECBackgroundScale(fCard, iSystematic);
  
  // Bin borders that are searched from the background scaler
  std::pair<double,double> centralityBinBorders;
  std::pair<double,double> jetPtBinBorders;
  double trackPtLowBorder;  // We only care about the lower border in track pT bins

  // Helper histograms to help that we do not need to load all histograms to memory
  TH1D* rawHistogram;
  TH1D* signalFakeHistogram1;
  TH1D* signalFakeHistogram2;
  TH1D* sameConeFakeFakeHistogram1;
  TH1D* sameConeFakeFakeHistogram2;
  TH1D* mixedConeFakeFakeHistogram;
  TH1D* trueFakeBackground;

  bool isPbPb = fSystemAndEnergy.Contains("PbPb");
  
  // Loop over the energy-energy correlator histograms from the input file
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    
    // Only subtract background from selected types of histograms
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    // Loop over selected bin range
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      centralityBinBorders = fCard->GetBinBordersCentrality(iCentrality);
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        trackPtLowBorder = fCard->GetLowBinBorderTrackPtEEC(iTrackPt);
        jetPtBinBorders.first = fCard->GetLowBinBorderJetPtEEC(0);
        jetPtBinBorders.second = fCard->GetHighBinBorderJetPtEEC(fCard->GetNJetPtBinsEEC()-1);

        // =================================== //
        // Histograms without jet pT selection //
        // =================================== //
        
        // First find the correct histogram and normalize it to one
        rawHistogram = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::knSubeventCombinations);
        fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorNormalized] = (TH1D*) rawHistogram->Clone(Form("%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorNormalized], iCentrality, iTrackPt));
        
        normalizationFactor = fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorNormalized]->Integral("width");
        fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorNormalized]->Scale(1/normalizationFactor);

        // Mixed and perpendicular cone background subtraction methods. In these method, jet cone is placed on the location of the original jet in a minimum bias events, or the jet cone is rotated 90 degrees in phi. Particles from original jet cone are paired with particles in the mixed/perpendicular cone. For this pairing, the fake+fake pairs will be mistreated since local correlations there are killed, and pairs are double counted. To correct for this, we find another mixed event/perpendicular cone, subtract the pairing between two different mixed/perpendicular cones, and add back pairings from a single mixed/perpendicular cone. This gives the most accurate estimation of the background.
        if(iMethod < 2){

          signalFakeHistogram1 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kSignalMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
          signalFakeHistogram2 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kSignalSecondMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
          sameConeFakeFakeHistogram1 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
          sameConeFakeFakeHistogram2 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kSecondMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
          mixedConeFakeFakeHistogram = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kMixedMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);

          // Since pairing signal particles with mixed cone particles should give the same results regardless if we use the first or second mixed cones to do the pairing, average these two in order to suppress fluctuations
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground] = (TH1D*) signalFakeHistogram1->Clone(Form("%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorBackground], iCentrality, iTrackPt));
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground]->Add(signalFakeHistogram2);
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground]->Scale(0.5);

          // Since mixed/perpendicular cone pairs in two used mixed events/perpendicular cones both have the same information, get the average of them to suppress fluctuations in the background
          trueFakeBackground = (TH1D*) sameConeFakeFakeHistogram1->Clone(Form("trueFakeBackground_%d%d%d", iEnergyEnergyCorrelatorType, iCentrality, iTrackPt));
          trueFakeBackground->Add(sameConeFakeFakeHistogram2);
          trueFakeBackground->Scale(0.5);

          // Add the true fake background and subtract the fake fake background
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground]->Add(trueFakeBackground);
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground]->Add(mixedConeFakeFakeHistogram, -1);

          // Do not normalize the background distribution. This allows to do different normalizations in the end
          //fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground]->Scale(1/normalizationFactor);
          
        } 
         
        // Reflected cone background subtraction method. In this method, background is estimated by placing another jet cone in the event with reflecting the eta-value. Around midrapidity, the jet cone is shifted by 2R instead of reflecting to avoid overlap of the cone. Particles from signal cone are paired with the reflected eta cone. We know from simulation that the shape of the background distribution is correct in case where fake+fake background is negligible. We can correct for the normalization of the background by MC-based scaling factor. This gives a good background subtraction in a kinematic region where fake+fake background is negligible.
        else if (iMethod == 2){

          signalFakeHistogram1 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kSignalReflectedConePair, EECHistograms::knSubeventCombinations);

          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground] = (TH1D*) signalFakeHistogram1->Clone(Form("%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorBackground], iCentrality, iTrackPt));

          // Do not normalize here. This allows to do different normalizations later
          //fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground]->Scale(1/normalizationFactor);

          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground]->Scale(1/scaleProvider->GetEECBackgroundScale(centralityBinBorders, jetPtBinBorders, trackPtLowBorder, isPbPb));

        }
        
        // Subtracts the background from the raw distribution
        if(iMethod < 4){
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorSignal] = (TH1D*) rawHistogram->Clone(Form("%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorSignal], iCentrality, iTrackPt));
        } else {
          // For Monte Carlo, we can directly read the truth level distributions
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorSignal] =  (TH1D*) GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaPythia);
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorSignal]->SetName(Form("%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorSignal], iCentrality, iTrackPt));
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground] = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaHydjet);
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground]->SetName(Form("%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorBackground], iCentrality, iTrackPt));
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground]->Add(GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, fnJetPtBinsEEC, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kHydjetHydjet));
        }

        // Subtract the background unless specifically defined not to do that
        if(iMethod < 3){
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorSignal]->Add(fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][kEnergyEnergyCorrelatorBackground], -1);
        }
        
        // =================================== //
        // Histograms in different jet pT bins //
        // =================================== //
        
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
          jetPtBinBorders = fCard->GetBinBordersJetPtEEC(iJetPt);
          
          // First find the correct histogram and normalize it to one
          rawHistogram = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::knSubeventCombinations);
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorNormalized] = (TH1D*) rawHistogram->Clone(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorNormalized], iCentrality, iTrackPt, iJetPt));
          
          normalizationFactor = fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorNormalized]->Integral("width");
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorNormalized]->Scale(1/normalizationFactor);
          
          // Mixed and perpendicular cone background subtraction methods. In these method, jet cone is placed on the location of the original jet in a minimum bias events, or the jet cone is rotated 90 degrees in phi. Particles from original jet cone are paired with particles in the mixed/perpendicular cone. For this pairing, the fake+fake pairs will be mistreated since local correlations there are killed, and pairs are double counted. To correct for this, we find another mixed event/perpendicular cone, subtract the pairing between two different mixed/perpendicular cones, and add back pairings from a single mixed/perpendicular cone. This gives the most accurate estimation of the background.
          if(iMethod < 2){

            signalFakeHistogram1 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
            signalFakeHistogram2 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalSecondMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
            sameConeFakeFakeHistogram1 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
            sameConeFakeFakeHistogram2 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSecondMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
            mixedConeFakeFakeHistogram = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kMixedMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);

            // Since pairing signal particles with mixed cone particles should give the same results regardless if we use the first or second mixed cones to do the pairing, average these two in order to suppress fluctuations
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground] = (TH1D*) signalFakeHistogram1->Clone(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorBackground], iCentrality, iTrackPt, iJetPt));
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground]->Add(signalFakeHistogram2);
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground]->Scale(0.5);

            // Since mixed/perpendicular cone pairs in two used mixed events/perpendicular cones both have the same information, get the average of them to suppress fluctuations in the background
            trueFakeBackground = (TH1D*) sameConeFakeFakeHistogram1->Clone(Form("trueFakeBackground_%d%d%d%d", iEnergyEnergyCorrelatorType, iCentrality, iTrackPt, iJetPt));
            trueFakeBackground->Add(sameConeFakeFakeHistogram2);
            trueFakeBackground->Scale(0.5);

            // Add the true fake background and subtract the fake fake background
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground]->Add(trueFakeBackground);
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground]->Add(mixedConeFakeFakeHistogram, -1);

            // Do not normalize background. This allows to do different normalizations later
            //fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground]->Scale(1/normalizationFactor);
          } 
         
          // Reflected cone background subtraction method. In this method, background is estimated by placing another jet cone in the event with reflecting the eta-value. Around midrapidity, the jet cone is shifted by 2R instead of reflecting to avoid overlap of the cone. Particles from signal cone are paired with the reflected eta cone. We know from simulation that the shape of the background distribution is correct in case where fake+fake background is negligible. We can correct for the normalization of the background by MC-based scaling factor. This gives a good background subtraction in a kinematic region where fake+fake background is negligible.
          else if (iMethod == 2){

            signalFakeHistogram1 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalReflectedConePair, EECHistograms::knSubeventCombinations);

            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground] = (TH1D*) signalFakeHistogram1->Clone(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorBackground], iCentrality, iTrackPt, iJetPt));

            // Do not normalize background. This allows to do different normalizations later
            //fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground]->Scale(1/normalizationFactor);

            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground]->Scale(1/scaleProvider->GetEECBackgroundScale(centralityBinBorders, jetPtBinBorders, trackPtLowBorder, isPbPb));

          }
          
          // Now that the background is properly normalized, it can be subtracted from the total distribution to get the signal
          if(iMethod < 4){
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorSignal] = (TH1D*) rawHistogram->Clone(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorSignal], iCentrality, iTrackPt, iJetPt));
          } else {
            // For Monte Carlo, we can directly read the truth level distributions
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorSignal] =  GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaPythia);
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorSignal]->SetName(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorSignal], iCentrality, iTrackPt, iJetPt));
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground] = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaHydjet);
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground]->SetName(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorBackground], iCentrality, iTrackPt, iJetPt));
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground]->Add(GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kHydjetHydjet));
          }

          // Subtract the background unless specifically instructed not to do that
          if(iMethod < 3){
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorSignal]->Add(fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackground], -1);
          }
          
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
}

/*
 * Subtract the background from unfolded histograms
 *
 *  int iMethod = Background subtraction method
 *                0: Use mixed event background subtraction with two mixed events
 *                1: Use ALICE style perpendicular cone subtraction
 *                2: Use reflected cone background subtraction with MC based scaling
 *  const int iSystematic = Index for systematic uncertainty estimation for background subtraction.
 *                          0: Nominal results, no systematic uncertainty estimation
 *                          1: Systematic uncertainty derived from 2% centrality shifted simulation
 *                          2: Systematic uncertainty derived from 6% centrality shifted simulation
 *                          3: Systematic uncertainty by lower estimate for unfolded signal-to-background ratio scaling
 *                          4: Systematic uncertainty by higher estimate for unfolded signal-to-background ratio scaling
 */
void EECHistogramManager::SubtractBackgroundFromUnfolded(int iMethod, const int iSystematic){
  
  // Sanity check for selected background subtraction method
  if(iMethod < 0) iMethod = 0;
  if(iMethod > 2) iMethod = 2;

  double scalingFactor;
  EECBackgroundScale* scaleProvider = new EECBackgroundScale(true, iSystematic, fCard->GetWeightExponent()); // Always use gen level correction for unfolded histograms
  EECSignalToBackgroundUnfoldingScale* signalToBackgroundScaleProvider = new EECSignalToBackgroundUnfoldingScale(fCard->GetWeightExponent());
  bool isPbPb = fSystemAndEnergy.Contains("PbPb");
  
  // Bin borders that are searched from the background scaler
  std::pair<double,double> centralityBinBorders;
  std::pair<double,double> jetPtBinBorders;
  double trackPtLowBorder;  // We only care about the lower border in track pT bins
  TH1D* trueFakeBackground; // Properly estimated fake-fake background

  // Helper histograms to help that we do not need to load all histograms to memory
  TH1D* rawHistogram;
  TH1D* unfoldedHistogram;
  TH1D* signalFakeHistogram1;
  TH1D* signalFakeHistogram2;
  TH1D* sameConeFakeFakeHistogram1;
  TH1D* sameConeFakeFakeHistogram2;
  TH1D* mixedConeFakeFakeHistogram;

  // Loop over the energy-energy correlator histograms from the input file
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    
    // Only subtract background from selected types of histograms
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    // Loop over selected bin range
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      centralityBinBorders = fCard->GetBinBordersCentrality(iCentrality);
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        trackPtLowBorder = fCard->GetLowBinBorderTrackPtEEC(iTrackPt);
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
          jetPtBinBorders = fCard->GetBinBordersJetPtEEC(iJetPt);
          
          // First, clone the unfolded distribution
          unfoldedHistogram = GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, kEnergyEnergyCorrelatorUnfolded);
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorUnfoldedSignal] = (TH1D*) unfoldedHistogram->Clone(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorUnfoldedSignal], iCentrality, iTrackPt, iJetPt));
          
          // Next, find the scaling factor for background from integrals of the unfolded, and the non-unfolded distributions
          rawHistogram = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::knSubeventCombinations);
          scalingFactor = unfoldedHistogram->Integral("width") / rawHistogram->Integral("width");
          
          // Mixed and perpendicular cone background subtraction methods. In these method, jet cone is placed on the location of the original jet in a minimum bias events, or the jet cone is rotated 90 degrees in phi. Particles from original jet cone are paired with particles in the mixed/perpendicular cone. For this pairing, the fake+fake pairs will be mistreated since local correlations there are killed, and pairs are double counted. To correct for this, we find another mixed event/perpendicular cone, subtract the pairing between two different mixed/perpendicular cones, and add back pairings from a single mixed/perpendicular cone. This gives the most accurate estimation of the background.
          if(iMethod < 2){

            // Load the necessary histograms to memory
            signalFakeHistogram1 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
            signalFakeHistogram2 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalSecondMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
            sameConeFakeFakeHistogram1 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
            sameConeFakeFakeHistogram2 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSecondMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);
            mixedConeFakeFakeHistogram = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kMixedMixedConePair + 5*iMethod, EECHistograms::knSubeventCombinations);

            // Since pairing signal particles with mixed cone particles should give the same results regardless if we use the first or second mixed cones to do the pairing, average these two in order to suppress fluctuations
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding] = (TH1D*) signalFakeHistogram1->Clone(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorBackgroundAfterUnfolding], iCentrality, iTrackPt, iJetPt));
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding]->Add(signalFakeHistogram2);
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding]->Scale(0.5);

            // Since mixed/perpendicular cone pairs in two used mixed events/perpendicular cones both have the same information, get the average of them to suppress fluctuations in the background
            trueFakeBackground = (TH1D*) sameConeFakeFakeHistogram1->Clone(Form("trueFakeBackground_%d%d%d%d", iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt));
            trueFakeBackground->Add(sameConeFakeFakeHistogram2);
            trueFakeBackground->Scale(0.5);

            // Add the true fake background and subtract the fake fake background
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding]->Add(trueFakeBackground);
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding]->Add(mixedConeFakeFakeHistogram, -1);
          }

          // Reflected cone background subtraction method. In this method, background is estimated by placing another jet cone in the event with reflecting the eta-value. Around midrapidity, the jet cone is shifted by 2R instead of reflecting to avoid overlap of the cone. Particles from signal cone are paired with the reflected eta cone. We know from simulation that the shape of the background distribution is correct in case where fake+fake background is negligible. We can correct for the normalization of the background by MC-based scaling factor. This gives a good background subtraction in a kinematic region where fake+fake background is negligible.
          else if (iMethod == 2){

            signalFakeHistogram1 = GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalReflectedConePair, EECHistograms::knSubeventCombinations);

            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding] = (TH1D*) signalFakeHistogram1->Clone(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorBackgroundAfterUnfolding], iCentrality, iTrackPt, iJetPt));
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding]->Scale(1/scaleProvider->GetEECBackgroundScale(centralityBinBorders, jetPtBinBorders, trackPtLowBorder, isPbPb));
          }

          // Scale the reflected cone background estimate with the difference in yields before and after unfolding to keep to total background level correct
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding]->Scale(scalingFactor);

          // Subtract the background unless specifically instructed not to do that
          if(iMethod < 3){

            // We also need to scale the reflected cone background estimate with the difference in signal to background
            // ratios before and after unfolding in order to not oversubtract the background
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding]->Scale(signalToBackgroundScaleProvider->GetEECSignalToBackgroundUnfoldingScale(centralityBinBorders, jetPtBinBorders, trackPtLowBorder, isPbPb, iSystematic));
          
            // Now that the background is properly normalized, it can be subtracted from the unfolded distribution to get the signal
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorUnfoldedSignal]->Add(fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorBackgroundAfterUnfolding], -1);
          }
          
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
}

/*
 * Stabilize the background histograms by combining all jet pT bins for histograms that do not depend on jet pT
 */
void EECHistogramManager::StabilizeBackground(){

  // Helper variables
  double targetIntegral;
  int lowIntegralBin, highIntegralBin;

  // Loop over the energy-energy correlator histograms from the input file
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){

    // Only subtract background from selected types of histograms
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    // Loop over selected bin range
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        // Do the stabilization only for background+background contributions, which are for sure not dependent of jet pT
        for(int iPairingType : *fIncludedPairingTypes){
          if(iPairingType < EECHistograms::kReflectedConePair) continue;
          if(iPairingType == EECHistograms::kSignalMixedConePair) continue;
          if(iPairingType == EECHistograms::kSignalSecondMixedConePair) continue;
          if(iPairingType >= EECHistograms::kSignalPerpendicularConePair) continue; // No stabilization for perpendicular cone for now
          for(int iSubevent = 0; iSubevent <= EECHistograms::knSubeventCombinations; iSubevent++){

            // Subevent binning is only relevant for PbPb and pPb MC. Skip this for all other systems
            if(!fSystemAndEnergy.Contains("Pb MC") && iSubevent < EECHistograms::knSubeventCombinations) continue;

            // After the summation has been performed, put a properly scaled sum histogram in each jet pT bin as the new background histogram
            lowIntegralBin = fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent]->GetXaxis()->FindBin(0.008);
            highIntegralBin = fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent]->GetXaxis()->FindBin(0.39);
            for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
              targetIntegral = fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent]->Integral(lowIntegralBin, highIntegralBin, "width");

              // Only do the sacling if there is some content in the original histogram
              if(targetIntegral > 0){
                fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = (TH1D*) fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent]->Clone(Form("stabilizedBackground%d%d%d%d%d%d", iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, iPairingType, iSubevent));
                fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent]->Scale(targetIntegral / fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent]->Integral(lowIntegralBin, highIntegralBin, "width"));
              }
            } // Jet pT loop
          } // Subevent loop
        } // Pairing type loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
}

/*
 * Combine mixed cone histograms from other histogram manager to those in this histogram manager
 *
 *  Arguments:
 *   EECHistogramManager* anotherManager = Histogram manager from which the combined mixed cone histograms are found
 *
 */
void EECHistogramManager::CombineMixedConeBackgrounds(EECHistogramManager* anotherManager){

  // Matching bins from the other histogram manager
  int matchingCentralityBin;
  int matchingTrackPtBin;
  int matchingJetPtBin;

  // We need EECCard related to the other histogram manager to be sure we match the kinematic bins
  EECCard *anotherCard = anotherManager->GetCard();

  // Loop over the energy-energy correlator histograms
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){

    // Only combine the selected types of histograms
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;

    // Loop over all the kinematic bins and find matching bins from the other histogram manager
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Find matching centrality bin. Throw an error if none is found
      matchingCentralityBin = anotherCard->FindBinIndexCentrality(fCard->GetBinBordersCentrality(iCentrality));
      if(matchingCentralityBin == -1){
        throw std::invalid_argument(Form("EECHistogramManager::ERROR When combining mixed cone histograms, did not find matching centrality bin for %.0f-%.0f!", fCard->GetBinBordersCentrality(iCentrality).first, fCard->GetBinBordersCentrality(iCentrality).second));
      } 

      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

        // Find matching centrality bin. Throw an error if none is found
        matchingTrackPtBin = anotherCard->GetBinIndexTrackPtEEC(fCard->GetLowBinBorderTrackPtEEC(iTrackPt));
        if(matchingTrackPtBin == -1){
          throw std::invalid_argument(Form("EECHistogramManager::ERROR When combining mixed cone histograms, did not find matching track bin for %.1f!", fCard->GetLowBinBorderTrackPtEEC(iTrackPt)));
        } 

        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){

          // Find matching jet pT bin. Throw an error if none is found.
          matchingJetPtBin = anotherCard->FindBinIndexJetPtEEC(fCard->GetBinBordersJetPtEEC(iJetPt));
          if(matchingJetPtBin == -1){
            throw std::invalid_argument(Form("EECHistogramManager::ERROR When combining mixed cone histograms, did not find matching jet pT bin for %.0f-%.0f!", fCard->GetBinBordersJetPtEEC(iJetPt).first, fCard->GetBinBordersJetPtEEC(iJetPt).second));
          }

          // Loop over all the mixed event related pairing types
          for(int iPairingType = EECHistograms::kSignalMixedConePair; iPairingType < EECHistograms::knPairingTypes; iPairingType++){

            // Combining histograms with this method assumes that both histograms have similar number of events
            // More complicated weighting system is needed if the histograms have different numbers of events.
            fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations]->Add(anotherManager->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, matchingCentralityBin, matchingJetPtBin, matchingTrackPtBin, iPairingType));
            fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations]->Scale(0.5);

          } // Pairing type loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop

}

/*
 *  Normalize each bin in a histogram histogram to DeltaR bin area
 *
 *   Arguments: TH1D* histogramInNeedOfNormalization = Histogram for which the normalization is performed
 */
void EECHistogramManager::NormalizeToDeltaRBinArea(TH1D* histogramInNeedOfNormalization){
 
  // Variables needed to calculate the density from particle number
  double binLowBorder, binHighBorder;
  double binLowArea, binHighArea, binArea;
  double newBinContent, newBinError;
  
  // Loop over the bins and normalize each bin to the bin area
  for(int iBin = 1; iBin <= histogramInNeedOfNormalization->GetNbinsX(); iBin++){
    
    // Calculate the bin area
    binLowBorder = histogramInNeedOfNormalization->GetXaxis()->GetBinLowEdge(iBin);
    binHighBorder = histogramInNeedOfNormalization->GetXaxis()->GetBinUpEdge(iBin);
    binLowArea = TMath::Pi() * binLowBorder * binLowBorder;
    binHighArea = TMath::Pi() * binHighBorder * binHighBorder;
    binArea = binHighArea - binLowArea;
    
    // Calculate the new bin content and error for particle density
    newBinContent = histogramInNeedOfNormalization->GetBinContent(iBin) / binArea;
    newBinError = histogramInNeedOfNormalization->GetBinError(iBin) / binArea;
    histogramInNeedOfNormalization->SetBinContent(iBin, newBinContent);
    histogramInNeedOfNormalization->SetBinError(iBin, newBinError);
    
  } // Bin loop for normalizing particle densities
  
}

/*
 * Load all the selected histograms from the inputfile
 */
void EECHistogramManager::LoadHistograms(){

  // First, do a sanity check for the weight exponent
  CheckWeightExponent();
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhTriggers = (TH1D*) fInputFile->Get("triggers");                      // Trigger selection information
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");                    // Number of tracks surviving different track cuts
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events
    
    // Load HF energy histograms
    LoadHFEnergyHistograms();

    // Load multiplicity histograms
    LoadMultiplicityHistograms();
  }
  
  // Load jet histograms
  LoadJetHistograms();
  
  // Load track histograms
  LoadTrackHistograms();
  
  // Load multiplicity within the jet cone histograms
  LoadMultiplicityInJetConeHistograms();
  
  // Load the particle density around the jet axis histograms
  LoadParticleDensityHistograms();
  
  // Load the maximum particle pT inside the jet cone histograms
  LoadMaxParticlePtInJetConeHistograms();
  
  // Load energy-energy correlator histograms
  LoadEnergyEnergyCorrelatorHistograms();

  // Stabilize the background energy-energy correlator histograms
  StabilizeBackground();

  // Load reflected cone QA histograms
  LoadReflectedConeQAHistograms();
  
  // Load jet pT response matrices
  LoadJetPtResponseMatrix();

  // Load jet pT closure histograms
  LoadJetPtClosureHistograms();

  // Load the jet pT unfolding histograms
  LoadJetPtUnfoldingHistograms();

  // Load the covariance histogram that is used as an input for unfolding algorithm
  LoadJetPtUnfoldingCovariance();

  // Load the one dimensional jet pT unfolding histograms
  LoadJetPtOneDimensionalUnfoldingHistograms();

  // Load the track/particle matching histograms
  LoadTrackParticleMatchingHistograms();

  // Update the weight exponent to the card if necessary
  UpdateWeightExponent();
  
}

/*
 * Loader for HF energy histograms
 *
 * THnSparse for HF energy:
 *
 *     Histogram name        Axis index            Content of axis
 * --------------------------------------------------------------------------
 *        hfEnergy             Axis 0               Energy in HF plus
 *        hfEnergy             Axis 1              Energy in HF minus
 *        hfEnergy             Axis 2            Energy in HF calorimeters
 *        hfEnergy             Axis 3                Leading jet pT
 *        hfEnergy             Axis 4             Number of high-pT jets
 */
void EECHistogramManager::LoadHFEnergyHistograms(){
    
  // Define helper variables
  int duplicateRemoverJetPt = -1;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;
  THnSparseD* histogramArray;

  // Open the multidimensional histogram from which the histograms are projected
  histogramArray = (THnSparseD*) fInputFile->Get("hfEnergy");

  // Get the number of jet number bins from the file
  TH1D* hBinner;
  hBinner = FindHistogram(histogramArray,4);
  int nJetBins = hBinner->GetNbinsX();
  if(nJetBins > knMaxJetsBins) nJetBins = knMaxJetsBins;

  // Do not project NJet bins that have zero content.
  std::vector<bool> zeroSuppression;
  double thisBinContent;
  for(int iNJetBin = 1; iNJetBin <= nJetBins; iNJetBin++){
    thisBinContent = hBinner->GetBinContent(iNJetBin);
    zeroSuppression.push_back(thisBinContent == 0);
  }

  // Load the histograms without jet pT or nJets selection
  fhHFPlus[fnJetPtBinsEEC][knMaxJetsBins] = FindHistogram(histogramArray, 0);
  fhHFMinus[fnJetPtBinsEEC][knMaxJetsBins] = FindHistogram(histogramArray, 1);
  fhHFSum[fnJetPtBinsEEC][knMaxJetsBins] = FindHistogram(histogramArray, 2);
  fhHFPlusVsHFMinus[fnJetPtBinsEEC][knMaxJetsBins] = FindHistogram2D(histogramArray, 0, 1);

  // Loop over jet pT bins
  for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){

    // Reset the ranges for all the axes in the histogram array
    for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
      histogramArray->GetAxis(iAxis)->SetRange(0, 0);
    }
    
    // Select the jet pT bin indices
    lowerJetPtBin = fJetPtIndicesEEC[iJetPt];
    higherJetPtBin = fJetPtIndicesEEC[iJetPt+1]+duplicateRemoverJetPt;
          
    // Load the histograms with jet pT selection
    fhHFPlus[iJetPt][knMaxJetsBins] = FindHistogram(histogramArray, 0, 3, lowerJetPtBin, higherJetPtBin);
    fhHFMinus[iJetPt][knMaxJetsBins] = FindHistogram(histogramArray, 1, 3, lowerJetPtBin, higherJetPtBin);
    fhHFSum[iJetPt][knMaxJetsBins] = FindHistogram(histogramArray, 2, 3, lowerJetPtBin, higherJetPtBin);
    fhHFPlusVsHFMinus[iJetPt][knMaxJetsBins] = FindHistogram2D(histogramArray, 0, 1, 3, lowerJetPtBin, higherJetPtBin);

    // Loop over number of high-pT jets in the event
    for(int iNJet = 0; iNJet < nJetBins; iNJet++){

      if(zeroSuppression.at(iNJet)) continue; // Suppress zero content bins

      // Load the histograms with jet pT and number of jets selection
      fhHFPlus[iJetPt][iNJet] = FindHistogram(histogramArray, 0, 3, lowerJetPtBin, higherJetPtBin, 4, iNJet+1, iNJet+1);
      fhHFMinus[iJetPt][iNJet] = FindHistogram(histogramArray, 1, 3, lowerJetPtBin, higherJetPtBin, 4, iNJet+1, iNJet+1);
      fhHFSum[iJetPt][iNJet] = FindHistogram(histogramArray, 2, 3, lowerJetPtBin, higherJetPtBin, 4, iNJet+1, iNJet+1);
      fhHFPlusVsHFMinus[iJetPt][iNJet] = FindHistogram2D(histogramArray, 0, 1, 3, lowerJetPtBin, higherJetPtBin, 4, iNJet+1, iNJet+1);

    }
    
  } // Jet pT loop

  // Reset the ranges for all the axes in the histogram array
  for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
    histogramArray->GetAxis(iAxis)->SetRange(0, 0);
  }

  // Load the histograms without jet pT selection, but with number of jet selection:
  for(int iNJet = 0; iNJet < nJetBins; iNJet++){

    if(zeroSuppression.at(iNJet)) continue; // Suppress zero content bins

    // Load the histograms with jet pT and number of jets selection
    fhHFPlus[fnJetPtBinsEEC][iNJet] = FindHistogram(histogramArray, 0, 4, iNJet+1, iNJet+1);
    fhHFMinus[fnJetPtBinsEEC][iNJet] = FindHistogram(histogramArray, 1, 4, iNJet+1, iNJet+1);
    fhHFSum[fnJetPtBinsEEC][iNJet] = FindHistogram(histogramArray, 2, 4, iNJet+1, iNJet+1);
    fhHFPlusVsHFMinus[fnJetPtBinsEEC][iNJet] = FindHistogram2D(histogramArray, 0, 1, 4, iNJet+1, iNJet+1);

  }

}

/*
 * Loader for multiplicity histograms
 *
 * THnSparse for multiplicity:
 *
 *     Histogram name        Axis index            Content of axis
 * --------------------------------------------------------------------------
 *      multiplicity           Axis 0             Track multiplicity
 *      multiplicity           Axis 1       Efficiency weighted multiplicity
 *      multiplicity           Axis 2                 Centrality
 */
void EECHistogramManager::LoadMultiplicityHistograms(){
    
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  THnSparseD* histogramArray;

  // Open the multidimensional histogram from which the histograms are projected
  histogramArray = (THnSparseD*) fInputFile->Get("multiplicity");
  
  // Multiplicity histograms have more denser centrality binning than other histograms, so need to determine them separately
  TH1D* hBinner;
  hBinner = FindHistogram(histogramArray,2,0,0,0);
  double* centralityBinIndicesMultiplicity = new double[fnCentralityBins+1];
  for(int iBin = 0; iBin < fnCentralityBins+1; iBin++){
    centralityBinIndicesMultiplicity[iBin] = hBinner->GetXaxis()->FindBin(fCentralityBinBorders[iBin]);
  }
  
  for(int iCentralityBin = 0; iCentralityBin <= fnCentralityBins; iCentralityBin++){
    
    // Select the centrality bin indices
    if(iCentralityBin < fnCentralityBins){
      lowerCentralityBin = centralityBinIndicesMultiplicity[iCentralityBin];
      higherCentralityBin = centralityBinIndicesMultiplicity[iCentralityBin+1]+duplicateRemoverCentrality;
    } else {
      // Load the histograms without centrality selection to the index fnCentralityBins
      lowerCentralityBin = centralityBinIndicesMultiplicity[0];
      higherCentralityBin = centralityBinIndicesMultiplicity[fnCentralityBins]+duplicateRemoverCentrality;
    }
    
    fhMultiplicity[iCentralityBin] = FindHistogram(histogramArray, 0, 2, lowerCentralityBin, higherCentralityBin);
    fhMultiplicityWeighted[iCentralityBin] = FindHistogram(histogramArray, 1, 2, lowerCentralityBin, higherCentralityBin);
    
  } // Centrality loop
    
  // Reset the ranges for all the axes in the histogram array
  for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
    histogramArray->GetAxis(iAxis)->SetRange(0, 0);
  }

  fhMultiplicityMap = FindHistogram2D(histogramArray, 0, 2, 1, 0, 0, 0);
  fhMultiplicityMapWeighted = FindHistogram2D(histogramArray, 1, 2, 0, 0, 0, 0);

  // Remove the created double array
  delete[] centralityBinIndicesMultiplicity;
}


/*
 * Loader for jet histograms
 *
 * THnSparse for jets:
 *
 *   Histogram name: inclusiveJet
 *
 *     Axis index           Content of axis
 * --------------------------------------------------------
 *       Axis 0                 Jet pT
 *       Axis 1                 Jet phi
 *       Axis 2                 Jet eta
 *       Axis 3               Jet eta_CM
 *       Axis 4               Centrality
 *       Axis 5   Jet flavor (for MC) 1 = Quark, 2 = Gluon
 *
 */
void EECHistogramManager::LoadJetHistograms(){
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  THnSparseD* histogramArray;
  
  // Define arrays to help find the histograms
  int axisIndices[2] = {0};
  int lowLimits[2] = {0};
  int highLimits[2] = {0};
  
  int nAxes = 1;           // Number of constraining axes for this iteration
  
  // Open the multidimensional histogram from which the histograms are projected
  histogramArray = (THnSparseD*) fInputFile->Get(fJetHistogramName);

  // Backwards compatibility: older files only have 4 axes in jet histogram
  int backwardCompatibilityGiver = 0;
  if(histogramArray->GetNdimensions() == 4) backwardCompatibilityGiver = -1;
  
  for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
    higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
    
    axisIndices[0] = 4+backwardCompatibilityGiver; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;  // Centrality
    
    // If jet flavor is specified, only read jets of specific flavor
    if(fJetFlavor == 1 || fJetFlavor == 2){
      nAxes++;
      axisIndices[1] = 5+backwardCompatibilityGiver; lowLimits[1] = fJetFlavor; highLimits[1] = fJetFlavor;  // Jet flavor
    }
    
    // Always load jet pT histograms
    fhJetPt[iCentralityBin] = FindHistogram(histogramArray,0,nAxes,axisIndices,lowLimits,highLimits);
    
    if(!fLoadJets) continue;  // Only load the remaining jet histograms if selected
    
    fhJetPhi[iCentralityBin] = FindHistogram(histogramArray,1,nAxes,axisIndices,lowLimits,highLimits);
    fhJetEta[iCentralityBin] = FindHistogram(histogramArray,2,nAxes,axisIndices,lowLimits,highLimits);
    if(backwardCompatibilityGiver == 0) fhJetEtaCM[iCentralityBin] = FindHistogram(histogramArray,3,nAxes,axisIndices,lowLimits,highLimits);
    if(fLoad2DHistograms) fhJetEtaPhi[iCentralityBin] = FindHistogram2D(histogramArray,1,2,nAxes,axisIndices,lowLimits,highLimits);
    
  } // Loop over centrality bins
}

/*
 * Loader for track histograms
 *
 * THnSparse for tracks:
 *
 *   Histogram name: track/trackUncorrected
 *
 *     Axis index       Content of axis
 * ----------------------------------------
 *       Axis 0            Track pT
 *       Axis 1            Track phi
 *       Axis 2            Track eta
 *       Axis 3            Centrality
 */
void EECHistogramManager::LoadTrackHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[2] = {0};
  int lowLimits[2] = {0};
  int highLimits[2] = {0};
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  THnSparseD* histogramArray;
  
  // Loop over all track histograms
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only load the selected track types

    // Open the multidimensional histogram from which the histograms are projected
    histogramArray = (THnSparseD*) fInputFile->Get(fTrackHistogramNames[iTrackType]);
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){

      // Reset the ranges for all the axes in the histogram array
      for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
        histogramArray->GetAxis(iAxis)->SetRange(0, 0);
      }
      
      // Select the bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
      higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
      
      // Setup axes with restrictions, (3 = centrality)
      axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
      
      fhTrackPt[iTrackType][iCentralityBin] = FindHistogram(histogramArray,0,1,axisIndices,lowLimits,highLimits);
      fhTrackPhi[iTrackType][iCentralityBin][fnTrackPtBins] = FindHistogram(histogramArray,1,1,axisIndices,lowLimits,highLimits);
      fhTrackEta[iTrackType][iCentralityBin][fnTrackPtBins] = FindHistogram(histogramArray,2,1,axisIndices,lowLimits,highLimits);
      if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][fnTrackPtBins] = FindHistogram2D(histogramArray,1,2,1,axisIndices,lowLimits,highLimits);
      
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
        
        // Select the bin indices for track pT
        lowerTrackPtBin = fTrackPtBinIndices[iTrackPtBin];
        higherTrackPtBin = fTrackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
        
        // Add restriction for pT axis (0)
        axisIndices[1] = 0; lowLimits[1] = lowerTrackPtBin; highLimits[1] = higherTrackPtBin;
        
        // Read the angle histograms in track pT bins
        fhTrackPhi[iTrackType][iCentralityBin][iTrackPtBin] = FindHistogram(histogramArray,1,2,axisIndices,lowLimits,highLimits);
        fhTrackEta[iTrackType][iCentralityBin][iTrackPtBin] = FindHistogram(histogramArray,2,2,axisIndices,lowLimits,highLimits);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][iTrackPtBin] = FindHistogram2D(histogramArray,1,2,2,axisIndices,lowLimits,highLimits);
        
      } // Track pT loop
    } // Centrality loop
  } // Track category loop
}

/*
 * Loader for multiplicity histograms within the jet cone
 *
 * THnSparse for multiplicity within the jet cone:
 *
 *   Histogram name: multiplicityInJets
 *
 *     Axis index       Content of axis                       Note
 * ---------------------------------------------------------------------------------
 *       Axis 0           Multiplicity
 *       Axis 1              Jet pT
 *       Axis 2           Track pT cut
 *       Axis 3            Centrality
 *       Axis 4             Subevent          (0 = Pythia, 1 = Hydjet, 2 = Combined)
 */
void EECHistogramManager::LoadMultiplicityInJetConeHistograms(){
  
  if(!fLoadMultiplicityInJetHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[4] = {0};
  int lowLimits[4] = {0};
  int highLimits[4] = {0};
  
  // Define helper variables
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  THnSparseD* histogramArray;
  
  // Loop over multiplicity types
  for(int iMultiplicityType = 0; iMultiplicityType < knMultiplicityInJetConeTypes; iMultiplicityType++){

    // Open the multidimensional histogram from which the histograms are projected
    histogramArray = (THnSparseD*) fInputFile->Get(fMultiplicityInJetsHistogramNames[iMultiplicityType]);
    
    // Loop over centrality bins
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Select the centrality bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentrality];
      higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemover;
      
      // Setup axes with restrictions (3 = centrality)
      axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        
        // Reset the ranges for all the axes in the histogram array
        for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
          histogramArray->GetAxis(iAxis)->SetRange(0, 0);
        }

        // Select the track pT bin indices. Each bin has the total multiplicity above the lower threshold
        lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];
        higherTrackPtBin = fTrackPtIndicesEEC[iTrackPt+1]+duplicateRemover;
        
        // Add restriction for track pT axis (2 = track pT)
        axisIndices[1] = 2; lowLimits[1] = lowerTrackPtBin; highLimits[1] = higherTrackPtBin;
        
        // Loop over subevent types
        for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
          
          // Only do subevents for Pythia+Hydjet and Pythia+EPOS simulations
          if(!fSystemAndEnergy.Contains("Pb MC") && iSubevent < EECHistograms::knSubeventTypes) continue;
          
          // Add a restriction for the subevent axis (4 = subevent)
          axisIndices[2] = 4; lowLimits[2] = iSubevent+1; highLimits[2] = iSubevent+1;
          
          // Read the multiplicity histograms within the jet cone without jet pT restrictions
          fhMultiplicityInJetCone[iCentrality][fnJetPtBinsEEC][iTrackPt][iMultiplicityType][iSubevent] = FindHistogram(histogramArray, 0, 3, axisIndices, lowLimits, highLimits);
          
        } // Subevent loop
        
        // Loop over jet pT bins
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
          
          // Select the jet pT bin indices
          lowerJetPtBin = fJetPtIndicesEEC[iJetPt];
          higherJetPtBin = fJetPtIndicesEEC[iJetPt+1]+duplicateRemover;
          
          // Add restriction for jet pT axis (1 = jet pT)
          axisIndices[2] = 1; lowLimits[2] = lowerJetPtBin; highLimits[2] = higherJetPtBin;
          
          // Loop over subevent types
          for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
            
            // Only do subevents for Pythia+Hydjet and Pythia+EPOS simulations
            if(!fSystemAndEnergy.Contains("Pb MC") && iSubevent < EECHistograms::knSubeventTypes) continue;
            
            // Add a restriction for the subevent axis (4 = subevent)
            axisIndices[3] = 4; lowLimits[3] = iSubevent+1; highLimits[3] = iSubevent+1;
            
            // Read the multiplicity histograms within the jet cone without jet pT restrictions
            fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubevent] = FindHistogram(histogramArray, 0, 4, axisIndices, lowLimits, highLimits);
            
          } // Subevent loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Multiplicity type loop
}

/*
 * Loader for particle density histograms around the jet cone
 *
 * THnSparse for particle densities around the jet cone:
 *
 *   Histogram name: particleDensity/particlePtDensity
 *
 *     Axis index       Content of axis               Note
 * -----------------------------------------------------------------
 *       Axis 0              DeltaR             Constant bin size
 *       Axis 1              Jet pT
 *       Axis 2           Track pT cut
 *       Axis 3            Centrality
 *       Axis 4          Jet cone type       Signal cone/reflected cone
 *       Axis 5             Subevent           0 = Pythia, 1 = Hydjet
 *       Axis 6              DeltaR              Logarithmic bins
 */
void EECHistogramManager::LoadParticleDensityHistograms(){
    
  // Define arrays to help find the histograms
  int axisIndices[5] = {0};
  int lowLimits[5] = {0};
  int highLimits[5] = {0};
  
  // Define helper variables
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  THnSparseD* histogramArray;
  
  // Loop over particle density types
  for(int iParticleDensityType = 0; iParticleDensityType < knParticleDensityAroundJetAxisTypes; iParticleDensityType++){

    // Only load the selected particle density types
    if(!fLoadParticleDensityAroundJetsHistograms[iParticleDensityType]) continue;
    
    // For track pT bins, we are looking at all the tracks above the lower threshold
    histogramArray = (THnSparseD*) fInputFile->Get(fParticleDensityAroundJetsHistogramNames[iParticleDensityType]);
    higherTrackPtBin = histogramArray->GetAxis(2)->GetNbins()+1;

      // Loop over pairing types
    for(int iJetConeType = 0; iJetConeType < EECHistograms::knJetConeTypes; iJetConeType++){

        // Setup axes with restrictions, (4 = jet cone type)
      axisIndices[0] = 4; lowLimits[0] = iJetConeType+1; highLimits[0] = iJetConeType+1;

        // Loop over centrality bins
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

          // Select the centrality bin indices
        lowerCentralityBin = fCentralityBinIndices[iCentrality];
        higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemover;

          // Setup axes with restrictions, (3 = centrality)
        axisIndices[1] = 3; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin;

          // Loop over track pT bins
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

          // Reset the ranges for all the axes in the histogram array
          for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
            histogramArray->GetAxis(iAxis)->SetRange(0, 0);
          }

          // Select the track pT bin indices. Notice that we do not change the higher bin index
          lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];

          // We want to also get the histograms in finer pT bins
          if(iParticleDensityType == kParticleDensityAroundJetAxisPtBinned || iParticleDensityType == kParticlePtDensityAroundJetAxisPtBinned){
            higherTrackPtBin = fTrackPtIndicesEEC[iTrackPt+1]+duplicateRemover;
          }

          // Add restriction for track pT axis (2 = track pT)
          axisIndices[2] = 2; lowLimits[2] = lowerTrackPtBin; highLimits[2] = higherTrackPtBin;

          // Read the particle density histograms without jet pT restrictions
          fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes] = FindHistogram(histogramArray, 0, 3, axisIndices, lowLimits, highLimits, false);

          // After the histograms are read, normalize each bin to the bin area to make the contents particle density
          NormalizeToDeltaRBinArea(fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes]);

          // For PbPb and pPb MC, read the particle density histograms without jet pT restrictions in subevent bins
          if(fSystemAndEnergy.Contains("Pb MC")) {
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){

              // Add a restriction for the subevent axis (5 = subevent)
              axisIndices[3] = 5;
              lowLimits[3] = iSubevent+1;
              highLimits[3] = iSubevent+1;

              // Read the particle density histograms in subevent bins
              fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][iSubevent] = FindHistogram(histogramArray, 0, 4, axisIndices, lowLimits, highLimits, false);

              // After the histograms are read, normalize each bin to the bin area to make the contents particle density
              NormalizeToDeltaRBinArea(fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][iSubevent]);

            }  // Subevent loop

            // Reset the range of the subevent axis before proceeding
            histogramArray->GetAxis(5)->SetRange(0,0);

          }    // PbPb or pPb MC requirement

            // Loop over jet pT bins
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){

            // Select the jet pT bin indices
            lowerJetPtBin = fJetPtIndicesEEC[iJetPt];
            higherJetPtBin = fJetPtIndicesEEC[iJetPt+1]+duplicateRemover;
            
            // Add restriction for jet pT axis (1 = jet pT)
            axisIndices[3] = 1; lowLimits[3] = lowerJetPtBin; highLimits[3] = higherJetPtBin;
            
            // Read the particle density histograms
            fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes] = FindHistogram(histogramArray, 0, 4, axisIndices, lowLimits, highLimits, false);
            
            // After the histograms are read, normalize each bin to the bin area to make the contents particle density
            NormalizeToDeltaRBinArea(fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes]);
            
            // For PbPb and pPb MC, loop over subevent types
            if(fSystemAndEnergy.Contains("Pb MC")){
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){

                // Add a restriction for the subevent axis (5 = subevent)
                axisIndices[4] = 5; lowLimits[4] = iSubevent+1; highLimits[4] = iSubevent+1;
                
                // Read the particle density histograms
                fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubevent] = FindHistogram(histogramArray, 0, 5, axisIndices, lowLimits, highLimits, false);
                
                // After the histograms are read, normalize each bin to the bin area to make the contents particle density
                NormalizeToDeltaRBinArea(fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubevent]);
                
              } // Subevent loop

              // Reset the range of the subevent axis before proceeding
              histogramArray->GetAxis(5)->SetRange(0,0);

            } // PbPb or pPb MC requirement
            
          } // Jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // Jet cone type loop (signal cone/reflected cone)
  } // Particle density type loop (regular/pT weighted)
  
}

/*
 * Loader for maximum particle pT within the jet cone
 *
 * THnSparse for maximum particle pT within the jet cone:
 *
 *   Histogram name: maxParticlePtInJet
 *
 *     Axis index         Content of axis                      Note
 * -------------------------------------------------------------------------
 *       Axis 0               Jet pT
 *       Axis 1        Maximum signal particle pT
 *       Axis 2      Maximum background particle pT     Only relevant for MC
 *       Axis 3             Centrality

 */
void EECHistogramManager::LoadMaxParticlePtInJetConeHistograms(){
    
  if(!fLoadMaxParticlePtWithinJetConeHistograms) return;
  
  // Define arrays to help find the histograms
  int axisIndices[3] = {0};
  int lowLimits[3] = {0};
  int highLimits[3] = {0};
  
  // Define helper variables
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  int highestTrackPtBin = 0;
  THnSparseD* histogramArray;
  
  // Determine the highest track pT bin from the file
  histogramArray = (THnSparseD*) fInputFile->Get(fMaxParticlePtInJetConeHistogramName);
  highestTrackPtBin = histogramArray->GetAxis(1)->GetNbins()+1;
  
  // Loop over maximum particle pT within jet cone types
  for(int iMaxParticlePtType = 0; iMaxParticlePtType < knMaxParticlePtWithinJetConeTypes; iMaxParticlePtType++){
    
    // Loop over centrality bins
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Select the centrality bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentrality];
      higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemover;
      
      // Setup axes with restrictions, (3 = centrality)
      axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
      
      // Loop over jet pT bins
      for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
        
        // Reset the ranges for all the axes in the histogram array
        for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
          histogramArray->GetAxis(iAxis)->SetRange(0,0);
        }

        // Select the jet pT bin indices
        lowerJetPtBin = fJetPtIndicesEEC[iJetPt];
        higherJetPtBin = fJetPtIndicesEEC[iJetPt+1]+duplicateRemover;
        
        // Add restriction for jet pT axis (0 = jet pT)
        axisIndices[1] = 0; lowLimits[1] = lowerJetPtBin; highLimits[1] = higherJetPtBin;
        
        // Maximum particle pT within the jet cone histograms without any particle pT restrictions
        fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][knProjectedMaxParticlePtBins] = FindHistogram(histogramArray, 1+iMaxParticlePtType, 2, axisIndices, lowLimits, highLimits);
        
        // Loop over track pT bins
        for(int iTrackPt = 0; iTrackPt < knProjectedMaxParticlePtBins; iTrackPt++){
          
          // Find the bin indices for the bin boundaries
          lowerTrackPtBin = histogramArray->GetAxis(1)->FindBin(fProjectedMaxParticlePtBinBorders[iTrackPt]+0.01);
          higherTrackPtBin = histogramArray->GetAxis(1)->FindBin(fProjectedMaxParticlePtBinBorders[iTrackPt+1]-0.01);
          
          // Add restriction for particle pT axis (1 = signal particle pT, 2 = background particle pT)
          axisIndices[2] = 2-iMaxParticlePtType; lowLimits[2] = lowerTrackPtBin; highLimits[2] = higherTrackPtBin;
          fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt] = FindHistogram(histogramArray, 1+iMaxParticlePtType, 3, axisIndices, lowLimits, highLimits);
          
          // Add a cut for particle pT axis (1 = signal particle pT, 2 = background particle pT)
          axisIndices[2] = 2-iMaxParticlePtType; lowLimits[2] = lowerTrackPtBin; highLimits[2] = highestTrackPtBin;
          fhMaxParticlePtInJetConePtCut[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt] = FindHistogram(histogramArray, 1+iMaxParticlePtType, 3, axisIndices, lowLimits, highLimits);
          
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Maximum particle pT within jet cone type loop
}

/*
 * Loader for energy-energy correlator histograms
 *
 * THnSparse for energy-energy correlators:
 *
 *   Histogram name: energyEnergyCorrelator/energyEnergyCorrelatorEfficiencyVariationPlus/energyEnergyCorrelatorEfficiencyVariationMinus/
 *                   energyEnergyCorrelatorPairEfficiencyVariationPlus/energyEnergyCorrelatorPairEfficiencyVariationMinus
 *
 *     Axis index       Content of axis                   Note
 * ----------------------------------------------------------------------
 *       Axis 0              DeltaR
 *       Axis 1              Jet pT
 *       Axis 2           Track pT cut
 *       Axis 3            Centrality
 *       Axis 4           Pairing type           Same jet/reflected cone
 *       Axis 5       Subevent combination        Only relevant for MC
 *       Axis 6        Energy weight index      Only present in new files
 */
void EECHistogramManager::LoadEnergyEnergyCorrelatorHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[6] = {0};
  int lowLimits[6] = {0};
  int highLimits[6] = {0};
  
  // Define helper variables
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  int weightExponentBin = 0;
  int weightRestricted = 0;
  THnSparseD* histogramArray;
  
  // Loop over all different energy-energy correlator histograms
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;  // Only load the selected energy-energy correlators
    
    // For track pT bins, we are looking at all the tracks above the lower threshold
    histogramArray = (THnSparseD*) fInputFile->Get(fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
    higherTrackPtBin = histogramArray->GetAxis(2)->GetNbins()+1;

    // Check the number of dimensions in the histogram array. It will be 6 in older files and 7 in newer files
    // We need this information in newer files to project the desired energy weight from the THnSparse
    // We also only need to restrict the axis if there are more than one bin.
    if(histogramArray->GetNdimensions() >= 7 && fCard->GetNWeightExponents() > 1){

      // First, we need to find the index in histogram axis that the defined jet exponent corresponds to
      weightExponentBin = fCard->FindWeightExponentIndex(fLoadedWeightExponent);

      // After we have determined a bin index for the desired weight exponent, add it as a constraint to the energy weight axis
      axisIndices[0] = 6; lowLimits[0] = weightExponentBin; highLimits[0] = weightExponentBin;

      // If we are using weight exponent axis as restriction, all other restricton axis indices must be shifted by one
      weightRestricted = 1;
    } else {
      weightRestricted = 0;
    }
    
    // Loop over all pairing types that are included in the input file
    for(int iPairingType : *fIncludedPairingTypes){
      
      // Setup axes with restrictions, (4 = pairing type)
      axisIndices[weightRestricted] = 4; 
      lowLimits[weightRestricted] = iPairingType+1; 
      highLimits[weightRestricted] = iPairingType+1;
      
      // Loop over centrality bins
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Select the centrality bin indices
        lowerCentralityBin = fCentralityBinIndices[iCentrality];
        higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemover;
        
        // Setup axes with restrictions, (3 = centrality)
        axisIndices[weightRestricted+1] = 3; 
        lowLimits[weightRestricted+1] = lowerCentralityBin; 
        highLimits[weightRestricted+1] = higherCentralityBin;
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
          
          // Reset the ranges for all the axes in the histogram array
          for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
            histogramArray->GetAxis(iAxis)->SetRange(0,0);
          }

          // Select the track pT bin indices. Notice that we do not change the higher bin index
          lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];
          
          // Add restriction for track pT axis (2 = track pT)
          axisIndices[weightRestricted+2] = 2; 
          lowLimits[weightRestricted+2] = lowerTrackPtBin; 
          highLimits[weightRestricted+2] = higherTrackPtBin;
          
          // Read the energy-energy correlator histograms without jet pT restrictions
          fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations] = FindHistogram(histogramArray, 0, 3+weightRestricted, axisIndices, lowLimits, highLimits);
          
          // For PbPb and pPb MC, read the energy-energy correlator histograms without jet pT restrictions in subevent bins
          if(fSystemAndEnergy.Contains("Pb MC")){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
              
              // Add a restriction for the subevent axis (5 = subevent)
              
              // If we are not doing signal-reflected cone pairing, there is no meaningful differentiation between Pythia+Hydjet and Hydjet+Pythia
              if(iPairingType != EECHistograms::kSignalReflectedConePair && (iSubevent == EECHistograms::kPythiaHydjet || iSubevent == EECHistograms::kHydjetPythia)){
                axisIndices[weightRestricted+3] = 5; 
                lowLimits[weightRestricted+3] = EECHistograms::kPythiaHydjet+1; 
                highLimits[weightRestricted+3] = EECHistograms::kHydjetPythia+1;
              } else {
                axisIndices[weightRestricted+3] = 5; 
                lowLimits[weightRestricted+3] = iSubevent+1; 
                highLimits[weightRestricted+3] = iSubevent+1;
              }
              
              // Read the energy-energy correlator histograms
              fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent] = FindHistogram(histogramArray, 0, 4+weightRestricted, axisIndices, lowLimits, highLimits);
              
            } // Subevent loop

            // Reset the range of the subevent axis before proceeding
            histogramArray->GetAxis(5)->SetRange(0,0);

          } // PbPb or pPb MC requirement
          
          // Loop over jet pT bins
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
            
            // Select the jet pT bin indices
            lowerJetPtBin = fJetPtIndicesEEC[iJetPt];
            higherJetPtBin = fJetPtIndicesEEC[iJetPt+1]+duplicateRemover;
            
            // Add restriction for jet pT axis (1 = jet pT)
            axisIndices[weightRestricted+3] = 1; 
            lowLimits[weightRestricted+3] = lowerJetPtBin; 
            highLimits[weightRestricted+3] = higherJetPtBin;
            
            // Read the energy-energy correlator histograms
            fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations] = FindHistogram(histogramArray, 0, 4+weightRestricted, axisIndices, lowLimits, highLimits);
            
            // For PbPb and pPb MC, loop over subevent types
            if(fSystemAndEnergy.Contains("Pb MC")){
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
                
                // Add a restriction for the subevent axis (5 = subevent)
                
                // If we are not doing signal-reflected cone pairing, there is no meaningful differentiation between Pythia+Hydjet and Hydjet+Pythia
                if(iPairingType != EECHistograms::kSignalReflectedConePair && (iSubevent == EECHistograms::kPythiaHydjet || iSubevent == EECHistograms::kHydjetPythia)){
                  axisIndices[weightRestricted+4] = 5; 
                  lowLimits[weightRestricted+4] = EECHistograms::kPythiaHydjet+1; 
                  highLimits[weightRestricted+4] = EECHistograms::kHydjetPythia+1;
                } else {
                  axisIndices[weightRestricted+4] = 5; 
                  lowLimits[weightRestricted+4] = iSubevent+1; 
                  highLimits[weightRestricted+4] = iSubevent+1;
                }
                
                // Read the energy-energy correlator histograms
                fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = FindHistogram(histogramArray, 0, 5+weightRestricted, axisIndices, lowLimits, highLimits);

              } // Subevent loop

              // Reset the range of the subevent axis before proceeding
              histogramArray->GetAxis(5)->SetRange(0,0);

            } // PbPb or pPb MC requirement
            
          } // Jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // Pairing type loop (same jet/reflected cone)
  } // Energy-energy correlator type type loop
}

/*
 * Loader for reflected cone QA histograms
 *
 * THnSparse for jets above 25 GeV within the reflected cone:
 *
 *   Histogram name: jetNumberInReflectedCone
 *
 *     Axis index               Content of axis
 * -----------------------------------------------------------
 *       Axis 0        Number of jets within reflected cone
 *       Axis 1                   Centrality
 *
 * 
 *  THnSparse for jet pT within reflected cone:
 *
 *   Histogram name: jetPtInReflectedCone
 *
 *     Axis index               Content of axis
 * -----------------------------------------------------------
 *       Axis 0         pT of jets within reflected cone
 *       Axis 1                   Centrality
 */
void EECHistogramManager::LoadReflectedConeQAHistograms(){
  
  if(!fLoadReflectedConeQAHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  THnSparseD* histogramArray;

  // Load the number of jets within reflected cone histograms
  histogramArray = (THnSparseD*)fInputFile->Get(fNumberOfJetsWithinReflectedConeName);
  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

    // Select the centrality bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentrality];
    higherCentralityBin = fCentralityBinIndices[iCentrality+1] + duplicateRemoverCentrality;

    fhNumberOfJetsWithinReflectedCone[iCentrality] = FindHistogram(histogramArray,0,1,lowerCentralityBin,higherCentralityBin,0,0,0,false);

  } // Centrality loop

  // Load the jets pT within reflected cone histograms
  histogramArray = (THnSparseD*)fInputFile->Get(fJetPtWithinReflectedConeName);
  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

    // Select the centrality bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentrality];
    higherCentralityBin = fCentralityBinIndices[iCentrality+1] + duplicateRemoverCentrality;

    fhJetPtWithinReflectedCone[iCentrality] = FindHistogram(histogramArray,0,1,lowerCentralityBin,higherCentralityBin,0,0,0,false);

  } // Centrality loop
}

/*
 * Loader for jet pT response matrix
 *
 * THnSparse for closure histograms is used for this:
 *
 *   Histogram name: jetPtClosure
 *
 *     Axis index                  Content of axis
 * -----------------------------------------------------------
 *       Axis 0              Matched generator level jet pT
 *       Axis 1               Matched reconstructed jet pT
 *       Axis 2                         Jet eta
 *       Axis 3                       Centrality
 *       Axis 4                      Quark / gluon
 *       Axis 5             Matched reco to gen jet pT ratio
 */
void EECHistogramManager::LoadJetPtResponseMatrix(){
  
  if(!fLoadJetPtResponseMatrix) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[1] = {0};
  int lowLimits[1] = {0};
  int highLimits[1] = {0};
  int nRestrictionAxes = 1;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;

  // Find the histogram array from which the projections are made
  THnSparseD* histogramArray = (THnSparseD*)fInputFile->Get("jetPtClosure");
  AlgorithmLibrary* twoDeeRebinner = new AlgorithmLibrary();
  const int nBinsPt = 16;
  const double binBorderPt[nBinsPt+1] = {50,80,100,120,140,160,180,200,220,240,260,280,300,350,400,450,500};

  // Load all the histograms from the file
  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
    // Select the bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentrality];
    higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemoverCentrality;
        
    // Setup centrality axis restrictions
    axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin; // Centrality

    // Project the response matrix from the closure histogram
    fhJetPtResponseMatrix[iCentrality] = FindHistogram2D(histogramArray,1,0,nRestrictionAxes,axisIndices,lowLimits,highLimits,false);
        
    // Rebin the reconstructed and generator level jet pT axes
    fhJetPtResponseMatrix[iCentrality] = twoDeeRebinner->RebinHistogram(fhJetPtResponseMatrix[iCentrality], nBinsPt, binBorderPt, nBinsPt, binBorderPt, false, false);
        
  } // Centrality loop
}

/*
 * Loader for jet pT closure histograms
 *
 * THnSparse for closure histograms:
 *
 *   Histogram name: jetPtClosure
 *
 *     Axis index                  Content of axis
 * -----------------------------------------------------------
 *       Axis 0              Matched generator level jet pT
 *       Axis 1               Matched reconstructed jet pT
 *       Axis 2                         Jet eta
 *       Axis 3                       Centrality
 *       Axis 4                      Quark / gluon
 *       Axis 5             Matched reco to gen jet pT ratio
 *       Axis 6                         Jet phi
 */
void EECHistogramManager::LoadJetPtClosureHistograms(){
  
  if(!fLoadJetPtClosureHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[4] = {0};
  int lowLimits[4] = {0};
  int highLimits[4] = {0};
  int nRestrictionAxes = 3;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  THnSparseD* histogramArray = (THnSparseD*)fInputFile->Get("jetPtClosure");
  
  // Load all the histograms from the file
  for(int iGenJetPt = 0; iGenJetPt < knGenJetPtBins; iGenJetPt++){
    for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Reset the ranges for all the axes in the histogram array
        for(int iAxis = 0; iAxis < histogramArray->GetNdimensions(); iAxis++){
          histogramArray->GetAxis(iAxis)->SetRange(0,0);
        }

        // Select the bin indices
        lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
        higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
        
        // Setup the axes with restrictions
        nRestrictionAxes = 3;
        axisIndices[0] = 0; lowLimits[0] = iGenJetPt+1;    highLimits[0] = iGenJetPt+1;             // Gen jet pT
        axisIndices[1] = 3; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin; // Centrality
        axisIndices[2] = 4; lowLimits[2] = iClosureParticle+1; highLimits[2] = iClosureParticle+1;  // Quark/gluon
        
        // For the last closure particle bin no restrictions for quark/gluon jets
        if(iClosureParticle == EECHistograms::knClosureParticleTypes){
          
          // Remove the closure particle requirement from the restriction axes
          nRestrictionAxes--;
        }
        
        fhJetPtClosure[iGenJetPt][knJetEtaBins][knJetPhiBins][iCentralityBin][iClosureParticle] = FindHistogram(histogramArray,5,nRestrictionAxes,axisIndices,lowLimits,highLimits);

        // Reset the range for the generator level jet pT axis
        histogramArray->GetAxis(0)->SetRange(0,0);
        
        // Eta binning for the closure histogram
        for(int iJetEta = 0; iJetEta < knJetEtaBins; iJetEta++){
          
          // For the last closure particle bin no restrictions for quark/gluon jets
          /*if(iClosureParticle == EECHistograms::knClosureParticleTypes){
           nRestrictionAxes = 3;
           axisIndices[2] = 2; lowLimits[2] = iJetEta+1; highLimits[2] = iJetEta+1; // Jet eta
           } else {
           nRestrictionAxes = 4;
           axisIndices[3] = 2; lowLimits[3] = iJetEta+1; highLimits[3] = iJetEta+1; // Jet eta
           }
           
           fhJetPtClosure[iClosureType][iGenJetPt][iJetEta][knJetPhiBins][iCentralityBin][iClosureParticle] = FindHistogram(fInputFile,"jetPtClosure",5,nRestrictionAxes,axisIndices,lowLimits,highLimits);*/
          
          // Fill the pT integrated eta slices only once
          if(iGenJetPt == 0){
            
            // Setup the axes with restrictions
            nRestrictionAxes = 3;
            axisIndices[0] = 2; lowLimits[0] = iJetEta+1;    highLimits[0] = iJetEta+1;                 // Jet eta
            axisIndices[1] = 3; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin; // Centrality
            axisIndices[2] = 4; lowLimits[2] = iClosureParticle+1; highLimits[2] = iClosureParticle+1;  // Quark/gluon
            
            // For the last closure particle bin no restrictions for quark/gluon jets
            if(iClosureParticle == EECHistograms::knClosureParticleTypes){
              
              // Remove the last set array bin
              nRestrictionAxes--;
            }
            
            fhJetPtClosure[knGenJetPtBins][iJetEta][knJetPhiBins][iCentralityBin][iClosureParticle] = FindHistogram(histogramArray,5,nRestrictionAxes,axisIndices,lowLimits,highLimits);
          }
          
        } // Jet eta bin loop

        // Reset the range for the jet eta axis
        histogramArray->GetAxis(2)->SetRange(0,0);

        // Phi binning for the closure histogram
        for(int iJetPhi = 0; iJetPhi < knJetPhiBins; iJetPhi++){
          
          // Fill the pT integrated phi slices only once
          if(iGenJetPt == 0){
            
            // Setup the axes with restrictions
            nRestrictionAxes = 3;
            axisIndices[0] = 6; lowLimits[0] = iJetPhi+1;    highLimits[0] = iJetPhi+1;                 // Jet phi
            axisIndices[1] = 3; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin; // Centrality
            axisIndices[2] = 4; lowLimits[2] = iClosureParticle+1; highLimits[2] = iClosureParticle+1;  // Quark/gluon
            
            // For the last closure particle bin no restrictions for quark/gluon jets
            if(iClosureParticle == EECHistograms::knClosureParticleTypes){
              
              // Remove the last set array bin
              nRestrictionAxes--;
            }
            
            fhJetPtClosure[knGenJetPtBins][knJetEtaBins][iJetPhi][iCentralityBin][iClosureParticle] = FindHistogram(histogramArray,5,nRestrictionAxes,axisIndices,lowLimits,highLimits);
          }
          
        } // Jet eta bin loop
      } // Centrality loop
    } // Closure particle loop
  } // Gen jet pT loop
}

/*
 * Loader for jet pT unfolding histograms
 *
 * THnSparse for unfolding distributions:
 *
 *   Histogram name: jetPtUnfoldingMeasured/jetPtUnfoldingTruth
 *
 *     Axis index               Content of axis                       Note
 * -----------------------------------------------------------------------------------
 *       Axis 0              DeltaR in jet pT bins
 *       Axis 1                    Track pT
 *       Axis 2                   Centrality
 *       Axis 3               Energy weight index           Only present in new files
 *
 * 
 *  THnSparse for unfolding response matrix:
 *
 *   Histogram name: jetPtUnfoldingResponse
 *
 *     Axis index               Content of axis                       Note
 * ------------------------------------------------------------------------------------
 *       Axis 0         DeltaR in reconstructed jet pT bins
 *       Axis 1        DeltaR in generator level jet pT bins
 *       Axis 2                    Track pT
 *       Axis 3                   Centrality
 *       Axis 4               Energy weight index           Only present in new files
 */
void EECHistogramManager::LoadJetPtUnfoldingHistograms(){
  
  if(!fLoadJetPtUnfoldingHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[3] = {0};
  int lowLimits[3] = {0};
  int highLimits[3] = {0};
  int nRestrictionAxes = 2;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  int weightExponentBin = 0;
  int weightRestricted = 0;
  THnSparseD* histogramArray;

  // Load the jet pT unfolding distribution
  for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
    histogramArray = (THnSparseD*)fInputFile->Get(fJetPtUnfoldingDistributionName[iUnfoldType]);
    higherTrackPtBin = histogramArray->GetAxis(1)->GetNbins()+1;

    // Check the number of dimensions in the histogram array. It will be 3 in older files and 4 in newer files
    // We need this information in newer files to project the desired energy weight from the THnSparse
    // We also only need to restrict the axis if there are more than one bin.
    if(histogramArray->GetNdimensions() >= 4 && fCard->GetNWeightExponents() > 1){

      // First, we need to find the index in histogram axis that the defined jet exponent corresponds to
      weightExponentBin = fCard->FindWeightExponentIndex(fLoadedWeightExponent);

      // After we have determined a bin index for the desired weight exponent, add it as a constraint to the energy weight axis
      axisIndices[0] = 3; lowLimits[0] = weightExponentBin; highLimits[0] = weightExponentBin;

      // If we are using weight exponent axis as restriction, all other restricton axis indices must be shifted by one
      weightRestricted = 1;
    } else {
      weightRestricted = 0;
    }

    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Select the centrality bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentrality];
      higherCentralityBin = fCentralityBinIndices[iCentrality+1] + duplicateRemoverCentrality;

      // Add restriction for centrality axis (2 = centrality)
      axisIndices[weightRestricted] = 2; 
      lowLimits[weightRestricted] = lowerCentralityBin; 
      highLimits[weightRestricted] = higherCentralityBin;

      // Loop over track pT bins
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

        // Select the track pT bin indices. Notice that we do not change the higher bin index
        lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];

        // Add restriction for track pT axis (1 = track pT)
        axisIndices[weightRestricted+1] = 1; 
        lowLimits[weightRestricted+1] = lowerTrackPtBin; 
        highLimits[weightRestricted+1] = higherTrackPtBin;

        fhJetPtUnfoldingDistribution[iUnfoldType][iCentrality][iTrackPt] = FindHistogram(histogramArray, 0, nRestrictionAxes+weightRestricted, axisIndices, lowLimits, highLimits, false);
      }

    } // Centrality loop
  } // Unfold distribution type loop

  // Load the jet pT unfolding response
  histogramArray = (THnSparseD*)fInputFile->Get(fJetPtResponseMatrixName);
  higherTrackPtBin = histogramArray->GetAxis(2)->GetNbins()+1;

  // Check the number of dimensions in the histogram array. It will be 4 in older files and 5 in newer files
  // We need this information in newer files to project the desired energy weight from the THnSparse
  // We also only need to restrict the axis if there are more than one bin.
  if(histogramArray->GetNdimensions() >= 5 && fCard->GetNWeightExponents() > 1){

    // First, we need to find the index in histogram axis that the defined jet exponent corresponds to
    weightExponentBin = fCard->FindWeightExponentIndex(fLoadedWeightExponent);

    // After we have determined a bin index for the desired weight exponent, add it as a constraint to the energy weight axis
    axisIndices[0] = 4; lowLimits[0] = weightExponentBin; highLimits[0] = weightExponentBin;

    // If we are using weight exponent axis as restriction, all other restricton axis indices must be shifted by one
    weightRestricted = 1;
  } else {
    weightRestricted = 0;
  }

  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

    // Select the centrality bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentrality];
    higherCentralityBin = fCentralityBinIndices[iCentrality + 1] + duplicateRemoverCentrality;

    // Add restriction for centrality axis (3 = centrality)
    axisIndices[weightRestricted] = 3; 
    lowLimits[weightRestricted] = lowerCentralityBin; 
    highLimits[weightRestricted] = higherCentralityBin;

    // Loop over track pT bins
    for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

      // Select the track pT bin indices. Notice that we do not change the higher bin index
      lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];

      // Add restriction for track pT axis (2 = track pT)
      axisIndices[weightRestricted+1] = 2; 
      lowLimits[weightRestricted+1] = lowerTrackPtBin; 
      highLimits[weightRestricted+1] = higherTrackPtBin;

      fhJetPtUnfoldingResponse[iCentrality][iTrackPt] = FindHistogram2D(histogramArray, 0, 1, nRestrictionAxes+weightRestricted, axisIndices, lowLimits, highLimits, false);
    }

  }  // Centrality loop
}

/*
 * Loader for covariance histograms used in jet pT unfolding algorithm
 * 
 *  THnSparse for unfolding response matrix:
 *
 *   Histogram name: jetPtUnfoldingCovariance
 *
 *     Axis index               Content of axis                       Note
 * -----------------------------------------------------------------------------------
 *       Axis 0              DeltaR in jet pT bins
 *       Axis 1              DeltaR in jet pT bins
 *       Axis 2                    Track pT
 *       Axis 3                   Centrality
 *       Axis 4               Energy weight index           Only present in new files
 */
void EECHistogramManager::LoadJetPtUnfoldingCovariance(){
  
  if(!fLoadJetPtUnfoldingCovariance) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[3] = {0};
  int lowLimits[3] = {0};
  int highLimits[3] = {0};
  int nRestrictionAxes = 2;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  int weightExponentBin = 0;
  int weightRestricted = 0;
  THnSparseD* histogramArray;

  // Load the jet pT unfolding distribution
  histogramArray = (THnSparseD*)fInputFile->Get(fJetPtCovarianceMatrixName[kCovarianceMatrixMeasured]);
  higherTrackPtBin = histogramArray->GetAxis(2)->GetNbins()+1;

  // Check the number of dimensions in the histogram array. It will be 4 in older files and 5 in newer files
  // We need this information in newer files to project the desired energy weight from the THnSparse
  // We also only need to restrict the axis if there are more than one bin.
  if(histogramArray->GetNdimensions() >= 5 && fCard->GetNWeightExponents() > 1){

    // First, we need to find the index in histogram axis that the defined jet exponent corresponds to
    weightExponentBin = fCard->FindWeightExponentIndex(fLoadedWeightExponent);

    // After we have determined a bin index for the desired weight exponent, add it as a constraint to the energy weight axis
    axisIndices[0] = 4; lowLimits[0] = weightExponentBin; highLimits[0] = weightExponentBin;

    // If we are using weight exponent axis as restriction, all other restricton axis indices must be shifted by one
    weightRestricted = 1;
  } else {
    weightRestricted = 0;
  }

  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

    // Select the centrality bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentrality];
    higherCentralityBin = fCentralityBinIndices[iCentrality+1] + duplicateRemoverCentrality;

    // Add restriction for centrality axis (3 = centrality)
    axisIndices[weightRestricted] = 3; 
    lowLimits[weightRestricted] = lowerCentralityBin; 
    highLimits[weightRestricted] = higherCentralityBin;

    // Loop over track pT bins
    for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

      // Select the track pT bin indices. Notice that we do not change the higher bin index
      lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];

      // Add restriction for track pT axis (2 = track pT)
      axisIndices[weightRestricted+1] = 2; 
      lowLimits[weightRestricted+1] = lowerTrackPtBin; 
      highLimits[weightRestricted+1] = higherTrackPtBin;

      fhJetPtUnfoldingCovariance[kCovarianceMatrixMeasured][iCentrality][iTrackPt] = FindHistogram2D(histogramArray, 0, 1, nRestrictionAxes+weightRestricted, axisIndices, lowLimits, highLimits, false);

    } // Track pT loop

  } // Centrality loop
}

/*
 * Loader for one-dimensional jet pT unfolding histograms
 *
 * THnSparse for one-dimensional unfolding distributions:
 *
 *   Histogram name: oneDimensionalJetPtUnfoldingMeasured/oneDimensionalJetPtUnfoldingTruth
 *
 *     Axis index               Content of axis
 * -----------------------------------------------------------
 *       Axis 0              Measured/true jet pT
 *       Axis 1                   Centrality
 *
 * 
 *  THnSparse for one-dimensional unfolding response matrix:
 *
 *   Histogram name: oneDimensionalJetPtUnfoldingResponse
 *
 *     Axis index               Content of axis
 * -----------------------------------------------------------
 *       Axis 0                Measured jet pT
 *       Axis 1                  True jet pT
 *       Axis 2                   Centrality
 */
void EECHistogramManager::LoadJetPtOneDimensionalUnfoldingHistograms(){
  
  if(!fLoadJetPtOneDimensionalUnfoldingHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  THnSparseD* histogramArray;

  // Load the jet pT distributions for one-dimensional unfolding
  for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
    histogramArray = (THnSparseD*)fInputFile->Get(fJetPtOneDimensionalUnfoldingDistributionName[iUnfoldType]);

    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Select the centrality bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentrality];
      higherCentralityBin = fCentralityBinIndices[iCentrality+1] + duplicateRemoverCentrality;

      fhOneDimensionalJetPtUnfoldingDistribution[iUnfoldType][iCentrality] = FindHistogram(histogramArray,0,1,lowerCentralityBin,higherCentralityBin, 0, 0, 0, false);

    } // Centrality loop
  } // Unfold distribution type loop

  // Load the jet pT unfolding response
  histogramArray = (THnSparseD*)fInputFile->Get(fJetPtOneDimensionalResponseMatrixName);

  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

    // Select the centrality bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentrality];
    higherCentralityBin = fCentralityBinIndices[iCentrality+1] + duplicateRemoverCentrality;

    fhOneDimensionalJetPtUnfoldingResponse[iCentrality] = FindHistogram2D(histogramArray, 0, 1, 2, lowerCentralityBin, higherCentralityBin, 0, 0, 0, false);

  }  // Centrality loop
}

/*
 * Loader for track/particle matching histograms
 *
 * THnSparse for particle number close to tracks and has matching particle flag histograms:
 *
 *   Histogram name: particlesCloseToTracks/tracksWithMatchedParticle
 *
 *     Axis index                            Content of axis
 * ----------------------------------------------------------------------------------
 *       Axis 0        Number of particles close to a track / Matching particle flag
 *       Axis 1                                  Jet pT
 *       Axis 2                                 Track pT
 *       Axis 3                                Centrality
 *
 * 
 *  THnSparse for track/particle pair DeltaR response matrix:
 *
 *   Histogram name: particleDeltaRResponseMatrix
 *
 *     Axis index                Content of axis
 * ---------------------------------------------------------------
 *       Axis 0          DeltaR for reconstructed track pairs
 *       Axis 1        DeltaR for generator level particle pairs
 *       Axis 2                       Jet pT
 *       Axis 3                      Track pT
 *       Axis 4                     Centrality
 *
 *
 *  THnSparse for track/particle pair pT response matrix:
 *
 *   Histogram name: particlePtResponseMatrix
 *
 *     Axis index                  Content of axis                               Note
 * ----------------------------------------------------------------------------------------------
 *       Axis 0          pT1*pT2 for reconstructed track pairs
 *       Axis 1        pT1*pT2 for generator level particle pairs
 *       Axis 2                        Jet pT
 *       Axis 3                       Track pT
 *       Axis 4                      Centrality
 *       Axis 5     reconstructed pT1*pT2 / generator level pT1*pT2
 *       Axis 6                 Energy weight index                    Only present in new files
 */
void EECHistogramManager::LoadTrackParticleMatchingHistograms(){

  if(!fLoadTrackParticleMatchingHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[3] = {0};
  int lowLimits[3] = {0};
  int highLimits[3] = {0};
  int nRestrictionAxes = 3;
  
  // Define helper variables
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;
  int weightExponentBin = 0;
  int weightRestricted = 0;
  THnSparseD* histogramArray;

  // First, load the track/particle matching QA histograms
  for(int iTrackParticleQAHistogram = 0; iTrackParticleQAHistogram < knTrackParticleMatchingQAHistograms; iTrackParticleQAHistogram++){

    histogramArray = (THnSparseD*) fInputFile->Get(fTrackParticleMatchingQAName[iTrackParticleQAHistogram]);
    higherTrackPtBin = histogramArray->GetAxis(2)->GetNbins() + 1;

    // Loop over centrality bins
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Select the centrality bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentrality];
      higherCentralityBin = fCentralityBinIndices[iCentrality + 1] + duplicateRemover;

      // Setup axes with restrictions, (3 = centrality)
      axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;

      // Loop over track pT bins
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

        // Select the track pT bin indices. Notice that we do not change the higher bin index
        lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];

        // Add restriction for track pT axis (2 = track pT)
        axisIndices[1] = 2; lowLimits[1] = lowerTrackPtBin; highLimits[1] = higherTrackPtBin;

        // Loop over jet pT bins
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){

          // Select the jet pT bin indices
          lowerJetPtBin = fJetPtIndicesEEC[iJetPt];
          higherJetPtBin = fJetPtIndicesEEC[iJetPt+1] + duplicateRemover;

          // Add restriction for jet pT axis (1 = jet pT)
          axisIndices[2] = 1; lowLimits[2] = lowerJetPtBin; highLimits[2] = higherJetPtBin;

          // Read the track/particle matching QA histograms
          fhTrackParticleMatchQA[iTrackParticleQAHistogram][iCentrality][iJetPt][iTrackPt] = FindHistogram(histogramArray, 0, nRestrictionAxes, axisIndices, lowLimits, highLimits, false);

        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Track/particle QA histogram loop

  // Then, load the track/particle pair DeltaR and pT1*pT2 response matrices
  for(int iTrackParticleResponseHistogram = 0; iTrackParticleResponseHistogram < knTrackParticleMatchingResponseTypes; iTrackParticleResponseHistogram++){

    histogramArray = (THnSparseD*) fInputFile->Get(fTrackParticleMatchingResponseName[iTrackParticleResponseHistogram]);
    higherTrackPtBin = histogramArray->GetAxis(3)->GetNbins()+1;

    // Loop over centrality bins
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Select the centrality bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentrality];
      higherCentralityBin = fCentralityBinIndices[iCentrality + 1] + duplicateRemover;

      // Setup axes with restrictions, (4 = centrality)
      axisIndices[0] = 4; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;

      // Loop over track pT bins
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

        // Select the track pT bin indices. Notice that we do not change the higher bin index
        lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];

        // Add restriction for track pT axis (3 = track pT)
        axisIndices[1] = 3; lowLimits[1] = lowerTrackPtBin; highLimits[1] = higherTrackPtBin;

        // Loop over jet pT bins
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){

          // Select the jet pT bin indices
          lowerJetPtBin = fJetPtIndicesEEC[iJetPt];
          higherJetPtBin = fJetPtIndicesEEC[iJetPt+1] + duplicateRemover;

          // Add restriction for jet pT axis (2 = jet pT)
          axisIndices[2] = 2; lowLimits[2] = lowerJetPtBin; highLimits[2] = higherJetPtBin;

          // Read the track/particle response matrices
          fhTrackParticleResponse[iTrackParticleResponseHistogram][iCentrality][iJetPt][iTrackPt] = FindHistogram2D(histogramArray, 0, 1, nRestrictionAxes, axisIndices, lowLimits, highLimits, false);

        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Track/particle response matrix type loop

  // Finally, load the track/particle pair pT1*pT2 closure histograms
  histogramArray = (THnSparseD*)fInputFile->Get("particlePtResponseMatrix");
  higherTrackPtBin = histogramArray->GetAxis(3)->GetNbins()+1;

  // Check the number of dimensions in the histogram array. It will be 6 in older files and 7 in newer files
  // We need this information in newer files to project the desired energy weight from the THnSparse
  // We also only need to restrict the axis if there are more than one bin.
  if(histogramArray->GetNdimensions() >= 7 && fCard->GetNWeightExponents() > 1){

    // First, we need to find the index in histogram axis that the defined jet exponent corresponds to
    weightExponentBin = fCard->FindWeightExponentIndex(fLoadedWeightExponent);

    // After we have determined a bin index for the desired weight exponent, add it as a constraint to the energy weight axis
    axisIndices[0] = 6; lowLimits[0] = weightExponentBin; highLimits[0] = weightExponentBin;

    // If we are using weight exponent axis as restriction, all other restricton axis indices must be shifted by one
    weightRestricted = 1;
  } else {
    weightRestricted = 0;
  }

  // Loop over centrality bins
  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

    // Select the centrality bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentrality];
    higherCentralityBin = fCentralityBinIndices[iCentrality+1] + duplicateRemover;

    // Setup axes with restrictions, (4 = centrality)
    axisIndices[weightRestricted] = 4; 
    lowLimits[weightRestricted] = lowerCentralityBin; 
    highLimits[weightRestricted] = higherCentralityBin;

    // Loop over track pT bins
    for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

      // Select the track pT bin indices. Notice that we do not change the higher bin index
      lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];

      // Add restriction for track pT axis (3 = track pT)
      axisIndices[weightRestricted+1] = 3; 
      lowLimits[weightRestricted+1] = lowerTrackPtBin; 
      highLimits[weightRestricted+1] = higherTrackPtBin;

      // Loop over jet pT bins
      for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){

        // Select the jet pT bin indices
        lowerJetPtBin = fJetPtIndicesEEC[iJetPt];
        higherJetPtBin = fJetPtIndicesEEC[iJetPt+1] + duplicateRemover;

        // Add restriction for jet pT axis (2 = jet pT)
        axisIndices[weightRestricted+2] = 2; 
        lowLimits[weightRestricted+2] = lowerJetPtBin; 
        highLimits[weightRestricted+2] = higherJetPtBin;

        // Read the track/particle pair pT1*pT2 closure histograms
        fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt] = FindHistogram(histogramArray, 5, nRestrictionAxes+weightRestricted, axisIndices, lowLimits, highLimits, false);

      }  // Jet pT loop
    }    // Track pT loop
  }      // Centrality loop
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Inputfile containing the THnSparse to be read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* EECHistogramManager::FindHistogram2D(THnSparseD* histogramArray, int xAxis, int yAxis, const bool normalizeToBinWidth){
  
  // Create a unique name for eeach histogram that is read from the file
  TString newName = histogramArray->GetName();
  newName.Append("_inclusive");
  
  // Project out the histogram and give it the created unique name
  TH2D* projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName.Data());
  
  // Apply bin width normalization to the projected histogram
  if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Inputfile containing the THnSparse to be read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int nAxes = Number of axes that are restained for the projection
 *   int* axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int* lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int* highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* EECHistogramManager::FindHistogram2D(THnSparseD* histogramArray, int xAxis, int yAxis, int nAxes, int* axisNumber, int* lowBinIndex, int* highBinIndex, const bool normalizeToBinWidth){
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH2D* projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName.Data());
  
  // Apply bin width normalization to the projected histogram
  if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Inputfile containing the THnSparse to be read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH2D* EECHistogramManager::FindHistogram2D(THnSparseD* histogramArray, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2, const bool normalizeToBinWidth){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram2D(histogramArray,xAxis,yAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices,normalizeToBinWidth);
}

/*
 * Extract a histogram without restrictions on other axes in THnSparse
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Histogram array from which the desired histograms are projected
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* EECHistogramManager::FindHistogram(THnSparseD* histogramArray, int xAxis, const bool normalizeToBinWidth){
  
  // Create a unique name for each histogram that is read from the file
  TString newName = histogramArray->GetName();
  newName.Append("_inclusive");
  
  // Project out the histogram and give it the created unique name
  TH1D* projectedHistogram = NULL;
  
  // Check that we are not trying to project a non-existing axis
  if(xAxis < histogramArray->GetNdimensions()){
    projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName.Data());
  
    // Apply bin width normalization to the projected histogram
    if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  }
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a histogram with given restrictions on other axes in THnSparse
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Histogram array from which the desired histograms are projected
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int* axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int* lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int* highBinIndex = Indices of the highest considered bins in the restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* EECHistogramManager::FindHistogram(THnSparseD* histogramArray, int xAxis, int nAxes, int* axisNumber, int* lowBinIndex, int* highBinIndex, const bool normalizeToBinWidth){
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for each histogram that is read from the file
  TString newName = histogramArray->GetName();
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    newName.Append(Form("_%d=%d-%d",axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]));
  }
  
  // Project out the histogram and give it the created unique name
  TH1D* projectedHistogram = NULL;
  
  // Check that we are not trying to project a non-existing axis
  if(xAxis < histogramArray->GetNdimensions()){
    projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName.Data());
  
    // Apply bin width normalization to the projected histogram
    if(normalizeToBinWidth) projectedHistogram->Scale(1.0,"width");
  }
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   THnSparseD* histogramArray = Histogram array from which the desired histograms are projected
 *   int xAxis = Index for the axis in THnSparse that is projected to TH1D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 *   const bool normalizeToBinWidth = Flag for normalizing the projected histogram to the bin width
 */
TH1D* EECHistogramManager::FindHistogram(THnSparseD* histogramArray, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2, const bool normalizeToBinWidth){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram(histogramArray,xAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices,normalizeToBinWidth);
}

/*
 * Write all the loaded histograms into a file
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void EECHistogramManager::Write(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile* outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  TString histogramNamer;
  
  // Write the event information histograms to the output file
  if(fLoadEventInformation){
    fhEvents->Write("",TObject::kOverwrite);             // Number of events surviving different event cuts
    fhTriggers->Write("",TObject::kOverwrite);           // Trigger selection information
    fhVertexZ->Write("",TObject::kOverwrite);            // Vertex z position
    fhVertexZWeighted->Write("",TObject::kOverwrite);    // MC weighted vertex z position
    fhTrackCuts->Write("",TObject::kOverwrite);          // Number of tracks surviving different track cuts
    fhCentrality->Write("",TObject::kOverwrite);         // Centrality in all events
    fhCentralityWeighted->Write("",TObject::kOverwrite); // MC weighted centrality in all events
    fhPtHat->Write("",TObject::kOverwrite);              // pT hat for MC events (only meaningful for MC)
    fhPtHatWeighted->Write("",TObject::kOverwrite);      // Weighted pT hat distribution (only meaningful for MC)

    // Create a directory for the HF energy histograms if it does not already exist
    if(!gDirectory->GetDirectory("HFEnergy")) gDirectory->mkdir("HFEnergy");
    gDirectory->cd("HFEnergy");

    // Write the histograms without jet pT selection
    fhHFPlus[fnJetPtBinsEEC][knMaxJetsBins]->Write("HFPlus",TObject::kOverwrite);                   // Energy in HF plus calorimeters
    fhHFMinus[fnJetPtBinsEEC][knMaxJetsBins]->Write("HFMinus",TObject::kOverwrite);                 // Energy in HF minus calorimeters
    fhHFSum[fnJetPtBinsEEC][knMaxJetsBins]->Write("HFSum",TObject::kOverwrite);                     // Energy in HF calorimeters
    fhHFPlusVsHFMinus[fnJetPtBinsEEC][knMaxJetsBins]->Write("HFPlusVsHFMinus",TObject::kOverwrite); // 2D-map between HF plus and HF minus

    // Loop over jet pT
    for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
      if(fhHFPlus[iJetPt][knMaxJetsBins]) fhHFPlus[iJetPt][knMaxJetsBins]->Write(Form("HFPlus_J%d", iJetPt), TObject::kOverwrite);                       // Energy in HF plus calorimeters
      if(fhHFMinus[iJetPt][knMaxJetsBins]) fhHFMinus[iJetPt][knMaxJetsBins]->Write(Form("HFMinus_J%d", iJetPt), TObject::kOverwrite);                       // Energy in HF minus calorimeters
      if(fhHFSum[iJetPt][knMaxJetsBins]) fhHFSum[iJetPt][knMaxJetsBins]->Write(Form("HFSum_J%d", iJetPt), TObject::kOverwrite);                       // Energy in HF calorimeters
      if(fhHFPlusVsHFMinus[iJetPt][knMaxJetsBins]) fhHFPlusVsHFMinus[iJetPt][knMaxJetsBins]->Write(Form("HFPlusVsHFMinus_J%d", iJetPt), TObject::kOverwrite); // 2D-map between HF plus and HF minus

      // Loop over number of jets
      for(int iNJet = 0; iNJet < knMaxJetsBins; iNJet++){
        if(fhHFPlus[iJetPt][iNJet]) fhHFPlus[iJetPt][iNJet]->Write(Form("HFPlus_J%d_N%d", iJetPt, iNJet), TObject::kOverwrite);                       // Energy in HF plus calorimeters
        if(fhHFMinus[iJetPt][iNJet]) fhHFMinus[iJetPt][iNJet]->Write(Form("HFMinus_J%d_N%d", iJetPt, iNJet), TObject::kOverwrite);                       // Energy in HF minus calorimeters
        if(fhHFSum[iJetPt][iNJet]) fhHFSum[iJetPt][iNJet]->Write(Form("HFSum_J%d_N%d", iJetPt, iNJet), TObject::kOverwrite);                       // Energy in HF calorimeters
        if(fhHFPlusVsHFMinus[iJetPt][iNJet]) fhHFPlusVsHFMinus[iJetPt][iNJet]->Write(Form("HFPlusVsHFMinus_J%d_N%d", iJetPt, iNJet), TObject::kOverwrite);  // 2D-map between HF plus and HF minus
      }
    }

    // Loop over number of jets
    for(int iNJet = 0; iNJet < knMaxJetsBins; iNJet++){
      if(fhHFPlus[fnJetPtBinsEEC][iNJet]) fhHFPlus[fnJetPtBinsEEC][iNJet]->Write(Form("HFPlus_N%d", iNJet), TObject::kOverwrite);                       // Energy in HF plus calorimeters
      if(fhHFMinus[fnJetPtBinsEEC][iNJet]) fhHFMinus[fnJetPtBinsEEC][iNJet]->Write(Form("HFMinus_N%d", iNJet), TObject::kOverwrite);                       // Energy in HF minus calorimeters
      if(fhHFSum[fnJetPtBinsEEC][iNJet]) fhHFSum[fnJetPtBinsEEC][iNJet]->Write(Form("HFSum_N%d", iNJet), TObject::kOverwrite);                       // Energy in HF calorimeters
      if(fhHFPlusVsHFMinus[fnJetPtBinsEEC][iNJet]) fhHFPlusVsHFMinus[fnJetPtBinsEEC][iNJet]->Write(Form("HFPlusVsHFMinus_N%d", iNJet), TObject::kOverwrite);  // 2D-map between HF plus and HF minus
    }

    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the multiplicity histograms if it does not already exist
    if(!gDirectory->GetDirectory("multiplicity")) gDirectory->mkdir("multiplicity");
    gDirectory->cd("multiplicity");
    
    // Loop over centrality
    for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
      histogramNamer = Form("multiplicity_C%d", iCentrality);
      if(fhMultiplicity[iCentrality]) fhMultiplicity[iCentrality]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      histogramNamer = Form("multiplicityWeighted_C%d", iCentrality);
      if(fhMultiplicityWeighted[iCentrality]) fhMultiplicityWeighted[iCentrality]->Write(histogramNamer.Data(), TObject::kOverwrite);
    }
    
    // Histograms without cerntality binning
    if(fhMultiplicity[fnCentralityBins]) fhMultiplicity[fnCentralityBins]->Write("multiplicity", TObject::kOverwrite);
    
    if(fhMultiplicityWeighted[fnCentralityBins]) fhMultiplicityWeighted[fnCentralityBins]->Write("multiplicityWeighted", TObject::kOverwrite);
    
    // Multiplicity vs. centrality maps
    if(fhMultiplicityMap) fhMultiplicityMap->Write("multiplicityMap", TObject::kOverwrite);
    
    if(fhMultiplicityMapWeighted) fhMultiplicityMapWeighted->Write("multiplicityMapWeighted", TObject::kOverwrite);
    
    // Return back to main directory
    gDirectory->cd("../");
  }
 
  // Write the jet histograms to the output file
  WriteJetHistograms();
  
  // Write the track histograms to the output file
  WriteTrackHistograms();
  
  // Write the multiplicity histograms within the jet cone to the output file
  WriteMultiplicityInJetConeHistograms();
  
  // Write the particle density histograms around the jet axis to the output file
  WriteParticleDensityAroundJetsHistograms();
  
  // Write the maximum particle pT within the jet cone histograms to the output file
  WriteMaxParticlePtWithinJetConeHistograms();
  
  // Write the energy-energy correlator histograms to the output file
  WriteEnergyEnergyCorrelatorHistograms();

  // Write the reflected cone QA histogram to the output file
  WriteReflectedConeQAHistograms();
  
  // Write the jet pT response matrices
  WriteJetPtResponseMatrix(); 

  // Write the jet pT closure histograms to the output file
  WriteClosureHistograms();

  // Write the jet pT unfolding histograms to the output file
  WriteJetPtUnfoldingHistograms();

  // Write the covariance histograms used in jet pT unfolding
  WriteJetPtUnfoldingCovariance();

  // Write the jet pT one-dimensional unfolding histograms to the output file
  WriteJetPtOneDimensionalUnfoldingHistograms();

  // Write the track/particle matching histograms to the output file
  WriteTrackParticleMatchingHistograms();
  
  // Write the card to the output file if it is not already written
  if(!gDirectory->GetDirectory("JCard")) fCard->Write(outputFile);
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Write the jet histograms to the file that is currently open
 */
void EECHistogramManager::WriteJetHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Write the jet histograms to the output file
  if(!fLoadJets) return;  // Only write the jet histograms if they are loaded
  
  // Create a directory for the histograms if it does not already exist
  if(!gDirectory->GetDirectory(fJetHistogramName)) gDirectory->mkdir(fJetHistogramName);
  gDirectory->cd(fJetHistogramName);
  
  for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
    
    // Check that the histograms are actually there before trying to save them.
    if(fhJetPt[iCentralityBin] == NULL) {
      cout << "Could not find histograms of type " << fJetHistogramName << " to write. Will skip writing these." << endl;
      continue;
    }
    
    // Jet pT
    histogramNamer = Form("%sPt_C%d",fJetHistogramName,iCentralityBin);
    if(fhJetPt[iCentralityBin]) fhJetPt[iCentralityBin]->Write(histogramNamer.Data(), TObject::kOverwrite);
    
    // Jet phi
    histogramNamer = Form("%sPhi_C%d",fJetHistogramName,iCentralityBin);
    if(fhJetPhi[iCentralityBin]) fhJetPhi[iCentralityBin]->Write(histogramNamer.Data(), TObject::kOverwrite);
    
    // Jet eta
    histogramNamer = Form("%sEta_C%d",fJetHistogramName,iCentralityBin);
    if(fhJetEta[iCentralityBin]) fhJetEta[iCentralityBin]->Write(histogramNamer.Data(), TObject::kOverwrite);

    // Jet eta_CM
    histogramNamer = Form("%sEtaCM_C%d",fJetHistogramName,iCentralityBin);
    if(fhJetEtaCM[iCentralityBin]) fhJetEtaCM[iCentralityBin]->Write(histogramNamer.Data(), TObject::kOverwrite);
    
    // Jet eta-phi
    histogramNamer = Form("%sEtaPhi_C%d",fJetHistogramName,iCentralityBin);
    if(fLoad2DHistograms && fhJetEtaPhi[iCentralityBin]) fhJetEtaPhi[iCentralityBin]->Write(histogramNamer.Data(), TObject::kOverwrite);
    
  } // Loop over centrality bins
  
  // Return back to main directory
  gDirectory->cd("../");
  
}

/*
 * Write the track histograms to the file that is currently open
 */
void EECHistogramManager::WriteTrackHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Write the track histograms to the output file
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only write the loaded track types
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fTrackHistogramNames[iTrackType])) gDirectory->mkdir(fTrackHistogramNames[iTrackType]);
    gDirectory->cd(fTrackHistogramNames[iTrackType]);
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Track pT
      histogramNamer = Form("%sPt_C%d", fTrackHistogramNames[iTrackType], iCentralityBin);
      fhTrackPt[iTrackType][iCentralityBin]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      // pT integrated track phi
      histogramNamer = Form("%sPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, fnTrackPtBins);
      fhTrackPhi[iTrackType][iCentralityBin][fnTrackPtBins]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      // pT integrated track eta
      histogramNamer = Form("%sEta_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, fnTrackPtBins);
      fhTrackEta[iTrackType][iCentralityBin][fnTrackPtBins]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      // pT integrated track eta-phi
      histogramNamer = Form("%sEtaPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, fnTrackPtBins);
      if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][fnTrackPtBins]->Write(histogramNamer.Data(), TObject::kOverwrite);
      
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
        
        // Track phi in track pT bins
        histogramNamer = Form("%sPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, iTrackPtBin);
        fhTrackPhi[iTrackType][iCentralityBin][iTrackPtBin]->Write(histogramNamer.Data(), TObject::kOverwrite);
        
        // Track eta in track pT bins
        histogramNamer = Form("%sEta_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, iTrackPtBin);
        fhTrackEta[iTrackType][iCentralityBin][iTrackPtBin]->Write(histogramNamer.Data(), TObject::kOverwrite);
        
        // Track eta-phi in track pT bins
        histogramNamer = Form("%sEtaPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, iTrackPtBin);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][iTrackPtBin]->Write(histogramNamer.Data(), TObject::kOverwrite);
        
      } // Track pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Track category loop

}

/*
 * Write the multiplicity histograms within the jet cone to the file that is currently open
 */
void EECHistogramManager::WriteMultiplicityInJetConeHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Only write the multiplicity within jet cone histograms if they are previously loaded
  if(fLoadMultiplicityInJetHistograms){
    
    for(int iMultiplicityType = 0; iMultiplicityType < knMultiplicityInJetConeTypes; iMultiplicityType++){
      
      // Create a directory for the histograms if it does not already exist
      if(!gDirectory->GetDirectory(fMultiplicityInJetsHistogramNames[iMultiplicityType])) gDirectory->mkdir(fMultiplicityInJetsHistogramNames[iMultiplicityType]);
      gDirectory->cd(fMultiplicityInJetsHistogramNames[iMultiplicityType]);
      
      // Loop over centrality
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Loop over track pT
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
          
          // Write histograms without jet pT binning
          histogramNamer = Form("%s_C%dT%d", fMultiplicityInJetsHistogramNames[iMultiplicityType], iCentrality, iTrackPt);
          if(fhMultiplicityInJetCone[iCentrality][fnJetPtBinsEEC][iTrackPt][iMultiplicityType][EECHistograms::knSubeventTypes]) fhMultiplicityInJetCone[iCentrality][fnJetPtBinsEEC][iTrackPt][iMultiplicityType][EECHistograms::knSubeventTypes]->Write(histogramNamer.Data(), TObject::kOverwrite);
          
          // For PbPb or pPb MC, write histograms without jet pT and with subevent type binning
          if(fSystemAndEnergy.Contains("Pb MC")){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){
              
              // Write the energy-energy correlator histograms with subevent binning
              histogramNamer = Form("%s_C%dT%dS%d", fMultiplicityInJetsHistogramNames[iMultiplicityType], iCentrality, iTrackPt, iSubevent);
              if(fhMultiplicityInJetCone[iCentrality][fnJetPtBinsEEC][iTrackPt][iMultiplicityType][iSubevent]) fhMultiplicityInJetCone[iCentrality][fnJetPtBinsEEC][iTrackPt][iMultiplicityType][iSubevent]->Write(histogramNamer.Data(), TObject::kOverwrite);
              
            } // Subevent type loop
          } // Data is PbPb or pPb MC
          
          // Loop over jet pT
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
            
            // Write the energy-energy correlator histograms
            histogramNamer = Form("%s_C%dT%dJ%d", fMultiplicityInJetsHistogramNames[iMultiplicityType], iCentrality, iTrackPt, iJetPt);
            if(fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][EECHistograms::knSubeventTypes]) fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][EECHistograms::knSubeventTypes]->Write(histogramNamer.Data(), TObject::kOverwrite);
            
            // For PbPb or pPb MC, write histograms without jet pT and with subevent type binning
            if(fSystemAndEnergy.Contains("Pb MC")){
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){
                
                // Write the energy-energy correlator histograms with subevent binning
                histogramNamer = Form("%s_C%dT%dJ%dS%d", fMultiplicityInJetsHistogramNames[iMultiplicityType], iCentrality, iTrackPt, iJetPt, iSubevent);
                if(fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubevent]) fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubevent]->Write(histogramNamer.Data(), TObject::kOverwrite);
                
              } // Subevent type loop
            } // Data is PbPb or pPb MC
            
          } // Loop over jet pT bins
        } // Loop over track pT bins
      } // Loop over centrality bins
      
      // Return back to main directory
      gDirectory->cd("../");
    } // Multiplicity type loop (regular/uncorrected)
    
  } // Check the multiplicity histograms are loaded
  
}

/*
 * Write the particle density histograms around the jet axes to the file that is currently open
 */
void EECHistogramManager::WriteParticleDensityAroundJetsHistograms(){

  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Loop over particle density types
  for(int iParticleDensityType = 0; iParticleDensityType < knParticleDensityAroundJetAxisTypes; iParticleDensityType++){

    // Only write the histograms that have been loaded
    if(!fLoadParticleDensityAroundJetsHistograms[iParticleDensityType]) continue;
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fParticleDensityAroundJetsSaveNames[iParticleDensityType])) gDirectory->mkdir(fParticleDensityAroundJetsSaveNames[iParticleDensityType]);
    gDirectory->cd(fParticleDensityAroundJetsSaveNames[iParticleDensityType]);

    // Loop over jet cone type
    for(int iJetConeType = 0; iJetConeType < EECHistograms::knJetConeTypes; iJetConeType++){

      // Loop over centrality
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

        // Loop over track pT
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

          // Write histograms without jet pT binning
          histogramNamer = Form("%s%s_C%dT%d", fParticleDensityAroundJetsSaveNames[iParticleDensityType], fJetConeTypeSaveName[iJetConeType], iCentrality, iTrackPt);
          if(fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes]) fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes]->Write(histogramNamer.Data(), TObject::kOverwrite);

          // For PbPb and pPb MC, write histograms without jet pT and with subevent type binning
          if(fSystemAndEnergy.Contains("Pb MC")){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){

              // Write the energy-energy correlator histograms with subevent binning
              histogramNamer = Form("%s%s_C%dT%dS%d", fParticleDensityAroundJetsSaveNames[iParticleDensityType], fJetConeTypeSaveName[iJetConeType], iCentrality, iTrackPt, iSubevent);
              if(fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][iSubevent]) fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][iSubevent]->Write(histogramNamer.Data(), TObject::kOverwrite);

            } // Subevent type loop
          } // Data is PbPb or pPb MC

          // Loop over jet pT
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){

            // Write the energy-energy correlator histograms
            histogramNamer = Form("%s%s_C%dT%dJ%d", fParticleDensityAroundJetsSaveNames[iParticleDensityType], fJetConeTypeSaveName[iJetConeType], iCentrality, iTrackPt, iJetPt);
            if(fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes]) fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes]->Write(histogramNamer.Data(), TObject::kOverwrite);

            // For PbPb and pPb MC, loop over subevent types
            if(fSystemAndEnergy.Contains("Pb MC")){
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){

              // Write the energy-energy correlator histograms with subevent binning
                histogramNamer = Form("%s%s_C%dT%dJ%dS%d", fParticleDensityAroundJetsSaveNames[iParticleDensityType], fJetConeTypeSaveName[iJetConeType], iCentrality, iTrackPt, iJetPt, iSubevent);
                if(fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubevent]) fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubevent]->Write(histogramNamer.Data(), TObject::kOverwrite);

              } // Subevent type loop
            } // Data is PbPb or pPb MC
            
          } // Loop over jet pT bins
        } // Loop over track pT bins
      } // Loop over centrality bins
    }  // Loop over jet cone types

    // Return back to main directory
    gDirectory->cd("../");

  } // Loop over different particle density types
  
}

/*
 * Write the maximum particle pT within the jet cone histograms to the file that is currently open
 */
void EECHistogramManager::WriteMaxParticlePtWithinJetConeHistograms(){
  
  // Can only write the histograms if they are loaded
  if(!fLoadMaxParticlePtWithinJetConeHistograms) return;
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  for(int iMaxParticlePtType = 0; iMaxParticlePtType < knMaxParticlePtWithinJetConeTypes; iMaxParticlePtType++){
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fMaxParticlePtInJetConeSaveName[iMaxParticlePtType])) gDirectory->mkdir(fMaxParticlePtInJetConeSaveName[iMaxParticlePtType]);
    gDirectory->cd(fMaxParticlePtInJetConeSaveName[iMaxParticlePtType]);
    
    // Loop over centrality
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Loop over jet pT
      for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
        
        // Write the histograms without track pT selection
        histogramNamer = Form("%s_C%dJ%d",fMaxParticlePtInJetConeSaveName[iMaxParticlePtType], iCentrality, iJetPt);
        if(fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][knProjectedMaxParticlePtBins]) fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][knProjectedMaxParticlePtBins]->Write(histogramNamer.Data(), TObject::kOverwrite);
        
        // Loop over track pT
        for(int iTrackPt = 0; iTrackPt < knProjectedMaxParticlePtBins; iTrackPt++){
          
          // Write track pT binned histograms
          histogramNamer = Form("%s_C%dJ%dT%d",fMaxParticlePtInJetConeSaveName[iMaxParticlePtType], iCentrality, iJetPt, iTrackPt);
          if(fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt]) fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
          
          // Write track pT cutted histograms
          histogramNamer = Form("%sTrackPtCut_C%dJ%dT%d",fMaxParticlePtInJetConeSaveName[iMaxParticlePtType], iCentrality, iJetPt, iTrackPt);
          if(fhMaxParticlePtInJetConePtCut[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt]) fhMaxParticlePtInJetConePtCut[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
          
        } // Loop over track pT bins
      } // Loop over jet pT bins
    } // Loop over centrality bins
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Maximum particle pT within jet cone type loop

}

/*
 * Write the energy-energy correlator histograms to the file that is currently open
 */
void EECHistogramManager::WriteEnergyEnergyCorrelatorHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Loop over energy-energy correlator histogram types
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    
    // Only write the histograms that have been loaded
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType])) gDirectory->mkdir(fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
    gDirectory->cd(fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
    
    // Loop over all pair types that are available in the file
    for(int iPairingType : *fIncludedPairingTypes){

      // Create a subdirectories for different pairing types.
      if(!gDirectory->GetDirectory(fPairingTypeSaveName[iPairingType])) gDirectory->mkdir(fPairingTypeSaveName[iPairingType]);
      gDirectory->cd(fPairingTypeSaveName[iPairingType]);
      
      // Loop over centrality
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Loop over track pT
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
          
          // Write histograms without jet pT binning
          histogramNamer = Form("%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt);
          if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations]) fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations]->Write(histogramNamer.Data(), TObject::kOverwrite);
          
          // For PbPb and pPb MC, write histograms without jet pT and with subevent type binning
          if(fSystemAndEnergy.Contains("Pb MC")){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
              
              // Write the energy-energy correlator histograms with subevent binning
              histogramNamer = Form("%s%s_C%dT%dS%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iSubevent);
              if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent]) fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent]->Write(histogramNamer.Data(), TObject::kOverwrite);
              
            } // Subevent type loop
          } // Data is PbPb or pPb MC
          
          // Loop over jet pT
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
            
            // Write the energy-energy correlator histograms
            histogramNamer = Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iJetPt);
            if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations]) fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations]->Write(histogramNamer.Data(), TObject::kOverwrite);
            
            // For PbPb and pPb MC, loop over subevent types
            if(fSystemAndEnergy.Contains("Pb MC")){
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
                
                // Write the energy-energy correlator histograms with subevent binning
                histogramNamer = Form("%s%s_C%dT%dJ%dS%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iJetPt, iSubevent);
                if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent]) fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent]->Write(histogramNamer.Data(), TObject::kOverwrite);
                
              } // Subevent type loop
            } // Data is PbPb or pPb MC
            
          } // Loop over jet pT bins
        } // Loop over track pT bins
      } // Loop over centrality bins

      // Leave the pairing type directory
      gDirectory->cd("../");

    } // Pairing type loop (same jet/reflected cone)
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Loop over different energy-energy correlator types
  
}

/*
 * Write the jet pT response matrices to the file that is currently open
 */
void EECHistogramManager::WriteJetPtResponseMatrix(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Only write the jet pT closure histograms if they are previously loaded
  if(fLoadJetPtResponseMatrix){
    
    // Create a directory for the histograms if it does not already exist
    histogramNamer = "jetPtResponseMatrix";
    if(!gDirectory->GetDirectory(histogramNamer.Data())) gDirectory->mkdir(histogramNamer.Data());
    gDirectory->cd(histogramNamer.Data());

    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Only write histograms that are non-NULL
      if(fhJetPtResponseMatrix[iCentrality]) {
        histogramNamer = Form("jetPtResponseMatrix_C%d", iCentrality);
        fhJetPtResponseMatrix[iCentrality]->Write(histogramNamer.Data(), TObject::kOverwrite);
      }

    }  // Centrality loop

    // Return back to main directory
    gDirectory->cd("../");
    
  } // Writing jet pT response matrices
  
}

/*
 * Write the reflected cone QA histograms to the file that is currently open
 */
void EECHistogramManager::WriteReflectedConeQAHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Write the jet histograms to the output file
  if(!fLoadReflectedConeQAHistograms) return;  // Only write the reflected cone QA histograms if they are loaded
  
  // Create a directory for the histograms if it does not already exist
  if(!gDirectory->GetDirectory(fReflectedConeQAFolderName)) gDirectory->mkdir(fReflectedConeQAFolderName);
  gDirectory->cd(fReflectedConeQAFolderName);
  
  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
    
    // Number of jets within the reflected cone
    histogramNamer = Form("%s_C%d",fNumberOfJetsWithinReflectedConeName,iCentrality);
    if(fhNumberOfJetsWithinReflectedCone[iCentrality]) fhNumberOfJetsWithinReflectedCone[iCentrality]->Write(histogramNamer.Data(), TObject::kOverwrite);
    
    // Jet pT within the reflected cone
    histogramNamer = Form("%s_C%d",fJetPtWithinReflectedConeName,iCentrality);
    if(fhJetPtWithinReflectedCone[iCentrality]) fhJetPtWithinReflectedCone[iCentrality]->Write(histogramNamer.Data(), TObject::kOverwrite);
    
  } // Loop over centrality bins
  
  // Return back to main directory
  gDirectory->cd("../");
  
}

/*
 * Write the closure histograms to the file that is currently open
 */
void EECHistogramManager::WriteClosureHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Only write the jet pT closure histograms if they are previously loaded
  if(fLoadJetPtClosureHistograms){
    
    // Create a directory for the histograms if it does not already exist
    histogramNamer = Form("jetPtClosure_%s",fJetHistogramName);
    if(!gDirectory->GetDirectory(histogramNamer.Data())) gDirectory->mkdir(histogramNamer.Data());
    gDirectory->cd(histogramNamer.Data());
    
    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Loop over closure particles (quark/gluon/no selection)
      for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
        
        // Loop over generator level jet pT bins
        for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
          
          // Loop over jet eta bins
          for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){

            // Loop over jet phi bins
            for(int iJetPhi = 0; iJetPhi <= knJetPhiBins; iJetPhi++){
            
              // Only write histogram that are non-NULL
              if(fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iClosureParticle]){
                histogramNamer = Form("jetPtClosure_%s%s_C%d", fJetHistogramName, fClosureParticleName[iClosureParticle], iCentrality);
                if(iGenJetPt < knGenJetPtBins) histogramNamer.Append(Form("T%d",iGenJetPt));
                if(iJetEta < knJetEtaBins) histogramNamer.Append(Form("E%d",iJetEta));
                if(iJetPhi < knJetPhiBins) histogramNamer.Append(Form("P%d",iJetPhi));
                fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iClosureParticle]->Write(histogramNamer.Data(), TObject::kOverwrite);
              }
            
            } // JEt phi bin loop
          } // Jet eta bin loop
        } // Generator level jet pT loop
      } // Closure particle type (quark/gluon) loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Writing jet pT closure histograms
  
}

/*
 * Write the jet pT unfolding histograms to the file that is currently open
 */
void EECHistogramManager::WriteJetPtUnfoldingHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Only write the jet pT unfolding histograms if they are previously loaded
  if(fLoadJetPtUnfoldingHistograms){
    
    for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
      if(!gDirectory->GetDirectory(fJetPtUnfoldingDistributionName[iUnfoldType])) gDirectory->mkdir(fJetPtUnfoldingDistributionName[iUnfoldType]);
      gDirectory->cd(fJetPtUnfoldingDistributionName[iUnfoldType]);

      // Centrality loop
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

        // Track pT loop
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

          // Only write histograms that are non-NULL
          if(fhJetPtUnfoldingDistribution[iUnfoldType][iCentrality][iTrackPt]){
            histogramNamer = Form("%s_C%dT%d", fJetPtUnfoldingDistributionName[iUnfoldType], iCentrality, iTrackPt);
            fhJetPtUnfoldingDistribution[iUnfoldType][iCentrality][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
          }
        } // Track pT loop
      } // Centrality loop

      // Return back to main directory
      gDirectory->cd("../");
    } // Unfold type loop

    if(!gDirectory->GetDirectory(fJetPtResponseMatrixName)) gDirectory->mkdir(fJetPtResponseMatrixName);
    gDirectory->cd(fJetPtResponseMatrixName);

    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Track pT loop
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

        // Only write histograms that are non-NULL
        if(fhJetPtUnfoldingResponse[iCentrality][iTrackPt]){
          histogramNamer = Form("%s_C%dT%d", fJetPtResponseMatrixName, iCentrality, iTrackPt);
          fhJetPtUnfoldingResponse[iCentrality][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
        }

      } // Track pT loop
    } // Centrality loop

    // Return back to main directory
    gDirectory->cd("../");
  }
}

/*
 * Write the covariance histograms used in jet pT unfolding
 */
void EECHistogramManager::WriteJetPtUnfoldingCovariance(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Only write the covariance histograms used in jet pT unfolding if they are previously loaded
  if(fLoadJetPtUnfoldingCovariance){

    // Covariance matrix type loop
    for(int iCoverianceMatrixType = 0; iCoverianceMatrixType < knCovarianceMatrixTypes; iCoverianceMatrixType++){

      if(!gDirectory->GetDirectory(fJetPtCovarianceMatrixName[iCoverianceMatrixType])) gDirectory->mkdir(fJetPtCovarianceMatrixName[iCoverianceMatrixType]);
      gDirectory->cd(fJetPtCovarianceMatrixName[iCoverianceMatrixType]);

      // Centrality loop
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

        // Track pT loop
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

          // Only write histograms that are non-NULL
          if(fhJetPtUnfoldingCovariance[iCoverianceMatrixType][iCentrality][iTrackPt]){
            histogramNamer = Form("%s_C%dT%d", fJetPtCovarianceMatrixName[iCoverianceMatrixType], iCentrality, iTrackPt);
            fhJetPtUnfoldingCovariance[iCoverianceMatrixType][iCentrality][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
          }

        } // Track pT loop
      } // Centrality loop
      
      // Return back to main directory
      gDirectory->cd("../");

    } // Covariance matrix type loop  
  }
}

/*
 * Write the jet pT one dimensional unfolding histograms to the file that is currently open
 */
void EECHistogramManager::WriteJetPtOneDimensionalUnfoldingHistograms(){
  
  // Helper variable for histogram naming
  TString histogramNamer;
  
  // Only write the jet pT unfolding histograms if they are previously loaded
  if(fLoadJetPtOneDimensionalUnfoldingHistograms){
    
    for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
      if(!gDirectory->GetDirectory(fJetPtOneDimensionalUnfoldingDistributionName[iUnfoldType])) gDirectory->mkdir(fJetPtOneDimensionalUnfoldingDistributionName[iUnfoldType]);
      gDirectory->cd(fJetPtOneDimensionalUnfoldingDistributionName[iUnfoldType]);

      // Centrality loop
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

        // Only write histograms that are non-NULL
        if(fhOneDimensionalJetPtUnfoldingDistribution[iUnfoldType][iCentrality]){
          histogramNamer = Form("%s_C%d", fJetPtOneDimensionalUnfoldingDistributionName[iUnfoldType], iCentrality);
          fhOneDimensionalJetPtUnfoldingDistribution[iUnfoldType][iCentrality]->Write(histogramNamer.Data(), TObject::kOverwrite);
        }
      } // Centrality loop

      // Return back to main directory
      gDirectory->cd("../");
    } // Unfold type loop

    if(!gDirectory->GetDirectory(fJetPtOneDimensionalResponseMatrixName)) gDirectory->mkdir(fJetPtOneDimensionalResponseMatrixName);
    gDirectory->cd(fJetPtOneDimensionalResponseMatrixName);

    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Only write histograms that are non-NULL
      if(fhOneDimensionalJetPtUnfoldingResponse[iCentrality]){
        histogramNamer = Form("%s_C%d", fJetPtOneDimensionalResponseMatrixName, iCentrality);
        fhOneDimensionalJetPtUnfoldingResponse[iCentrality]->Write(histogramNamer.Data(), TObject::kOverwrite);
      }

    } // Centrality loop

    // Return back to main directory
    gDirectory->cd("../");
  }
}

/*
 * Write the track/particle matching study histograms to the file that is currently open
 */
void EECHistogramManager::WriteTrackParticleMatchingHistograms(){
  
  // Only write the track/particle matching histograms if they are loaded
  if(!fLoadTrackParticleMatchingHistograms) return;  

  // Helper variable for histogram naming
  TString histogramNamer;
  
  // First, write the track/particle matching QA histograms
  for(int iTrackParticleQAHistogram = 0; iTrackParticleQAHistogram < knTrackParticleMatchingQAHistograms; iTrackParticleQAHistogram++){

    // Change to a folder where the track/particle matching QA histograms are saved
    if(!gDirectory->GetDirectory(fTrackParticleMatchingQAName[iTrackParticleQAHistogram])) gDirectory->mkdir(fTrackParticleMatchingQAName[iTrackParticleQAHistogram]);
    gDirectory->cd(fTrackParticleMatchingQAName[iTrackParticleQAHistogram]);

    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Jet pT loop
      for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
        
        // Track pT loop
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

          // Only write histograms that are non-NULL
          if(fhTrackParticleMatchQA[iTrackParticleQAHistogram][iCentrality][iJetPt][iTrackPt]){
            histogramNamer = Form("%s_C%dJ%dT%d", fTrackParticleMatchingQAName[iTrackParticleQAHistogram], iCentrality, iJetPt, iTrackPt);
            fhTrackParticleMatchQA[iTrackParticleQAHistogram][iCentrality][iJetPt][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
          }
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

    // Return back to main directory
    gDirectory->cd("../");

  } // Track/particle matching QA type loop

  // Then, write the track/particle matching response matrices
  for(int iTrackParticleResponseHistogram = 0; iTrackParticleResponseHistogram < knTrackParticleMatchingResponseTypes; iTrackParticleResponseHistogram++){

    // Change to a folder where the track/particle matching response histograms are saved
    if(!gDirectory->GetDirectory(fTrackParticleMatchingResponseName[iTrackParticleResponseHistogram])) gDirectory->mkdir(fTrackParticleMatchingResponseName[iTrackParticleResponseHistogram]);
    gDirectory->cd(fTrackParticleMatchingResponseName[iTrackParticleResponseHistogram]);

    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      // Jet pT loop
      for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
        
        // Track pT loop
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

          // Only write histograms that are non-NULL
          if(fhTrackParticleResponse[iTrackParticleResponseHistogram][iCentrality][iJetPt][iTrackPt]){
            histogramNamer = Form("%s_C%dJ%dT%d", fTrackParticleMatchingResponseName[iTrackParticleResponseHistogram], iCentrality, iJetPt, iTrackPt);
            fhTrackParticleResponse[iTrackParticleResponseHistogram][iCentrality][iJetPt][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
          }
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

    // Return back to main directory
    gDirectory->cd("../");

  } // Track/particle matching QA type loop

  // Finally, write the track/particle matching pair pT closure histograms

  // Change to a folder where the track/particle matching pT closure histograms are saved
  if(!gDirectory->GetDirectory(fTrackParticlePtClosureSaveName)) gDirectory->mkdir(fTrackParticlePtClosureSaveName);
  gDirectory->cd(fTrackParticlePtClosureSaveName);
  
  // Centrality loop
  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
    // Jet pT loop
    for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
      // Track pT loop
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        // Only write histograms that are non-NULL
        if(fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt]){
          histogramNamer = Form("%s_C%dJ%dT%d", fTrackParticlePtClosureSaveName, iCentrality, iJetPt, iTrackPt);
          fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
        }
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Return back to main directory
  gDirectory->cd("../");
}


/*
 * Write the processed histograms into a file
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void EECHistogramManager::WriteProcessed(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile* outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  TString folderNamer;
  
  // Go through the processed histograms and write them to file
  
  // Loop over energy-energy correlator histogram types
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    
    // Only write the histograms that have been loaded
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    // Create a directory for the histograms if it does not already exist
    folderNamer = Form("%sProcessed", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
    if(!gDirectory->GetDirectory(folderNamer)) gDirectory->mkdir(folderNamer);
    gDirectory->cd(folderNamer);
    
    // Loop over selected bin range
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        
        // =================================== //
        // Histograms without jet pT selection //
        // =================================== //
        
        // Loop over the different processing steps
        for(int iProcessingLevel = 0; iProcessingLevel < kEnergyEnergyCorrelatorUnfolded; iProcessingLevel++){
          
          // Write processed histograms
          if(fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iProcessingLevel]) fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iProcessingLevel]->Write(nullptr, TObject::kOverwrite);
          
        } // Loop over different processing levels
        
        // =================================== //
        // Histograms in different jet pT bins //
        // =================================== //
        
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
          
          // Loop over the different processing steps
          for(int iProcessingLevel = 0; iProcessingLevel < kEnergyEnergyCorrelatorUnfolded; iProcessingLevel++){
            
            // Write processed histograms
            if(fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel]) fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel]->Write(nullptr, TObject::kOverwrite);
            
          } // Loop over different processing levels
        } // Loop over jet pT bins
      } // Loop over track pT bins
    } // Loop over centrality bins
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Loop over different energy-energy correlator types
  
  // If there a JCard directory, update the git hash for histogram processing. If not, write a JCard directory with all the information
  if(!gDirectory->GetDirectory("JCard")){
    fCard->Write(outputFile);
  } else {
    fCard->WriteProcessHash(outputFile);
    fCard->WriteBackground(outputFile);
  }
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Write the mixed cone histograms after combining them
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void EECHistogramManager::WriteCombinedMixedConeHistograms(const char* fileName, const char* fileOption){

  // Create the output file
  TFile* outputFile = new TFile(fileName,fileOption);

  // Write the energy-energy correlator histograms to the output file
  WriteEnergyEnergyCorrelatorHistograms();

  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;

}

/*
 * Write the unfolded energy-energy correlators to a file
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void EECHistogramManager::WriteUnfoldedEnergyEnergyCorrelators(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile* outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  TString folderNamer;
  
  // Go through the unfolded histograms and write them to file
  
  // Loop over energy-energy correlator histogram types
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    
    // Use this flag in the other macro to indicate which types of correlators have been unfolded
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    // Create a directory for the histograms if it does not already exist
    folderNamer = Form("%sProcessed", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
    if(!gDirectory->GetDirectory(folderNamer)) gDirectory->mkdir(folderNamer);
    gDirectory->cd(folderNamer);
    
    // Loop over selected bin range
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
            
          // Write unfolded histograms
          if(fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorUnfolded]) fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorUnfolded]->Write(nullptr, TObject::kOverwrite);
            
        } // Loop over jet pT bins
      } // Loop over track pT bins
    } // Loop over centrality bins
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Loop over different energy-energy correlator types
  
  // If a JCard does not exist in the target file, create one 
  if(!gDirectory->GetDirectory("JCard")){
    fCard->Write(outputFile);
  } else {
    fCard->WriteUnfoldInfo(outputFile);
  }
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Write the processed energy-energy correlators after unfolding into a file
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void EECHistogramManager::WriteProcessedAfterUnfolding(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile* outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  TString folderNamer;
  
  // Go through the processed histograms and write them to file
  
  // Loop over energy-energy correlator histogram types
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    
    // Only write the histograms that have been loaded
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    // Create a directory for the histograms if it does not already exist
    folderNamer = Form("%sProcessed", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
    if(!gDirectory->GetDirectory(folderNamer)) gDirectory->mkdir(folderNamer);
    gDirectory->cd(folderNamer);
    
    // Loop over selected bin range
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
          
          // Loop over the different processing steps
          for(int iProcessingLevel = kEnergyEnergyCorrelatorBackgroundAfterUnfolding; iProcessingLevel < knEnergyEnergyCorrelatorProcessingLevels; iProcessingLevel++){
            
            // Write processed histograms
            if(fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel]) fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel]->Write(nullptr, TObject::kOverwrite);
            
          } // Loop over different processing levels
        } // Loop over jet pT bins
      } // Loop over track pT bins
    } // Loop over centrality bins
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Loop over different energy-energy correlator types
  
  // If there a JCard directory, update the git hash for histogram processing. If not, write a JCard directory with all the information
  if(!gDirectory->GetDirectory("JCard")){
    fCard->Write(outputFile);
  } else {
    fCard->WriteProcessHash(outputFile);
    fCard->WriteBackground(outputFile);
  }
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Write the covariance matrices after unfolding
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void EECHistogramManager::WriteCovarianceMatrixAfterUnfolding(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile* outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  TString histogramNamer;
  
  // If there is not a directory for unfolded covariance matrices yet, create one
  if(!gDirectory->GetDirectory(fJetPtCovarianceMatrixName[kCovarianceMatrixUnfolded])) gDirectory->mkdir(fJetPtCovarianceMatrixName[kCovarianceMatrixUnfolded]);
  gDirectory->cd(fJetPtCovarianceMatrixName[kCovarianceMatrixUnfolded]);

  // Centrality loop
  for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

    // Track pT loop
    for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

      // Only write histograms that are non-NULL
      if(fhJetPtUnfoldingCovariance[kCovarianceMatrixUnfolded][iCentrality][iTrackPt]){
        histogramNamer = Form("%s_C%dT%d", fJetPtCovarianceMatrixName[kCovarianceMatrixUnfolded], iCentrality, iTrackPt);
        fhJetPtUnfoldingCovariance[kCovarianceMatrixUnfolded][iCentrality][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
      }

    } // Track pT loop
  } // Centrality loop

  // Return back to main directory
  gDirectory->cd("../");  
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Load the selected histograms from a file containing readily processed histograms
 *
 * OBSOLETE! This code should not be used anymore. It is left here for backwards compatibility.
 * The new histogram getters load histograms when needed, so pre-loading is not needed anymore.
 * This makes the code use significantly less memory.
 */
[[deprecated("EECHistogramManager::LoadProcessedHistograms is obsolete. It should not be anymore to improve memory usage!")]]
void EECHistogramManager::LoadProcessedHistograms(){

  cout << endl;
  cout << "EECHistogramManager::LoadProcessedHistograms is obsolete. It should not be anymore to improve memory usage!" << endl;
  cout << endl;
  
  // Helper variable for finding names of loaded histograms
  TString histogramNamer;
  TString folderNamer;
  TH1D* testHistogram;
  bool legacyEnergyEnergyCorrelatorMode = false; // Older files have different subevent combination indexing for energy-energy correlators. Take that into account here
  bool preMixingCompatibilityMode = false; // Before event mixing for energy-energy correlators were added, there was a different folder structure to hold energy-energy correlator histograms. Take that into account here 
  int subeventIndex;
  const char* currentPairingType0;
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhTriggers = (TH1D*) fInputFile->Get("triggers");                      // Trigger selection information
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");                    // Number of tracks surviving different track cuts
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events

    fhHFPlus[fnJetPtBinsEEC][knMaxJetsBins] = (TH1D*) fInputFile->Get("hfEnergy/HFPlus");    // Energy in HF plus calorimeters
    fhHFMinus[fnJetPtBinsEEC][knMaxJetsBins] = (TH1D*) fInputFile->Get("hfEnergy/HFMinus");  // Energy in HF minus calorimeters
    fhHFSum[fnJetPtBinsEEC][knMaxJetsBins] = (TH1D*) fInputFile->Get("hfEnergy/HFSum");      // Energy in HF calorimeters
    
    for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
      
      histogramNamer = Form("multiplicity/multiplicity_C%d", iCentrality);
      fhMultiplicity[iCentrality] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      histogramNamer = Form("multiplicity/multiplicityWeighted_C%d", iCentrality);
      fhMultiplicityWeighted[iCentrality] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
    }
    
    // Multiplicity histograms without centrality bins
    fhMultiplicity[fnCentralityBins] = (TH1D*) fInputFile->Get("multiplicity/multiplicity");
    
    fhMultiplicityWeighted[fnCentralityBins] = (TH1D*) fInputFile->Get("multiplicity/multiplicityWeighted");
    
    // Multiplicity vs. centrality maps
    fhMultiplicityMap = (TH2D*) fInputFile->Get("multiplicity/multiplicityMap");
    
    fhMultiplicityMapWeighted = (TH2D*) fInputFile->Get("multiplicity/multiplicityMapWeighted");
    
  }
  
  // Load the jet histograms from the input file
  
  for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
    
    // Always load jet pT histograms
    histogramNamer = Form("%s/%sPt_C%d", fJetHistogramName, fJetHistogramName, iCentralityBin);
    fhJetPt[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer.Data());
    
    if(!fLoadJets) continue;  // Only load the loaded the selected histograms
    
    // Jet phi
    histogramNamer = Form("%s/%sPhi_C%d", fJetHistogramName, fJetHistogramName, iCentralityBin);
    fhJetPhi[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer.Data());
    
    // Jet eta
    histogramNamer = Form("%s/%sEta_C%d", fJetHistogramName, fJetHistogramName, iCentralityBin);
    fhJetEta[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer.Data());

    // Jet eta_CM
    histogramNamer = Form("%s/%sEtaCM_C%d", fJetHistogramName, fJetHistogramName, iCentralityBin);
    fhJetEtaCM[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer.Data());
    
    // Jet eta-phi
    histogramNamer = Form("%s/%sEtaPhi_C%d", fJetHistogramName, fJetHistogramName, iCentralityBin);
    if(fLoad2DHistograms) fhJetEtaPhi[iCentralityBin] = (TH2D*) fInputFile->Get(histogramNamer.Data());
  } // Loop over centrality bins
  
  
  // Load the track histograms from the input file
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only load the selected track types
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Track pT
      histogramNamer = Form("%s/%sPt_C%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin);
      fhTrackPt[iTrackType][iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      // pT integrated track phi
      histogramNamer = Form("%s/%sPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,fnTrackPtBins);
      fhTrackPhi[iTrackType][iCentralityBin][fnTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      // pT integrated track eta
      histogramNamer = Form("%s/%sEta_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,fnTrackPtBins);
      fhTrackEta[iTrackType][iCentralityBin][fnTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      
      // pT integrated track eta-phi
      histogramNamer = Form("%s/%sEtaPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,fnTrackPtBins);
      if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][fnTrackPtBins] = (TH2D*) fInputFile->Get(histogramNamer.Data());
      
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
        
        // Track phi in track pT bins
        histogramNamer = Form("%s/%sPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,iTrackPtBin);
        fhTrackPhi[iTrackType][iCentralityBin][iTrackPtBin] = (TH1D*) fInputFile->Get(histogramNamer.Data());
        
        // Track eta in track pT bins
        histogramNamer = Form("%s/%sEta_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,iTrackPtBin);
        fhTrackEta[iTrackType][iCentralityBin][iTrackPtBin] = (TH1D*) fInputFile->Get(histogramNamer.Data());
        
        // Track eta-phi in track pT bins
        histogramNamer = Form("%s/%sEtaPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,iTrackPtBin);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][iTrackPtBin] = (TH2D*) fInputFile->Get(histogramNamer.Data());
        
      } // Track pT loop
    } // Centrality loop
  } // Track category loop
  
  // Load the multiplicity histograms within the jet cone
  if(fLoadMultiplicityInJetHistograms){
    for(int iMultiplicityType = 0; iMultiplicityType < knMultiplicityInJetConeTypes; iMultiplicityType++){
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
          
          // Load the histograms without jet pT binning
          histogramNamer = Form("%s/%s_C%dT%d", fMultiplicityInJetsHistogramNames[iMultiplicityType], fMultiplicityInJetsHistogramNames[iMultiplicityType], iCentrality, iTrackPt);
          fhMultiplicityInJetCone[iCentrality][fnJetPtBinsEEC][iTrackPt][iMultiplicityType][EECHistograms::knSubeventTypes] = (TH1D*) fInputFile->Get(histogramNamer.Data());
          
          // For PbPb and pPb MC, load histograms without jet pT and with subevent type binning
          if(fSystemAndEnergy.Contains("Pb MC")){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){

              // Load the particle density histograms with subevent binning
              histogramNamer = Form("%s/%s_C%dT%dS%d", fMultiplicityInJetsHistogramNames[iMultiplicityType], fMultiplicityInJetsHistogramNames[iMultiplicityType], iCentrality, iTrackPt, iSubevent);
              fhMultiplicityInJetCone[iCentrality][fnJetPtBinsEEC][iTrackPt][iMultiplicityType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());

            } // Subevent type loop
          } // Data is PbPb or pPb MC
          
          // Loop over jet pT
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
            
            // Make sure that jet pT integrated histogram is not overwritten by null
            if(fLastLoadedJetPtBinEEC >= fnJetPtBinsEEC) continue;
            
            // Load the multiplicity within the jet cone histograms
            histogramNamer = Form("%s/%s_C%dT%dJ%d", fMultiplicityInJetsHistogramNames[iMultiplicityType], fMultiplicityInJetsHistogramNames[iMultiplicityType], iCentrality, iTrackPt, iJetPt);
            fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][EECHistograms::knSubeventTypes]  = (TH1D*) fInputFile->Get(histogramNamer.Data());
            
            // For PbPb and pPb MC, load histograms without jet pT and with subevent type binning
            if(fSystemAndEnergy.Contains("Pb MC")){
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){

                // Load the particle density histograms with subevent binning
                histogramNamer = Form("%s/%s_C%dT%dJ%dS%d", fMultiplicityInJetsHistogramNames[iMultiplicityType], fMultiplicityInJetsHistogramNames[iMultiplicityType], iCentrality, iTrackPt, iJetPt, iSubevent);
                fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());

              } // Subevent type loop
            } // Data is PbPb or pPb MC
            
          } // Jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // Multiplicity type loop
  }
  
  // Load the particle density around the jet axis histograms
  for(int iParticleDensityType = 0; iParticleDensityType < knParticleDensityAroundJetAxisTypes; iParticleDensityType++){

    // Only load the selected types of histograms
    if(!fLoadParticleDensityAroundJetsHistograms[iParticleDensityType]) continue;
    
    for(int iJetConeType = 0; iJetConeType < EECHistograms::knJetConeTypes; iJetConeType++){
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

          // Load the histograms without jet pT binning
          histogramNamer = Form("%s/%s%s_C%dT%d", fParticleDensityAroundJetsSaveNames[iParticleDensityType], fParticleDensityAroundJetsSaveNames[iParticleDensityType], fJetConeTypeSaveName[iJetConeType], iCentrality, iTrackPt);
          fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes] = (TH1D*) fInputFile->Get(histogramNamer.Data());

          // For PbPb and pPb MC, load histograms without jet pT and with subevent type binning
          if(fSystemAndEnergy.Contains("Pb MC")){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){

              // Load the particle density histograms with subevent binning
              histogramNamer = Form("%s/%s%s_C%dT%dS%d", fParticleDensityAroundJetsSaveNames[iParticleDensityType], fParticleDensityAroundJetsSaveNames[iParticleDensityType], fJetConeTypeSaveName[iJetConeType], iCentrality, iTrackPt, iSubevent);
              fhParticleDensityAroundJetAxis[iCentrality][fnJetPtBinsEEC][iTrackPt][iJetConeType][iParticleDensityType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());

            } // Subevent type loop
          } // Data is PbPb or pPb MC

          // Loop over jet pT
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){

            // Make sure that jet pT integrated histogram is not overwritten by null
            if(fLastLoadedJetPtBinEEC >= fnJetPtBinsEEC) continue;

            // Load the particle density histograms
            histogramNamer = Form("%s/%s%s_C%dT%dJ%d", fParticleDensityAroundJetsSaveNames[iParticleDensityType], fParticleDensityAroundJetsSaveNames[iParticleDensityType], fJetConeTypeSaveName[iJetConeType], iCentrality, iTrackPt, iJetPt);
            fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][EECHistograms::knSubeventTypes] = (TH1D*) fInputFile->Get(histogramNamer.Data());

            // For PbPb and pPb MC, loop over subevent types
            if(fSystemAndEnergy.Contains("Pb MC")){
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){

                // Load the particle density histograms with subevent binning
                histogramNamer = Form("%s/%s%s_C%dT%dJ%dS%d", fParticleDensityAroundJetsSaveNames[iParticleDensityType], fParticleDensityAroundJetsSaveNames[iParticleDensityType], fJetConeTypeSaveName[iJetConeType], iCentrality, iTrackPt, iJetPt, iSubevent);
                fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());

              } // Subevent type loop
            } // Data is PbPb or pPb MC

          } // Jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // Jet cone type loop
  } // Particle density type loop
  
  // Load the maximum particle pT within the jet cone histograms
  if(fLoadMaxParticlePtWithinJetConeHistograms){
    
    // Loop over maximum particle pT type types
    for(int iMaxParticlePtType = 0; iMaxParticlePtType < knMaxParticlePtWithinJetConeTypes; iMaxParticlePtType++){
  
      // Loop over centrality
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Loop over jet pT
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
          
          // Load the histograms without track pT selection
          histogramNamer = Form("%s/%s_C%dJ%d",fMaxParticlePtInJetConeSaveName[iMaxParticlePtType], fMaxParticlePtInJetConeSaveName[iMaxParticlePtType], iCentrality, iJetPt);
          fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][knProjectedMaxParticlePtBins] = (TH1D*) fInputFile->Get(histogramNamer.Data());
          
          // Loop over track pT
          for(int iTrackPt = 0; iTrackPt < knProjectedMaxParticlePtBins; iTrackPt++){
            
            // Load track pT binned histograms
            histogramNamer = Form("%s/%s_C%dJ%dT%d",fMaxParticlePtInJetConeSaveName[iMaxParticlePtType], fMaxParticlePtInJetConeSaveName[iMaxParticlePtType], iCentrality, iJetPt, iTrackPt);
            fhMaxParticlePtInJetConePtBin[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());
            
            // Load track pT cutted histograms
            histogramNamer = Form("%s/%sTrackPtCut_C%dJ%dT%d",fMaxParticlePtInJetConeSaveName[iMaxParticlePtType], fMaxParticlePtInJetConeSaveName[iMaxParticlePtType], iCentrality, iJetPt, iTrackPt);
            fhMaxParticlePtInJetConePtCut[iMaxParticlePtType][iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());
            
          } // Loop over track pT bins
        } // Loop over jet pT bins
      } // Loop over centrality bins
    } // Maximum particle pT within jet cone type loop
    
  } // Load the maximum particle pT within the jet cone histograms
  
  // Load the energy-energy correlator histograms from the input file
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    
    // Only load the selected types of histograms
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    // New folder structure was added together with event mixing for reflected cones. We need to figure out if the files are produced before or after the folder structure change to find the correct paths for histograms
    histogramNamer = Form("%s/%s/%s%s_C0T0", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[0], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[0]);
    testHistogram = (TH1D*) fInputFile->Get(histogramNamer.Data());
    if(testHistogram == NULL) preMixingCompatibilityMode = true;

    // Old file compatibility mode. There are less subevent combinations in old files. Take this into account when loading older files.
    if(preMixingCompatibilityMode){
      histogramNamer = Form("%s/%s_C0T0S3", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
      testHistogram = (TH1D*) fInputFile->Get(histogramNamer.Data());
      if(testHistogram == NULL) legacyEnergyEnergyCorrelatorMode = true;
    }

    // The name for the first pairing type was also changed after event mixing was added. Take this into account here
    if(preMixingCompatibilityMode){
      currentPairingType0 = fPairingTypeSaveName[0];
      fPairingTypeSaveName[0] = "";
    }
    
    for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){

      if(preMixingCompatibilityMode){
        folderNamer = Form("%s/", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
      } else {
        folderNamer = Form("%s/%s/", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType]);
      }

      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
          
          // Load the histograms without jet pT binning
          histogramNamer = Form("%s%s%s_C%dT%d", folderNamer.Data(), fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt);
          fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations] = (TH1D*) fInputFile->Get(histogramNamer.Data());
          
          // For PbPb and pPb MC, load histograms without jet pT and with subevent type binning
          if(fSystemAndEnergy.Contains("Pb MC")){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
              
              // Take into account that one subevent index is missing from the older files
              subeventIndex = iSubevent;
              if(legacyEnergyEnergyCorrelatorMode && iSubevent > EECHistograms::kPythiaHydjet) subeventIndex = iSubevent - 1;
              
              // Load the energy-energy correlator histograms with subevent binning
              histogramNamer = Form("%s%s%s_C%dT%dS%d", folderNamer.Data(), fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, subeventIndex);
              fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());
              
            } // Subevent type loop
          } // Data is PbPb or pPb MC
          
          // Loop over jet pT
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
            
            // Make sure that jet pT integrated histogram is not overwritten by null
            if(fLastLoadedJetPtBinEEC >= fnJetPtBinsEEC) continue;
            
            // Load the energy-energy correlator histograms
            histogramNamer = Form("%s%s%s_C%dT%dJ%d", folderNamer.Data(), fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iJetPt);
            fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][EECHistograms::knSubeventCombinations] = (TH1D*) fInputFile->Get(histogramNamer.Data());
            
            // For PbPb and pPb MC, loop over subevent types
            if(fSystemAndEnergy.Contains("Pb MC")){
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
                
                // Take into account that one subevent index is missing from the older files
                subeventIndex = iSubevent;
                if(legacyEnergyEnergyCorrelatorMode && iSubevent > EECHistograms::kPythiaHydjet) subeventIndex = iSubevent - 1;
                
                // Load the energy-energy correlator histograms with subevent binning
                histogramNamer = Form("%s%s%s_C%dT%dJ%dS%d", folderNamer.Data(), fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iJetPt, subeventIndex);
                fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());
                
              } // Subevent type loop
            } // Data is PbPb or pPb MC
            
          } // Jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // Pairing type loop

    if(preMixingCompatibilityMode){
      fPairingTypeSaveName[0] = currentPairingType0;
    }
    
    // Load also the processed energy-energy correlator histograms
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
          
          // Loop over the different processing steps
          for(int iProcessingLevel = 0; iProcessingLevel < knEnergyEnergyCorrelatorProcessingLevels; iProcessingLevel++){
            
            // Load processed energy-energt correlator histograms
            histogramNamer = Form("%sProcessed/%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[iProcessingLevel], iCentrality, iTrackPt, iJetPt);
            fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel] = (TH1D*) fInputFile->Get(histogramNamer.Data());
            
          } // Loop over different processing levels
        } // Loop over jet pT bins
        
        // =================================== //
        // Histograms without jet pT selection //
        // =================================== //
        
        // Loop over the different processing steps
        for(int iProcessingLevel = 0; iProcessingLevel < knEnergyEnergyCorrelatorProcessingLevels; iProcessingLevel++){
          
          // Load processed energy-energt correlator histograms
          histogramNamer = Form("%sProcessed/%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[iProcessingLevel], iCentrality, iTrackPt);
          fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iProcessingLevel] = (TH1D*) fInputFile->Get(histogramNamer.Data());
          
        } // Loop over different processing levels
        
      } // Loop over track pT bins
    } // Loop over centrality bins
    
  } // Energy-energy correlator type loop

  // Load the reflected cone QA histograms
  if(fLoadReflectedConeQAHistograms){

    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      histogramNamer = Form("%s/%s_C%d", fReflectedConeQAFolderName, fNumberOfJetsWithinReflectedConeName, iCentrality);
      fhNumberOfJetsWithinReflectedCone[iCentrality] = (TH1D*)fInputFile->Get(histogramNamer.Data());

      histogramNamer = Form("%s/%s_C%d", fReflectedConeQAFolderName, fJetPtWithinReflectedConeName, iCentrality);
      fhJetPtWithinReflectedCone[iCentrality] = (TH1D*)fInputFile->Get(histogramNamer.Data());

    }  // Centrality loop
  }
  
  // Load the jet pT response matrices from a processed file
  if(fLoadJetPtResponseMatrix){

    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){

      histogramNamer = Form("jetPtResponseMatrix/jetPtResponseMatrix_C%d", iCentrality);
      fhJetPtResponseMatrix[iCentrality] = (TH2D*)fInputFile->Get(histogramNamer.Data());

    }  // Centrality loop
  }

  // Load the jet pT closure histograms from a processed file
  if(fLoadJetPtClosureHistograms){
    
    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Loop over closure particles (quark/gluon/no selection)
      for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
        
        // Loop over generator level jet pT bins
        for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
          
          // Loop over jet eta bins
          for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){

            // Loop over jet phi bins
            for(int iJetPhi = 0; iJetPhi <= knJetPhiBins; iJetPhi++){
            
              // Legacy naming convention
              histogramNamer = Form("jetPtClosure_%s/jetPtClosure_%s%s_C%dT%dE%d", fJetHistogramName, fJetHistogramName, fClosureParticleName[iClosureParticle], iCentrality, iGenJetPt, iJetEta);
              fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iClosureParticle] = (TH1D*) fInputFile->Get(histogramNamer.Data());

              // If we do not find a histogram, try the updated naming convention
              if(fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iClosureParticle] == NULL){
                histogramNamer = Form("jetPtClosure_%s/jetPtClosure_%s%s_C%d", fJetHistogramName, fJetHistogramName, fClosureParticleName[iClosureParticle], iCentrality);
                if(iGenJetPt < knGenJetPtBins) histogramNamer.Append(Form("T%d",iGenJetPt));
                if(iJetEta < knJetEtaBins) histogramNamer.Append(Form("E%d",iJetEta));
                if(iJetPhi < knJetPhiBins) histogramNamer.Append(Form("P%d",iJetPhi));
                fhJetPtClosure[iGenJetPt][iJetEta][iJetPhi][iCentrality][iClosureParticle] = (TH1D*) fInputFile->Get(histogramNamer.Data());
              }
            
            } // Jet phi bin loop
          } // Jet eta bin loop
        } // Generator level jet pT loop
      } // Closure particle type (quark/gluon) loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Opening jet pT closure histograms

  // Load the jet pT unfolding histograms from a processed file
  if(fLoadJetPtUnfoldingHistograms){

    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
        for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
          histogramNamer = Form("%s/%s_C%dT%d", fJetPtUnfoldingDistributionName[iUnfoldType], fJetPtUnfoldingDistributionName[iUnfoldType], iCentrality, iTrackPt);
          fhJetPtUnfoldingDistribution[iUnfoldType][iCentrality][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());
        }
        histogramNamer = Form("%s/%s_C%dT%d", fJetPtResponseMatrixName, fJetPtResponseMatrixName, iCentrality, iTrackPt);
        fhJetPtUnfoldingResponse[iCentrality][iTrackPt] = (TH2D*) fInputFile->Get(histogramNamer.Data());
      } // Track pT loop
    } // Centrality loop


  } // Loading jet pT unfolding histograms

  // Load the covariance matrices used for jet pT unfolding
  if(fLoadJetPtUnfoldingCovariance){

    for(int iCoverianceMatrixType = 0; iCoverianceMatrixType < knCovarianceMatrixTypes; iCoverianceMatrixType++){
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
          histogramNamer = Form("%s/%s_C%dT%d", fJetPtCovarianceMatrixName[iCoverianceMatrixType], fJetPtCovarianceMatrixName[iCoverianceMatrixType], iCentrality, iTrackPt);
          fhJetPtUnfoldingCovariance[iCoverianceMatrixType][iCentrality][iTrackPt] = (TH2D*) fInputFile->Get(histogramNamer.Data());
        } // Track pT loop
      } // Centrality loop
    } // Covariance matrix type loop


  } // Loading jet pT unfolding histograms

  // Load the jet pT one dimensional unfolding histograms from a processed file
  if(fLoadJetPtOneDimensionalUnfoldingHistograms){

    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iUnfoldType = 0; iUnfoldType < knUnfoldingDistributionTypes; iUnfoldType++){
        histogramNamer = Form("%s/%s_C%d", fJetPtOneDimensionalUnfoldingDistributionName[iUnfoldType], fJetPtOneDimensionalUnfoldingDistributionName[iUnfoldType], iCentrality);
        fhOneDimensionalJetPtUnfoldingDistribution[iUnfoldType][iCentrality] = (TH1D*) fInputFile->Get(histogramNamer.Data());
      }
      histogramNamer = Form("%s/%s_C%d", fJetPtOneDimensionalResponseMatrixName, fJetPtOneDimensionalResponseMatrixName, iCentrality);
      fhOneDimensionalJetPtUnfoldingResponse[iCentrality] = (TH2D*) fInputFile->Get(histogramNamer.Data());
    } // Centrality loop


  } // Loading jet pT one dimensional unfolding histograms

  // Load the track/particle matching histograms from a processed file
  if(fLoadTrackParticleMatchingHistograms){
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){

          // Load the track/particle pair pT closure histograms
          histogramNamer = Form("%s/%s_C%dJ%dT%d", fTrackParticlePtClosureSaveName, fTrackParticlePtClosureSaveName, iCentrality, iJetPt, iTrackPt);
          fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());

          // Load the track/particle matching QA histograms
          for(int iTrackParticleQAHistogram = 0; iTrackParticleQAHistogram < knTrackParticleMatchingQAHistograms; iTrackParticleQAHistogram++){
            histogramNamer = Form("%s/%s_C%dJ%dT%d", fTrackParticleMatchingQAName[iTrackParticleQAHistogram], fTrackParticleMatchingQAName[iTrackParticleQAHistogram], iCentrality, iJetPt, iTrackPt);
            fhTrackParticleMatchQA[iTrackParticleQAHistogram][iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());
          }

          // Load the track/particle matching response matrices
          for(int iTrackParticleResponseHistogram = 0; iTrackParticleResponseHistogram < knTrackParticleMatchingResponseTypes; iTrackParticleResponseHistogram++){
            histogramNamer = Form("%s/%s_C%dJ%dT%d", fTrackParticleMatchingResponseName[iTrackParticleResponseHistogram], fTrackParticleMatchingResponseName[iTrackParticleResponseHistogram], iCentrality, iJetPt, iTrackPt);
            fhTrackParticleResponse[iTrackParticleResponseHistogram][iCentrality][iJetPt][iTrackPt] = (TH2D*) fInputFile->Get(histogramNamer.Data());
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Load track/particle matching histograms 
}

/*
 * Load the event histogram
 */
void EECHistogramManager::LoadEventsHistogram(){
  fhEvents = (TH1D*) fInputFile->Get("nEvents");
}

/*
 * Load the jet pT histogram
 */
void EECHistogramManager::LoadJetPtHistogram(const int iCentrality){
  fhJetPt[iCentrality] = (TH1D*) fInputFile->Get(Form("%s/%sPt_C%d", fJetHistogramName, fJetHistogramName, iCentrality));
}

/*
 * Check that the weight exponent requested is present in the input data.
 */
void EECHistogramManager::CheckWeightExponent(){

  // If we have only one weight exponent, we are golden and can just use that. Nothing needs to be done here.
  if(fCard->GetNWeightExponents() == 1) return;

  // We start by checking if the defined weight exponent exists in the card
  int weightExponentIndex = fCard->FindWeightExponentIndex(fLoadedWeightExponent);

  // If there are several weight exponents, we need to make sure that the one we are looking for exists in the card.
  // If it does not, we should throw an exception
  if(weightExponentIndex < 0){
    std::cout << "EECHistogramManager::ERROR! The weight exponent " << fLoadedWeightExponent << " is not present in the file " << fInputFile->GetName() << std::endl;
    std::cout << "The exponents available in this file are ";
    for(int iIndex = 1; iIndex <= fCard->GetNWeightExponents(); iIndex++){
      std::cout << fCard->GetWeightExponent(iIndex) << " ";
    }
    std::cout << std::endl;
    std::cout << "Enjoy your exception!" << std::endl;
    throw std::invalid_argument("Defined weight exponent in EECHistogramManager is not present in input file!");
  }

}

/*
 * If one weight exponent from many is projected, update the information in card
 */
void EECHistogramManager::UpdateWeightExponent(){

  // Start by finding the index for the weight exponent
  int weightExponentIndex = fCard->FindWeightExponentIndex(fLoadedWeightExponent);

  // In case index is 0, weight exponent is undefined in the card. This happens in older files before this was added to the card.
  // In this case we know that wegiht exponent must be 1. Add it to the card.
  if(weightExponentIndex == 0){
    fCard->AddOneDimensionalVector(EECCard::kWeightExponent, 1);
    return;
  }

  // If we have only one weight exponent, nothing needs to be updated
  if(fCard->GetNWeightExponents() == 1) return;

  // If the weight exponent exists, update the card such that we only leave single weight exponent value there
  fCard->AddOneDimensionalVector(EECCard::kWeightExponent, fLoadedWeightExponent);

}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   int* binIndices = Array of integers to be filled with bin index information read from the file
 *   const double* binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 */
void EECHistogramManager::SetBinIndices(const char* histogramName, const int nBins, int* binIndices, const double* binBorders, const int iAxis){
  THnSparseD* histogramArray = (THnSparseD*) fInputFile->Get(histogramName);
  TH1D* hBinner = FindHistogram(histogramArray,iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   double* copyBinBorders = Array to which a copy of bin borders is made
 *   int* binIndices = Array of integers to be filled with bin index information read from the file
 *   const double* binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 *   const bool setIndices = true: Set both bin indices and bin borders, false: Set only bin borders
 */
void EECHistogramManager::SetBinBordersAndIndices(const char* histogramName, const int nBins, double* copyBinBorders, int* binIndices, const double* binBorders, const int iAxis, const bool setIndices){
  TH1D* hBinner;
  THnSparseD* histogramArray;
  if(setIndices) {
    histogramArray = (THnSparseD*) fInputFile->Get(histogramName);
    hBinner = FindHistogram(histogramArray,iAxis,0,0,0);
  }
  for(int iBin = 0; iBin < nBins+1; iBin++){
    copyBinBorders[iBin] = binBorders[iBin];
    if(setIndices) binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Set up generic bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const char* histogramName = Name of the histogram from which the bin indices are searched
 *  const int iAxis = Axis from which the set indices can be found
 *  int nSetBins = Number of bins that is set
 *  double* setBinBorders = Bin borders that are set
 *  int* setBinIndices = Bin indices that are set
 *  const int nBins = New number of bins that is given
 *  const double* binBorders = New bin borders that are given
 *  const char* errorMessage = Type of the set bins to be printed in possible error message
 *  const int maxBins = Maximum number of allowed bins of this type
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double* binBorders, const char* errorMessage, const int maxBins, const bool setIndices){
  
  // If bins are read from file, do not use the given bin borders
  if(readBinsFromFile){
    if(setIndices) SetBinIndices(histogramName, nSetBins, setBinIndices, setBinBorders, iAxis);
  } else { // If the given bin borders are use, update the number of bins and bin borders according to input
    if(nBins <= maxBins){
      nSetBins = nBins;
      SetBinBordersAndIndices(histogramName, nSetBins, setBinBorders, setBinIndices, binBorders, iAxis, setIndices);
    } else {
      cout << "Error! Too many " << errorMessage << " bins given. Maximum number is " << maxBins << ". Will not set bins." << endl;
    }
  }
}

/*
 * Set up centrality bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given centrality bins
 *  const double* binBorders = New bin borders for centrality
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetCentralityBins(const bool readBinsFromFile, const int nBins, const double* binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, fEnergyEnergyCorrelatorHistogramNames[kEnergyEnergyCorrelator], 3, fnCentralityBins, fCentralityBinBorders, fCentralityBinIndices, nBins, binBorders, "centrality", kMaxCentralityBins, setIndices);
  
}

/*
 * Set up track pT bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double* binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetTrackPtBins(const bool readBinsFromFile, const int nBins, const double* binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, "track", 0, fnTrackPtBins, fTrackPtBinBorders, fTrackPtBinIndices, nBins, binBorders, "track pT", kMaxTrackPtBins, setIndices);
  
}

/*
 * Set up jet pT bin indices for energy-energy correlator according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double* binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetJetPtBinsEEC(const bool readBinsFromFile, const int nBins, const double* binBorders, bool setIndices){
  
  SetGenericBins(readBinsFromFile, fEnergyEnergyCorrelatorHistogramNames[kEnergyEnergyCorrelator], 1, fnJetPtBinsEEC, fJetPtBinBordersEEC, fJetPtIndicesEEC, nBins, binBorders, "jet pT EEC", kMaxJetPtBinsEEC, setIndices);
  
}

/*
 * Set up track pT bin indices for energy-e nergy correlator according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double* binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetTrackPtBinsEEC(const bool readBinsFromFile, const int nBins, const double* binBorders, bool setIndices ){
 
  SetGenericBins(readBinsFromFile, fEnergyEnergyCorrelatorHistogramNames[kEnergyEnergyCorrelator], 2, fnTrackPtBinsEEC, fTrackPtBinBordersEEC, fTrackPtIndicesEEC, nBins, binBorders, "track pT EEC", kMaxTrackPtBinsEEC, setIndices);
  
}

// Setter for loading event information
void EECHistogramManager::SetLoadEventInformation(const bool loadOrNot){
  fLoadEventInformation = loadOrNot;
}

// Setter for loading jet histograms
void EECHistogramManager::SetLoadJetHistograms(const bool loadOrNot){
  fLoadJets = loadOrNot;
}

// Setter for loading tracks
void EECHistogramManager::SetLoadTracks(const bool loadOrNot){
  fLoadTracks[kTrack] = loadOrNot;
}

// Setter for loading uncorrected tracks
void EECHistogramManager::SetLoadTracksUncorrected(const bool loadOrNot){
  fLoadTracks[kUncorrectedTrack] = loadOrNot;
}

// Setter for loading all track histograms
void EECHistogramManager::SetLoadAllTracks(const bool loadTracks, const bool loadUncorrected){
  SetLoadTracks(loadTracks);
  SetLoadTracksUncorrected(loadUncorrected);
}

// Setter for loading track multiplicity within the jet cone
void EECHistogramManager::SetLoadMultiplicityInJets(const bool loadOrNot){
  fLoadMultiplicityInJetHistograms = loadOrNot;
}

// Setter for loading particle density histograms around the jet axis
void EECHistogramManager::SetLoadParticleDensityAroundJets(const bool loadOrNot){
  fLoadParticleDensityAroundJetsHistograms[kParticleDensityAroundJetAxis] = loadOrNot;
}

// Setter for loading particle pT density histograms around the jet axis
void EECHistogramManager::SetLoadParticlePtDensityAroundJets(const bool loadOrNot){
  fLoadParticleDensityAroundJetsHistograms[kParticlePtDensityAroundJetAxis] = loadOrNot;
}

// Setter for loading pT binned particle density histograms around the jet axis
void EECHistogramManager::SetLoadParticleDensityAroundJetsPtBinned(const bool loadOrNot){
  fLoadParticleDensityAroundJetsHistograms[kParticleDensityAroundJetAxisPtBinned] = loadOrNot;
}

// Setter for loading pT binned particle pT density histograms around the jet axis
void EECHistogramManager::SetLoadParticlePtDensityAroundJetsPtBinned(const bool loadOrNot){
  fLoadParticleDensityAroundJetsHistograms[kParticlePtDensityAroundJetAxisPtBinned] = loadOrNot;
}

// Setter for loading particle pT density histograms around the jet axis
void EECHistogramManager::SetLoadAllParticleDensitiesAroundJets(const bool loadRegular, const bool loadPtWeighted, const bool loadPtBinned, const bool loadPtBinnedPtWeighted){
  SetLoadParticleDensityAroundJets(loadRegular);
  SetLoadParticlePtDensityAroundJets(loadPtWeighted);
  SetLoadParticleDensityAroundJetsPtBinned(loadPtBinned);
  SetLoadParticlePtDensityAroundJetsPtBinned(loadPtBinnedPtWeighted);
}

// Setter for loading maximum particle pT within the jet cone histograms
void EECHistogramManager::SetLoadMaxParticlePtWithinJetCone(const bool loadOrNot){
  fLoadMaxParticlePtWithinJetConeHistograms = loadOrNot;
}

// Setter for loading energy-energy correlators
void EECHistogramManager::SetLoadEnergyEnergyCorrelators(const bool loadOrNot){
  fLoadEnergyEnergyCorrelatorHistograms[kEnergyEnergyCorrelator] = loadOrNot;
}

// Setter for loading energy-energy correlators with positive track efficiency variation
void EECHistogramManager::SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(const bool loadOrNot){
  fLoadEnergyEnergyCorrelatorHistograms[kEnergyEnergyCorrelatorEfficiencyVariationPlus] = loadOrNot;
}

// Setter for loading energy-energy correlators with negative track efficiency variation
void EECHistogramManager::SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(const bool loadOrNot){
  fLoadEnergyEnergyCorrelatorHistograms[kEnergyEnergyCorrelatorEfficiencyVariationMinus] = loadOrNot;
}

// Setter for loading energy-energy correlators with positive track pair efficiency variation 
void EECHistogramManager::SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(const bool loadOrNot){
  fLoadEnergyEnergyCorrelatorHistograms[kEnergyEnergyCorrelatorPairEfficiencyVariationPlus] = loadOrNot;
}

// Setter for loading energy-energy correlators with negative track pair efficiency variation  
void EECHistogramManager::SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(const bool loadOrNot){
  fLoadEnergyEnergyCorrelatorHistograms[kEnergyEnergyCorrelatorPairEfficiencyVariationMinus] = loadOrNot;
} 

// Setter for loading all energy-energy correlators
void EECHistogramManager::SetLoadAllEnergyEnergyCorrelators(const bool loadRegular, const bool loadEfficiencyVariationPlus, const bool loadEfficiencyVariationMinus, const bool loadPairEfficiencyVariationPlus, const bool loadPairEfficiencyVariationMinus){
  SetLoadEnergyEnergyCorrelators(loadRegular);
  SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(loadEfficiencyVariationPlus);
  SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(loadEfficiencyVariationMinus);
  SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(loadPairEfficiencyVariationPlus);
  SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(loadPairEfficiencyVariationMinus);
}

// Setter for loading reflected cone QA histograms
void EECHistogramManager::SetLoadReflectedConeQAHistograms(const bool loadOrNot){
  fLoadReflectedConeQAHistograms = loadOrNot;
}

// Setter for loading histograms needed in the jet pT unfolding study
void EECHistogramManager::SetLoadJetPtUnfoldingHistograms(const bool loadOrNot){
  fLoadJetPtUnfoldingHistograms = loadOrNot;
}

// Setter for loading the covariance histograms for input to jet pT unfolding
void EECHistogramManager::SetLoadJetPtUnfoldingCovariance(const bool loadOrNot){
  fLoadJetPtUnfoldingCovariance = loadOrNot;
}

// Setter for loading histograms needed in the one dimensional jet pT unfolding study
void EECHistogramManager::SetLoadJetPtOneDimensionalUnfoldingHistograms(const bool loadOrNot){
  fLoadJetPtOneDimensionalUnfoldingHistograms = loadOrNot;
}

// Setter for loading histograms needed in track/particle matching study
void EECHistogramManager::SetLoadTrackParticleMatchingHistograms(const bool loadOrNot){
  fLoadTrackParticleMatchingHistograms = loadOrNot;
}

// Setter for flavor selection for jets in jet histograms
void EECHistogramManager::SetJetFlavor(const int iFlavor){
  fJetFlavor = iFlavor;
}

// Define the weight exponent value that is searched from the file
void EECHistogramManager::SetLoadedWeightExponent(const double weightExponent){
  fLoadedWeightExponent = weightExponent;
}

 // Setter for loading two-dimensional histograms
void EECHistogramManager::SetLoad2DHistograms(const bool loadOrNot){
  fLoad2DHistograms = loadOrNot;
}

// Setter for loading jet pT closure histograms
void EECHistogramManager::SetLoadJetPtClosureHistograms(const bool loadOrNot){
  fLoadJetPtClosureHistograms = loadOrNot;
}

// Setter for loading jet pT response matrix
void EECHistogramManager::SetLoadJetPtResponseMatrix(const bool loadOrNot){
  fLoadJetPtResponseMatrix = loadOrNot;
}

// Setter for loaded centrality bins
void EECHistogramManager::SetCentralityBinRange(const int first, const int last){
  fFirstLoadedCentralityBin = first;
  fLastLoadedCentralityBin = last;
  
  // Sanity check for centrality bins
  BinSanityCheck(fnCentralityBins,fFirstLoadedCentralityBin,fLastLoadedCentralityBin);
}

// Setter for loaded track pT bins
void EECHistogramManager::SetTrackPtBinRange(const int first, const int last){
  fFirstLoadedTrackPtBin = first;
  fLastLoadedTrackPtBin = last;
  
  // Sanity check for track pT bins
  BinSanityCheck(fnTrackPtBins,fFirstLoadedTrackPtBin,fLastLoadedTrackPtBin);
}

// Setter for jet pT bin range in energy-energy correlator histograms
void EECHistogramManager::SetJetPtBinRangeEEC(const int first, const int last){
  fFirstLoadedJetPtBinEEC = first;
  fLastLoadedJetPtBinEEC = last;
  
  // Sanity check for jet pT bins in energy-energy correlator histograms
  BinSanityCheck(fnJetPtBinsEEC,fFirstLoadedJetPtBinEEC,fLastLoadedJetPtBinEEC);
}

// Setter for track pT bin range in energy-energy correlator histograms
void EECHistogramManager::SetTrackPtBinRangeEEC(const int first, const int last){
  fFirstLoadedTrackPtBinEEC = first;
  fLastLoadedTrackPtBinEEC = last;
  
  // Sanity check for track pT bins in energy-energy correlator histograms
  BinSanityCheck(fnTrackPtBinsEEC,fFirstLoadedTrackPtBinEEC,fLastLoadedTrackPtBinEEC);
}

// Unfolding is done in a separate macro. Thus provide setter for unfolded energy-energy correlators so they can be stored in the histogram manager
void EECHistogramManager::SetUnfoldedEnergyEnergyCorrelator(const TH1D* unfoldedEnergyEnergyCorrelator, const int iEnergyEnergyCorrelatorType, const int iCentrality, const int iJetPt, const int iTrackPt){

  // Do a sanity check for the input bin indices
  if(iEnergyEnergyCorrelatorType < 0 || iEnergyEnergyCorrelatorType >= knEnergyEnergyCorrelatorTypes){
    cout << "EECHistogramManager::ERROR in SetUnfoldedEnergyEnergyCorrelator" << endl;
    cout << "Energy-energy correlator type index " << iEnergyEnergyCorrelatorType << " is out of range 0-" << knEnergyEnergyCorrelatorTypes-1 << "!" << endl;
    cout << "Cannot set the unfolded energy-energy correlator. Please check your code." << endl;
    return;
  }

  if(iCentrality < 0 || iCentrality >= kMaxCentralityBins){
    cout << "EECHistogramManager::ERROR in SetUnfoldedEnergyEnergyCorrelator" << endl;
    cout << "Centrality index " << iCentrality << " is out of range 0-" << kMaxCentralityBins-1 << "!" << endl;
    cout << "Cannot set the unfolded energy-energy correlator. Please check your code." << endl;
    return;
  }

  if(iJetPt < 0 || iJetPt >= kMaxJetPtBinsEEC){
    cout << "EECHistogramManager::ERROR in SetUnfoldedEnergyEnergyCorrelator" << endl;
    cout << "Jet pT index " << iJetPt << " is out of range 0-" << kMaxJetPtBinsEEC-1 << "!" << endl;
    cout << "Cannot set the unfolded energy-energy correlator. Please check your code." << endl;
    return;
  }

  if(iTrackPt < 0 || iTrackPt >= kMaxTrackPtBinsEEC){
    cout << "EECHistogramManager::ERROR in SetUnfoldedEnergyEnergyCorrelator" << endl;
    cout << "Track pT index " << iTrackPt << " is out of range 0-" << kMaxTrackPtBinsEEC-1 << "!" << endl;
    cout << "Cannot set the unfolded energy-energy correlator. Please check your code." << endl;
    return;
  }

  // If we are not out of bounds from the histogram array dimensions, copy the energy-energy correlator histogram to the array
  fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][kEnergyEnergyCorrelatorUnfolded] = (TH1D*) unfoldedEnergyEnergyCorrelator->Clone(Form("%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[kEnergyEnergyCorrelatorUnfolded], iCentrality, iTrackPt, iJetPt));

}

// Unfolding is done in a separate macro. Thus provide setter for unfolded covariance matrices so they can be stored in the histogram manager
void EECHistogramManager::SetUnfoldedCoverianceMatrix(const TH2D* unfoldedCovarianceMatrix, const int iCentrality, const int iTrackPt){

  // Do a sanity check for the input bin indices
  if(iCentrality < 0 || iCentrality >= kMaxCentralityBins){
    cout << "EECHistogramManager::ERROR in SetUnfoldedCoverianceMatrix" << endl;
    cout << "Centrality index " << iCentrality << " is out of range 0-" << kMaxCentralityBins-1 << "!" << endl;
    cout << "Cannot set the unfolded covariance matrix. Please check your code." << endl;
    return;
  }

  if(iTrackPt < 0 || iTrackPt >= kMaxTrackPtBinsEEC){
    cout << "EECHistogramManager::ERROR in SetUnfoldedCoverianceMatrix" << endl;
    cout << "Track pT index " << iTrackPt << " is out of range 0-" << kMaxTrackPtBinsEEC-1 << "!" << endl;
    cout << "Cannot set the unfolded covariance matrix. Please check your code." << endl;
    return;
  }

  // If we are not out of bounds from the histogram array dimensions, copy the energy-energy correlator histogram to the array
  fhJetPtUnfoldingCovariance[kCovarianceMatrixUnfolded][iCentrality][iTrackPt] = (TH2D*) unfoldedCovarianceMatrix->Clone(Form("%s_C%dT%d", fJetPtCovarianceMatrixName[kCovarianceMatrixUnfolded], iCentrality, iTrackPt));

}

// Sanity check for set bins
void EECHistogramManager::BinSanityCheck(const int nBins, int& first, int& last){
  if(first < 0) first = 0;
  if(last < first) last = first;
  if(last > nBins-1) last = nBins-1;
}

// Sanity check for input bin index
int EECHistogramManager::BinIndexCheck(const int nBins, const int binIndex) const{
  if(binIndex < 0) return 0;
  if(binIndex > nBins-1) return nBins-1;
  return binIndex;
}

// Getter for the number of centrality bins
int EECHistogramManager::GetNCentralityBins() const{
  return fnCentralityBins;
}

// Getter for the number of track pT bins
int EECHistogramManager::GetNTrackPtBins() const{
  return fnTrackPtBins;
}

// Getter for the number of jet pT bins in energy-energy correlator histograms
int EECHistogramManager::GetNJetPtBinsEEC() const{
  return fnJetPtBinsEEC;
}

// Getter for the number of track pT bins in energy-energy correlator histograms
int EECHistogramManager::GetNTrackPtBinsEEC() const{
  return fnTrackPtBinsEEC;
}

// Getter for track histogram name
const char* EECHistogramManager::GetTrackHistogramName(int iTrackType) const{
  iTrackType = BinIndexCheck(knTrackCategories,iTrackType);
  return fTrackHistogramNames[iTrackType];
}

// Getter for name suitable for x-axis in a given track histogram
const char* EECHistogramManager::GetTrackAxisName(int iTrackType) const{
  iTrackType = BinIndexCheck(knTrackCategories,iTrackType);
  return fTrackAxisNames[iTrackType];
}

// Getter for the jet histogram name
const char* EECHistogramManager::GetJetHistogramName() const{
  return fJetHistogramName;
}

// Getter for name suitable for x-axis in a given jet histogram
const char* EECHistogramManager::GetJetAxisName() const{
  return fJetAxisName;
}

// Getter for multiplicity in jet cone histogram name
const char* EECHistogramManager::GetMultiplicityInJetConeHistogramName(int iMultiplicityType) const{
  return fMultiplicityInJetsHistogramNames[iMultiplicityType];
}

// Getter for multiplicity in jet cone axis name
const char* EECHistogramManager::GetMultiplicityInJetConeAxisName(int iMultiplicityType) const{
  return fMultiplicityInJetsAxisNames[iMultiplicityType];
}

// Getter for the particle density around jet axis histogram name
const char* EECHistogramManager::GetParticleDensityAroundJetAxisHistogramName(int iParticleDensityType) const{
  return fParticleDensityAroundJetsHistogramNames[iParticleDensityType];
}

// Getter for the particle density around jet axis save name
const char* EECHistogramManager::GetParticleDensityAroundJetAxisSaveName(int iParticleDensityType) const{
  return fParticleDensityAroundJetsSaveNames[iParticleDensityType];
}

// Getter for the particle density around jet axis axis name
const char* EECHistogramManager::GetParticleDensityAroundJetAxisAxisName(int iParticleDensityType) const{
  return fParticleDensityAroundJetsAxisNames[iParticleDensityType];
}

// Getter for energy-energy correlator histogram name
const char* EECHistogramManager::GetEnergyEnergyCorrelatorHistogramName(int iEnergyEnergyCorrelatorType) const{
  iEnergyEnergyCorrelatorType = BinIndexCheck(knEnergyEnergyCorrelatorTypes, iEnergyEnergyCorrelatorType);
  return fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType];
}

// Getter for energy-energy correlator axis name
const char* EECHistogramManager::GetEnergyEnergyCorrelatorAxisName(int iEnergyEnergyCorrelatorType) const{
  iEnergyEnergyCorrelatorType = BinIndexCheck(knEnergyEnergyCorrelatorTypes, iEnergyEnergyCorrelatorType);
  return fEnergyEnergyCorrelatorAxisNames[iEnergyEnergyCorrelatorType];
}

// Getter for energy-energy correlator processing save name
const char* EECHistogramManager::GetEnergyEnergyCorrelatorProcessSaveName(int iProcessingLevel) const{
  iProcessingLevel = BinIndexCheck(knEnergyEnergyCorrelatorProcessingLevels, iProcessingLevel);
  return fEnergyEnergyCorrelatorProcessedSaveString[iProcessingLevel];
}

// Getter for subevent type
const char* EECHistogramManager::GetSubeventType(const int iSubeventType) const{
  if(iSubeventType < 0) return "";
  if(iSubeventType >= EECHistograms::knSubeventTypes) return "";
  return fSubeventTypeName[iSubeventType];
}

// Getter for subevent combination
const char* EECHistogramManager::GetSubeventCombination(const int iSubeventType) const{
  if(iSubeventType < 0) return "All combinations";
  if(iSubeventType >= EECHistograms::knSubeventCombinations) return "All combinations";
  return fSubeventCombinationName[iSubeventType];
}

// Getter for a well thought save name for subevent combination
TString EECHistogramManager::GetSubeventCombinationSaveName(const int iSubeventType) const{
  if(iSubeventType < 0) return "";
  if(iSubeventType >= EECHistograms::knSubeventCombinations) return "";
  TString niceSaveName = fSubeventCombinationName[iSubeventType];
  niceSaveName.ReplaceAll("-","");
  return niceSaveName;
}

// Getter for pairing type save names
const char* EECHistogramManager::GetPairingTypeSaveName(const int iPairingType) const{
  if(iPairingType < 0) return "NonsensicalIndex";
  if(iPairingType >= EECHistograms::knPairingTypes) return "NonsensicalIndex";
  return fPairingTypeSaveName[iPairingType];
}

// Getter for jet cone type save name
const char* EECHistogramManager::GetJetConeTypeSaveName(const int iJetConeType) const{
  return fJetConeTypeSaveName[iJetConeType];
}

// Getter for maximum particle pT within jet cone save name
const char* EECHistogramManager::GetMaxParticlePtWithinJetConeSaveName(const int iMaxParticlePtWithinJetConeType) const{
  return fMaxParticlePtInJetConeSaveName[iMaxParticlePtWithinJetConeType];
}

// Getter for maximum particle pT within jet cone axis name
const char* EECHistogramManager::GetMaxParticlePtWithinJetConeAxisName(const int iMaxParticlePtWithinJetConeType) const{
  return fMaxParticlePtInJetConeAxisName[iMaxParticlePtWithinJetConeType];
}

// Getter for the save name for covariance matrices
const char* EECHistogramManager::GetJetPtUnfoldingCovarianceSaveName(const int iCovarianceType) const{
  return fJetPtCovarianceMatrixName[iCovarianceType];
}

// Getter for collision system
TString EECHistogramManager::GetSystem() const{
  return fCard->GetDataType();
}

// Getter for i:th centrality bin border
double EECHistogramManager::GetCentralityBinBorder(const int iCentrality) const{
  return fCentralityBinBorders[iCentrality];
}

// Getter for i:th track pT bin border
double EECHistogramManager::GetTrackPtBinBorder(const int iTrackPt) const{
  return fTrackPtBinBorders[iTrackPt];
}

// Getter for i:th jet pT bin border in energy-energy correlator histograms
double EECHistogramManager::GetJetPtBinBorderEEC(const int iJetPt) const{
  return fJetPtBinBordersEEC[iJetPt];
}

// Getter for i:th track pT bin border in energy-energy correlator histograms
double EECHistogramManager::GetTrackPtBinBorderEEC(const int iTrackPt) const{
  return fTrackPtBinBordersEEC[iTrackPt];
}

// Getter for i:th track pT bin border in projections for maximum particle pT within the jet cone
double EECHistogramManager::GetMaxTrackPtWithinJetConeBinBorder(const int iTrackPt) const{
  return fProjectedMaxParticlePtBinBorders[iTrackPt];
}

// Getters for event information histograms

// Getter for z-vertex histogram
TH1D* EECHistogramManager::GetHistogramVertexZ(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhVertexZ == NULL) fhVertexZ = (TH1D*) fInputFile->Get("vertexZ"); 

  return fhVertexZ;
}

// Getter for z-vertex histogram
TH1D* EECHistogramManager::GetHistogramVertexZWeighted(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhVertexZWeighted == NULL) fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted"); 

  return fhVertexZWeighted;
}

// Getter for histogram for number of events surviving different event cuts
TH1D* EECHistogramManager::GetHistogramEvents(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhEvents == NULL) LoadEventsHistogram();

  return fhEvents;
}

// Getter for histogram for trigger selection
TH1D* EECHistogramManager::GetHistogramTriggers(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhTriggers == NULL) fhTriggers = (TH1D*) fInputFile->Get("triggers"); 

  return fhTriggers;
}

// Getter for histogram for number of tracks surviving different track cuts
TH1D* EECHistogramManager::GetHistogramTrackCuts(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhTrackCuts == NULL) fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts"); 

  return fhTrackCuts;
}

// Getter for centrality histogram in all events
TH1D* EECHistogramManager::GetHistogramCentrality(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhCentrality == NULL) fhCentrality = (TH1D*) fInputFile->Get("centrality"); 

  return fhCentrality;
}

// Getter for weighted centrality histogram in all events
TH1D* EECHistogramManager::GetHistogramCentralityWeighted(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhCentralityWeighted == NULL) fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted"); 

  return fhCentralityWeighted;
}

// Getter for pT hat histogram
TH1D* EECHistogramManager::GetHistogramPtHat(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhPtHat == NULL) fhPtHat = (TH1D*) fInputFile->Get("pthat");

  return fhPtHat;
}

// Getter for pT hat histogram
TH1D* EECHistogramManager::GetHistogramPtHatWeighted(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhPtHatWeighted == NULL) fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted"); 

  return fhPtHatWeighted;
}

// Getter for multiplicity histogram in all events
TH1D* EECHistogramManager::GetHistogramMultiplicity(int iCentrality){
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp

  // If the histogram is NULL, try to load the processed version of it
  if(fhMultiplicity[iCentrality] == NULL){
    TString histogramNamer = "multiplicity/multiplicity";
    if(iCentrality != fnCentralityBins) histogramNamer.Append(Form("_C%d", iCentrality));
    fhMultiplicity[iCentrality] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhMultiplicity[iCentrality];
}

// Getter for track efficiency weighted multiplicity histogram in all events
TH1D* EECHistogramManager::GetHistogramMultiplicityWeighted(int iCentrality){
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp

  // If the histogram is NULL, try to load the processed version of it
  if(fhMultiplicityWeighted[iCentrality] == NULL){
    TString histogramNamer = "multiplicity/multiplicityWeighted";
    if(iCentrality != fnCentralityBins) histogramNamer.Append(Form("_C%d", iCentrality));
    fhMultiplicity[iCentrality] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhMultiplicityWeighted[iCentrality];
}

// Getter for multiplicity vs. centrality map
TH2D* EECHistogramManager::GetHistogramMultiplicityMap(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhMultiplicityMap == NULL) fhMultiplicityMap = (TH2D*) fInputFile->Get("multiplicity/multiplicityMap");

  return fhMultiplicityMap;
}

// Getter for weighted multiplicity vs. centrality map
TH2D* EECHistogramManager::GetHistogramWeightedMultiplicityMap(){

  // If the histogram is NULL, try to load the processed version of it
  if(fhMultiplicityMapWeighted == NULL) fhMultiplicityMapWeighted = (TH2D*) fInputFile->Get("multiplicity/multiplicityMapWeighted");

  return fhMultiplicityMapWeighted;
}

// Getter for energy in HF plus calorimeters
TH1D* EECHistogramManager::GetHistogramHFPlus(int iJetPt, int nJets){

  // If the jet pT index is out of bounds, give the jet pT integrated histogram
  if(iJetPt < 0 || iJetPt > fnJetPtBinsEEC) iJetPt = fnJetPtBinsEEC;

  // If the nJets index in out of bounds, give the histogram without selection on number of jets
  if(nJets < 0 || nJets > knMaxJetsBins) nJets = knMaxJetsBins;

  // If the histogram is NULL, try to load the processed version of it
  if(fhHFPlus[iJetPt][nJets] == NULL) {
    TString histogramNamer = "HFEnergy/HFPlus";
    if(iJetPt != fnJetPtBinsEEC) histogramNamer.Append(Form("_J%d", iJetPt));
    if(nJets != knMaxJetsBins) histogramNamer.Append(Form("_N%d", nJets));
    fhHFPlus[iJetPt][nJets] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhHFPlus[iJetPt][nJets];
}

// Getter for energy in HF minus calorimeters          
TH1D* EECHistogramManager::GetHistogramHFMinus(int iJetPt, int nJets){

  // If the jet pT index is out of bounds, give the jet pT integrated histogram
  if(iJetPt < 0 || iJetPt > fnJetPtBinsEEC) iJetPt = fnJetPtBinsEEC;

  // If the nJets index in out of bounds, give the histogram without selection on number of jets
  if(nJets < 0 || nJets > knMaxJetsBins) nJets = knMaxJetsBins;

  // If the histogram is NULL, try to load the processed version of it
  if(fhHFMinus[iJetPt][nJets] == NULL) {
    TString histogramNamer = "HFEnergy/HFMinus";
    if(iJetPt != fnJetPtBinsEEC) histogramNamer.Append(Form("_J%d", iJetPt));
    if(nJets != knMaxJetsBins) histogramNamer.Append(Form("_N%d", nJets));
    fhHFMinus[iJetPt][nJets] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhHFMinus[iJetPt][nJets];
}

// Getter for energy in HF calorimeters            
TH1D* EECHistogramManager::GetHistogramHFSum(int iJetPt, int nJets){

  // If the jet pT index is out of bounds, give the jet pT integrated histogram
  if(iJetPt < 0 || iJetPt > fnJetPtBinsEEC) iJetPt = fnJetPtBinsEEC;

  // If the nJets index in out of bounds, give the histogram without selection on number of jets
  if(nJets < 0 || nJets > knMaxJetsBins) nJets = knMaxJetsBins;

  // If the histogram is NULL, try to load the processed version of it
  if(fhHFSum[iJetPt][nJets] == NULL) {
    TString histogramNamer = "HFEnergy/HFSum";
    if(iJetPt != fnJetPtBinsEEC) histogramNamer.Append(Form("_J%d", iJetPt));
    if(nJets != knMaxJetsBins) histogramNamer.Append(Form("_N%d", nJets));
    fhHFSum[iJetPt][nJets] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhHFSum[iJetPt][nJets];
}

// Getter for energy in HF plus calorimeters
TH2D* EECHistogramManager::GetHistogramHFPlusVsHFMinus(int iJetPt, int nJets){

  // If the jet pT index is out of bounds, give the jet pT integrated histogram
  if(iJetPt < 0 || iJetPt > fnJetPtBinsEEC) iJetPt = fnJetPtBinsEEC;

  // If the nJets index in out of bounds, give the histogram without selection on number of jets
  if(nJets < 0 || nJets > knMaxJetsBins) nJets = knMaxJetsBins;

  // If the histogram is NULL, try to load the processed version of it
  if(fhHFPlusVsHFMinus[iJetPt][nJets] == NULL) {
    TString histogramNamer = "HFEnergy/HFPlusVsHFMinus";
    if(iJetPt != fnJetPtBinsEEC) histogramNamer.Append(Form("_J%d", iJetPt));
    if(nJets != knMaxJetsBins) histogramNamer.Append(Form("_N%d", nJets));
    fhHFPlusVsHFMinus[iJetPt][nJets] = (TH2D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhHFPlusVsHFMinus[iJetPt][nJets];
}

// Getter for energy in HF plus calorimeters
TH1D* EECHistogramManager::GetHistogramHFPlus(std::pair<double,double> jetPtBin, int nJets){
  return GetHistogramHFPlus(fCard->FindBinIndexJetPtEEC(jetPtBin), nJets);
} 

 // Getter for energy in HF minus calorimeters
TH1D* EECHistogramManager::GetHistogramHFMinus(std::pair<double,double> jetPtBin, int nJets){
  return GetHistogramHFMinus(fCard->FindBinIndexJetPtEEC(jetPtBin), nJets);
}

// Getter for energy in HF calorimeters    
TH1D* EECHistogramManager::GetHistogramHFSum(std::pair<double,double> jetPtBin, int nJets){
  return GetHistogramHFSum(fCard->FindBinIndexJetPtEEC(jetPtBin), nJets);
}

// Getter for 2D map between HF plus and HF minus
TH2D* EECHistogramManager::GetHistogramHFPlusVsHFMinus(std::pair<double,double> jetPtBin, int nJets){
  return GetHistogramHFPlusVsHFMinus(fCard->FindBinIndexJetPtEEC(jetPtBin), nJets);
}

// Getters for jet histograms

// Getter for jet pT histograms
TH1D* EECHistogramManager::GetHistogramJetPt(int iCentrality){
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp

  // If the histogram is NULL, try to load the processed version of it
  if(fhJetPt[iCentrality] == NULL) LoadJetPtHistogram(iCentrality);

  return fhJetPt[iCentrality];
}

// Getter for jet phi histograms
TH1D* EECHistogramManager::GetHistogramJetPhi(int iCentrality){
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp

  // If the histogram is NULL, try to load the processed version of it
  if(fhJetPhi[iCentrality] == NULL) fhJetPhi[iCentrality] = (TH1D*) fInputFile->Get(Form("%s/%sPhi_C%d", fJetHistogramName, fJetHistogramName, iCentrality));

  return fhJetPhi[iCentrality];
}

// Getter for jet eta histograms
TH1D* EECHistogramManager::GetHistogramJetEta(int iCentrality){
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp

  // If the histogram is NULL, try to load the processed version of it
  if(fhJetEta[iCentrality] == NULL) fhJetEta[iCentrality] = (TH1D*) fInputFile->Get(Form("%s/%sEta_C%d", fJetHistogramName, fJetHistogramName, iCentrality));

  return fhJetEta[iCentrality];
}

// Getter for jet eta_CM histograms
TH1D* EECHistogramManager::GetHistogramJetEtaCM(int iCentrality){
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp

  // If the histogram is NULL, try to load the processed version of it
  if(fhJetEtaCM[iCentrality] == NULL) fhJetEtaCM[iCentrality] = (TH1D*) fInputFile->Get(Form("%s/%sEtaCM_C%d", fJetHistogramName, fJetHistogramName, iCentrality));

  return fhJetEtaCM[iCentrality];
}

// Getter for 2D eta-phi histogram for jets
TH2D* EECHistogramManager::GetHistogramJetEtaPhi(int iCentrality){
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp

  // If the histogram is NULL, try to load the processed version of it
  if(fhJetEtaPhi[iCentrality] == NULL) fhJetEtaPhi[iCentrality] = (TH2D*) fInputFile->Get(Form("%s/%sEtaPhi_C%d", fJetHistogramName, fJetHistogramName, iCentrality));

  return fhJetEtaPhi[iCentrality];
}

// Getters for histograms for tracks

// Getter for track pT histograms
TH1D* EECHistogramManager::GetHistogramTrackPt(const int iTrackType, const int iCentrality){

  // If the histogram is NULL, try to load the processed version of it
  if(fhTrackPt[iTrackType][iCentrality] == NULL) fhTrackPt[iTrackType][iCentrality] = (TH1D*) fInputFile->Get(Form("%s/%sPt_C%d",fTrackHistogramNames[iTrackType], fTrackHistogramNames[iTrackType], iCentrality));

  return fhTrackPt[iTrackType][iCentrality];
}

// Getter for track phi histograms
TH1D* EECHistogramManager::GetHistogramTrackPhi(const int iTrackType, const int iCentrality, const int iTrackPt){

  // If the histogram is NULL, try to load the processed version of it
  if(fhTrackPhi[iTrackType][iCentrality][iTrackPt] == NULL) fhTrackPhi[iTrackType][iCentrality][iTrackPt] = (TH1D*) fInputFile->Get(Form("%s/%sPhi_C%dT%d", fTrackHistogramNames[iTrackType], fTrackHistogramNames[iTrackType], iCentrality, iTrackPt));

  return fhTrackPhi[iTrackType][iCentrality][iTrackPt];
}

// Getter for track eta histograms
TH1D* EECHistogramManager::GetHistogramTrackEta(const int iTrackType, const int iCentrality, const int iTrackPt){

  // If the histogram is NULL, try to load the processed version of it
  if(fhTrackEta[iTrackType][iCentrality][iTrackPt] == NULL) fhTrackEta[iTrackType][iCentrality][iTrackPt] = (TH1D*) fInputFile->Get(Form("%s/%sEta_C%dT%d", fTrackHistogramNames[iTrackType], fTrackHistogramNames[iTrackType], iCentrality, iTrackPt));

  return fhTrackEta[iTrackType][iCentrality][iTrackPt];
}

// Getter for 2D eta-phi histogram for track
TH2D* EECHistogramManager::GetHistogramTrackEtaPhi(const int iTrackType, const int iCentrality, const int iTrackPt){

  // If the histogram is NULL, try to load the processed version of it
  if(fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] == NULL) fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] = (TH2D*) fInputFile->Get(Form("%s/%sEtaPhi_C%dT%d", fTrackHistogramNames[iTrackType], fTrackHistogramNames[iTrackType], iCentrality, iTrackPt));

  return fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt];
}

// Getter for multiplicity histogram within the jet cone
TH1D* EECHistogramManager::GetHistogramMultiplicityInJetCone(const int iCentrality, const int iJetPt, const int iTrackPt, const int iMultiplicityType, const int iSubevent){

  // If the histogram is NULL, try to load the processed version of it
  if(fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubevent] == NULL){
    TString histogramNamer = Form("%s/%s_C%dT%d", fMultiplicityInJetsHistogramNames[iMultiplicityType], fMultiplicityInJetsHistogramNames[iMultiplicityType], iCentrality, iTrackPt);
    if(iJetPt != fnJetPtBinsEEC) histogramNamer.Append(Form("J%d", iJetPt));
    if(iSubevent != EECHistograms::knSubeventCombinations) histogramNamer.Append(Form("S%d", iSubevent));
    fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhMultiplicityInJetCone[iCentrality][iJetPt][iTrackPt][iMultiplicityType][iSubevent];
}

// Getter for particle density histogram around the jet cone
TH1D* EECHistogramManager::GetHistogramParticleDensityAroundJetAxis(const int iCentrality, const int iJetPt, const int iTrackPt, const int iJetConeType, const int iParticleDensityType, const int iSubevent){

  // If the histogram is NULL, try to load the processed version of it
  if(fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubevent] == NULL){
    TString histogramNamer = Form("%s/%s%s_C%dT%d", fParticleDensityAroundJetsSaveNames[iParticleDensityType], fParticleDensityAroundJetsSaveNames[iParticleDensityType], fJetConeTypeSaveName[iJetConeType], iCentrality, iTrackPt);
    if(iJetPt != fnJetPtBinsEEC) histogramNamer.Append(Form("J%d", iJetPt));
    if(iSubevent != EECHistograms::knSubeventCombinations) histogramNamer.Append(Form("S%d", iSubevent));
    fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhParticleDensityAroundJetAxis[iCentrality][iJetPt][iTrackPt][iJetConeType][iParticleDensityType][iSubevent];
}

// Maximum particle pT in jet cone
TH1D* EECHistogramManager::GetHistogramMaxParticlePtInJetCone(const int iMaxParticlePtWithinJetConeType, const int iCentrality, const int iJetPt, const int iTrackPt){

  // If the histogram is NULL, try to load the processed version of it
  if(fhMaxParticlePtInJetConePtBin[iMaxParticlePtWithinJetConeType][iCentrality][iJetPt][iTrackPt] == NULL){
    TString histogramNamer = Form("%s/%s_C%dJ%d",fMaxParticlePtInJetConeSaveName[iMaxParticlePtWithinJetConeType], fMaxParticlePtInJetConeSaveName[iMaxParticlePtWithinJetConeType], iCentrality, iJetPt);
    if(iTrackPt != knProjectedMaxParticlePtBins) histogramNamer.Append(Form("T%d", iTrackPt));
    fhMaxParticlePtInJetConePtBin[iMaxParticlePtWithinJetConeType][iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhMaxParticlePtInJetConePtBin[iMaxParticlePtWithinJetConeType][iCentrality][iJetPt][iTrackPt];
}

// Maximum particle pT in jet cone with pT cut for background particles
TH1D* EECHistogramManager::GetHistogramMaxParticlePtInJetConePtCut(const int iMaxParticlePtWithinJetConeType, const int iCentrality, const int iJetPt, const int iTrackPt){

  // If the histogram is NULL, try to load the processed version of it
  if(fhMaxParticlePtInJetConePtCut[iMaxParticlePtWithinJetConeType][iCentrality][iJetPt][iTrackPt] == NULL){
    TString histogramNamer = Form("%s/%sTrackPtCut_C%dJ%dT%d",fMaxParticlePtInJetConeSaveName[iMaxParticlePtWithinJetConeType], fMaxParticlePtInJetConeSaveName[iMaxParticlePtWithinJetConeType], iCentrality, iJetPt, iTrackPt);
    fhMaxParticlePtInJetConePtCut[iMaxParticlePtWithinJetConeType][iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhMaxParticlePtInJetConePtCut[iMaxParticlePtWithinJetConeType][iCentrality][iJetPt][iTrackPt];
}

// Getter for maximum particle pT in jet cone
TH1D* EECHistogramManager::GetHistogramMaxSignalParticlePtInJetCone(const int iCentrality, const int iJetPt, const int iTrackPt){
  return GetHistogramMaxParticlePtInJetCone(kMaxSignalParticlePt,iCentrality,iJetPt,iTrackPt);
}

// Getter for maximum particle pT in jet cone with pT cut for background particles
TH1D* EECHistogramManager::GetHistogramMaxSignalParticlePtInJetConePtCut(const int iCentrality, const int iJetPt, const int iTrackPt){
  return GetHistogramMaxParticlePtInJetConePtCut(kMaxSignalParticlePt,iCentrality,iJetPt,iTrackPt);
}

// Getter for maximum background particle pT in jet cone
TH1D* EECHistogramManager::GetHistogramMaxBackgroundParticlePtInJetCone(const int iCentrality, const int iJetPt, const int iTrackPt){
  return GetHistogramMaxParticlePtInJetCone(kMaxBackgroundParticlePt,iCentrality,iJetPt,iTrackPt);
}

// Maximum background particle pT in jet cone with pT cut for signal particles
TH1D* EECHistogramManager::GetHistogramMaxBackgroundParticlePtInJetConePtCut(const int iCentrality, const int iJetPt, const int iTrackPt){
  return GetHistogramMaxParticlePtInJetConePtCut(kMaxBackgroundParticlePt,iCentrality,iJetPt,iTrackPt);
}

// Getter for energy-energy correlator histograms
TH1D* EECHistogramManager::GetHistogramEnergyEnergyCorrelator(const int iEnergyEnergyCorrelatorType, const int iCentrality, const int iJetPt, const int iTrackPt, const int iPairingType, const int iSubevent){

  // If the histogram is NULL, try to load it from the input file
  if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] == NULL){
    TString histogramNamer = Form("%s/%s/%s%s_C%dT%d", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt);

    if(iJetPt != fnJetPtBinsEEC) histogramNamer.Append(Form("J%d", iJetPt));
    if(iSubevent != EECHistograms::knSubeventCombinations) histogramNamer.Append(Form("S%d", iSubevent));

    fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());

  }

  // If the histogram is still NULL, we might be dealing with files that have old naming convention. Try that name instead.
  if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] == NULL){
    TString histogramNamer = Form("%s/%s", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);

    if(iPairingType != EECHistograms::kSameJetPair) histogramNamer.Append(fPairingTypeSaveName[iPairingType]);
    histogramNamer.Append(Form("_C%dT%d", iCentrality, iTrackPt));
    if(iJetPt != fnJetPtBinsEEC) histogramNamer.Append(Form("J%d", iJetPt));
    if(iSubevent != EECHistograms::knSubeventCombinations) histogramNamer.Append(Form("S%d", iSubevent));

    fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer.Data());

  }

  return fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent];
}

// Getter for energy-energy correlator histograms using bin boundaries instead of indices
TH1D* EECHistogramManager::GetHistogramEnergyEnergyCorrelator(const int iEnergyEnergyCorrelatorType, const std::pair<int,int> centralityBin, const std::pair<int,int> jetPtBin, const double trackPtCut, const int iPairingType, const int iSubevent){

  // PbPb and OO data/MC have centrality selections. pp and pPb use bin zero.
  int iCentrality = 0;
  const TString collisionSystem = fCard->GetDataType();
  if(collisionSystem.Contains("PbPb") || collisionSystem.Contains("OO")){
    iCentrality = fCard->FindBinIndexCentrality(centralityBin);
    if(iCentrality < 0) return nullptr;
  }

  // Always find the correct jet and track pT bins
  int iJetPt = fCard->FindBinIndexJetPtEEC(jetPtBin);
  int iTrackPt = fCard->GetBinIndexTrackPtEEC(trackPtCut);
  if(iJetPt < 0 || iTrackPt < 0) return nullptr;

  return GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelatorType, iCentrality, iJetPt, iTrackPt, iPairingType, iSubevent);

}

// Getter for processed energy-energy correlator histograms
TH1D* EECHistogramManager::GetHistogramEnergyEnergyCorrelatorProcessed(const int iEnergyEnergyCorrelatorType, const int iCentrality, const int iJetPt, const int iTrackPt, const int iProcessingLevel){

  // If the histogram is NULL, try to load it from the input file
  if(fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel] == NULL){
    TString histogramNamer = Form("%sProcessed/%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorProcessedSaveString[iProcessingLevel], iCentrality, iTrackPt);

    if(iJetPt != fnJetPtBinsEEC) histogramNamer.Append(Form("J%d", iJetPt));
    
    fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhEnergyEnergyCorrelatorProcessed[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iProcessingLevel];
}

// Getter for histograms showing number of jets above 25 GeV within the reflected cone
TH1D* EECHistogramManager::GetHistogramNumberOfJetsWithinReflectedCone(const int iCentrality){

  // If the histogram is NULL, try to load it from the input file
  if(fhNumberOfJetsWithinReflectedCone[iCentrality] == NULL) {
    fhNumberOfJetsWithinReflectedCone[iCentrality] = (TH1D*)fInputFile->Get(Form("%s/%s_C%d", fReflectedConeQAFolderName, fNumberOfJetsWithinReflectedConeName, iCentrality));
  }

  return fhNumberOfJetsWithinReflectedCone[iCentrality];
}

// Getter for pT of the jets that are found from the reflected cone
TH1D* EECHistogramManager::GetHistogramJetPtWithinReflectedCone(const int iCentrality){

  // If the histogram is NULL, try to load it from the input file
  if(fhJetPtWithinReflectedCone[iCentrality] == NULL){
    fhJetPtWithinReflectedCone[iCentrality] = (TH1D*)fInputFile->Get(Form("%s/%s_C%d", fReflectedConeQAFolderName, fJetPtWithinReflectedConeName, iCentrality));
  }
      
  return fhJetPtWithinReflectedCone[iCentrality];
}

// Getter for jet pT response matrix
TH2D* EECHistogramManager::GetHistogramJetPtResponseMatrix(const int iCentrality){

  // If the histogram is NULL, try to load it from the input file
  if(fhJetPtResponseMatrix[iCentrality] == NULL){
    fhJetPtResponseMatrix[iCentrality] = (TH2D*)fInputFile->Get(Form("jetPtResponseMatrix/jetPtResponseMatrix_C%d", iCentrality));
  }

  return fhJetPtResponseMatrix[iCentrality];
}

// Getter for jet pT closure histograms
TH1D* EECHistogramManager::GetHistogramJetPtClosure(const int iGenPtBin, const int iEtaBin, const int iPhiBin, const int iCentrality, const int iClosureParticle){

  // If the histogram is NULL, try to load it from the input file
  if(fhJetPtClosure[iGenPtBin][iEtaBin][iPhiBin][iCentrality][iClosureParticle] == NULL){
    TString histogramNamer = Form("jetPtClosure_%s/jetPtClosure_%s%s_C%d", fJetHistogramName, fJetHistogramName, fClosureParticleName[iClosureParticle], iCentrality);
    if(iGenPtBin < knGenJetPtBins) histogramNamer.Append(Form("T%d",iGenPtBin));
    if(iEtaBin < knJetEtaBins) histogramNamer.Append(Form("E%d",iEtaBin));
    if(iPhiBin < knJetPhiBins) histogramNamer.Append(Form("P%d",iPhiBin));
    fhJetPtClosure[iGenPtBin][iEtaBin][iPhiBin][iCentrality][iClosureParticle] = (TH1D*) fInputFile->Get(histogramNamer.Data());
  }

  return fhJetPtClosure[iGenPtBin][iEtaBin][iPhiBin][iCentrality][iClosureParticle];
}

// Getter for measured jet pT unfolding distribution
TH1D* EECHistogramManager::GetHistogramJetPtUnfoldingMeasured(const int iCentrality, const int iTrackPt){

  // If the histogram is NULL, try to load it from the input file
  if(fhJetPtUnfoldingDistribution[kUnfoldingMeasured][iCentrality][iTrackPt] == NULL){
    fhJetPtUnfoldingDistribution[kUnfoldingMeasured][iCentrality][iTrackPt] = (TH1D*) fInputFile->Get(Form("%s/%s_C%dT%d", fJetPtUnfoldingDistributionName[kUnfoldingMeasured], fJetPtUnfoldingDistributionName[kUnfoldingMeasured], iCentrality, iTrackPt));

  }

  return fhJetPtUnfoldingDistribution[kUnfoldingMeasured][iCentrality][iTrackPt];
}

// Getter for truth jet pT unfolding distribution
TH1D* EECHistogramManager::GetHistogramJetPtUnfoldingTruth(const int iCentrality, const int iTrackPt){

  // If the histogram is NULL, try to load it from the input file
  if(fhJetPtUnfoldingDistribution[kUnfoldingTruth][iCentrality][iTrackPt] == NULL){
    fhJetPtUnfoldingDistribution[kUnfoldingTruth][iCentrality][iTrackPt] = (TH1D*) fInputFile->Get(Form("%s/%s_C%dT%d", fJetPtUnfoldingDistributionName[kUnfoldingTruth], fJetPtUnfoldingDistributionName[kUnfoldingTruth], iCentrality, iTrackPt));

  }

  return fhJetPtUnfoldingDistribution[kUnfoldingTruth][iCentrality][iTrackPt];
}

// Getter for jet pT unfolding response
TH2D* EECHistogramManager::GetHistogramJetPtUnfoldingResponse(const int iCentrality, const int iTrackPt){

  // If the histogram is NULL, try to load it from the input file
  if(fhJetPtUnfoldingResponse[iCentrality][iTrackPt] == NULL){
    fhJetPtUnfoldingResponse[iCentrality][iTrackPt] = (TH2D*) fInputFile->Get(Form("%s/%s_C%dT%d", fJetPtResponseMatrixName, fJetPtResponseMatrixName, iCentrality, iTrackPt));
  }

  return fhJetPtUnfoldingResponse[iCentrality][iTrackPt];
}

// Getter for jet pT unfolding covariance
TH2D* EECHistogramManager::GetHistogramJetPtUnfoldingCovariance(const int iCovarianceMatrixType, const int iCentrality, const int iTrackPt){

  // If the histogram is NULL, try to load it from the input file
  if(fhJetPtUnfoldingCovariance[iCovarianceMatrixType][iCentrality][iTrackPt] == NULL){
    fhJetPtUnfoldingCovariance[iCovarianceMatrixType][iCentrality][iTrackPt] = (TH2D*) fInputFile->Get(Form("%s/%s_C%dT%d", fJetPtCovarianceMatrixName[iCovarianceMatrixType], fJetPtCovarianceMatrixName[iCovarianceMatrixType], iCentrality, iTrackPt));
  }

  return fhJetPtUnfoldingCovariance[iCovarianceMatrixType][iCentrality][iTrackPt];
}

// Getter for measured one dimensional jet pT unfolding distribution
TH1D* EECHistogramManager::GetHistogramJetPtOneDimensionalUnfoldingMeasured(const int iCentrality){

  // If the histogram is NULL, try to load it from the input file
  if(fhOneDimensionalJetPtUnfoldingDistribution[kUnfoldingMeasured][iCentrality] == NULL){
    fhOneDimensionalJetPtUnfoldingDistribution[kUnfoldingMeasured][iCentrality] = (TH1D*) fInputFile->Get(Form("%s/%s_C%d", fJetPtOneDimensionalUnfoldingDistributionName[kUnfoldingMeasured], fJetPtOneDimensionalUnfoldingDistributionName[kUnfoldingMeasured], iCentrality));
  }

  return fhOneDimensionalJetPtUnfoldingDistribution[kUnfoldingMeasured][iCentrality];
}

// Getter for truth one dimensional jet pT unfolding distribution
TH1D* EECHistogramManager::GetHistogramJetPtOneDimensionalUnfoldingTruth(const int iCentrality){

  // If the histogram is NULL, try to load it from the input file
  if(fhOneDimensionalJetPtUnfoldingDistribution[kUnfoldingTruth][iCentrality] == NULL){
    fhOneDimensionalJetPtUnfoldingDistribution[kUnfoldingTruth][iCentrality] = (TH1D*) fInputFile->Get(Form("%s/%s_C%d", fJetPtOneDimensionalUnfoldingDistributionName[kUnfoldingTruth], fJetPtOneDimensionalUnfoldingDistributionName[kUnfoldingTruth], iCentrality));
  }

  return fhOneDimensionalJetPtUnfoldingDistribution[kUnfoldingTruth][iCentrality];
}

// Getter for one dimensional jet pT unfolding response
TH2D* EECHistogramManager::GetHistogramJetPtOneDimensionalUnfoldingResponse(const int iCentrality){

  // If the histogram is NULL, try to load it from the input file
  if(fhOneDimensionalJetPtUnfoldingResponse[iCentrality] == NULL){
    fhOneDimensionalJetPtUnfoldingResponse[iCentrality] = (TH2D*) fInputFile->Get(Form("%s/%s_C%d", fJetPtOneDimensionalResponseMatrixName, fJetPtOneDimensionalResponseMatrixName, iCentrality));
  }

  return fhOneDimensionalJetPtUnfoldingResponse[iCentrality];
}


// Getter for number of particles close to tracks
TH1D* EECHistogramManager::GetHistogramParticlesCloseToTrack(const int iCentrality, const int iJetPt, const int iTrackPt){

  // If the histogram is NULL, try to load it from the input file
  if(fhTrackParticleMatchQA[kNumberOfParticlesCloseToTrack][iCentrality][iJetPt][iTrackPt] == NULL){
    fhTrackParticleMatchQA[kNumberOfParticlesCloseToTrack][iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(Form("%s/%s_C%dJ%dT%d", fTrackParticleMatchingQAName[kNumberOfParticlesCloseToTrack], fTrackParticleMatchingQAName[kNumberOfParticlesCloseToTrack], iCentrality, iJetPt, iTrackPt));
  }

  return fhTrackParticleMatchQA[kNumberOfParticlesCloseToTrack][iCentrality][iJetPt][iTrackPt];
}

// Getter for flag if a matching particle is found
TH1D* EECHistogramManager::GetHistogramHasMatchingParticle(const int iCentrality, const int iJetPt, const int iTrackPt){

  // If the histogram is NULL, try to load it from the input file
  if(fhTrackParticleMatchQA[kHasMatchingParticle][iCentrality][iJetPt][iTrackPt] == NULL){
    fhTrackParticleMatchQA[kHasMatchingParticle][iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(Form("%s/%s_C%dJ%dT%d", fTrackParticleMatchingQAName[kHasMatchingParticle], fTrackParticleMatchingQAName[kHasMatchingParticle], iCentrality, iJetPt, iTrackPt));
  }

  return fhTrackParticleMatchQA[kHasMatchingParticle][iCentrality][iJetPt][iTrackPt];
}

// Getter for deltaR response matrix between track pairs and matched particle pairs
TH2D* EECHistogramManager::GetHistogramTrackParticleDeltaRResponse(const int iCentrality, const int iJetPt, const int iTrackPt){

  // If the histogram is NULL, try to load it from the input file
  if(fhTrackParticleMatchQA[kTrackParticleMatchingDeltaRRresponse][iCentrality][iJetPt][iTrackPt] == NULL){
    fhTrackParticleResponse[kTrackParticleMatchingDeltaRRresponse][iCentrality][iJetPt][iTrackPt] = (TH2D*) fInputFile->Get(Form("%s/%s_C%dJ%dT%d", fTrackParticleMatchingResponseName[kTrackParticleMatchingDeltaRRresponse], fTrackParticleMatchingResponseName[kTrackParticleMatchingDeltaRRresponse], iCentrality, iJetPt, iTrackPt));
  }

  return fhTrackParticleResponse[kTrackParticleMatchingDeltaRRresponse][iCentrality][iJetPt][iTrackPt];
}

// Getter for pT response matrix between pT1*pT2 from track pairs and particle pairs
TH2D* EECHistogramManager::GetHistogramTrackParticlePtResponse(const int iCentrality, const int iJetPt, const int iTrackPt){

  // If the histogram is NULL, try to load it from the input file
  if(fhTrackParticleMatchQA[kTrackParticleMatchingPtResponse][iCentrality][iJetPt][iTrackPt] == NULL){
    fhTrackParticleResponse[kTrackParticleMatchingPtResponse][iCentrality][iJetPt][iTrackPt] = (TH2D*) fInputFile->Get(Form("%s/%s_C%dJ%dT%d", fTrackParticleMatchingResponseName[kTrackParticleMatchingPtResponse], fTrackParticleMatchingResponseName[kTrackParticleMatchingPtResponse], iCentrality, iJetPt, iTrackPt));
  }

  return fhTrackParticleResponse[kTrackParticleMatchingPtResponse][iCentrality][iJetPt][iTrackPt];
}

// Getter for track pT1*pT2 / particle pT1*pT2 histograms
TH1D* EECHistogramManager::GetHistogramTrackParticlePtClosure(const int iCentrality, const int iJetPt, const int iTrackPt){

  // If the histogram is NULL, try to load it from the input file
  if(fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt]){
    fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(Form("%s/%s_C%dJ%dT%d", fTrackParticlePtClosureSaveName, fTrackParticlePtClosureSaveName, iCentrality, iJetPt, iTrackPt));
  }

  return fhTrackParticlePtClosure[iCentrality][iJetPt][iTrackPt];
}

/*
 * Getter for any one-dimensional histogram based on input string
 *
 *  Arguments:
 *   TString name = Name corresponding to loaded histogram
 *   int bin1 = First bin index for the histogram
 *   int bin2 = Second bin index for the histogram
 *   int bin3 = Third bin index for the histogram
 *   int bin4 = Fourth bin index for the histogram
 *   int bin5 = Fifth bin index for the histogram
 *   int bin6 = Sixth bin index for the histogram
 *   int bin7 = Seventh bin index for the histogram
 *
 *   return: Histogram corresponding to given name and bins
 */
TH1D* EECHistogramManager::GetOneDimensionalHistogram(TString name, int bin1, int bin2, int bin3, int bin4, int bin5, int bin6, int bin7){
  if(name.EqualTo("vertexz",TString::kIgnoreCase) || name.EqualTo("fhvertexz",TString::kIgnoreCase)) return GetHistogramVertexZ();
  if(name.EqualTo("events",TString::kIgnoreCase) || name.EqualTo("fhevents",TString::kIgnoreCase)) return GetHistogramEvents();
  if(name.EqualTo("triggers",TString::kIgnoreCase) || name.EqualTo("fhtriggers",TString::kIgnoreCase)) return GetHistogramTriggers();
  if(name.EqualTo("trackcuts",TString::kIgnoreCase) || name.EqualTo("fhtrackcuts",TString::kIgnoreCase)) return GetHistogramTrackCuts();
  if(name.EqualTo("centrality",TString::kIgnoreCase) || name.EqualTo("fhcentrality",TString::kIgnoreCase)) return GetHistogramCentrality();
  if(name.EqualTo("centralityweighted",TString::kIgnoreCase) || name.EqualTo("fhcentralityweighted",TString::kIgnoreCase)) return GetHistogramCentralityWeighted();
  if(name.EqualTo("multiplicity",TString::kIgnoreCase) || name.EqualTo("fhmultiplicity",TString::kIgnoreCase)) return GetHistogramMultiplicity(bin1);
  if(name.EqualTo("multiplicityweighted",TString::kIgnoreCase) || name.EqualTo("fhmultiplicityweighted",TString::kIgnoreCase)) return GetHistogramMultiplicityWeighted(bin1);
  if(name.EqualTo("jetpt",TString::kIgnoreCase) || name.EqualTo("fhjetpt",TString::kIgnoreCase)) return GetHistogramJetPt(bin1);
  if(name.EqualTo("jetphi",TString::kIgnoreCase) || name.EqualTo("fhjetphi",TString::kIgnoreCase)) return GetHistogramJetPhi(bin1);
  if(name.EqualTo("jeteta",TString::kIgnoreCase) || name.EqualTo("fhjeteta",TString::kIgnoreCase)) return GetHistogramJetEta(bin1);
  if(name.EqualTo("jetetacm",TString::kIgnoreCase) || name.EqualTo("fhjetetacm",TString::kIgnoreCase)) return GetHistogramJetEtaCM(bin1);
  if(name.EqualTo("trackpt",TString::kIgnoreCase) || name.EqualTo("fhtrackpt",TString::kIgnoreCase)) return GetHistogramTrackPt(bin1,bin2);
  if(name.EqualTo("trackphi",TString::kIgnoreCase) || name.EqualTo("fhtrackphi",TString::kIgnoreCase)) return GetHistogramTrackPhi(bin1,bin2,bin3);
  if(name.EqualTo("tracketa",TString::kIgnoreCase) || name.EqualTo("fhtracketa",TString::kIgnoreCase)) return GetHistogramTrackEta(bin1,bin2,bin3);
  if(name.EqualTo("multiplicityinjetcone",TString::kIgnoreCase) || name.EqualTo("fhmultiplicityinjetcone",TString::kIgnoreCase)) return GetHistogramMultiplicityInJetCone(bin1,bin2,bin3,bin4,bin5);
  if(name.EqualTo("particledensityaroundjetaxis",TString::kIgnoreCase) || name.EqualTo("fhparticledensityaroundjetaxis",TString::kIgnoreCase)) return GetHistogramParticleDensityAroundJetAxis(bin1,bin2,bin3,bin4,bin5,bin6);
  if(name.EqualTo("energyenergycorrelator",TString::kIgnoreCase) || name.EqualTo("fhenergyenergycorrelator",TString::kIgnoreCase) || name.EqualTo("eec",TString::kIgnoreCase)) return GetHistogramEnergyEnergyCorrelator(bin1,bin2,bin3,bin4,bin5,bin6);
  if(name.EqualTo("energyenergycorrelatorprocessed",TString::kIgnoreCase) || name.EqualTo("fhenergyenergycorrelatorprocessed",TString::kIgnoreCase)) return GetHistogramEnergyEnergyCorrelatorProcessed(bin1,bin2,bin3,bin4,bin5);
  return NULL;
}

/*
 * Getter for any two-dimensional histogram based on input string
 *
 *  Arguments:
 *   TString name = Name corresponding to loaded histogram
 *   int bin1 = First bin index for the histogram
 *   int bin2 = Second bin index for the histogram
 *   int bin3 = Third bin index for the histogram
 *   int bin4 = Fourth bin index for the histogram
 *   int bin5 = Fifth bin index for the histogram
 *
 *   return: Histogram corresponding to given name and bins
 */
TH2D* EECHistogramManager::GetTwoDimensionalHistogram(TString name, int bin1, int bin2, int bin3, int bin4, int bin5){
  if(name.EqualTo("jetetaphi",TString::kIgnoreCase) || name.EqualTo("fhjetetaphi",TString::kIgnoreCase)) return GetHistogramJetEtaPhi(bin1);
  if(name.EqualTo("tracketaphi",TString::kIgnoreCase) || name.EqualTo("fhtracketaphi",TString::kIgnoreCase)) return GetHistogramTrackEtaPhi(bin1,bin2,bin3);
  if(name.EqualTo("multiplicitymap",TString::kIgnoreCase) || name.EqualTo("fhmultiplicitymap",TString::kIgnoreCase)) return GetHistogramMultiplicityMap();
  if(name.EqualTo("multiplicitymapweighted",TString::kIgnoreCase) || name.EqualTo("fhmultiplicitymapweighted",TString::kIgnoreCase)) return GetHistogramMultiplicityMap();
  return NULL;
}

// Get the first loaded centrality bin
int EECHistogramManager::GetFirstCentralityBin() const{
  return fFirstLoadedCentralityBin;
}

// Get the last loaded centrality bin
int EECHistogramManager::GetLastCentralityBin() const{
  return fLastLoadedCentralityBin;
}

// Get the first loaded track pT bin
int EECHistogramManager::GetFirstTrackPtBin() const{
  return fFirstLoadedTrackPtBin;
}

// Get the last loaded track pT bin
int EECHistogramManager::GetLastTrackPtBin() const{
  return fLastLoadedTrackPtBin;
}

// Get the first loaded energy-energy correlator jet pT bin
int EECHistogramManager::GetFirstJetPtBinEEC() const{
  return fFirstLoadedJetPtBinEEC;
}

// Get the last loaded energy-energy correlator jet pT bin
int EECHistogramManager::GetLastJetPtBinEEC() const{
  return fLastLoadedJetPtBinEEC;
}

// Get the first loaded energy-energy correlator track pT bin
int EECHistogramManager::GetFirstTrackPtBinEEC() const{
  return fFirstLoadedTrackPtBinEEC;
}

// Get the last loaded energy-energy correlator track pT bin
int EECHistogramManager::GetLastTrackPtBinEEC() const{
  return fLastLoadedTrackPtBinEEC;
}

// Getter for the number of events passing the cuts
int EECHistogramManager::GetNEvents(){

  // If the events histogram is not loaded, load it first
  if(fhEvents == NULL){
    LoadEventsHistogram();
  }

  return fhEvents->GetBinContent(fhEvents->FindBin(EECHistograms::kVzCut));
}

// Getter for the JCard
EECCard* EECHistogramManager::GetCard() const{
  return fCard;
}

// Getter for integral over inclusive jet pT. Include the overflow bin in the integral.
double EECHistogramManager::GetJetPtIntegral(const int iCentrality){

  // If the jet pT histogram is not loaded, load it first
  if(fhJetPt[iCentrality] == NULL){
    LoadJetPtHistogram(iCentrality);
  }

  return fhJetPt[iCentrality]->Integral(1,fhJetPt[iCentrality]->GetNbinsX()+1,"width");
}

/*
 * Getter for integral over inclusive jet pT over specified range
 *
 *  const int iCentrality = Centrality bin
 *  const double minPt = Lower pT range for integral calculation
 *  const double maxPt = Higher pT range for integral calculation
 */
double EECHistogramManager::GetJetPtIntegral(const int iCentrality, const double minPt, const double maxPt){

  // If the jet pT histogram is not loaded, load it first
  if(fhJetPt[iCentrality] == NULL){
    LoadJetPtHistogram(iCentrality);
  }

  return fhJetPt[iCentrality]->Integral(fhJetPt[iCentrality]->FindBin(minPt+0.001), fhJetPt[iCentrality]->FindBin(maxPt-0.001), "width");
}

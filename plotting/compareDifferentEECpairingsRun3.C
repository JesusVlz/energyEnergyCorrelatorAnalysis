#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"
#include "SystematicUncertaintyOrganizer.h"
#include "../src/EECHistograms.h"
#include "Run3ClosureTools.h"

/*
 * Macro for comparing final energy-energy correlators results
 */
// Copy of compareDifferentEECpairings.C at 792d76c0; original remains unchanged.
// plotKind: signal, background, signalWithTotal, or methods. The default is
// raw-scale closure (the legacy closure presets), NOT normalized pPb shapes.
// averageSameCone=false preserves the actual legacy C1C1 arithmetic. Set true
// explicitly for the two-cone-averaged definition; never mix definitions.
void compareDifferentEECpairingsRun3(
  TString fileName="merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root",
  TString plotKind="signal", bool correlatedRatio=false,
  double minJetPt=80., double maxJetPt=120.,
  double minCentrality=0., double maxCentrality=30., double trackThreshold=2.,
  TString outputDirectory="figures/run3LegacyClosure", bool averageSameCone=false,
  bool normalizeShapes=false, double ratioMinimum=0.8, double ratioMaximum=1.2,
  bool extendRatioRange=false){
  using Run3Closure::Require;
  Require(plotKind=="signal" || plotKind=="background" || plotKind=="signalWithTotal" || plotKind=="methods","Unsupported plotKind");
  Require(!correlatedRatio || ((plotKind=="signal" || plotKind=="signalWithTotal") && !normalizeShapes),
    "Shared-signal approximation is only defined for raw-scale SIGNAL closure, not background/normalized ratios");
  Require(minJetPt<maxJetPt && minCentrality<maxCentrality && ratioMinimum<ratioMaximum,"Invalid selection/range");
  Require(gSystem->mkdir(outputDirectory,true)==0 || !gSystem->AccessPathName(outputDirectory),"Cannot create output directory");

  
  // Predefined configuration for public figures
  bool truthToMixedConeBackgroundComparison = false;
  bool truthSignalToBackgroundSubtractedComparison = false;
  bool backgroundComponentsFromData = false;
  bool backgroundComponentsFromSimulation = false;

  // Sanity check for public figure drawing
  if(truthToMixedConeBackgroundComparison + truthSignalToBackgroundSubtractedComparison + backgroundComponentsFromData + backgroundComponentsFromSimulation > 1){
    cout << "ERROR! Only one predefined figure can be activated at a time." << endl;
    cout << "Please select only one figure to draw," << endl;
    return;
  } 

  // Explicit Run-3 file; public-figure presets must not override the input.
  TFile* inputFile=Run3Closure::Open(fileName);

  EECCard* card = new EECCard(inputFile);

  Run3Closure::CheckCard(*card);
  const bool hasCentrality=true; // OO bins, no PbPb MC centrality shift.
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();

  // Native OO selections are resolved by physical boundaries, never by PbPb indices.
  std::vector<std::pair<double,double>> comparedCentralityBin{{minCentrality,maxCentrality}};
  std::vector<std::pair<double,double>> comparedJetPtBin{{minJetPt,maxJetPt}};
  std::vector<double> comparedTrackPtBin{trackThreshold};
  int adjustCentralityPercentage=0;
  TString systemString="Pythia+Hijing OO 5.36 TeV";

  // Enumeration for different pairing types
  enum enumPairingType{kSameJetPair, kSignalReflectedConePair, kReflectedConePair, kSignalMixedConePair, kMixedConePair, kSignalSecondMixedConePair, kMixedMixedConePair, kSecondMixedConePair, kSignalPerpendicularConePair, kPerpendicularConePair, kSignalSecondPerpendicularConePair, kPerpendicularPerpendicularConePair, kSecondPerpendicularConePair, kCompiledBackgroundMixed, kCompiledBackgroundPerpendicularCone, kCompiledBackgroundTruthLevel, kCombinedSignalMixedConePair, kCombinedMixedConePair, kMixedBackgroundSubtracted, kPerpendicularBackgroundSubtracted, knPairingTypes};
  enum enumSubeventType{kPythiaPythia, kPythiaHydjet, kHydjetPythia, kHydjetHydjet, knSubeventCombinations, kAllBackground, knSubeventTypes};

  // Indices that go into pairing here:
  // First index = Jet cone pairing index
  //   EECHistograms::kSameJetPair
  //   EECHistograms::kSignalReflectedConePair
  //   EECHistograms::kReflectedConePair
  //   EECHistograms::kSignalMixedConePair
  //   EECHistograms::kReflectedMixedConePair
  //   EECHistograms::kMixedConePair
  //   EECHistograms::kSignalSecondMixedConePair
  //   EECHistograms::kReflectedSecondMixedConePair
  //   EECHistograms::kMixedMixedConePair
  //   EECHistograms::kSecondMixedConePair
  // Second index = Subevent pairing index
  //   EECHistograms::kPythiaPythia
  //   EECHistograms::kPythiaHydjet
  //   EECHistograms::kHydjetPythia
  //   EECHistograms::kHydjetHydjet
  //   EECHistograms::knSubeventCombinations (accept any subevent combination)
  std::vector<std::pair<int,int>> comparedEnergyEnergyCorrelatorPairings;
  comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, kPythiaPythia));
  comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kMixedBackgroundSubtracted, knSubeventCombinations));
  comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kPerpendicularBackgroundSubtracted, knSubeventCombinations));
  
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, kPythiaHydjet));
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, kHydjetHydjet));
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kMixedMixedConePair, knSubeventCombinations));
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kCombinedMixedConePair, knSubeventCombinations));
  
  // Choose legen position for the plots
  enum enumLegendPosition{kBottomLeft, kBottomRight, kTopLeft, kPaperClosure};
  int legendPosition = kPaperClosure;

  // Option to manually provide legend text for the compared distributions
  bool useManualLegend = true;
  std::vector<TString> manualLegend;
  manualLegend.push_back("True signal");
  manualLegend.push_back("Mixed cone signal");
  manualLegend.push_back("Perpendicular cone signal");
  //manualLegend.push_back("Pythia+EPOS");
  //manualLegend.push_back("EPOS+EPOS");

  // Flag for normalizing distributions to one over studied deltaR range
  bool normalizeDistributions = normalizeShapes;

  // Compare integrals of the distributions
  bool compareIntegrals = false;

  // It makes no sense to normalize distributions to 1 if you want to compare their integrals...
  if(compareIntegrals) normalizeDistributions = false;

  // OO uses the stored centrality directly; never apply the PbPb +4% shift.

  // Retain legacy local pairing enums and rendering; select the Run-3 closure.
  if(plotKind=="background") {
    comparedEnergyEnergyCorrelatorPairings={{kSameJetPair,kAllBackground},
      {kCompiledBackgroundMixed,knSubeventCombinations},{kCompiledBackgroundPerpendicularCone,knSubeventCombinations}};
    manualLegend={"True background","Mixed cone estimate","Perpendicular cone estimate"};
    legendPosition=kPaperClosure;
  } else if(plotKind=="methods") {
    comparedEnergyEnergyCorrelatorPairings={{kCompiledBackgroundPerpendicularCone,knSubeventCombinations},
      {kCompiledBackgroundMixed,knSubeventCombinations}};
    manualLegend={"Perpendicular cone background","Mixed cone background"};
    legendPosition=kPaperClosure;
  } else {
    legendPosition=kBottomLeft;
    if(plotKind=="signalWithTotal") {
      comparedEnergyEnergyCorrelatorPairings={{kSameJetPair,kPythiaPythia},
        {kSameJetPair,knSubeventCombinations},{kMixedBackgroundSubtracted,knSubeventCombinations},
        {kPerpendicularBackgroundSubtracted,knSubeventCombinations}};
      manualLegend={"True signal","Signal + background","Mixed cone signal","Perpendicular cone signal"};
    }
  }

  // ====================================================
  //                Drawing configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  TString saveComment = "_Run3_"+plotKind+(correlatedRatio?"_sharedSignalApprox":"_uncorrelated");   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  int weightExponent = card->GetWeightExponent();
  if(weightExponent == 2){
    saveComment.Prepend("_energyWeightSquared");
  } else {
    saveComment.Prepend("_nominalEnergyWeight");
  }

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(ratioMinimum, ratioMaximum);
  std::pair<double, double> eecZoom = std::make_pair(0.2, 30);
  bool logRatio = false;
  bool automaticZoom = true;

  // Configuration for predefined figures
  if(truthToMixedConeBackgroundComparison){
    ratioZoom = std::make_pair(0.91, 1.09);
    logRatio = false;
    automaticZoom = true;
  }

  if(truthSignalToBackgroundSubtractedComparison){
    //ratioZoom = std::make_pair(0, 14);
    ratioZoom = std::make_pair(0.8, 1.2);
    logRatio = false;
    automaticZoom = true;
  }

  // Configuration for predefined figures
  if(backgroundComponentsFromData || backgroundComponentsFromSimulation){
    ratioZoom = std::make_pair(0.01, 2);
    logRatio = true;
    automaticZoom = true;
  }

  // Sanity check for input. Ensure that all the selected bins actually exist in the input file.

  // Sanity check for centrality bins
  if(hasCentrality){
    for(auto centralityBin : comparedCentralityBin){
      if(card->FindBinIndexCentrality(centralityBin) < 0){
        cout << "ERROR! Centrality bin " << centralityBin.first << "-" << centralityBin.second << " does not exist in file " << fileName.Data() << endl;
        cout << "Please only choose centrality bins that are included in the input files." << endl;
        return;
      }
    }
  }

  // Sanity check for jet pT bins
  for(auto jetPtBin : comparedJetPtBin){
    if(card->FindBinIndexJetPtEEC(jetPtBin) < 0){
      cout << "ERROR! Jet pT bin " << jetPtBin.first << "-" << jetPtBin.second << " does not exist in file " << fileName.Data() << endl;
      cout << "Please only choose jet pT bins that are included in the input files." << endl;
      return;
    }
  }

  // Sanity check for track pT bins
  for(auto trackPtBin : comparedTrackPtBin){
    if(Run3Closure::Track(*card,trackPtBin) < 0){
      cout << "ERROR! Track pT cut > " << trackPtBin << " GeV does not exist in file " << fileName.Data() << endl;
      cout << "Please only choose track pT bins that are included in the input files." << endl;
      return;
    }
  } 
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms = new EECHistogramManager(inputFile, card);
  for(const auto& c:comparedCentralityBin) for(const auto& j:comparedJetPtBin) for(double t:comparedTrackPtBin)
    Run3Closure::Preflight(*histograms,card->FindBinIndexCentrality(c),card->FindBinIndexJetPtEEC(j),Run3Closure::Track(*card,t));


  // Load energy-energy correlators
  histograms->SetLoadEnergyEnergyCorrelators(true);

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][knPairingTypes][knSubeventTypes];
  TH1D* hEnergyEnergyCorrelatorRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][knPairingTypes][knSubeventTypes];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iPairing = 0; iPairing < knPairingTypes; iPairing++){
          for(int iSubevent = 0; iSubevent < knSubeventTypes; iSubevent++){
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][iPairing][iSubevent] = NULL;
            hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][iPairing][iSubevent] = NULL;
          } // Subevent type loop
        } // Pairing type loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // Helper variables
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality, iJetPt, iTrackPt;
  TH1D* helperHistogram;
  TH1D* backgroundHistogram;

  // Because we might define also derived histograms, we first need to load all possible pairing types and subevent combinations
  for(auto centralityBin : comparedCentralityBin){
    if(hasCentrality){
      iCentrality = card->FindBinIndexCentrality(centralityBin);
    } else {
      iCentrality = 0;
    }
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
        for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
          for(int iSubevent = 0; iSubevent <= EECHistograms::knSubeventCombinations; iSubevent++){

            // For data, there is no subevent division for the histograms
            if(!card->GetDataType().Contains("MC") && (iSubevent != EECHistograms::knSubeventCombinations)) continue;

            // Read the histogram corresponding to the defined bin
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = histograms->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType, iSubevent);
            if(auto* original=hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent])
              hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent]=
                Run3Closure::Clone(original,Form("run3_input_C%dJ%dT%dP%dS%d",iCentrality,iJetPt,iTrackPt,iPairingType,iSubevent));


          } // Subevent loop
        } // Pairing type loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // If we have derived histograms, derive them before any normalizations take place
  for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
  
    // Compile total background contribution from different parts in truth level
    if(pairingType.first == kCompiledBackgroundTruthLevel){
      for(auto centralityBin : comparedCentralityBin){
        if(hasCentrality){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-reflected cone pairs, which is the regular background estimate
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalReflectedConePair][knSubeventCombinations]->Clone(Form("compiledBackground%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Subtract from the background fake-fake contribution between cones, since this is badly modeled
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalReflectedConePair][kHydjetHydjet], -1);

            // Add the real fake+fake contribution back to the background estimate
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSameJetPair][kHydjetHydjet]);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Compiling total background from different components

    // Compile total background contribution from different parts using jet cones in mixed events
    if(pairingType.first == kCompiledBackgroundMixed || pairingType.first == kCompiledBackgroundPerpendicularCone){

      int perpendicularShifter = 0;
      if(pairingType.first == kCompiledBackgroundPerpendicularCone){
        perpendicularShifter = 5;
      }

      for(auto centralityBin : comparedCentralityBin){
        if(hasCentrality){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-mixed cone pairs, which is the basis of the background estimate
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalMixedConePair+perpendicularShifter][knSubeventCombinations]->Clone(Form("compiledBackground%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Add the other signal-mixed cone component and normalize to keep the correct scaling
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalSecondMixedConePair+perpendicularShifter][knSubeventCombinations]);
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Scale(0.5);

            // Subtract from the background fake-fake contribution between cones, since this is badly modeled
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedMixedConePair+perpendicularShifter][knSubeventCombinations], -1);

            // To improve statistics for the mixed cone only component, tae average of the two equivalent components
            helperHistogram = (TH1D*)hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]->Clone(Form("temporaryBackground%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));
            helperHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSecondMixedConePair+perpendicularShifter][knSubeventCombinations]);
            helperHistogram->Scale(0.5);

            // Add the real fake+fake contribution back to the background estimate
            // Default reproduces legacy arithmetic; optional averaged C1C1/C2C2 is explicit.
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(averageSameCone ? helperHistogram : hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Compiling total background from different components

    // Add together all background related subevents
    if(pairingType.second == kAllBackground){
      for(auto centralityBin : comparedCentralityBin){
        if(hasCentrality){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start from the signal+fake background
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][kPythiaHydjet]->Clone(Form("allBackground%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // If we are dealing with correlation between signal and refglected cones, add fake+signal component
            if(pairingType.first == kSignalReflectedConePair){
              hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][kHydjetPythia]);
            }

            // Next, add the fake+fake component to the total background
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][kHydjetHydjet]);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Adding together all background related subevents

    // Combine the two signal-mixed event components into one amalgamation
    if(pairingType.first == kCombinedSignalMixedConePair){
      for(auto centralityBin : comparedCentralityBin){
        if(hasCentrality){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-mixed cone pairs
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalMixedConePair][knSubeventCombinations]->Clone(Form("combineSignalMixed%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Add the signal-mixed cone pairs from the other mixed event
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalSecondMixedConePair][knSubeventCombinations]);

            // Take the average of these two to keep the scaling
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Scale(0.5);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Combining the two signal-mixed event components

    // Combine the two mixed event one components into one amalgamation
    if(pairingType.first == kCombinedMixedConePair){
      for(auto centralityBin : comparedCentralityBin){
        if(hasCentrality){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-mixed cone pairs
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair][knSubeventCombinations]->Clone(Form("combineMixedOnly%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Add the signal-mixed cone pairs from the other mixed event
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSecondMixedConePair][knSubeventCombinations]);

            // Take the average of these two to keep the scaling
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Scale(0.5);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Combining the two mixed event only components

    // Compile total background contribution from different parts using jet cones in mixed events. Then subtract it from the regular distribution
    if(pairingType.first == kMixedBackgroundSubtracted || pairingType.first == kPerpendicularBackgroundSubtracted){

      int perpendicularShifter = 0;
      if(pairingType.first == kPerpendicularBackgroundSubtracted){
        perpendicularShifter = 5;
      }

      for(auto centralityBin : comparedCentralityBin){
        if(hasCentrality){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-mixed cone pairs
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSameJetPair][knSubeventCombinations]->Clone(Form("backgroundSubtracted%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Start by taking the signal-mixed cone pairs, which is the basis of the background estimate
            backgroundHistogram = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalMixedConePair+perpendicularShifter][knSubeventCombinations]->Clone(Form("backgroundForSubtraction%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Add the other signal-mixed cone component and normalize to keep the correct scaling
            backgroundHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalSecondMixedConePair+perpendicularShifter][knSubeventCombinations]);
            backgroundHistogram->Scale(0.5);

            // Subtract from the background fake-fake contribution between cones, since this is badly modeled
            backgroundHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedMixedConePair+perpendicularShifter][knSubeventCombinations], -1);

            // To improve statistics for the mixed cone only component, tae average of the two equivalent components
            helperHistogram = (TH1D*)hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]->Clone(Form("temporaryBackgroundForSubtraction%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));
            helperHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSecondMixedConePair+perpendicularShifter][knSubeventCombinations]);
            helperHistogram->Scale(0.5);

            // Add the real fake+fake contribution back to the background estimate
            backgroundHistogram->Add(averageSameCone ? helperHistogram : hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]);

            // Finally subtract the background from the regular distribution
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(backgroundHistogram, -1);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Compiling total background from different components

  }

  // Go through all the histograms, normalize them, and calculate their ratios
  for(auto centralityBin : comparedCentralityBin){
    if(hasCentrality){
      iCentrality = card->FindBinIndexCentrality(centralityBin);
    } else {
      iCentrality = 0;
    }
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
        for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){

          // Normalize the histogram to one in the studied range
          if(normalizeDistributions){
            lowNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.second - epsilon);
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Scale(1 / hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Integral(lowNormalizationBin, highNormalizationBin, "width"));
          }

          // Calculate the ratio with respect to the first pairing type
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Clone(Form("pairingRatio%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));
          auto* denominator=hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][comparedEnergyEnergyCorrelatorPairings.at(0).first][comparedEnergyEnergyCorrelatorPairings.at(0).second];
          Run3Closure::PositiveDenominator(denominator,drawingRange.first,drawingRange.second);
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Divide(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][comparedEnergyEnergyCorrelatorPairings.at(0).first][comparedEnergyEnergyCorrelatorPairings.at(0).second]);
          if(correlatedRatio && pairingType!=comparedEnergyEnergyCorrelatorPairings.front()) {
            std::unique_ptr<TH1D> trueBackground(Run3Closure::TrueBackground(*histograms,iCentrality,iJetPt,iTrackPt));
            Run3Closure::ApplySharedSignal(
              hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second],
              hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second],denominator,
              hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSameJetPair][knSubeventCombinations],trueBackground.get());
          }


        } // Pairing type loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // ==========================================================================
  //               Draw all the pairing types in the same figure
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  drawer->SetLogX(true);

  TString compactJetPtString = "";
  TString compactTrackPtString = "";
  TString compactCentralityString = "";
  TString energyWeightString = (weightExponent == 1) ? "n=1" : "n=2";
  TString legendString;
  int markerStyle[5] = {kFullCircle, kOpenSquare, kOpenCross, kFullStar, kFullCross};
  int color[] = {kBlack,kRed,kBlue,kGreen+3,kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};
  int styleIndex, legendIndex;
  bool firstLoop, secondLoop;

  double legendX1 = 0.18;
  double legendY1 = 0.04;
  double legendX2 = 0.45;
  double legendY2 = 0.58;

  if(legendPosition == kBottomRight){
    legendX1 = 0.58;
    legendY1 = 0.02;
    legendX2 = 0.85;
    legendY2 = 0.52;
  } else if (legendPosition == kTopLeft){
    legendX1 = 0.16;
    legendY1 = 0.43;
    legendX2 = 0.43;
    legendY2 = 0.95;
  } else if (legendPosition == kPaperClosure){
    legendX1 = 0.53;
    legendY1 = 0.03;
    legendX2 = 0.80;
    legendY2 = 0.46;
  }

  // Automatic zoom helper veriables
  double binContent;
  int firstBin, lastBin;
  TLegend* legend;
  TLegend* tagLegend;

  for(auto centralityBin : comparedCentralityBin){
    if(hasCentrality){
      iCentrality = card->FindBinIndexCentrality(centralityBin);
      compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
    } else {
      iCentrality = 0;
      compactCentralityString = "";
    }
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
      compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
        compactTrackPtString = Form("_T>%.1f",trackPtBin);
        compactTrackPtString.ReplaceAll(".","v");
          
        // Optional range extension is explicit: the unmodified legacy zoom remains default.
        if(extendRatioRange) {
          double low=ratioMinimum,high=ratioMaximum;
          for(size_t curve=1;curve<comparedEnergyEnergyCorrelatorPairings.size();++curve) {
            auto pairing=comparedEnergyEnergyCorrelatorPairings[curve];
            auto* h=hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairing.first][pairing.second];
            for(int b=1;b<=h->GetNbinsX();++b) {
              if(h->GetXaxis()->GetBinUpEdge(b)<=drawingRange.first || h->GetXaxis()->GetBinLowEdge(b)>=drawingRange.second) continue;
              low=std::min(low,h->GetBinContent(b)-h->GetBinError(b));
              high=std::max(high,h->GetBinContent(b)+h->GetBinError(b));
            }
          }
          const double pad=.08*(high-low);ratioZoom={low-pad,high+pad};
        }
        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();

        // Logarithmic EEC axis
        drawer->SetLogY(true);

        // Create the legend and add jet and track pT information to it
        legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.048);legend->SetTextFont(62);

        if(hasCentrality){
          legend->AddEntry((TObject*)0,systemString,"");
          legend->AddEntry((TObject*)0,Form("Cent: %.0f-%.0f%%, %s",centralityBin.first,centralityBin.second,Run3Closure::Mode(*card).Data()),"");
        }
        legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second), "");
        legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV, %s", trackPtBin, energyWeightString.Data()), "");

        // Set drawing style for all histograms
        styleIndex = 0;
        for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetMarkerStyle(markerStyle[styleIndex]);
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetMarkerColor(color[styleIndex]);
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetLineColor(color[styleIndex]); 
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetMarkerStyle(markerStyle[styleIndex]);
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetMarkerColor(color[styleIndex]);
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetLineColor(color[styleIndex++]);
        } // File loop

        // Automatic zooming for the drawn histograms
        if(automaticZoom){
          eecZoom.first = 1e12;
          eecZoom.second = 0;
          for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
            firstBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.first);
            lastBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.second);
            for(int iBin = firstBin; iBin <= lastBin; iBin++){
              binContent = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetBinContent(iBin);
              if(binContent < eecZoom.first){
                if(binContent > 0){
                  eecZoom.first = binContent;
                }
              } 
              if(binContent > eecZoom.second) eecZoom.second = binContent;
            } // Bin loop
          } // File loop
          eecZoom.first = eecZoom.first / 2.0;
          eecZoom.second = eecZoom.second * 2.0;
        }

        firstLoop = true;
        legendIndex = 0;
        for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);

          if(firstLoop){
            drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second], "#Deltar", "EEC (A.U.)", " ");
            firstLoop = false;
          } else {
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Draw("same");
          }

          if(useManualLegend){
            legendString = manualLegend.at(legendIndex++);
          } else {
            legendString = histograms->GetPairingTypeSaveName(pairingType.first);
          }

          legend->AddEntry(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second], legendString.Data(), "p");
        }

        // Draw the legends to the upper pad
        legend->Draw();

        // Add CMS simulation tag to the canvas
        //if(legendPosition == kPaperClosure){
        //  tagLegend = new TLegend(0.16, 0.74, 0.43, 0.78);
        //  tagLegend->SetFillStyle(0);tagLegend->SetBorderSize(0);tagLegend->SetTextSize(0.07);tagLegend->SetTextFont(62);
        //  tagLegend->AddEntry((TObject*) 0, "CMS simulation", "");
        //  tagLegend->Draw();
        //}
          
        // Set linear scale for ratio, unless specifically asked to be logarithmic
        if(!logRatio) drawer->SetLogY(false);

        // Draw the histograms
        drawer->SetGridY(true);

        firstLoop = true;
        secondLoop = false;
        for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second); 
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
          if(firstLoop){
            firstLoop = false;
            secondLoop = true;
          } else if(secondLoop) {
            
            if(useManualLegend){
              if(truthToMixedConeBackgroundComparison){
                legendString = Form("#frac{Mixed cone bkg}{%s}", manualLegend.at(0).Data());
              } else {
                legendString = Form("#frac{Color}{%s}", manualLegend.at(0).Data());
              }
              
            } else {
              legendString = Form("#frac{Color}{%s}", histograms->GetPairingTypeSaveName(comparedEnergyEnergyCorrelatorPairings.at(0).first));
            }

            drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second], "#Deltar", legendString.Data(), " ");
            secondLoop = false;
          } else {
            hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Draw("same");
          }
        }

        drawer->SetGridY(false);
          
        // Save the figures to a file
        if(saveFigures){
          TCanvas* canvas=gPad->GetCanvas();
          Run3Closure::Header(canvas,correlatedRatio,averageSameCone,plotKind);
          TString stem=outputDirectory+"/"+Form("C%.0f-%.0f_J%.0f-%.0f_T%.1f_n%d_",centralityBin.first,centralityBin.second,jetPtBin.first,jetPtBin.second,trackPtBin,weightExponent)+plotKind;
          stem+=(averageSameCone?"_averagedCones":"_legacyC1");
          stem+=(normalizeDistributions?"_unitArea":"_raw");
          stem+=(correlatedRatio?"_sharedSignalApprox":"_uncorrelated");
          std::vector<TH1D*> curves,ratios;
          for(auto pairing:comparedEnergyEnergyCorrelatorPairings) {
            curves.push_back(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairing.first][pairing.second]);
            ratios.push_back(hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairing.first][pairing.second]);
          }
          Run3Closure::Save(canvas,stem,curves,ratios,manualLegend,fileName,correlatedRatio,averageSameCone,normalizeDistributions);
          std::cout<<"RUN3_CLOSURE_SAVED "<<stem<<std::endl;
        }

        if(compareIntegrals){
          if(hasCentrality){
            cout << "Bin: " << Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second) << " " << Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second) << " " << Form("%.1f < track p_{T}", trackPtBin) << endl;
          } else {
            cout << "Bin: " << Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second) << " " << Form("%.1f < track p_{T}", trackPtBin) << endl;
          }
          legendIndex = 0;
          for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
            lowNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.second - epsilon);
            cout << "Integral of " << manualLegend.at(legendIndex++) << " is " << hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Integral(lowNormalizationBin, highNormalizationBin, "width") << endl;
            cout << "Ratio to " << manualLegend.at(0) << " is " << hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Integral(lowNormalizationBin, highNormalizationBin, "width") / hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][comparedEnergyEnergyCorrelatorPairings.at(0).first][comparedEnergyEnergyCorrelatorPairings.at(0).second]->Integral(lowNormalizationBin, highNormalizationBin, "width") << endl;
          } // Loop over compared bins
        } // Compare integrals
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

}

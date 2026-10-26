#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECDrawer.h"

#include <bitset>

/*
 * Copy of plotEEChistograms.C for total same-jet OO MC EEC versus minimum track pT.
 * Select OO 0-30%, jet pT 40-80 GeV, and track thresholds 0.5, 1.0, 1.5, 2.0 GeV.
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 */
void plotEEChistogramsTrackPtComparison(TString inputFileName = "projected/OOMC_GenGen_withSubevents.root", double minCentrality = 0., double maxCentrality = 30.){

  // Print the file name to console
  cout << "Plotting MC EEC histograms from " << inputFileName.Data() << endl;

  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================

  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";

  // Choose which figure sets to draw
  bool drawEventInformation = false;
  bool drawJets = false;
  bool drawTracks = false;
  bool drawUncorrectedTracks = false;
  bool drawEnergyEnergyCorrelators = true;
  bool drawEnergyEnergyCorrelatorsEfficiencyVariationPlus = false;
  bool drawEnergyEnergyCorrelatorsEfficiencyVariationMinus = false;
  bool drawEnergyEnergyCorrelatorsPairEfficiencyVariationPlus = false;
  bool drawEnergyEnergyCorrelatorsPairEfficiencyVariationMinus = false;
  bool drawMultiplicityHistograms = false;
  bool drawParticleDensityAroundJets = false;
  bool drawParticlePtDensityAroundJets = false;
  bool drawParticleDensityAroundJetsPtBinned = false;
  bool drawParticlePtDensityAroundJetsPtBinned = false;
  bool drawMaxParticlePtWithinJetCone = false;
  bool drawMaxBackgroundParticlePtWithinJetCone = false;
  bool drawCovarianceMatrix = false;

  // Open the input file
  TFile* inputFile = TFile::Open(inputFileName);

  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // Load the card from the file and read the collision system
  EECCard* card = new EECCard(inputFile);
  TString collisionSystem = card->GetDataType();

  // ====================================================
  //  Binning configuration for the drawn histograms
  // ====================================================

  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nTrackPtBins = card->GetNTrackPtBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();

  // Updated OO MC card binning
  // centrality = {0,30,50,90}
  // track pT = {0.5,1,1.5,2,2.5,3,3.5,4,300}
  // jet pT in energy-energy correlators = {40,80,120,160,200,300,500,5020}
  // track pT in energy-energy correlators = {0.5,1,1.5,2,2.5,3,3.5,4,300}

  // Bin range to be drawn
  int firstDrawnCentralityBin = card->FindBinIndexCentrality(minCentrality,maxCentrality);
  if(firstDrawnCentralityBin < 0){
    cout << "Requested centrality bin is unavailable: " << minCentrality << "-" << maxCentrality << endl;
    return;
  }
  int lastDrawnCentralityBin = firstDrawnCentralityBin;

  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = 3;

  int firstDrawnJetPtBinEEC = 0;
  int lastDrawnJetPtBinEEC = 0; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin

  int firstDrawnTrackPtBinEEC = 0;
  int lastDrawnTrackPtBinEEC = 3; // Note: Tracks integrated over all pT ranges are in nTrackPtBinsEEC bin

  // Remove centrality selection from pp data
  if(collisionSystem.Contains("pp")){
    lastDrawnCentralityBin = 0;
  }

  // Select with multiplicity histograms to draw from all possible ones
  bool drawMultiplicityInJetCone = true && drawMultiplicityHistograms;
  bool drawMultiplicityInReflectedCone = false && drawMultiplicityHistograms;
  bool drawMultiplicityInJetConeUncorrected = false && drawMultiplicityHistograms;
  bool drawMultiplicityInReflectedConeUncorrected = false && drawMultiplicityHistograms;

  // Select track pairing type to be draw
  const bool drawSameJetEnergyEnergyCorrelator = true;       // Draw energy-energy correlator where tracks from the same jet are paired
  const bool drawSignalReflectedConeEnergyEnergyCorrelator = false; // Draw energy-energy correlator where tracks from jet cone are paired with tracks from reflected jet cone
  const bool drawReflectedConeOnlyEnergyEnergyCorrelator = false; // Draw energy-energy correlator where tracks from reflected jet cone are paired with tracks from reflected jet cone

  // Select which processed energy-energy correlators to draw
  const bool drawEnergyEnergyCorrelatorNormalized = false;    // Draw normalized energy-energy correlators
  const bool drawEnergyEnergyCorrelatorBackground = false;   // Draw normalized energy-energy correlator background estimate
  const bool drawEnergyEnergyCorrelatorSignal = false;       // Draw background subtracted energy-energy correlators

  // Figure saving
  const bool saveFigures = true;  // Save figures to file
  const char* figureFormat = "pdf";
  TString figureNameSuffix = Form("_OO_C%.0f-%.0f_GenGen_J40-80_T0p5-2p0",minCentrality,maxCentrality); // Suffix for figure names

  // Comment that can be added to figure legends
  TString legendComment = "";

  // Logarithmic scales for figures
  const bool logPt = true;          // pT distributions
  const bool logDeltaR = true;      // DeltaR axis for energy-energy correlators
  const bool logEEC = true;         // EEC axis for energy-energy correlators

  // Plotting style for 2D and 3D plots
  const int colorPalette = kLightTemperature;  // kRainBow kTemperatureMap kLightTemperature
  const char* style2D = "colz";
  const char* style3D = "surf1";

  // Select the style of histograms drawn for particle density histograms
  const bool drawIndividualParticleDensities = false;
  const bool drawParticleDensitiesForConstantJetPt = false;

  // Select the style of histograms drawn for energy-energy correlators
  const bool drawIndividualEnergyEnergyCorrelators = false;
  const bool drawEnergyEnergyCorrelatorsForConstantJetPt = true;
  const bool drawEnergyEnergyCorrelatorsForConstantTrackPt = false;
  bool drawEnergyEnergyCorrelatorsSubevent = false; // Draw energy-energy correlators with subevent decomposition (only for PbPb and pPb MC)

  // Select which subevents to draw
  bool drawAllSubevents = true;   // Draw histograms without subevent selection
  bool drawPythiaOnly = false;    // Draw only Pythia histograms in Pythia+Hydjet simulation
  bool drawHydjetOnly = false;    // Draw only Hydjet histograms in Pythia+Hydjet simulation

  bool drawAllSubeventPairs = true;  // Draw energy-energy correlators without subevent selection
  bool drawSignalOnly = false;        // Draw Pythia+Pythia correlations from MC
  bool drawSignalFake = false;        // Draw Pythia+Hydjet correlations from MC
  bool drawFakeFake = false;          // Draw Hydjet+Hydjet correlations from MC

  // If the collision system in not PbPb or pPb MC, we cannot draw subevent decomposition
  if(!(collisionSystem.Contains("Pb MC") || collisionSystem.Contains("OO MC"))){    drawEnergyEnergyCorrelatorsSubevent = false;

    drawPythiaOnly = false;
    drawHydjetOnly = false;
    drawAllSubevents = true;

    drawSignalOnly = false;
    drawSignalFake = false;
    drawFakeFake = false;
    drawAllSubeventPairs = true;
  }

  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================

  // ============================ //
  //     EECHistogramManager    //
  // ============================ //

  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager* histograms = new EECHistogramManager(inputFile,card);

  // Set the bin range we want to draw
  histograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  histograms->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  histograms->SetJetPtBinRangeEEC(firstDrawnJetPtBinEEC,lastDrawnJetPtBinEEC);
  histograms->SetTrackPtBinRangeEEC(firstDrawnTrackPtBinEEC,lastDrawnTrackPtBinEEC);

  // ============================ //
  //           EECDrawer          //
  // ============================ //

  // Create a new EECDrawer
  EECDrawer* resultDrawer = new EECDrawer(histograms);

  // Set which histograms to draw and the drawing style to use
  resultDrawer->SetDrawEventInformation(drawEventInformation);
  resultDrawer->SetDrawJetHistograms(drawJets);
  resultDrawer->SetDrawTracks(drawTracks);
  resultDrawer->SetDrawTracksUncorrected(drawUncorrectedTracks);
  resultDrawer->SetDrawEnergyEnergyCorrelor(drawEnergyEnergyCorrelators);
  resultDrawer->SetDrawEnergyEnergyCorrelorEfficiencyVariationPlus(drawEnergyEnergyCorrelatorsEfficiencyVariationPlus);
  resultDrawer->SetDrawEnergyEnergyCorrelorEfficiencyVariationMinus(drawEnergyEnergyCorrelatorsEfficiencyVariationMinus);
  resultDrawer->SetDrawEnergyEnergyCorrelorPairEfficiencyVariationPlus(drawEnergyEnergyCorrelatorsPairEfficiencyVariationPlus);
  resultDrawer->SetDrawEnergyEnergyCorrelorPairEfficiencyVariationMinus(drawEnergyEnergyCorrelatorsPairEfficiencyVariationMinus);

  resultDrawer->SetDrawMultiplicityInJetCone(drawMultiplicityInJetCone);
  resultDrawer->SetDrawMultiplicityInReflectedCone(drawMultiplicityInReflectedCone);
  resultDrawer->SetDrawMultiplicityInJetConeUncorrected(drawMultiplicityInJetConeUncorrected);
  resultDrawer->SetDrawMultiplicityInReflectedConeUncorrected(drawMultiplicityInReflectedConeUncorrected);

  resultDrawer->SetDrawParticleDensityAroundJetAxis(drawParticleDensityAroundJets);
  resultDrawer->SetDrawParticlePtDensityAroundJetAxis(drawParticlePtDensityAroundJets);
  resultDrawer->SetDrawParticleDensityAroundJetAxisPtBinned(drawParticleDensityAroundJetsPtBinned);
  resultDrawer->SetDrawParticlePtDensityAroundJetAxisPtBinned(drawParticlePtDensityAroundJetsPtBinned);

  resultDrawer->SetDrawSingleParticleDensityHistograms(drawIndividualParticleDensities);
  resultDrawer->SetDrawParticleDensityForConstantJetPt(drawParticleDensitiesForConstantJetPt);

  resultDrawer->SetDrawMaxParticlePtWithinJetCone(drawMaxParticlePtWithinJetCone);
  resultDrawer->SetDrawMaxBackgroundParticlePtWithinJetCone(drawMaxBackgroundParticlePtWithinJetCone);

  resultDrawer->SetDrawSingleEnergyEnergyCorrelators(drawIndividualEnergyEnergyCorrelators);
  resultDrawer->SetDrawEnergyEnergyCorrelatorsForConstantJetPt(drawEnergyEnergyCorrelatorsForConstantJetPt);
  resultDrawer->SetDrawEnergyEnergyCorrelatorsForConstantTrackPt(drawEnergyEnergyCorrelatorsForConstantTrackPt);
  resultDrawer->SetDrawEnergyEnergyCorrelatorsSubevent(drawEnergyEnergyCorrelatorsSubevent);

  resultDrawer->SetDrawSameJetEnergyEnergyCorrelators(drawSameJetEnergyEnergyCorrelator);
  resultDrawer->SetDrawSignalReflectedConeEnergyEnergyCorrelators(drawSignalReflectedConeEnergyEnergyCorrelator);
  resultDrawer->SetDrawReflectedConeOnlyEnergyEnergyCorrelators(drawReflectedConeOnlyEnergyEnergyCorrelator);

  resultDrawer->SetDrawEnergyEnergyCorrelatorNormalized(drawEnergyEnergyCorrelatorNormalized);
  resultDrawer->SetDrawEnergyEnergyCorrelatorBackground(drawEnergyEnergyCorrelatorBackground);
  resultDrawer->SetDrawEnergyEnergyCorrelatorSignal(drawEnergyEnergyCorrelatorSignal);

  resultDrawer->SetDrawCovarianceMatrices(drawCovarianceMatrix);

  resultDrawer->SetDrawAllSubeventTypes(drawAllSubevents, drawPythiaOnly, drawHydjetOnly);
  resultDrawer->SetDrawAllSubeventCombinations(drawAllSubeventPairs, drawSignalOnly, drawSignalFake, drawFakeFake);

  resultDrawer->SetSaveFigures(saveFigures,figureFormat,figureNameSuffix);
  resultDrawer->SetLogPt(logPt);
  resultDrawer->SetLogDeltaR(logDeltaR);
  resultDrawer->SetLogEEC(logEEC);
  resultDrawer->SetLegendComment(legendComment);
  resultDrawer->SetDrawingStyles(colorPalette,style2D,style3D);

  // Draw the selected histograms
  resultDrawer->DrawHistograms();

}

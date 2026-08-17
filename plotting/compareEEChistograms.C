#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECComparingDrawer.h"
#include "EECCard.h"

/*
 * Macro for configuring the EECComparingDrawer and defining which histograms are drawn
 */
void compareEEChistograms(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Define the used data files, and a comment describing the data in each file
  const int nDatasets = 5;
  TString inputFileName[] = {"/afs/cern.ch/user/j/jvelazqu/EECAnalysis/CMSSW_15_0_9_patch4/src/energyEnergyCorrelatorAnalysis/merged_R4_axisWTA_minbias_processed.root", 
    "/afs/cern.ch/user/j/jvelazqu/EECAnalysis/CMSSW_15_0_9_patch4/src/energyEnergyCorrelatorAnalysis/merged_R4_axisWTA_primaryvertex_processed.root", 
    "/afs/cern.ch/user/j/jvelazqu/EECAnalysis/CMSSW_15_0_9_patch4/src/energyEnergyCorrelatorAnalysis/merged_R4_axisWTA_zvertex_processed.root",
    "/afs/cern.ch/user/j/jvelazqu/EECAnalysis/CMSSW_15_0_9_patch4/src/energyEnergyCorrelatorAnalysis/merged_R4_axisWTA_pileupf_processed.root",
    "/afs/cern.ch/user/j/jvelazqu/EECAnalysis/CMSSW_15_0_9_patch4/src/energyEnergyCorrelatorAnalysis/merged_R4_axisWTA_hfcoinc_processed.root"
  };
  //TString inputFileName[] = {"data/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_processed_2025-02-19.root", "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_truthReference_processed_2024-01-10.root"};
  // eecAnalysis_akFlowJets_updatedMultiplicityAndDensity_eschemeAxis_preprocessed_2022-10-17.root
  // eecAnalysis_akFlowJets_updatedMultiplicityAndDensity_wtaAxis_preprocessed_2022-10-17.root
  // eecAnalysis_akFlowJets_removeBadAcceptance_wtaAxis_processed_2022-10-25.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_finalMcWeight_matchJets_processed_2023-03-06.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_finalMcWeight_processed_2023-03-08.root
  // data/MinBiasHydjet_RecoGen_eecAnalysis_akFlowJet_firstMinBiasScan_noTrigger_preprocessed_2022-10-10.root
  
  TString legendBase[] = {"Min Bias Selection", "+ Primary Vertex Filter", "+ Z Vertex Filter", "+ Pileup Filter", "+ HF Coincidence Filter"};
  TString legendComment[nDatasets];
  
  // Try to open the files
  TFile* inputFile[nDatasets];
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    // Open the file for the given dataset
    inputFile[iDataset] = TFile::Open(inputFileName[iDataset]);
    if(inputFile[iDataset] == NULL){
      cout << "Error! The file " << inputFileName[iDataset].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Please give a file that exists. Will not execute the code" << endl;
      return;
    }
  }
   
  // Read the EECCard from the first file. This is used to gain binning information, and assumed to be the same in other files too
  EECCard* card = new EECCard(inputFile[0]);
  
  // Choose which figure sets to draw
  bool drawEventInformation = false;
  bool drawJets = false;
  bool drawTracks = false;
  bool drawUncorrectedTracks = false;
  
  // Multiplicity in jet cone
  bool drawMultiplicityInJetCone = false;
  bool drawMultiplicityInReflectedCone = false;
  bool drawMultiplicityInJetConeUncorrected = false;
  bool drawMultiplicityInReflectedConeUncorrected = false;
  
  // Particle density around jets
  bool drawParticleDensityAroundJets = false;
  bool drawParticlePtDensityAroundJets = false;
  bool drawParticleDensityAroundJetsPtBinned = false;
  bool drawParticlePtDensityAroundJetsPtBinned = false;
  
  // Energy-energy correlators
  bool drawEnergyEnergyCorrelators = true;
  bool drawEnergyEnergyCorrelatorsEfficiencyVariationPlus = false;
  bool drawEnergyEnergyCorrelatorsEfficiencyVariationMinus = false;
  bool drawEnergyEnergyCorrelatorsPairEfficiencyVariationPlus = false;
  bool drawEnergyEnergyCorrelatorsPairEfficiencyVariationMinus = false;
  
  // Select which pairing types to draw
  const bool drawSameJetEnergyEnergyCorrelator = true;       // Draw energy-energy correlator where tracks from the same jet are paired
  const bool drawSignalReflectedConeEnergyEnergyCorrelator = false; // Draw energy-energy correlator where tracks from jet cone are paired with tracks from reflected jet cone
  const bool drawReflectedConeOnlyEnergyEnergyCorrelator = false; // Draw energy-energy correlator where tracks from reflected jet cone are paired with tracks from reflected jet cone
  
  // Select which processed energy-energy correlators to draw
  const bool drawEnergyEnergyCorrelatorNormalized = false;    // Draw normalized energy-energy correlators
  const bool drawEnergyEnergyCorrelatorBackground = false;   // Draw normalized energy-energy correlator background estimate
  const bool drawEnergyEnergyCorrelatorSignal = false;       // Draw background subtracted energy-energy correlators
  const bool drawEnergyEnergyCorrelatorUnfolded = false;            // Draw unfolded energy-energy correlators
  const bool drawEnergyEnergyCorrelatorUnfoldedBackground = false;  // Draw energy-energy correlator background estimate after unfolding
  const bool drawEnergyEnergyCorrelatorUnfoldedSignal = false;      // Draw background subtracted energy-energy correlators after unfolding
  
  // Select which subevents to draw
  bool drawAllSubevents = true;   // Draw histograms without subevent selection
  bool drawPythiaOnly = false;    // Draw only Pythia histograms in Pythia+Hydjet simulation
  bool drawHydjetOnly = false;    // Draw only Hydjet histograms in Pythia+Hydjet simulation
  
  bool drawAllSubeventPairs = true;  // Draw energy-energy correlators without subevent selection
  bool drawSignalOnly = false;        // Draw Pythia+Pythia correlations from MC
  bool drawSignalFake = false;        // Draw Pythia+Hydjet correlations from MC
  bool drawFakeFake = false;          // Draw Hydjet+Hydjet correlations from MC
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = true;
  const char* figureFormat = "pdf";
  const char* figureComment = "_jetAxisComparisonSquared";
  
  // Logarithmic scales for figures
  bool logPt = true;       // pT axis for jet
  bool logDeltaR = true;   // deltaR axis for energy-energy correlators
  bool logEEC = true;      // EEC axis for energy-energy correlators
  bool logParticleDensity = false;  // Logarithmic y-axis for particle densities
  
  // Plotting style for 2D and 3D plots
  int colorPalette = kRainBow;
  const char* style2D = "colz";
  const char* style3D = "surf1";
  
  // Drawing style for the histograms
  const int lineWidth = 1; // Line width for the drawn histograms
  
  // Settings for ratios
  bool useDifferenceInsteadOfRatio = false;
  double minZoom = 0.8;
  double maxZoom = 1.2;
  //TString ratioLabel = "#frac{MC}{Data}";
  TString ratioLabel = "Ratio to Min Bias";

  bool manualLegend = true; // Set this true if you want to set legend manually in EECComparingDrawer.cxx code instead of using automatic legend generation
  bool addSystemToLegend = true;  // Add the collision system from first file to legend. Useful if all files are from same system
  bool includeMCtype = false;      // Include MC type in the system
  bool addEnergyToLegend = true;  // Add the collision energy from the first file to legend. Useful if all files are from same system
  
  // Scaling for histograms
  int scaleHistograms = 1; // 0 = Do not scale histograms. 1 = Scale integral to one. 2 = Scale average to one
  bool rebinJetPt = false;
  
  // ====================================================
  //  Binning configuration for the drawn histograms
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nTrackPtBins = card->GetNTrackPtBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,2,3,4,6,300}
  
  // Bin range to be drawn
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = 0;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = 0;
  
  int firstDrawnJetPtBinEEC = 2;
  int lastDrawnJetPtBinEEC = 2; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin

  int firstDrawnTrackPtBinEEC =3;
  int lastDrawnTrackPtBinEEC = 3;
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // ==================================================================
  // ================= Setup EECHistogramManagers =====================
  // ==================================================================
  EECHistogramManager* histograms[nDatasets];
  
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    
    // Create a new histogram manager
    histograms[iDataset] = new EECHistogramManager(inputFile[iDataset]);

    // Set which histograms to draw from the input file
    /*histograms[iDataset]->SetLoadEventInformation(drawEventInformation);
    histograms[iDataset]->SetLoadJetHistograms(drawJets);
    histograms[iDataset]->SetLoadTracks(drawTracks);
    histograms[iDataset]->SetLoadTracksUncorrected(drawUncorrectedTracks);
    histograms[iDataset]->SetLoadMultiplicityInJets(drawMultiplicityInJetCone || drawMultiplicityInReflectedCone || drawMultiplicityInJetConeUncorrected || drawMultiplicityInReflectedConeUncorrected);
    histograms[iDataset]->SetLoadParticleDensityAroundJets(drawParticleDensityAroundJets);
    histograms[iDataset]->SetLoadParticlePtDensityAroundJets(drawParticlePtDensityAroundJets);
    histograms[iDataset]->SetLoadParticleDensityAroundJetsPtBinned(drawParticleDensityAroundJetsPtBinned);
    histograms[iDataset]->SetLoadParticlePtDensityAroundJetsPtBinned(drawParticlePtDensityAroundJetsPtBinned);
    histograms[iDataset]->SetLoadEnergyEnergyCorrelators(drawEnergyEnergyCorrelators);
    histograms[iDataset]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(drawEnergyEnergyCorrelatorsEfficiencyVariationPlus);
    histograms[iDataset]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(drawEnergyEnergyCorrelatorsEfficiencyVariationMinus);
    histograms[iDataset]->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(drawEnergyEnergyCorrelatorsPairEfficiencyVariationPlus);
    histograms[iDataset]->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(drawEnergyEnergyCorrelatorsPairEfficiencyVariationMinus);*/
    
    histograms[iDataset]->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
    histograms[iDataset]->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
    histograms[iDataset]->SetJetPtBinRangeEEC(firstDrawnJetPtBinEEC,lastDrawnJetPtBinEEC);
    histograms[iDataset]->SetTrackPtBinRangeEEC(firstDrawnTrackPtBinEEC,lastDrawnTrackPtBinEEC);

    // Build legend string with event count. GetNEvents() lee el histograma
    // Build legend string with the ENTRIES of the specific EEC projection
    // (not the raw event count). Bin combo must match firstDrawnCentralityBin /
    // firstDrawnJetPtBinEEC / firstDrawnTrackPtBinEEC set above (currently C2/J0/T4).
    // iSubevent = EECHistograms::knSubeventCombinations selects the inclusive/combined
    // subevent bin, correct for real data (no MC subevent split).
    TH1D* legendRefHisto = histograms[iDataset]->GetHistogramEnergyEnergyCorrelator(
        EECHistogramManager::kEnergyEnergyCorrelator,
        firstDrawnCentralityBin,
        firstDrawnJetPtBinEEC,
        firstDrawnTrackPtBinEEC,
        EECHistograms::kSameJetPair,
        EECHistograms::knSubeventCombinations
    );
    legendComment[iDataset] = Form("%s ( %.0f )", legendBase[iDataset].Data(), legendRefHisto->GetEntries());
    
    // Load the histograms from the file
    //histograms[iDataset]->LoadProcessedHistograms();

  } // Loop over datasets
  
  // ==================================================================
  // ==== Setup EECComparingDrawer and draw the comparison plots ======
  // ==================================================================
  
  EECComparingDrawer* drawer = new EECComparingDrawer(histograms[0]);
  drawer->AddLegendComment(legendComment[0]);
  for(int i = 1; i < nDatasets; i++){
    drawer->AddHistogramToDraw(histograms[i]);
    drawer->AddLegendComment(legendComment[i]);
  }
  
  drawer->SetDrawEventInformation(drawEventInformation);
  drawer->SetDrawJets(drawJets);
  drawer->SetDrawAllTracks(drawTracks,drawUncorrectedTracks);
  
  drawer->SetDrawMultiplicityInJetCone(drawMultiplicityInJetCone);
  drawer->SetDrawMultiplicityInReflectedCone(drawMultiplicityInReflectedCone);
  drawer->SetDrawMultiplicityInJetConeUncorrected(drawMultiplicityInJetConeUncorrected);
  drawer->SetDrawMultiplicityInReflectedConeUncorrected(drawMultiplicityInReflectedConeUncorrected);
  
  drawer->SetDrawParticleDensityAroundJetAxis(drawParticleDensityAroundJets);
  drawer->SetDrawParticlePtDensityAroundJetAxis(drawParticlePtDensityAroundJets);
  drawer->SetDrawParticleDensityAroundJetAxisPtBinned(drawParticleDensityAroundJetsPtBinned);
  drawer->SetDrawParticlePtDensityAroundJetAxisPtBinned(drawParticlePtDensityAroundJetsPtBinned);
  
  drawer->SetDrawEnergyEnergyCorrelator(drawEnergyEnergyCorrelators);
  drawer->SetDrawEnergyEnergyCorrelatorEfficiencyVariationPlus(drawEnergyEnergyCorrelatorsEfficiencyVariationPlus);
  drawer->SetDrawEnergyEnergyCorrelatorEfficiencyVariationMinus(drawEnergyEnergyCorrelatorsEfficiencyVariationMinus);
  drawer->SetDrawEnergyEnergyCorrelatorPairEfficiencyVariationPlus(drawEnergyEnergyCorrelatorsPairEfficiencyVariationPlus);
  drawer->SetDrawEnergyEnergyCorrelatorPairEfficiencyVariationMinus(drawEnergyEnergyCorrelatorsPairEfficiencyVariationMinus);
  
  drawer->SetDrawSameJetEnergyEnergyCorrelators(drawSameJetEnergyEnergyCorrelator);
  drawer->SetDrawSignalReflectedConeEnergyEnergyCorrelators(drawSignalReflectedConeEnergyEnergyCorrelator);
  drawer->SetDrawReflectedConeOnlyEnergyEnergyCorrelators(drawReflectedConeOnlyEnergyEnergyCorrelator);
  
  drawer->SetDrawEnergyEnergyCorrelatorNormalized(drawEnergyEnergyCorrelatorNormalized);
  drawer->SetDrawEnergyEnergyCorrelatorBackground(drawEnergyEnergyCorrelatorBackground);
  drawer->SetDrawEnergyEnergyCorrelatorSignal(drawEnergyEnergyCorrelatorSignal);
  drawer->SetDrawEnergyEnergyCorrelatorUnfolded(drawEnergyEnergyCorrelatorUnfolded);
  drawer->SetDrawEnergyEnergyCorrelatorUnfoldedBackground(drawEnergyEnergyCorrelatorUnfoldedBackground);
  drawer->SetDrawEnergyEnergyCorrelatorUnfoldedSignal(drawEnergyEnergyCorrelatorUnfoldedSignal);
  
  drawer->SetDrawAllSubeventTypes(drawAllSubevents, drawPythiaOnly, drawHydjetOnly);
  drawer->SetDrawAllSubeventCombinations(drawAllSubeventPairs, drawSignalOnly, drawSignalFake, drawFakeFake);
  
  drawer->SetSaveFigures(saveFigures,figureFormat,figureComment);
  drawer->SetLogAxes(logPt, logDeltaR, logEEC);
  drawer->SetLogParticleDensity(logParticleDensity);
  drawer->SetDrawingStyles(colorPalette,style2D,style3D);
  drawer->SetUseDifferenceInRatioPlot(useDifferenceInsteadOfRatio);
  drawer->SetRatioZoom(minZoom,maxZoom);
  drawer->SetRatioLabel(ratioLabel);
  drawer->SetApplyScaling(scaleHistograms);
  drawer->SetJetPtRebin(rebinJetPt);
  drawer->SetManualLegend(manualLegend);
  drawer->SetAddSystemToLegend(addSystemToLegend, includeMCtype);
  drawer->SetAddEnergyToLegend(addEnergyToLegend);
  
  drawer->SetLineWidth(lineWidth);

  // Draw the selected histograms
  drawer->DrawHistograms();
  
}

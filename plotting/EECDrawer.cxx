/*
 * Implementation of EECDrawer
 */

// Root includes
#include <TPad.h>

// Own includes
#include "EECDrawer.h"

/*
 * Constructor
 */
EECDrawer::EECDrawer(EECHistogramManager *inputHistograms) :
  fHistograms(inputHistograms),
  fFigureSaveNameAppend(""),
  fLegendComment(""),
  fDrawEventInformation(false),
  fDrawJets(false),
  fDrawCovarianceMatrix(false),
  fDrawIndividualParticleDensities(true),
  fDrawParticleDensitiesForConstantJetPt(false),
  fDrawIndividualEnergyEnergyCorrelators(true),
  fDrawEnergyEnergyCorrelatorsForConstantJetPt(false),
  fDrawEnergyEnergyCorrelatorsForConstantTrackPt(false),
  fDrawEnergyEnergyCorrelatorsSubevents(false),
  fSaveFigures(false),
  fFigureFormat("pdf"),
  fLogPt(true),
  fLogDeltaR(true),
  fLogEEC(true),
  fColorPalette(kRainBow),
  fStyle2D("colz"),
  fStyle3D("surf1")
{
  
  // Read card from inputfile and collision system from card
  //TString collisionSystem = fHistograms->GetSystem();
  TString collisionSystem = fHistograms->GetCard()->GetAlternativeDataType(false);
  
  // Make a string for collision system based on information on the card
  TString collisionEnergy = "5.02 TeV";
  if(collisionSystem.Contains("OO")) collisionEnergy = "5.36 TeV";
  fSystemAndEnergy = Form("%s %s", collisionSystem.Data(), collisionEnergy.Data());  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Create a new drawer
  fDrawer = new JDrawer();

  // Do not draw anything by default
  for(int iTrackType = 0; iTrackType < EECHistogramManager::knTrackCategories; iTrackType++){
    fDrawTracks[iTrackType] = false;
  }
  for(int iMultiplicityType = 0; iMultiplicityType < EECHistogramManager::knMultiplicityInJetConeTypes; iMultiplicityType++){
    fDrawMultiplicityInJetCone[iMultiplicityType] = false;
  }
  for(int iParticleDensityType = 0; iParticleDensityType < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensityType++){
    fDrawParticleDensityAroundJets[iParticleDensityType] = false;
  }
  for(int iMaxParticlePtWithinJetCone = 0; iMaxParticlePtWithinJetCone < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iMaxParticlePtWithinJetCone++){
    fDrawMaxParticlePtWithinJetCone[iMaxParticlePtWithinJetCone] = false;
  }
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    fDrawEnergyEnergyCorrelators[iEnergyEnergyCorrelator] = false;
  }
  for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
    fDrawPairingType[iPairingType] = false;
  }
  for(int iProcessLevel = 0; iProcessLevel < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels; iProcessLevel++){
    fDrawProcessingLevel[iProcessLevel] = false;
  }
  
  // By default, only draw the energy-energy correlators without subevent selection
  fDrawSubeventType[EECHistograms::knSubeventTypes] = true;
  for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){
    fDrawSubeventType[iSubevent] = false;
  }
  fDrawSubeventCombination[EECHistograms::knSubeventCombinations] = true;
  for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
    fDrawSubeventCombination[iSubevent] = false;
  }

  // Setup the centrality, track pT and energy-energy correlator jet and track pT bins to be drawn
  fFirstDrawnCentralityBin = fHistograms->GetFirstCentralityBin();
  fLastDrawnCentralityBin = fHistograms->GetLastCentralityBin();
  fFirstDrawnTrackPtBin = fHistograms->GetFirstTrackPtBin();
  fLastDrawnTrackPtBin = fHistograms->GetLastTrackPtBin();
  fFirstDrawnJetPtBinEEC = fHistograms->GetFirstJetPtBinEEC();
  fLastDrawnJetPtBinEEC = fHistograms->GetLastJetPtBinEEC();
  fFirstDrawnTrackPtBinEEC = fHistograms->GetFirstTrackPtBinEEC();
  fLastDrawnTrackPtBinEEC = fHistograms->GetLastTrackPtBinEEC();
  
}

/*
 * Destructor
 */
EECDrawer::~EECDrawer(){
  delete fDrawer;
}

/*
 * Draw all the selected histograms using JDrawer
 */
void EECDrawer::DrawHistograms(){
  
  // Draw the event information histograms
  DrawEventInformation();

  // Draw the jet histograms
  DrawJetHistograms();
  
  // Draw the track histograms
  DrawTrackHistograms();
  
  // Draw the mulplicity within the jet cone
  DrawMultiplicityInJetCone();
  
  // Draw the particle density around the jet axis
  DrawParticleDensityAroundJetAxis();
  
  // Draw the maximum particle pT within the jet cone histograms
  DrawMaxParticlePtWithinJetCone();
  
  // Draw the energy-energy correlation histograms
  DrawEnergyEnergyCorrelationHistograms();
  
  // Draw the processed energy-energy correlator histograms
  DrawProcessedEnergyEnergyCorrelators();

  // Draw the covariance matrices
  DrawCovarianceMatrices();
  
}

/*
 * Draw event information histograms
 */
void EECDrawer::DrawEventInformation(){
  
  if(!fDrawEventInformation) return;  // Only draw the event information histograms if they are selected for drawing
  
  // Helper variables for histogram drawing
  TH1D* drawnHistogram;
  TLegend* legend;
  
  // === Vertex z-position ===
  drawnHistogram = fHistograms->GetHistogramVertexZ();
  drawnHistogram->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(drawnHistogram,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
  legend = new TLegend(0.65,0.75,0.85,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("vz");
  
  // === Event cuts ===
  drawnHistogram = fHistograms->GetHistogramEvents();
  fDrawer->DrawHistogram(drawnHistogram," ","Number of events", " ");
  legend = new TLegend(0.17,0.22,0.37,0.37);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("eventCuts");
  
  // === Track cuts ===
  drawnHistogram = fHistograms->GetHistogramTrackCuts();
  fDrawer->DrawHistogram(drawnHistogram," ","Number of tracks", " ");
  legend = new TLegend(0.65,0.75,0.85,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("trackCuts");
  
  // === Centrality ===
  drawnHistogram = fHistograms->GetHistogramCentrality();
  drawnHistogram->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(drawnHistogram,"Centrality percentile","N"," ");
  legend = new TLegend(0.63,0.75,0.83,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("centrality");
  
}

/*
 * Draw jet histograms
 */
void EECDrawer::DrawJetHistograms(){
  
  // If jet histograms are not drawn, there is nothing to do here
  if(!fDrawJets) return;
  
  // Helper variables for histogram drawing
  TH1D* drawnHistogram;
  TH2D* drawnHistogram2D;
  TLegend* legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString namerX;
  TString namerY;
  
  // Loop over centrality
  for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
    compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
    
    // Select logarithmic drawing for pT
    fDrawer->SetLogY(fLogPt);
    
    // === Jet pT ===
    drawnHistogram = fHistograms->GetHistogramJetPt(iCentrality);
    namerX = Form("%s p_{T}  (GeV)",fHistograms->GetJetAxisName());
    fDrawer->DrawHistogram(drawnHistogram,namerX.Data(),"#frac{dN}{dp_{T}}  (1/GeV)"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    namerX = Form("%sPt",fHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // Set linear drawing
    fDrawer->SetLogY(false);
    
    // === Jet phi ===
    drawnHistogram = fHistograms->GetHistogramJetPhi(iCentrality);
    namerX = Form("%s #varphi",fHistograms->GetJetAxisName());
    fDrawer->DrawHistogram(drawnHistogram,namerX.Data(),"#frac{dN}{d#varphi}"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    namerX = Form("%sPhi",fHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // === Jet eta ===
    drawnHistogram = fHistograms->GetHistogramJetEta(iCentrality);
    namerX = Form("%s #eta",fHistograms->GetJetAxisName());
    fDrawer->DrawHistogram(drawnHistogram,namerX.Data(),"#frac{dN}{d#eta}"," ");
    legend = new TLegend(0.4,0.20,0.82,0.35);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    namerX = Form("%sEta",fHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // Change the right margin better suited for 2D-drawing
    fDrawer->SetRightMargin(0.1);
    
    // === Jet eta vs. phi ===
    drawnHistogram2D = fHistograms->GetHistogramJetEtaPhi(iCentrality);
    namerX = Form("%s #varphi",fHistograms->GetJetAxisName());
    namerY = Form("%s #eta",fHistograms->GetJetAxisName());
    fDrawer->DrawHistogram(drawnHistogram2D,namerX.Data(),namerY.Data()," ",fStyle2D);
    legend = new TLegend(0.17,0.78,0.37,0.93);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figures to file
    namerX = Form("%sEtaPhi",fHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // Change right margin back to 1D-drawing
    fDrawer->SetRightMargin(0.06);
    
  } // Centrality loop
}

/*
 * Draw the track histograms
 */
void EECDrawer::DrawTrackHistograms(){
  
  // Helper variables for histogram drawing
  TH1D* drawnHistogram;
  TH2D* drawnHistogram2D;
  TLegend* legend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString namerX;
  TString namerY;
  
  // Number of events for normalization
  int numberOfEvents = fHistograms->GetNEvents();  // Normalize with the number of all events for inclusive histograms. TODO: Fix the normalization, currently not reliable
  
  // Loop over track types
  for(int iTrackType = 0; iTrackType < EECHistogramManager::knTrackCategories; iTrackType++){
    if(!fDrawTracks[iTrackType]) continue;  // Only draw selected tracks
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCard()->GetLowBinBorderCentrality(iCentrality),fHistograms->GetCard()->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCard()->GetLowBinBorderCentrality(iCentrality),fHistograms->GetCard()->GetHighBinBorderCentrality(iCentrality));

      // Select logarithmic drawing for pT
      fDrawer->SetLogY(fLogPt);
      
      // === Track pT ===
      drawnHistogram = fHistograms->GetHistogramTrackPt(iTrackType,iCentrality);
      drawnHistogram->Scale(1.0/numberOfEvents);  // Normalize with the number of events
      namerX = Form("%s p_{T}  (GeV)",fHistograms->GetTrackAxisName(iTrackType));
      fDrawer->DrawHistogram(drawnHistogram,namerX.Data(),"#frac{1}{N_{event}} #frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      namerX = Form("%sPt",fHistograms->GetTrackHistogramName(iTrackType));
      SaveFigure(namerX,compactCentralityString);
      
      // Select linear drawing
      fDrawer->SetLogY(false);
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fHistograms->GetNTrackPtBins(); iTrackPt++){
        
        // Draw the selected track pT bins and the special bin containing integrated distributions
        if(iTrackPt > fLastDrawnTrackPtBin && iTrackPt != fHistograms->GetNTrackPtBins()) continue;
        
        // Set the correct track pT bins
        trackPtString = Form("Track pT: %.1f-%.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString = Form("_pT=%.1f-%.1f",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString.ReplaceAll(".","v");
        
        // No pT selection for integrated distributions
        if(iTrackPt == fHistograms->GetNTrackPtBins()){
          trackPtString = "";
          compactTrackPtString = "";
        }
        
        // === Track phi ===
        drawnHistogram = fHistograms->GetHistogramTrackPhi(iTrackType,iCentrality,iTrackPt);
        drawnHistogram->Scale(1.0/numberOfEvents);  // Normalize with the number of events
        namerX = Form("%s #varphi",fHistograms->GetTrackAxisName(iTrackType));
        fDrawer->DrawHistogram(drawnHistogram,namerX.Data(),"#frac{1}{N_{event}} #frac{dN}{d#varphi}"," ");
        legend = new TLegend(0.17,0.20,0.37,0.35);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        namerX = Form("%sPhi",fHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
        // === Track eta ===
        drawnHistogram = fHistograms->GetHistogramTrackEta(iTrackType,iCentrality,iTrackPt);
        drawnHistogram->Scale(1.0/numberOfEvents);  // Normalize with the number of events
        
        // Set nice position for the legend
        legendY1 = 0.20; legendY2 = legendY1+0.15;
        if(iTrackPt == fHistograms->GetNTrackPtBins()){
          legendX1 = 0.4; legendX2 = legendX1+0.2;
        } else if (iTrackPt == fHistograms->GetNTrackPtBins() - 1){
          legendX1 = 0.32; legendX2 = legendX1+0.2;
        } else {
          legendX1 = 0.34; legendX2 = legendX1+0.2;
        }
        
        namerX = Form("%s #eta",fHistograms->GetTrackAxisName(iTrackType));
        fDrawer->DrawHistogram(drawnHistogram,namerX.Data(),"#frac{1}{N_{event}} #frac{dN}{d#eta}"," ");
        legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        namerX = Form("%sEta",fHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
        // Change the right margin better suited for 2D-drawing
        fDrawer->SetRightMargin(0.13);
        fDrawer->SetLeftMargin(0.12);
        fDrawer->SetTitleOffsetY(1.0);
        
        // === Track eta-phi ===
        drawnHistogram2D = fHistograms->GetHistogramTrackEtaPhi(iTrackType,iCentrality,iTrackPt);
        drawnHistogram2D->Scale(1.0/numberOfEvents);  // Normalize with the number of events
        namerX = Form("%s #varphi",fHistograms->GetTrackAxisName(iTrackType));
        namerY = Form("%s #eta",fHistograms->GetTrackAxisName(iTrackType));
        fDrawer->DrawHistogram(drawnHistogram2D,namerX.Data(),namerY.Data()," ",fStyle2D);
        legend = new TLegend(0.17,0.78,0.37,0.93);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        namerX = Form("%sEtaPhi",fHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
        // Change right margin back to 1D-drawing
        fDrawer->SetRightMargin(0.06);
        fDrawer->SetLeftMargin(0.15);
        fDrawer->SetTitleOffsetY(1.1);
        
      } // Track pT loop
    } // Centrality loop
  } // Track type loop
}

/*
 * Draw multiplicity in the jet cone histograms
 */
void EECDrawer::DrawMultiplicityInJetCone(){

  // Helper variables for histogram drawing
  TH1D* drawnHistogram;
  TLegend* legend;

  // Helper variables for naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString jetPtString;
  TString compactJetPtString;
  TString subeventString;
  TString namerX;

  // Loop over multiplicity types
  for(int iMultiplicityType = 0; iMultiplicityType < EECHistogramManager::knMultiplicityInJetConeTypes; iMultiplicityType++){

    // Only draw the selected multiplicity types
    if(!fDrawMultiplicityInJetCone[iMultiplicityType]) continue;

    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){

      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));

      // Loop over different subevent combinations
      for(int iSubevent = 0; iSubevent <= EECHistograms::knSubeventTypes; iSubevent++){
        
        // Only draw selected subevent combinations
        if(!fDrawSubeventType[iSubevent]) continue;
        
        subeventString = fHistograms->GetSubeventType(iSubevent);
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Loop over jet pT bins
          for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
            
            // Set the jet pT information for legends and figure saving
            if(iJetPt == fHistograms->GetNJetPtBinsEEC()){
              jetPtString = Form("jet p_{T} > %.0f GeV", fHistograms->GetCard()->GetJetPtCut());
              compactJetPtString = "";
            } else {
              jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
              compactJetPtString = Form("_J=%.0f-%.0f", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
            }
            
            // === Multiplicity within the jet ===
            drawnHistogram = fHistograms->GetHistogramMultiplicityInJetCone(iCentrality, iJetPt, iTrackPt, iMultiplicityType, iSubevent);
            namerX = Form("%s",fHistograms->GetMultiplicityInJetConeAxisName(iMultiplicityType));
            fDrawer->DrawHistogram(drawnHistogram,namerX.Data(),"Counts"," ");
            legend = new TLegend(0.62,0.75,0.82,0.9);
            SetupLegend(legend,centralityString,subeventString,jetPtString,trackPtString);
            legend->Draw();
            
            // Save the figure to a file
            namerX = Form("%s%s",fHistograms->GetMultiplicityInJetConeHistogramName(iMultiplicityType), subeventString.Data());
            SaveFigure(namerX,compactCentralityString, compactJetPtString, compactTrackPtString);
            
          } // Jet pT loop
        } // Track pT loop
      } // Subevent loop
    } // Centrality loop
  } // Multiplicity type loop
}

/*
 * Draw histograms from particle density around the jet axis
 */
void EECDrawer::DrawParticleDensityAroundJetAxis(){
  
  // Helper variables for histogram drawing
  TH1D* drawnHistogram;
  TLegend* legend;

  // Helper variables for naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString jetPtString;
  TString compactJetPtString;
  TString subeventString;
  TString jetConeTypeString;
  TString namerY;
  
  // Helper variables for histogram normalization
  double normalizationFactor = 1;
  int color[10] = {kBlack, kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure-1, kOrange-1, kGray};

  // Loop over particle density around the jet axis types
  for(int iParticleDensityType = 0; iParticleDensityType < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensityType++){

    // Only draw the selected particle density around the jet axis types
    if(!fDrawParticleDensityAroundJets[iParticleDensityType]) continue;

    for(int iJetConeType = 0; iJetConeType < EECHistograms::knJetConeTypes; iJetConeType++){
      
      jetConeTypeString = fHistograms->GetJetConeTypeSaveName(iJetConeType);
      
      // Loop over centrality
      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
        
        centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        
        // Loop over different subevent combinations
        for(int iSubevent = 0; iSubevent <= EECHistograms::knSubeventTypes; iSubevent++){
          
          // Only draw selected subevent combinations
          if(!fDrawSubeventType[iSubevent]) continue;
          
          subeventString = fHistograms->GetSubeventType(iSubevent);
          
          if(fDrawIndividualParticleDensities){
            
            // Loop over track pT bins
            for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
              
              // Track pT binning is different depending on the particle density type
              if(iParticleDensityType == EECHistogramManager::kParticleDensityAroundJetAxisPtBinned || iParticleDensityType == EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned){
                trackPtString = Form("%.1f < p_{T}^{ch} < %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt), fHistograms->GetTrackPtBinBorderEEC(iTrackPt+1));
                compactTrackPtString = Form("_T%.1f-%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt), fHistograms->GetTrackPtBinBorderEEC(iTrackPt+1));
                compactTrackPtString.ReplaceAll(".","v");
              } else {
                trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
                compactTrackPtString = Form("_T>%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
                compactTrackPtString.ReplaceAll(".","v");
              }
              
              // Loop over jet pT bins
              for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
                
                // Set the jet pT information for legends and figure saving
                if(iJetPt == fHistograms->GetNJetPtBinsEEC()){
                  jetPtString = Form("jet p_{T} > %.0f GeV", fHistograms->GetCard()->GetJetPtCut());
                  compactJetPtString = "";
                  normalizationFactor = fHistograms->GetJetPtIntegral(iCentrality);
                } else {
                  jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
                  compactJetPtString = Form("_J=%.0f-%.0f", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
                  normalizationFactor = fHistograms->GetJetPtIntegral(iCentrality, fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
                }
                
                // === Particle density around the jet axis ===
                drawnHistogram = fHistograms->GetHistogramParticleDensityAroundJetAxis(iCentrality, iJetPt, iTrackPt, iJetConeType, iParticleDensityType, iSubevent);
                drawnHistogram->Scale(1/normalizationFactor);
                namerY = Form("#frac{1}{N_{jets}} %s",fHistograms->GetParticleDensityAroundJetAxisAxisName(iParticleDensityType));
                fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY.Data()," ");
                legend = new TLegend(0.62,0.75,0.82,0.9);
                SetupLegend(legend,centralityString,jetConeTypeString,subeventString,jetPtString,trackPtString);
                legend->Draw();
                
                // Save the figure to a file
                namerY = Form("%s%s%s",fHistograms->GetParticleDensityAroundJetAxisHistogramName(iParticleDensityType), jetConeTypeString.Data(), subeventString.Data());
                SaveFigure(namerY,compactCentralityString, compactJetPtString, compactTrackPtString);
                
              } // Jet pT loop
            } // Track pT loop
          } // Draw individual particle density histograms
          
          // Draw all track pT cuts in the same plot with the jet pT integrated histogram TODO: Automatic scaling for y-axes, add jet pT bins
          if(fDrawParticleDensitiesForConstantJetPt){
            
            // Loop over jet pT bins
            for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
              
              // Set the jet pT information for legends and figure saving
              if(iJetPt == fHistograms->GetNJetPtBinsEEC()){
                jetPtString = Form("jet p_{T} > %.0f GeV", fHistograms->GetCard()->GetJetPtCut());
                compactJetPtString = Form("%.0f", fHistograms->GetCard()->GetJetPtCut());
              } else {
                jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
                compactJetPtString = Form("_J=%.0f-%.0f", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
              }
              
              // Only one legend for the plot
              legend = new TLegend(0.62,0.35,0.82,0.9);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
              if(iSubevent < EECHistograms::knSubeventTypes) legend->AddEntry((TObject*) 0, subeventString.Data(), "");
              if(jetConeTypeString != "") legend->AddEntry((TObject*) 0, jetConeTypeString.Data(), "");
              legend->AddEntry((TObject*) 0, centralityString.Data(), "");
              legend->AddEntry((TObject*) 0, jetPtString.Data(), "");
              
              
              // Loop over track pT bins
              for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
                
                // Track pT binning is different depending on the particle density type
                if(iParticleDensityType == EECHistogramManager::kParticleDensityAroundJetAxisPtBinned || iParticleDensityType == EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned){
                  trackPtString = Form("%.1f < p_{T}^{ch} < %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt), fHistograms->GetTrackPtBinBorderEEC(iTrackPt+1));
                } else {
                  trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
                }
                
                drawnHistogram = fHistograms->GetHistogramParticleDensityAroundJetAxis(iCentrality, iJetPt, iTrackPt, iJetConeType, iParticleDensityType, iSubevent);
                drawnHistogram->Scale(1/drawnHistogram->Integral(1, drawnHistogram->FindBin(0.3999), "width")); // To compare shapes, just normalize everything to one within 0 < DeltaR < 0.4
                drawnHistogram->SetLineColor(color[iTrackPt]);
                
                if(iTrackPt == fFirstDrawnTrackPtBinEEC){
                  namerY = Form("#frac{1}{N_{jets}} %s",fHistograms->GetParticleDensityAroundJetAxisAxisName(iParticleDensityType));
                  fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY.Data()," ");
                } else {
                  drawnHistogram->Draw("same");
                }
                
                legend->AddEntry(drawnHistogram, trackPtString.Data(), "l");
                
              } // Track pT loop
              
              legend->Draw();
              
              // Save the figure to a file
              namerY = Form("%s%s%sConstantJetPt%s", fHistograms->GetParticleDensityAroundJetAxisHistogramName(iParticleDensityType), jetConeTypeString.Data(), subeventString.Data(), compactJetPtString.Data());
              SaveFigure(namerY, compactCentralityString);
              
            } // Jet pT loop
          } // Draw all track pT cuts in the same plot with a constant jet pT histogram
          
        } // Subevent loop
      } // Centrality loop
    } // Jet cone tpye loop
  } // Particle density type loop
}

/*
 * Draw the maximum particle pT within the jet cone histograms
 */
void EECDrawer::DrawMaxParticlePtWithinJetCone(){
  
  // Helper variables for histogram drawing
  TH1D* drawnHistogram;
  TLegend* legend;

  // Helper variables for naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString trackPtCutString;
  TString compactTrackPtCutString;
  TString jetPtString;
  TString compactJetPtString;
  TString namer;
  TString trackTypeString[] = {"background", "track"};
  
  // Loop over maximum particle pT within the jet cone types
  for(int iMaxParticlePtWithinJetCone = 0; iMaxParticlePtWithinJetCone < EECHistogramManager::knMaxParticlePtWithinJetConeTypes; iMaxParticlePtWithinJetCone++){
    
    // Only draw the selected maximum particle pT within the jet cone types
    if(!fDrawMaxParticlePtWithinJetCone[iMaxParticlePtWithinJetCone]) continue;
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        
        // Loop over jet pT bins
        for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == fHistograms->GetNJetPtBinsEEC()){
            jetPtString = Form("jet p_{T} > %.0f GeV", fHistograms->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          // Draw the histogram without any track pT selections
          drawnHistogram = fHistograms->GetHistogramMaxParticlePtInJetCone(iMaxParticlePtWithinJetCone, iCentrality, iJetPt);
          drawnHistogram->Scale(1/drawnHistogram->Integral("width")); // Normalize to the number of jets
          fDrawer->DrawHistogram(drawnHistogram,fHistograms->GetMaxParticlePtWithinJetConeAxisName(iMaxParticlePtWithinJetCone),"#frac{1}{N_{jets}} counts"," ");
          legend = new TLegend(0.45,0.7,0.82,0.9);
          SetupLegend(legend,centralityString,jetPtString);
          legend->Draw();
          
          // Save the figure to a file
          namer = fHistograms->GetMaxParticlePtWithinJetConeSaveName(iMaxParticlePtWithinJetCone);
          SaveFigure(namer,compactCentralityString, compactJetPtString);
          
          // Loop over track pT bins
          for(int iTrackPt = 0; iTrackPt < EECHistogramManager::knProjectedMaxParticlePtBins; iTrackPt++){
            
            // Two different track pT binnings are projected from the histograms
            trackPtString = Form("%.1f < %s p_{T} < %.1f GeV",fHistograms->GetMaxTrackPtWithinJetConeBinBorder(iTrackPt), trackTypeString[iMaxParticlePtWithinJetCone].Data(), fHistograms->GetMaxTrackPtWithinJetConeBinBorder(iTrackPt+1));
            compactTrackPtString = Form("_T%.1f-%.1f",fHistograms->GetMaxTrackPtWithinJetConeBinBorder(iTrackPt), fHistograms->GetMaxTrackPtWithinJetConeBinBorder(iTrackPt+1));
            compactTrackPtString.ReplaceAll(".","v");
            trackPtCutString = Form("%.1f GeV < %s p_{T}",fHistograms->GetMaxTrackPtWithinJetConeBinBorder(iTrackPt), trackTypeString[iMaxParticlePtWithinJetCone].Data());
            compactTrackPtCutString = Form("_T>%.1f",fHistograms->GetMaxTrackPtWithinJetConeBinBorder(iTrackPt));
            compactTrackPtCutString.ReplaceAll(".","v");
            
            // === Track pT in bins ===
            drawnHistogram = fHistograms->GetHistogramMaxParticlePtInJetCone(iMaxParticlePtWithinJetCone, iCentrality, iJetPt, iTrackPt);
            drawnHistogram->Scale(1/drawnHistogram->Integral("width"));
            fDrawer->DrawHistogram(drawnHistogram,fHistograms->GetMaxParticlePtWithinJetConeAxisName(iMaxParticlePtWithinJetCone),"#frac{1}{N_{jets}} counts"," ");
            legend = new TLegend(0.45,0.65,0.82,0.9);
            SetupLegend(legend,centralityString,jetPtString,trackPtString);
            legend->Draw();
            
            // Save the figure to a file
            namer = fHistograms->GetMaxParticlePtWithinJetConeSaveName(iMaxParticlePtWithinJetCone);
            SaveFigure(namer, compactCentralityString, compactJetPtString, compactTrackPtString);
            
            // === Track pT in cuts ===
            drawnHistogram = fHistograms->GetHistogramMaxParticlePtInJetConePtCut(iMaxParticlePtWithinJetCone, iCentrality, iJetPt, iTrackPt);
            drawnHistogram->Scale(1/drawnHistogram->Integral("width"));
            fDrawer->DrawHistogram(drawnHistogram,fHistograms->GetMaxParticlePtWithinJetConeAxisName(iMaxParticlePtWithinJetCone),"#frac{1}{N_{jets}} counts"," ");
            legend = new TLegend(0.45,0.65,0.82,0.9);
            SetupLegend(legend,centralityString,jetPtString,trackPtCutString);
            legend->Draw();
            
            // Save the figure to a file
            namer = fHistograms->GetMaxParticlePtWithinJetConeSaveName(iMaxParticlePtWithinJetCone);
            SaveFigure(namer, compactCentralityString, compactJetPtString, compactTrackPtCutString);
            
          } // Track pT loop
        } // Jet pT loop
    } // Centrality loop
  } // Maximum particle pT within the jet cone type loop
}

/*
 * Draw energy-energy correlator histograms
 */
void EECDrawer::DrawEnergyEnergyCorrelationHistograms(){
  
  // Helper variables for histogram drawing
  TH1D* drawnHistogram;
  TLegend* legend;
  TLegend* anotherLegend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString jetPtString;
  TString compactJetPtString;
  TString subeventString;
  TString compactSubeventString;
  TString namerY;

  double legendY1;
  
  int colorAdder = 1;
  
  double normalizationFactor;
  int color[10] = {kBlack, kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure-1, kOrange-1, kGray};
  int style[10] = {kFullSquare, kFullCircle, kFullDiamond, kFullCross, kOpenFourTrianglesPlus, kOpenCrossX, kOpenStar, kOpenTriangleUp, kOpenTriangleDown, kOpenDoubleDiamond};
  
  // Set logarithmic drawing for deltaR and EEC
  fDrawer->SetLogX(fLogDeltaR);
  fDrawer->SetLogY(fLogEEC);
    
  // Loop over energy-energy correlator types
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only draw the selected histograms
    if(!fDrawEnergyEnergyCorrelators[iEnergyEnergyCorrelator]) continue;
    
    // Loop over different pairing type (same jet/reflected cone)
    for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
      
      // Only draw selected pairing types
      if(!fDrawPairingType[iPairingType]) continue;
      
      // Loop over centrality
      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
        
        centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        //centralityString = "Cent: 0-10%";
        compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        
        // Draw individual energy-energy correlator histograms
        if(fDrawIndividualEnergyEnergyCorrelators){
          
          // Loop over different subevent combinations
          for(int iSubevent = 0; iSubevent <= EECHistograms::knSubeventCombinations; iSubevent++){
            
            // Only draw selected subevent combinations
            if(!fDrawSubeventCombination[iSubevent]) continue;
            
            subeventString = fHistograms->GetSubeventCombination(iSubevent);
            compactSubeventString = fHistograms->GetSubeventCombinationSaveName(iSubevent);

            // Loop over track pT bins
            for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
              
              trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString = Form("_T>%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString.ReplaceAll(".","v");
              
              // Loop over jet pT bins
              for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
                
                jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
                compactJetPtString = Form("_J=%.0f-%.0f",fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
                
                // === Energy-energy correlator ===
                drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType, iSubevent);
                drawnHistogram->Scale(1/drawnHistogram->Integral("width")); // For now, just normalize the integral to one
                namerY = Form("%s %s", fHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), fHistograms->GetPairingTypeSaveName(iPairingType));
                
                // For logarithmic x-axis, cannot go all the way to zero
                if(fLogDeltaR){
                  double jetR = fHistograms->GetCard()->GetJetRadius();
                  if(jetR > 0.6){                                              // R = 0.8
                    drawnHistogram->GetXaxis()->SetRangeUser(0.008, 1.6);
                  } else {                                                     // R = 0.4 (default)
                    drawnHistogram->GetXaxis()->SetRangeUser(0.008, 0.8);
                  }
                }                
                fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY.Data()," ");
                legend = new TLegend(0.62,0.7,0.82,0.9);
                legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
                legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
                if(iSubevent < EECHistograms::knSubeventCombinations) legend->AddEntry((TObject*) 0, subeventString.Data(), "");
                legend->AddEntry((TObject*) 0, centralityString.Data(), "");
                legend->AddEntry((TObject*) 0, jetPtString.Data(), "");
                legend->AddEntry((TObject*) 0, trackPtString.Data(), "");
                legend->Draw();
                
                // Save the figure to a file
                namerY = Form("%s%s%s", fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), fHistograms->GetPairingTypeSaveName(iPairingType), compactSubeventString.Data());
                SaveFigure(namerY, compactCentralityString, compactTrackPtString, compactJetPtString);
                
              } // Jet pT loop
            } // Track pT loop
          } // Subevent type loop
        } // Draw individual energy-energy correlator histograms
        
        // Draw all track pT cuts in the same plot with the jet pT integrated histogram TODO: Automatic scaling for y-axes, add jet pT bins
        if(fDrawEnergyEnergyCorrelatorsForConstantJetPt){
          
          jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetCard()->GetLowBinBorderJetPtEEC(fFirstDrawnJetPtBinEEC), fHistograms->GetCard()->GetHighBinBorderJetPtEEC(fFirstDrawnJetPtBinEEC));
          
          // Loop over different subevent combinations
          for(int iSubevent = 0; iSubevent <= EECHistograms::knSubeventCombinations; iSubevent++){
            
            // Only draw selected subevent combinations
            if(!fDrawSubeventCombination[iSubevent]) continue;
            
            subeventString = fHistograms->GetSubeventCombination(iSubevent);
            compactSubeventString = fHistograms->GetSubeventCombinationSaveName(iSubevent);
            
            // Only one legend for the plot
            legend = new TLegend(0.18,0.20,0.50,0.70);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
            if(iSubevent < EECHistograms::knSubeventCombinations) legend->AddEntry((TObject*) 0, subeventString.Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(), "");
            //legend->AddEntry((TObject*) 0, "Cent: 0-10%", "");
            legend->AddEntry((TObject*) 0, Form("anti-k_{T} R = %.1f", fHistograms->GetCard()->GetJetRadius()), "");
            legend->AddEntry((TObject*) 0, jetPtString.Data(), "");
            
            
            // Loop over track pT bins
            for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
              
              trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
              
              drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, fFirstDrawnJetPtBinEEC, iTrackPt, iPairingType, iSubevent);
              drawnHistogram->Scale(1/drawnHistogram->Integral("width"));
              drawnHistogram->SetLineColor(color[iTrackPt-fFirstDrawnTrackPtBinEEC]);
              drawnHistogram->SetMarkerColor(color[iTrackPt-fFirstDrawnTrackPtBinEEC]);
              drawnHistogram->SetMarkerStyle(style[iTrackPt-fFirstDrawnTrackPtBinEEC]);
              
              // For logarithmic x-axis, cannot go all the way to zero
              if(fLogDeltaR){
                double jetR = fHistograms->GetCard()->GetJetRadius();
                if(jetR > 0.6){                                              // R = 0.8
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 1.6);
                } else {                                                     // R = 0.4 (default)
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 0.8);
                }
              }


              if(iTrackPt == fFirstDrawnTrackPtBinEEC){
                namerY = "EEC";
                drawnHistogram->GetYaxis()->SetRangeUser(2e-5, drawnHistogram->GetMaximum()*5);
                fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY.Data()," ","p");
              } else {
                drawnHistogram->Draw("p,same");
              }
              
              legend->AddEntry(drawnHistogram, trackPtString.Data(), "p");
              
            } // Track pT loop
            
            legend->Draw();
            
            // Save the figure to a file
            namerY = Form("%s%s%sConstantJetPt%.0f", fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), fHistograms->GetPairingTypeSaveName(iPairingType), compactSubeventString.Data(), fHistograms->GetCard()->GetJetPtCut());
            SaveFigure(namerY, compactCentralityString);
            
          } // Subevent loop
          
        } // Draw all track pT cuts in the same plot with the jet pT integrated histogram
        
        // Draw all jet pT selections for a single track cut TODO: Automatic scaling for y-axes
        if(fDrawEnergyEnergyCorrelatorsForConstantTrackPt){
          
          // Loop over track pT bins
          for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
            
            trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Loop over different subevent combinations
            for(int iSubevent = 0; iSubevent <= EECHistograms::knSubeventCombinations; iSubevent++){
              
              // Only draw selected subevent combinations
              if(!fDrawSubeventCombination[iSubevent]) continue;
              
              subeventString = fHistograms->GetSubeventCombination(iSubevent);
              compactSubeventString = fHistograms->GetSubeventCombinationSaveName(iSubevent);
              
              // Only one legend for the plot
              legend = new TLegend(0.18,0.20,0.50,0.70);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
              if(iSubevent < EECHistograms::knSubeventCombinations) legend->AddEntry((TObject*) 0, subeventString.Data(), "");
              legend->AddEntry((TObject*) 0, centralityString.Data(),"");
              //legend->AddEntry((TObject*) 0, "Cent: 0-10%","");
              legend->AddEntry((TObject*) 0, Form("anti-k_{T} R = %.1f", fHistograms->GetCard()->GetJetRadius()), "");
              legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
              
              // Loop over jet pT bins
              for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
                
                jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
                
                drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType, iSubevent);
                drawnHistogram->Scale(1/drawnHistogram->Integral("width"));
                drawnHistogram->SetLineColor(color[iJetPt-fFirstDrawnJetPtBinEEC]);
                drawnHistogram->SetMarkerColor(color[iJetPt-fFirstDrawnJetPtBinEEC]);
                drawnHistogram->SetMarkerStyle(style[iJetPt-fFirstDrawnJetPtBinEEC]);
                
                // For logarithmic x-axis, cannot go all the way to zero
                if(fLogDeltaR){
                double jetR = fHistograms->GetCard()->GetJetRadius();
                if(jetR > 0.6){                                              // R = 0.8
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 1.6);
                } else {                                                     // R = 0.4 (default)
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 0.8);
                }
                }
                
                if(iJetPt == fFirstDrawnJetPtBinEEC){
                  //namerY = Form("%s %s", fHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), fHistograms->GetPairingTypeSaveName(iPairingType));
                  namerY = "EEC";
                  drawnHistogram->GetYaxis()->SetRangeUser(2e-4, drawnHistogram->GetMaximum()*5);
                  fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY.Data()," ","p");
                } else {
                  drawnHistogram->Draw("p,same");
                }
                
                legend->AddEntry(drawnHistogram, jetPtString.Data(), "p");
                
              } // Jet pT loop
              
              legend->Draw();
              
              // Save the figure to a file
              namerY = Form("%s%s%sConstantTrackPt%s",fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), fHistograms->GetPairingTypeSaveName(iPairingType), compactSubeventString.Data(), compactTrackPtString.Data());
              SaveFigure(namerY, compactCentralityString);
              
            } // Subevent loop
            
          } // Track pT loop
        } // Draw all jet pT selections for a single track cut
        
        // Draw subevent decomposition for the energy-energy correlators
        if(fDrawEnergyEnergyCorrelatorsSubevents){
          
          // Loop over track pT bins
          for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
            
            trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Loop over jet pT bins
            for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
              
              if(iJetPt == fHistograms->GetNJetPtBinsEEC()){
                jetPtString = Form("jet p_{T} > %.0f GeV", fHistograms->GetCard()->GetJetPtCut());
                compactJetPtString = "";
              } else {
                jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
                compactJetPtString = Form("_J=%.0f-%.0f",fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
              }
              
              // Only one legend for the plot
              legendY1 = (fLegendComment != "") ? 0.5 : 0.55;
              legend = new TLegend(0.22,legendY1,0.42,0.8);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              TString sysLabel = fHistograms->GetCard()->GetAlternativeDataType(false);
              TString eLabel   = sysLabel.Contains("OO") ? "5.36 TeV" : "5.02 TeV";
              legend->AddEntry((TObject*) 0, Form("%s %s", sysLabel.Data(), eLabel.Data()), "");
              if(fLegendComment != "") legend->AddEntry((TObject*) 0, fLegendComment.Data(), "");
              legend->AddEntry((TObject*) 0, centralityString.Data(),"");
              //legend->AddEntry((TObject*) 0, "Cent: 0-10%","");
              legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
              legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
              

              anotherLegend = new TLegend(0.62,0.2,0.82,0.45);
              anotherLegend->SetFillStyle(0);anotherLegend->SetBorderSize(0);anotherLegend->SetTextSize(0.05);anotherLegend->SetTextFont(62);
              
              // First, draw the total distribution to the canvas
              drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType);
              normalizationFactor = 1.0/drawnHistogram->Integral("width");
              drawnHistogram->Scale(normalizationFactor);
              drawnHistogram->SetLineColor(color[0]);
              drawnHistogram->SetMarkerColor(color[0]);
              drawnHistogram->SetMarkerStyle(style[0]);
              anotherLegend->AddEntry(drawnHistogram, "All combinations", "p");
              
              // For logarithmic x-axis, cannot go all the way to zero. Drawing range based on weight exponent.
              if(fHistograms->GetCard()->GetWeightExponent() == 1){
                drawnHistogram->GetYaxis()->SetRangeUser(0.001, 10);
              } else {
                drawnHistogram->GetYaxis()->SetRangeUser(0.00005, 60);
              }
              if(fLogDeltaR){
                double jetR = fHistograms->GetCard()->GetJetRadius();
                if(jetR > 0.6){                                              // R = 0.8
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 1.6);
                } else {                                                     // R = 0.4 (default)
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 0.8);
                }
              }
              //namerY = Form("%s %s", fHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), fHistograms->GetPairingTypeSaveName(iPairingType));
              namerY = "EEC";
              fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY.Data()," ");
              
              // Draw the different subevent contributions to the same canvas
              colorAdder = 1;
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
                if(iPairingType == EECHistograms::kSignalCone && iSubevent == EECHistograms::kHydjetPythia) {
                  colorAdder = 0;
                  continue;
                }
                drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType, iSubevent);
                drawnHistogram->Scale(normalizationFactor);
                drawnHistogram->SetLineColor(color[iSubevent+colorAdder]);
                drawnHistogram->SetMarkerColor(color[iSubevent+colorAdder]);
                drawnHistogram->SetMarkerStyle(style[iSubevent+colorAdder]);
                drawnHistogram->Draw("same");
                
                anotherLegend->AddEntry(drawnHistogram, fHistograms->GetSubeventCombination(iSubevent), "p");
              } // Subevent type loop
              
              // Draw the legend
              legend->Draw();
              anotherLegend->Draw();
              
              // Save the figure to a file
              namerY = Form("%s%sSubeventDecomposition",fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), fHistograms->GetPairingTypeSaveName(iPairingType));
              SaveFigure(namerY, compactCentralityString, compactTrackPtString, compactJetPtString);
              
            } // Jet pT loop
          } // Track pT bin loop
          
        } // Draw subevent decomposition for the energy-energy correlators
        
      } // Centrality loop
    } // Pairing type loop (same jet/reflected cone)
  } // Energy-energy correlator type loop
  
  // Reset logarithmic drawing flags
  fDrawer->SetLogX(false);
  fDrawer->SetLogY(false);
}

/*
 * Draw processed energy-energy correlators
 */
void EECDrawer::DrawProcessedEnergyEnergyCorrelators(){
  
  // Helper variables for histogram drawing
  TH1D* drawnHistogram;
  TLegend* legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString jetPtString;
  TString compactJetPtString;
  TString namerY;
  
  int color[10] = {kBlack, kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure-1, kOrange-1, kGray};
  
  // Set logarithmic drawing for deltaR and EEC
  fDrawer->SetLogX(fLogDeltaR);
  fDrawer->SetLogY(fLogEEC);
    
  // Loop over energy-energy correlator types
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only draw the selected histograms
    if(!fDrawEnergyEnergyCorrelators[iEnergyEnergyCorrelator]) continue;
    
    // Loop over different processing levels (normalized, background, signal)
    for(int iProcessLevel = 0; iProcessLevel < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels; iProcessLevel++){
            
      // Only draw selected processing levels
      if(!fDrawProcessingLevel[iProcessLevel]) continue;
            
      // Loop over centrality
      for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
        
        centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
        
        // Draw individual energy-energy correlator histograms
        if(fDrawIndividualEnergyEnergyCorrelators){
                    
          // Loop over track pT bins
          for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
            
            trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("_T>%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Loop over jet pT bins
            for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
              
              jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
              compactJetPtString = Form("_J=%.0f-%.0f",fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
              
              // === Energy-energy correlator ===
              drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iProcessLevel);
              namerY = Form("%s %s", fHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), fHistograms->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel));
              
              // For logarithmic x-axis, cannot go all the way to zero
              if(fLogDeltaR){
                double jetR = fHistograms->GetCard()->GetJetRadius();
                if(jetR > 0.6){                                              // R = 0.8
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 1.6);
                } else {                                                     // R = 0.4 (default)
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 0.8);
                }
              }
              
              fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY.Data()," ");
              legend = new TLegend(0.62,0.7,0.82,0.9);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
              legend->AddEntry((TObject*) 0, centralityString.Data(), "");
              legend->AddEntry((TObject*) 0, jetPtString.Data(), "");
              legend->AddEntry((TObject*) 0, trackPtString.Data(), "");
              legend->Draw();
              
              // Save the figure to a file
              namerY = Form("%s%s", fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), fHistograms->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel));
              SaveFigure(namerY, compactCentralityString, compactTrackPtString, compactJetPtString);
              
            } // Jet pT loop
          } // Track pT loop
        } // Draw individual energy-energy correlator histograms
        
        // Draw all track pT cuts in the same plot with the jet pT integrated histogram TODO: Automatic scaling for y-axes, add jet pT bins
        if(fDrawEnergyEnergyCorrelatorsForConstantJetPt){
          
          jetPtString = Form("jet p_{T} > %.0f GeV", fHistograms->GetCard()->GetJetPtCut());
          
          // Only one legend for the plot
          legend = new TLegend(0.62,0.35,0.82,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
          legend->AddEntry((TObject*) 0, centralityString.Data(), "");
          legend->AddEntry((TObject*) 0, jetPtString.Data(), "");
          
          // Loop over track pT bins
          for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
            
            trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
            
            drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, iCentrality, fHistograms->GetNJetPtBinsEEC(), iTrackPt, iProcessLevel);
            drawnHistogram->SetLineColor(color[iTrackPt]);
            
            // For logarithmic x-axis, cannot go all the way to zero
            if(fLogDeltaR){
              double jetR = fHistograms->GetCard()->GetJetRadius();
              if(jetR > 0.6){                                              // R = 0.8
                drawnHistogram->GetXaxis()->SetRangeUser(0.008, 1.6);
              } else {                                                     // R = 0.4 (default)
                drawnHistogram->GetXaxis()->SetRangeUser(0.008, 0.8);
              }
            }
            
            if(iTrackPt == fFirstDrawnTrackPtBinEEC){
              namerY = Form("%s %s", fHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), fHistograms->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel));
              fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY.Data()," ");
            } else {
              drawnHistogram->Draw("same");
            }
            
            legend->AddEntry(drawnHistogram, trackPtString.Data(), "l");
            
          } // Track pT loop
          
          legend->Draw();
          
          // Save the figure to a file
          namerY = Form("%s%sConstantJetPt%.0f", fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), fHistograms->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel), fHistograms->GetCard()->GetJetPtCut());
          SaveFigure(namerY, compactCentralityString);
          
        } // Draw all track pT cuts in the same plot with the jet pT integrated histogram
        
        // Draw all jet pT selections for a single track cut TODO: Automatic scaling for y-axes
        if(fDrawEnergyEnergyCorrelatorsForConstantTrackPt){
          
          // Loop over track pT bins
          for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
            
            trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Only one legend for the plot
            legend = new TLegend(0.18,0.20,0.50,0.70);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            
            // Loop over jet pT bins
            for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
              
              jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
              
              drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iProcessLevel);
              drawnHistogram->SetLineColor(color[iJetPt]);
              
              // For logarithmic x-axis, cannot go all the way to zero
              if(fLogDeltaR){
                double jetR = fHistograms->GetCard()->GetJetRadius();
                if(jetR > 0.6){                                              // R = 0.8
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 1.6);
                } else {                                                     // R = 0.4 (default)
                  drawnHistogram->GetXaxis()->SetRangeUser(0.008, 0.8);
                }
              }
              
              if(iJetPt == fFirstDrawnJetPtBinEEC){
                namerY = Form("%s %s", fHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), fHistograms->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel));
                fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY.Data()," ");
              } else {
                drawnHistogram->Draw("same");
              }
              
              legend->AddEntry(drawnHistogram, jetPtString.Data(), "l");
              
            } // Jet pT loop
            
            legend->Draw();
            
            // Save the figure to a file
            namerY = Form("%s%sConstantTrackPt%s",fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), fHistograms->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel), compactTrackPtString.Data());
            SaveFigure(namerY, compactCentralityString);
            
          } // Track pT loop
        } // Draw all jet pT selections for a single track cut
        
      } // Centrality loop
    } // Processing level (normalized/background/signal)
  } // Energy-energy correlator type loop
  
  // Reset logarithmic drawing flags
  fDrawer->SetLogX(false);
  fDrawer->SetLogY(false);
}

/*
 * Draw covariance matrices
 */
void EECDrawer::DrawCovarianceMatrices(){
  
  // If jet histograms are not drawn, there is nothing to do here
  if(!fDrawCovarianceMatrix) return;
  
  // Helper variables for histogram drawing
  TH2D* drawnHistogram;
  TLegend* legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString covarianceString;
  TString namerX;
  TString namerY;

  // Change the right margin better suited for 2D-drawing
  fDrawer->SetRightMargin(0.12);
  fDrawer->SetLeftMargin(0.12);
  fDrawer->SetTitleOffsetY(0.8);
  fDrawer->SetTitleOffsetX(1.2);
  
  // Loop over covariance matrix types
  for(int iCovarianceMatrixType = 0; iCovarianceMatrixType < EECHistogramManager::knCovarianceMatrixTypes; iCovarianceMatrixType++){
    covarianceString = fHistograms->GetJetPtUnfoldingCovarianceSaveName(iCovarianceMatrixType);

    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));

      for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){

        trackPtString = Form("p_{T}^{ch} > %.1f GeV",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
        compactTrackPtString = Form("%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".","v");
    
        // Select logarithmic z-axis scale
        fDrawer->SetLogZ(true);
    
        // === Covariance matrix ===
        drawnHistogram = fHistograms->GetHistogramJetPtUnfoldingCovariance(iCovarianceMatrixType, iCentrality, iTrackPt);
        if(drawnHistogram == NULL) continue; // Cannot draw the histogram if it does not exist
        namerX = "#Deltar #otimes jet p_{T}";
        namerY = "#Deltar #otimes jet p_{T}";
        fDrawer->DrawHistogram(drawnHistogram,namerX.Data(),namerY.Data()," ",fStyle2D);
        legend = new TLegend(0.17,0.7,0.37,0.9);
        SetupLegend(legend,covarianceString,centralityString,trackPtString);
        legend->Draw();
    
        // Save the figures to file
        SaveFigure(covarianceString,compactCentralityString,compactTrackPtString);
    
      } // Track pT loop
    
    } // Centrality loop
  }

  // Change right margin back to 1D-drawing
  fDrawer->SetTitleOffsetY(1.1);
  fDrawer->SetTitleOffsetX(1.1);
  fDrawer->SetLeftMargin(0.15);
  fDrawer->SetRightMargin(0.06);

  // Reset logarithmic drawing setting
  fDrawer->SetLogZ(false);
}

/*
 * Common legend style setup for figures
 *
 *  TLegend *legend = Pointer to legend that needs setup
 *  TString centralityString = Collision centrality
 *  TString jetString = Jet pT information
 *  TString trackString = Track pT information
 *  TString extraString = Additional line to be put into the legend
 *  TString anotherString = Another string to be put to the legend
 */
void EECDrawer::SetupLegend(TLegend *legend, TString centralityString, TString jetString, TString trackString, TString extraString, TString anotherString){
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62); // Size: 0.05
  legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
  if(fSystemAndEnergy.Contains("PbPb") || fSystemAndEnergy.Contains("Hydjet") || fSystemAndEnergy.Contains("OO")) legend->AddEntry((TObject*) 0,centralityString.Data(),""); 
  if(jetString != "") legend->AddEntry((TObject*) 0,jetString.Data(),"");
  if(trackString != "") legend->AddEntry((TObject*) 0,trackString.Data(),"");
  if(extraString != "") legend->AddEntry((TObject*) 0,extraString.Data(),"");
  if(anotherString != "") legend->AddEntry((TObject*) 0,anotherString.Data(),"");
}

/*
 * Save the figure in current canvas to a file
 *
 *  TString figureName = Name for the saved figures
 *  TString centralityString = Information about collision centrality
 *  TString trackPtString = Information about track pT
 *  TString correlationTypeString = Information about correlation type (same/mixed event)
 *  TString deltaPhiString = Information about deltaPhi binning
 */
void EECDrawer::SaveFigure(TString figureName, TString centralityString, TString trackPtString, TString correlationTypeString, TString deltaPhiString){
  
  // Only save the figures if flag is set
  if(!fSaveFigures) return;
  
  // Write the figure to a file
  TString figName = Form("figures/%s_%s",figureName.Data(),fCompactSystemAndEnergy.Data());
  if(fCompactSystemAndEnergy.Contains("PbPb") || fCompactSystemAndEnergy.Contains("Hydjet")) figName.Append(centralityString);
  figName.Append(trackPtString);
  figName.Append(correlationTypeString);
  figName.Append(deltaPhiString);
  figName.Append(fFigureSaveNameAppend);
  gPad->GetCanvas()->SaveAs(Form("%s.%s",figName.Data(),fFigureFormat));
  
}

// Setter for drawing event information
void EECDrawer::SetDrawEventInformation(const bool drawOrNot){
  fDrawEventInformation = drawOrNot;
}

// Setter for drawing all jet histograms
void EECDrawer::SetDrawJetHistograms(const bool drawOrNot){
  fDrawJets = drawOrNot;
}

// Setter for drawing tracks
void EECDrawer::SetDrawTracks(const bool drawOrNot){
  fDrawTracks[EECHistogramManager::kTrack] = drawOrNot;
}

// Setter for drawing uncorrected tracks
void EECDrawer::SetDrawTracksUncorrected(const bool drawOrNot){
  fDrawTracks[EECHistogramManager::kUncorrectedTrack] = drawOrNot;
}

// Setter for drawing track histograms
void EECDrawer::SetDrawAllTracks(const bool drawTracks, const bool drawUncorrected){
  SetDrawTracks(drawTracks);
  SetDrawTracksUncorrected(drawUncorrected);
}

// Setter for drawing multiplicity within the jet cone
void EECDrawer::SetDrawMultiplicityInJetCone(const bool drawOrNot){
  fDrawMultiplicityInJetCone[EECHistogramManager::kMultiplicityInJetCone] = drawOrNot;
}

// Setter for drawing multiplicity within the reflected cone
void EECDrawer::SetDrawMultiplicityInReflectedCone(const bool drawOrNot){
  fDrawMultiplicityInJetCone[EECHistogramManager::kMultiplicityInReflectedCone] = drawOrNot;
}

// Setter for drawing uncorrected multiplicity within the jet cone
void EECDrawer::SetDrawMultiplicityInJetConeUncorrected(const bool drawOrNot){
  fDrawMultiplicityInJetCone[EECHistogramManager::kMultiplicityInJetConeUncorrected] = drawOrNot;
}

// Setter for drawing uncorrected multiplicity within the reflected cone
void EECDrawer::SetDrawMultiplicityInReflectedConeUncorrected(const bool drawOrNot){
  fDrawMultiplicityInJetCone[EECHistogramManager::kMultiplicityInReflectedConeUncorrected] = drawOrNot;
}

// Setter for drawing all multiplicity histograms within the jet cone
void EECDrawer::SetDrawAllMultiplicitiesInJetCone(const bool regular, const bool reflectedCone, const bool regularUncorrected, const bool reflectedConeUncorrected){
  SetDrawMultiplicityInJetCone(regular);
  SetDrawMultiplicityInReflectedCone(reflectedCone);
  SetDrawMultiplicityInJetConeUncorrected(regularUncorrected);
  SetDrawMultiplicityInReflectedConeUncorrected(reflectedConeUncorrected);
}

// Setter for drawing particle density around jet axis
void EECDrawer::SetDrawParticleDensityAroundJetAxis(const bool drawOrNot){
  fDrawParticleDensityAroundJets[EECHistogramManager::kParticleDensityAroundJetAxis] = drawOrNot;
}

// Setter for drawing particle pT density around jet axis
void EECDrawer::SetDrawParticlePtDensityAroundJetAxis(const bool drawOrNot){
  fDrawParticleDensityAroundJets[EECHistogramManager::kParticlePtDensityAroundJetAxis] = drawOrNot;
}

// Setter for drawing pT binned particle density around jet axis
void EECDrawer::SetDrawParticleDensityAroundJetAxisPtBinned(const bool drawOrNot){
  fDrawParticleDensityAroundJets[EECHistogramManager::kParticleDensityAroundJetAxisPtBinned] = drawOrNot;
}

// Setter for drawing pT binned particle pT density around jet axis
void EECDrawer::SetDrawParticlePtDensityAroundJetAxisPtBinned(const bool drawOrNot){
  fDrawParticleDensityAroundJets[EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned] = drawOrNot;
}

// Setter for drawing all particle densities around jet axis
void EECDrawer::SetDrawAllParticleDensitiesAroundJetAxis(const bool drawRegular, const bool drawPt, const bool drawPtBinned, const bool drawPtWeightedPtBinned){
  SetDrawParticleDensityAroundJetAxis(drawRegular);
  SetDrawParticlePtDensityAroundJetAxis(drawPt);
  SetDrawParticleDensityAroundJetAxisPtBinned(drawPtBinned);
  SetDrawParticlePtDensityAroundJetAxisPtBinned(drawPtWeightedPtBinned);
}

// Setter for drawing the individual particle density histograms
void EECDrawer::SetDrawSingleParticleDensityHistograms(const bool drawOrNot){
  fDrawIndividualParticleDensities = drawOrNot;
}

// Setter for drawing all track pT cuts to the same figure for constant jet pT selection
void EECDrawer::SetDrawParticleDensityForConstantJetPt(const bool drawOrNot){
  fDrawParticleDensitiesForConstantJetPt = drawOrNot;
}

// Setter for drawing the maximum particle pT within the jet cone histograms
void EECDrawer::SetDrawMaxParticlePtWithinJetCone(const bool drawOrNot){
  fDrawMaxParticlePtWithinJetCone[EECHistogramManager::kMaxSignalParticlePt] = drawOrNot;
}

// Setter for drawing the maximum background particle pT within the jet cone histograms
void EECDrawer::SetDrawMaxBackgroundParticlePtWithinJetCone(const bool drawOrNot){
  fDrawMaxParticlePtWithinJetCone[EECHistogramManager::kMaxBackgroundParticlePt] = drawOrNot;
}

// Setter for drawing energy-energy correlator
void EECDrawer::SetDrawEnergyEnergyCorrelor(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelator] = drawOrNot;
}

// Setter for drawing energy-energy correlators without single track efficiency corrections
void EECDrawer::SetDrawEnergyEnergyCorrelorEfficiencyVariationPlus(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus] = drawOrNot;
}

// Setter for drawing energy-energy correlators without single and pair track efficiency corrections
void EECDrawer::SetDrawEnergyEnergyCorrelorEfficiencyVariationMinus(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus] = drawOrNot;
}

// Setter for drawing energy-energy correlators with positive track pair efficiency variation
void EECDrawer::SetDrawEnergyEnergyCorrelorPairEfficiencyVariationPlus(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus] = drawOrNot;
}

// Setter for drawing energy-energy correlators with negative track pair efficiency variation
void EECDrawer::SetDrawEnergyEnergyCorrelorPairEfficiencyVariationMinus(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus] = drawOrNot;
}

// Setter for drawing all energy-energy correlators
void EECDrawer::SetDrawAllEnergyEnergyCorrelors(const bool drawRegular, const bool drawEfficiencyVariationPlus, const bool drawEfficiencyVariationMinus, const bool drawPairEfficiencyVariationPlus, const bool drawPairEfficiencyVariationMinus){
  SetDrawEnergyEnergyCorrelor(drawRegular);
  SetDrawEnergyEnergyCorrelorEfficiencyVariationPlus(drawEfficiencyVariationPlus);
  SetDrawEnergyEnergyCorrelorEfficiencyVariationMinus(drawEfficiencyVariationMinus);
  SetDrawEnergyEnergyCorrelorPairEfficiencyVariationPlus(drawPairEfficiencyVariationPlus);
  SetDrawEnergyEnergyCorrelorPairEfficiencyVariationMinus(drawPairEfficiencyVariationMinus);
}

// Setter for drawing the individual energy-energy correlator histograms
void EECDrawer::SetDrawSingleEnergyEnergyCorrelators(const bool drawOrNot){
  fDrawIndividualEnergyEnergyCorrelators = drawOrNot;
}

// Setter for drawing all track pT cuts to the same figure for constant jet pT selection
void EECDrawer::SetDrawEnergyEnergyCorrelatorsForConstantJetPt(const bool drawOrNot){
  fDrawEnergyEnergyCorrelatorsForConstantJetPt = drawOrNot;
}

// Setter for drawing all jet pT selections to the same figure for constant track pT cut
void EECDrawer::SetDrawEnergyEnergyCorrelatorsForConstantTrackPt(const bool drawOrNot){
  fDrawEnergyEnergyCorrelatorsForConstantTrackPt = drawOrNot;
}

// Setter for drawing subevent decomposition for energy-energy correlators
void EECDrawer::SetDrawEnergyEnergyCorrelatorsSubevent(const bool drawOrNot){
  fDrawEnergyEnergyCorrelatorsSubevents = drawOrNot;
}

// Setter for drawing same jet energy-energy correlators
void EECDrawer::SetDrawSameJetEnergyEnergyCorrelators(const bool drawOrNot){
  fDrawPairingType[EECHistograms::kSameJetPair] = drawOrNot;
}

// Setter for drawing signal-reflected cone energy-energy correlators
void EECDrawer::SetDrawSignalReflectedConeEnergyEnergyCorrelators(const bool drawOrNot){
  fDrawPairingType[EECHistograms::kSignalReflectedConePair] = drawOrNot;
}

// Setter for drawing reflected cone-reflected cone energy-energy correlators
void EECDrawer::SetDrawReflectedConeOnlyEnergyEnergyCorrelators(const bool drawOrNot){
  fDrawPairingType[EECHistograms::kReflectedConePair] = drawOrNot;
}

// Setter for drawing all different energy-energy correlator pairing types
void EECDrawer::SetDrawAllEnergyEnergyCorrelatorPairingTypes(const bool drawSameJet, const bool drawSignalReflectedCone, const bool drawReflectedConeOnly){
  SetDrawSameJetEnergyEnergyCorrelators(drawSameJet);
  SetDrawSignalReflectedConeEnergyEnergyCorrelators(drawSignalReflectedCone);
  SetDrawReflectedConeOnlyEnergyEnergyCorrelators(drawReflectedConeOnly);
}

// Setter for drawing all jet histograms
void EECDrawer::SetDrawCovarianceMatrices(const bool drawOrNot){
  fDrawCovarianceMatrix = drawOrNot;
}

// Setter for drawing histograms without subevent selection
void EECDrawer::SetDrawAllSubevents(const bool drawOrNot){
  fDrawSubeventType[EECHistograms::knSubeventTypes] = drawOrNot;
}

// Setter for drawing only Pythia histograms
void EECDrawer::SetDrawPythiaOnly(const bool drawOrNot){
  fDrawSubeventType[EECHistograms::kPythia] = drawOrNot;
}

// Setter for drawing only Hydjet histograms
void EECDrawer::SetDrawHydjetOnly(const bool drawOrNot){
  fDrawSubeventType[EECHistograms::kHydjet] = drawOrNot;
}

// Setter for drawing all subevent types
void EECDrawer::SetDrawAllSubeventTypes(const bool drawAll, const bool drawPythia, const bool drawHydjet){
  SetDrawAllSubevents(drawAll);
  SetDrawPythiaOnly(drawPythia);
  SetDrawHydjetOnly(drawHydjet);
}

// Setter for drawing all pairing combinations
void EECDrawer::SetDrawAllCombinations(const bool drawOrNot){
  fDrawSubeventCombination[EECHistograms::knSubeventCombinations] = drawOrNot;
}

// Setter for drawing Pythia+Pythia correlations for simulation
void EECDrawer::SetDrawSignalOnly(const bool drawOrNot){
  fDrawSubeventCombination[EECHistograms::kPythiaPythia] = drawOrNot;
}

// Setter for drawing Pythia+Hydjet correlations for simulation
void EECDrawer::SetDrawSignalFake(const bool drawOrNot){
  fDrawSubeventCombination[EECHistograms::kPythiaHydjet] = drawOrNot;
}

// Setter for drawing Hydjet+Hydjet correlations for simulation
void EECDrawer::SetDrawFakeFake(const bool drawOrNot){
  fDrawSubeventCombination[EECHistograms::kHydjetHydjet] = drawOrNot;
}

// Setter for drawing all subevent combinations for energy-energy correlators
void EECDrawer::SetDrawAllSubeventCombinations(const bool drawAll, const bool drawSignal, const bool drawSignalFake, const bool drawFakeFake){
  SetDrawAllCombinations(drawAll);
  SetDrawSignalOnly(drawSignal);
  SetDrawSignalFake(drawSignalFake);
  SetDrawFakeFake(drawFakeFake);
}

// Setter for drawing normalized energy-energy correlators
void EECDrawer::SetDrawEnergyEnergyCorrelatorNormalized(const bool drawOrNot){
  fDrawProcessingLevel[EECHistogramManager::kEnergyEnergyCorrelatorNormalized] = drawOrNot;
}

// Setter for drawing the normalized background estimate for energy-energy correlators
void EECDrawer::SetDrawEnergyEnergyCorrelatorBackground(const bool drawOrNot){
  fDrawProcessingLevel[EECHistogramManager::kEnergyEnergyCorrelatorBackground] = drawOrNot;
}

// Setter for drawing the background subtracted energy-energy correlators
void EECDrawer::SetDrawEnergyEnergyCorrelatorSignal(const bool drawOrNot){
  fDrawProcessingLevel[EECHistogramManager::kEnergyEnergyCorrelatorSignal] = drawOrNot;
}

// Setter for saving the figures to a file
void EECDrawer::SetSaveFigures(const bool saveOrNot, const char *format, const TString suffix){
  fSaveFigures = saveOrNot;
  fFigureFormat = format;
  fFigureSaveNameAppend = suffix;
}

// Setter for logarithmic pT axis
void EECDrawer::SetLogPt(const bool isLog){
  fLogPt = isLog;
}

// Setter for logarithmic deltaR axis in energy-energy correlators
void EECDrawer::SetLogDeltaR(const bool isLog){
  fLogDeltaR = isLog;
}

// Setter for logarithmic EEC axis in energy-energy correlators
void EECDrawer::SetLogEEC(const bool isLog){
  fLogEEC = isLog;
}

// Setter for legend comment
void EECDrawer::SetLegendComment(const TString newComment){
  fLegendComment = newComment;
}

// Setter for color palette
void EECDrawer::SetColorPalette(const int color){
  fColorPalette = color;
  gStyle->SetPalette(color);
}

// Setter for 2D drawing style
void EECDrawer::SetDrawingStyle2D(const char* style){
  fStyle2D = style;
}

// Setter for 3D drawing style
void EECDrawer::SetDrawingStyle3D(const char* style){
  fStyle3D = style;
}

// Setter for 2D drawing style
void EECDrawer::SetDrawingStyles(const int color, const char* style2D, const char* style3D){
  SetColorPalette(color);
  SetDrawingStyle2D(style2D);
  SetDrawingStyle3D(style3D);
}

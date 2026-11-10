#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "Run3ClosureTools.h"

/*
 * Macro for validating the background subtraction algorithm. In Pythia+Hydjet simulation, we can directly compare the extracted
 * signal to the true signal. It is done in this macro. Note that the input file must contain both regular and processed histograms.
 */
// Copy of validateBackgroundSubtraction.C at 792d76c0; original untouched.
// normalizeShapes=true retains the original normalized-shape comparison.
// For raw closure set normalizeShapes=false; only then may shared-signal
// approximation be requested. backgroundMethod: 0 mixed, 1 perpendicular,
// -1 read an already-processed signal (no shared-signal mode without provenance).
void validateBackgroundSubtractionRun3(
  TString inputFileName="merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root",
  int backgroundMethod=0, bool correlatedRatio=false, bool normalizeShapes=true,
  double minJetPt=80.,double maxJetPt=120.,double minCentrality=0.,double maxCentrality=30.,
  double trackThreshold=2.,TString outputDirectory="figures/run3LegacyValidation",bool averageSameCone=false){
  using Run3Closure::Require;
  Require(backgroundMethod>=-1 && backgroundMethod<=1,"Invalid background method");
  Require(!correlatedRatio || (!normalizeShapes && backgroundMethod>=0),"Shared-signal mode needs raw scale and explicitly reconstructed background");
  Require(gSystem->mkdir(outputDirectory,true)==0 || !gSystem->AccessPathName(outputDirectory),"Cannot create output directory");


  TFile* inputFile=Run3Closure::Open(inputFileName);

  EECCard* card = new EECCard(inputFile);
  Run3Closure::CheckCard(*card);
  const int weightExponent = card->GetWeightExponent();
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,1.5,2,2.5,3,3.5,4,300}
  
  int firstDrawnCentralityBin=card->FindBinIndexCentrality(minCentrality,maxCentrality);
  int lastDrawnCentralityBin=firstDrawnCentralityBin;
  int firstDrawnJetPtBinEEC=card->FindBinIndexJetPtEEC(minJetPt,maxJetPt);
  int lastDrawnJetPtBinEEC=firstDrawnJetPtBinEEC;
  int firstDrawnTrackPtBinEEC=Run3Closure::Track(*card,trackThreshold);
  int lastDrawnTrackPtBinEEC=firstDrawnTrackPtBinEEC;
  Require(firstDrawnCentralityBin>=0 && firstDrawnJetPtBinEEC>=0,"Requested OO bin is missing");

  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus] = false;
  
  // Find index of one energy-energy correlator that is drawn
  int studiedEnergyEnergyCorrelatorType = -1;
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    if(studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]){
      studiedEnergyEnergyCorrelatorType = iEnergyEnergyCorrelator;
      break;
    }
  }
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Axis zooming
  std::pair<double,double> ratioZoom = std::make_pair(0.9, 1.1);

  // Extra tag for the plots
  const bool simulationTag = false;
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "_energyWeightSquared_mixedEventBackgroundGenMatched";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager *histograms;
  histograms = new EECHistogramManager(inputFile,card);
  if(backgroundMethod>=0) Run3Closure::Preflight(*histograms,firstDrawnCentralityBin,firstDrawnJetPtBinEEC,firstDrawnTrackPtBinEEC,backgroundMethod);
    
  // Choose the energy-energy correlator types to load
  /*histograms->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus]);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus]);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus]);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus]);
    
  // Choose the bin ranges
  histograms->SetCentralityBinRange(0,nCentralityBins-1);
  histograms->SetJetPtBinRangeEEC(0,nJetPtBinsEEC);
  histograms->SetTrackPtBinRangeEEC(0,nTrackPtBinsEEC-1);
    
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();*/

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelatorSignal[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][2]; // Last bin, true signal/extracted signal
  TH1D* hEnergyEnergyCorrelatorSignalRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC]; // Ratio between true and extracted signal
  TH1D* hEnergyEnergyCorrelatorBackground[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][2]; // Last bin, true background/extracted background
  TH1D* hEnergyEnergyCorrelatorBackgroundRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC]; // Ratio between true and extracted background
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
          for(int iSignalType = 0; iSignalType < 2; iSignalType++){
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType] = NULL;
            hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType] = NULL;
          } // Signal type loop
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorBackgroundRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
  
  // Helper histograms
  TH1D* helperHistogram;
  double normalizationFactor;
  std::pair<double,double> drawingRange = std::make_pair(0.006, 0.39);
  double epsilon = 0.0000001;
  int lowSignalBin, highSignalBin;
  
  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only load the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
          
          // Read the true energy-energy correlator signal
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaPythia);
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]=Run3Closure::Clone(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0],"run3_validation_truth");

          // The extracted energy-energy correlator signal
          if(backgroundMethod<0) {
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1] = histograms->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorSignal);
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]=Run3Closure::Clone(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1],"run3_validation_extracted");
          } else {
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]=Run3Closure::Read(*histograms,iCentrality,iJetPt,iTrackPt,EECHistograms::kSameJetPair,EECHistograms::knSubeventCombinations,"run3_validation_extracted");
            std::unique_ptr<TH1D> background(Run3Closure::Background(*histograms,iCentrality,iJetPt,iTrackPt,backgroundMethod,averageSameCone));
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->Add(background.get(),-1.);
          }
          
          // Normalize the signal distributions to one in the region 0 to 0.4
          lowSignalBin = hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->FindBin(drawingRange.first+epsilon);
          highSignalBin = hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->FindBin(drawingRange.second-epsilon);
          if(normalizeShapes) for(int iSignalType = 0; iSignalType < 2; iSignalType++){
            Require(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType]->Integral(lowSignalBin,highSignalBin,"width")>0.,"Nonpositive shape normalization");
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType]->Scale(1/hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType]->Integral(lowSignalBin,highSignalBin,"width"));
          }
          
          // Calculate the extracted to true signal ratio
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->Clone(Form("signalRatio%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iTrackPt, iJetPt));
          
          Run3Closure::PositiveDenominator(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0],drawingRange.first,drawingRange.second);
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]);
          if(correlatedRatio) {
            std::unique_ptr<TH1D> total(Run3Closure::Read(*histograms,iCentrality,iJetPt,iTrackPt,EECHistograms::kSameJetPair,EECHistograms::knSubeventCombinations,"run3_validation_total"));
            std::unique_ptr<TH1D> truthBackground(Run3Closure::TrueBackground(*histograms,iCentrality,iJetPt,iTrackPt));
            Run3Closure::ApplySharedSignal(hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt],hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1],hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0],total.get(),truthBackground.get());
          }
          
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop

  // ==========================================================================
  //        Draw the selected energy-energy correlator signal ratios
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  if(logDeltaR) drawer->SetLogX(true);

  TLegend* legend;
  TLegend* simulationLegend;
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  double legendX1, legendX2;
  
  // Loop over all selected histograms
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only draw the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    // Loop over centrality bins
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      // Set the centrality information for legends and figure saving
      centralityString = Form("Cent: %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
      //centralityString = "Cent: 0-10%";
      compactCentralityString = Form("_C=%.0f-%.0f", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over jet pT bins
      for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
        
        // Set the jet pT information for legends and figure saving
        if(iJetPt == histograms->GetNJetPtBinsEEC()){
          jetPtString = Form("Jet p_{T} > %.0f", histograms->GetCard()->GetJetPtCut());
          compactJetPtString = "";
        } else {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
          compactJetPtString = Form("_J=%.0f-%.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
        }
        
        // Loop over track pT bins
        for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",histograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Create a legend for the figure
          legend = new TLegend(0.2,0.04,0.45,0.44);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, Form("Pythia+Hijing OO 5.36 TeV, %s",Run3Closure::Mode(*card).Data()), "");
          //legend->AddEntry((TObject*) 0, "Energy weight squared","");
          legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
          legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
          
          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();
          
          // Logarithmic EEC axis
          if(logEEC) drawer->SetLogY(true);
          
          // Set the x-axis drawing range
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
          // Draw the histograms to the upper canvas
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->SetLineColor(kBlack);
          drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0], "#Deltar", "EEC Signal", " ");
          legend->AddEntry(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0], "True signal", "l");
          
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->SetLineColor(kRed);
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->Draw("same");
          legend->AddEntry(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1], "Extracted signal", "l");
          
          // Draw the legends to the upper pad
          legend->Draw();

          // Add simulation tag to the figures
          if(simulationTag){
            legendX1 = (weightExponent == 1) ? 0.23 : 0.18;
            legendX2 = (weightExponent == 1) ? 0.48 : 0.43;
            simulationLegend = new TLegend(legendX1,0.5,legendX2,0.65);
            simulationLegend->SetFillStyle(0);simulationLegend->SetBorderSize(0);simulationLegend->SetTextSize(0.05);simulationLegend->SetTextFont(62);
            simulationLegend->AddEntry((TObject*) 0, "CMS simulation +","");
            simulationLegend->AddEntry((TObject*) 0, "   Private work","");
            simulationLegend->Draw();
          }
          
          // Linear scale for the ratio
          drawer->SetLogY(false);
          
          // Set the x-axis drawing range
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
          drawer->SetGridY(true);
          drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{Extracted}{True}", " ");
          drawer->SetGridY(false);
          
          // Save the figures to a file
          if(saveFigures){
            TCanvas* canvas=gPad->GetCanvas();
            Run3Closure::Header(canvas,correlatedRatio,averageSameCone,normalizeShapes?"Unit-area shape closure":"Raw-scale signal closure");
            TString stem=outputDirectory+"/"+Form("C%.0f-%.0f_J%.0f-%.0f_T%.1f_n%d_method%d",minCentrality,maxCentrality,minJetPt,maxJetPt,trackThreshold,weightExponent,backgroundMethod);
            stem+=(normalizeShapes?"_unitArea":"_raw");stem+=(averageSameCone?"_averagedCones":"_legacyC1");
            stem+=(correlatedRatio?"_sharedSignalApprox":"_uncorrelated");
            auto* truth=hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0];
            auto* extracted=hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1];
            auto* ratio=hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt];
            std::unique_ptr<TH1D> unity(Run3Closure::Clone(truth,"truth_reference_ratio")); unity->Divide(truth);
            Run3Closure::Save(canvas,stem,{truth,extracted},{unity.get(),ratio},{"True signal","Extracted signal"},inputFileName,correlatedRatio,averageSameCone,normalizeShapes);
          }
          
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
}

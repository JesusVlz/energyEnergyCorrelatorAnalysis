#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"
#include "AlgorithmLibrary.h"
#include "SystematicUncertaintyOrganizer.h"

/*
 * Macro for comparing final energy-energy correlators results
 */
void compareEECinDefinedBins(){
  
  // Files for comparison
  std::vector<TString> fileName;

  //fileName.push_back("data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_mixedConeSubtractedHFShift25_processed_2025-06-30.root");
  //fileName.push_back("data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_mixedConeSubtractedHFShift28_processed_2025-07-02.root");
  //fileName.push_back("data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_mixedConeSubtractedHFShift22_processed_2025-07-02.root");
  //fileName.push_back("data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_mixedConeNoHFShift_processed_2025-07-03.root");
  //fileName.push_back("data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_perpendicularConeSubtracted_processed_2025-06-30.root");
  //fileName.push_back("data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_noBackgroundSubtraction_processed_2025-06-30.root");

  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixedConeHFshift28_processed_2025-06-30.root");
  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixedConeHFshift31_smallStats_processed_2025-07-08.root");
  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixedConeHFshift25_smallStats_processed_2025-07-03.root");
  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixedEventSubtracted_processed_2025-06-05.root");
  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_perpendicularConeSubtracted_processed_2025-06-05.root");
  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfCsJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_perpendicularConeSubtracted_processed_2025-07-14.root");
  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_noSubtraction_processed_2025-06-05.root");

  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixHFPlusAndMinus_processed_2025-06-13.root");
  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixedEventSubtracted_processed_2025-06-05.root");
  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixHFPlus_processed_2025-06-13.root");
  //fileName.push_back("data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixHFMinus_processed_2025-06-13.root");

  //fileName.push_back("pythiaEposPilot_GenGen_truthLevelBackground_processed.root");
  //fileName.push_back("pythiaEposPilot_GenGen_mixingWithNominalHF_processed.root");
  //fileName.push_back("pythiaEposPilot_GenGen_mixingWithShiftedHF_processed.root");
  fileName.push_back("myRootFile_merged_processed.root");
    

  const int nComparisonFiles = fileName.size();

  TString uncertaintyFileName = "systematicUncertainties/systematicUncertainties_PbPb_nominalEnergyWeight_includeMCnonClosure_2024-03-08.root";
  // systematicUncertainties_PbPb_nominalEnergyWeight_includeMCnonClosure_2024-03-08.root
  // systematicUncertainties_PbPb_energyWeightSquared_includeMCnonClosure_2024-03-08.root


  std::vector<TString> fileDescription;
  //fileDescription.push_back("pPb 5.02 TeV, match HF p and pPb going");
  //fileDescription.push_back("pPb 5.02 TeV, match HF sum");
  //fileDescription.push_back("pPb 5.02 TeV, match HF Pb going");
  //fileDescription.push_back("pPb 5.02 TeV, match HF p going");
  //fileDescription.push_back("pp 5.02 TeV |#eta_{CM}| < 1.135");
  //fileDescription.push_back("pPb 8.16 TeV (p #rightarrow -#eta)");
  //fileDescription.push_back("pPb 8.16 TeV (p #rightarrow +#eta)");
  //fileDescription.push_back("pPb 5.02 TeV (p #rightarrow -#eta)");
  //fileDescription.push_back("pp 5 TeV, mixed cone");
  //fileDescription.push_back("pp 5 TeV, perpendicular cone");
  //fileDescription.push_back("pp 5 TeV, unsubtracted");
  //fileDescription.push_back("pPb 5 TeV, HF shift 28");
  //fileDescription.push_back("pPb 5 TeV, HF shift 31");
  //fileDescription.push_back("pPb 5 TeV, HF shift 25");
  //fileDescription.push_back("pPb 5 TeV, no HF shift");
  //fileDescription.push_back("pPb 5 TeV, unsubtracted");
  //fileDescription.push_back("pPb 5 TeV, perp cone sub");
  //fileDescription.push_back("pPb 5 TeV, perp cone and CS sub");
  //fileDescription.push_back("Pythia+EPOS, truth signal");
  //fileDescription.push_back("Mixed cone, no HF shift");
  //fileDescription.push_back("Mixed cone, with HF shift");
//  fileDescription.push_back("Perpendicular cone");
  fileDescription.push_back("OO 5.02 TeV");


  // Check that a description exists for each file
  if(fileDescription.size() < fileName.size()){
    cout << "ERROR! Not enough file descriptions given. Please give a description for all your files!" << endl;
    return;
  }
  
  // Open the files and check that they exist
  TFile* inputFile[nComparisonFiles];
  EECCard* card[nComparisonFiles];
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    
    inputFile[iFile] = TFile::Open(fileName.at(iFile));
    
    if(inputFile[iFile] == NULL){
      cout << "Error! The file " << fileName.at(iFile).Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    card[iFile] = new EECCard(inputFile[iFile]);
  }

  TFile* uncertaintyFile;
  EECCard* uncertaintyCard;
  const bool includeRelativeUncertainty = false; // Draw relative uncertainties in the ratio plot

  if(includeRelativeUncertainty){
    uncertaintyFile = TFile::Open(uncertaintyFileName);
  
    if(uncertaintyFile == NULL && includeRelativeUncertainty){
      cout << "Error! The file " << uncertaintyFileName.Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    uncertaintyCard = new EECCard(uncertaintyFile);
  }
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card[0]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[0]->GetNTrackPtBinsEEC();
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  //comparedCentralityBin.push_back(std::make_pair(0,10));
  comparedCentralityBin.push_back(std::make_pair(0,30));
  comparedCentralityBin.push_back(std::make_pair(30,50));
  comparedCentralityBin.push_back(std::make_pair(50,90));
  bool individualCentrality = true; // True = make different figure for each bin. False = plot all centrality bin to the same figure.

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(100,120));
  //comparedJetPtBin.push_back(std::make_pair(40,50));
  //comparedJetPtBin.push_back(std::make_pair(50,60));
  //comparedJetPtBin.push_back(std::make_pair(60,80));
  bool individualJetPt = true; // True = make different figure for each bin. False = plot all jet pT bin to the same figure.

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(0.7);  // T0 = lowest track pT bin
  //comparedTrackPtBin.push_back(1.5);
  //comparedTrackPtBin.push_back(2.0);
  //comparedTrackPtBin.push_back(2.5);
  //comparedTrackPtBin.push_back(3.0);
  bool individualTrackPt = true; // True = make different figure for each bin. False = plot all track pT bin to the same figure.

  // Different normalization options
  enum enumNormalization{kNoNormalization, kNormalizeToPairs, kNormalizeToJets};
  const int normalizeDistributions = kNormalizeToPairs;

  // Choose the type of draw energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorNormalized = Normalized energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorBackground = Estimated background
  // EECHistogramManager::kEnergyEnergyCorrelatorSignal = Background subtracted, but not unfolded energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorUnfolded = Unfolded energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorBackgroundAfterUnfolding = Estimated background after unfolding
  // EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal = Unfolded energy-energy correlator signal
  // EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels = Raw energy-energy correlator
  int drawnEnergyEnergyCorrelator = EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels;

  // Choose the pairing type if raw energy-energy correlator is drawn
  // EECHistograms::kSameJetPair;
  // EECHistograms::kSignalReflectedConePair;
  // EECHistograms::kReflectedConePair;
  // EECHistograms::kSignalMixedConePair;
  // EECHistograms::kReflectedMixedConePair;
  // EECHistograms::kMixedConePair;
  // EECHistograms::kSignalSecondMixedConePair; 
  // EECHistograms::kReflectedSecondMixedConePair; 
  // EECHistograms::kMixedMixedConePair; 
  // EECHistograms::kSecondMixedConePair;
  int iPairingType = EECHistograms::kSameJetPair;

  // If we are dealing with pp data, reset the centrality vector
  if(card[0]->GetDataType().Contains("pp") || card[0]->GetDataType().Contains("pPb")){
    comparedCentralityBin.clear();
    comparedCentralityBin.push_back(std::make_pair(-1,100));
  }

  // If we are dealing with MC, shift the centrality by 4% as is done in order to match background energy density
  //if(card[0]->GetDataType().Contains("MC")){
  //  for(auto& centralityBin : comparedCentralityBin){
  //    centralityBin.first += 4;
  //    centralityBin.second += 4;
  //  }
  //}

  // ====================================================
  //                Drawing configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = true;
  const char* saveComment = "_OO_raw";
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.81, 1.19);
  std::pair<double, double> eecZoom = std::make_pair(0.2, 30);
  const bool automaticZoom = true;

  // Sanity checks for input. Ensure that all the selected bins actually exist in the input files.
  // This check is only needed for unfolded bins, so skip it if only raw distribution is drawn.
  if(drawnEnergyEnergyCorrelator > EECHistogramManager::kEnergyEnergyCorrelatorSignal && drawnEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){

    for(int iFile = 0; iFile < nComparisonFiles; iFile++){

      // Sanity check for centrality bins
      for(auto centralityBin : comparedCentralityBin){

        // For MC, shift the centrality bin borders
        if(card[iFile]->GetDataType().Contains("MC")){
          centralityBin.first += 4;
          centralityBin.second += 4;
        }

        if(card[iFile]->FindBinIndexCentrality(centralityBin) < card[iFile]->GetFirstUnfoldedCentralityBin() || card[iFile]->FindBinIndexCentrality(centralityBin) > card[iFile]->GetLastUnfoldedCentralityBin()){
          cout << "ERROR! Centrality bin " << centralityBin.first << "-" << centralityBin.second << " does not exist in file " << fileName[iFile].Data() << endl;
          cout << "Please only choose centrality bins that are included in the input files." << endl;
          return;
        } 
      }

      // Sanity check for jet pT bins
      for(auto jetPtBin : comparedJetPtBin){
        if(card[iFile]->FindBinIndexJetPtEEC(jetPtBin) < card[iFile]->GetFirstUnfoldedJetPtBin() || card[iFile]->FindBinIndexJetPtEEC(jetPtBin) > card[iFile]->GetLastUnfoldedJetPtBin()){
          cout << "ERROR! Jet pT bin " << jetPtBin.first << "-" << jetPtBin.second << " does not exist in file " << fileName[iFile].Data() << endl;
          cout << "Please only choose jet pT bins that are included in the input files." << endl;
          return;
        }
      }

      // Sanity check for track pT bins
      for(auto trackPtBin : comparedTrackPtBin){
        if(card[iFile]->GetBinIndexTrackPtEEC(trackPtBin) < card[iFile]->GetFirstUnfoldedTrackPtBin() || card[iFile]->GetBinIndexTrackPtEEC(trackPtBin) > card[iFile]->GetLastUnfoldedTrackPtBin()){
          cout << "ERROR! Track pT cut > " << trackPtBin << " GeV does not exist in file " << fileName[iFile].Data() << endl;
          cout << "Please only choose track pT bins that are included in the input files." << endl;
          return;
        }
      } 
    } // File loop for input sanity check
  } // Unfolded distributions

  // Only allow one variable for which all bins are plotted to the same figure
  if(individualCentrality + individualJetPt + individualTrackPt < 2){
    cout << "You are tring to plot too many bins to the same figure!" << endl;
    cout << "This macro can only plot all the bins from one variable." << endl;
    cout << "Please check your configuration!" << endl; 
    return;
  }

  // Do not allow several bins plotted to same figure if several files are compared
  if(nComparisonFiles > 1 && (!individualTrackPt || !individualJetPt || !individualCentrality)){
    cout << "If you compare files, you cannot draw several bins from one variable to single figure." << endl;
    cout << "Please check your configuration!" << endl; 
    return;
  }

  // Based on the input from user, make a construct that can be handled later in the code to plot the selected bins
  enum enumTupleDecoder{kCentrality, kJetPt, kTrackPt}; // Components of the n-tuple in binning vector
  std::vector<std::tuple<std::vector<std::pair<double,double>>, std::vector<std::pair<double,double>>, std::vector<double>>> binningInformation;

  // Separate plots for each centrality bin
  std::vector<std::pair<double,double>> temporaryCentralityVector;
  std::vector<std::pair<double,double>> temporaryJetPtVector;
  std::vector<double> temporaryTrackPtVector;

  if(individualCentrality){
    for(auto centralityBin : comparedCentralityBin){
      temporaryCentralityVector.clear();
      temporaryCentralityVector.push_back(centralityBin);

      // Separate plots for each jet pT bin
      if(individualJetPt){

        for(auto jetPtBin : comparedJetPtBin){
          temporaryJetPtVector.clear();
          temporaryJetPtVector.push_back(jetPtBin);

          // Separate track pT bins for each bin
          if(individualTrackPt){

            for(auto trackPtBin : comparedTrackPtBin){
              temporaryTrackPtVector.clear();
              temporaryTrackPtVector.push_back(trackPtBin);

              binningInformation.push_back(std::make_tuple(temporaryCentralityVector, temporaryJetPtVector, temporaryTrackPtVector));

            }

          // All defined track pT bins in a single plot
          } else {

            binningInformation.push_back(std::make_tuple(temporaryCentralityVector, temporaryJetPtVector, comparedTrackPtBin));

          } // Track pT binning if-else

        } // Loop over ket pT bins

      // All defined jet pT bins in a single plot
      } else {

        // Separate track pT bins for each bin
        if(individualTrackPt){

          for(auto trackPtBin : comparedTrackPtBin){
            temporaryTrackPtVector.clear();
            temporaryTrackPtVector.push_back(trackPtBin);

            binningInformation.push_back(std::make_tuple(temporaryCentralityVector, comparedJetPtBin, temporaryTrackPtVector));

          }

          // All defined track pT bins in a single plot
          } else {

            binningInformation.push_back(std::make_tuple(temporaryCentralityVector, comparedJetPtBin, comparedTrackPtBin));

          } // Track pT binning if-else

      } // Jet pT binning if-else

    } // Loop over centrality bins

  // All defined centrality bins in a single plot 
  } else {

    // Separate plots for each jet pT bin
    if(individualJetPt){

      for(auto jetPtBin : comparedJetPtBin){
        temporaryJetPtVector.clear();
        temporaryJetPtVector.push_back(jetPtBin);

        // Separate track pT bins for each bin
        if(individualTrackPt){

          for(auto trackPtBin : comparedTrackPtBin){
            temporaryTrackPtVector.clear();
            temporaryTrackPtVector.push_back(trackPtBin);

            binningInformation.push_back(std::make_tuple(comparedCentralityBin, temporaryJetPtVector, temporaryTrackPtVector));

          }

        // All defined track pT bins in a single plot
        } else {

          binningInformation.push_back(std::make_tuple(comparedCentralityBin, temporaryJetPtVector, comparedTrackPtBin));

        } // Track pT binning if-else

      } // Loop over ket pT bins

    // All defined jet pT bins in a single plot
    } else {

      // Separate track pT bins for each bin
      if(individualTrackPt){

        for(auto trackPtBin : comparedTrackPtBin){
          temporaryTrackPtVector.clear();
          temporaryTrackPtVector.push_back(trackPtBin);

          binningInformation.push_back(std::make_tuple(comparedCentralityBin, comparedJetPtBin, temporaryTrackPtVector));

        }

        // All defined track pT bins in a single plot
        } else {

          binningInformation.push_back(std::make_tuple(comparedCentralityBin, comparedJetPtBin, comparedTrackPtBin));

        } // Track pT binning if-else

    } // Jet pT binning if-else
    
  } // Centrality binning if-else
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[nComparisonFiles];
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    histograms[iFile] = new EECHistogramManager(inputFile[iFile], card[iFile]);
  }

  // Create systematic uncertainty organizer to illustrate to draw the pp systematic uncertainties to plots
  SystematicUncertaintyOrganizer* uncertaintyOrganizer;
  if(includeRelativeUncertainty) uncertaintyOrganizer = new SystematicUncertaintyOrganizer(uncertaintyFile);

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorRatio[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hRelativeUncertainty[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        hRelativeUncertainty[iCentrality][iJetPt][iTrackPt] = NULL;
        for(int iFile = 0; iFile < nComparisonFiles; iFile++){
          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
        } // File loop
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop
  
  
  // Helper histograms
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality, iCentralityReference, iCentralityUncertainty;
  int iTrackPt, iTrackPtReference, iTrackPtUncertainty;
  int iJetPt, iJetPtReference, iJetPtUncertainty;
  std::pair<double,double> referenceCentralityBin;


  // Transformer to transform absolute uncertainties to relative ones
  AlgorithmLibrary* optimusPrimeTheTransformer = new AlgorithmLibrary();

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(auto jetPtBin : comparedJetPtBin){
      for(auto trackPtBin : comparedTrackPtBin){
        for(auto centralityBin : comparedCentralityBin){

          // For MC, shift the centrality bin borders
          referenceCentralityBin = centralityBin;
          if(card[0]->GetDataType().Contains("MC")){
            referenceCentralityBin.first += 4;
            referenceCentralityBin.second += 4;
          }
          if(card[iFile]->GetDataType().Contains("MC")){
            centralityBin.first += 4;
            centralityBin.second += 4;
          }

          // Find the proper binning and express it in term of the bins in the first file
          iJetPt = card[iFile]->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPt = card[iFile]->GetBinIndexTrackPtEEC(trackPtBin);
          iCentrality = card[iFile]->FindBinIndexCentrality(centralityBin);
          iJetPtReference = card[0]->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPtReference = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
          iCentralityReference = card[0]->FindBinIndexCentrality(referenceCentralityBin);

          if(card[0]->GetDataType().Contains("pp") || card[0]->GetDataType().Contains("pPb")){
            iCentrality = 0;
            iCentralityReference = 0;
          }

          // Load the selected energy-energy correlator histogram
          if(drawnEnergyEnergyCorrelator == EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
            hEnergyEnergyCorrelator[iFile][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms[iFile]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType);
          } else {
            hEnergyEnergyCorrelator[iFile][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, drawnEnergyEnergyCorrelator);
          }

          // Normalize the distributions to one in the drawingRange
          if(normalizeDistributions == kNormalizeToPairs){
            lowNormalizationBin = hEnergyEnergyCorrelator[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highNormalizationBin = hEnergyEnergyCorrelator[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.second - epsilon);

            hEnergyEnergyCorrelator[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1.0 / hEnergyEnergyCorrelator[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Integral(lowNormalizationBin, highNormalizationBin, "width"));
          } else if(normalizeDistributions == kNormalizeToJets){
            hEnergyEnergyCorrelator[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1.0 / histograms[iFile]->GetJetPtIntegral(iCentrality, jetPtBin.first, jetPtBin.second));
          }

          // Uncertainty histograms
          if(iFile == 0 && includeRelativeUncertainty){
            iCentralityUncertainty = uncertaintyCard->FindBinIndexCentrality(centralityBin);
            iJetPtUncertainty = uncertaintyCard->FindBinIndexJetPtEEC(jetPtBin);
            iTrackPtUncertainty = uncertaintyCard->GetBinIndexTrackPtEEC(trackPtBin);

            hRelativeUncertainty[iCentrality][iJetPt][iTrackPt] = (TH1D*) uncertaintyOrganizer->GetSystematicUncertainty(iCentralityUncertainty, iJetPtUncertainty, iTrackPtUncertainty)->Clone(Form("relativeUncertainty%d%d%d", iCentrality, iJetPt, iTrackPt));
            optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertainty[iCentrality][iJetPt][iTrackPt], true);
          }

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // File loop

  // After all the histograms have been read, calculate the ratios
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);
      iJetPtReference = individualJetPt ? iJetPt : card[0]->FindBinIndexJetPtEEC(comparedJetPtBin.at(0));
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
        iTrackPtReference = individualTrackPt ? iTrackPt : card[0]->GetBinIndexTrackPtEEC(comparedTrackPtBin.at(0));
        for(auto centralityBin : comparedCentralityBin){

          if(card[0]->GetDataType().Contains("MC")){
            centralityBin.first += 4;
            centralityBin.second += 4;
          }

          referenceCentralityBin = comparedCentralityBin.at(0);
          if(card[0]->GetDataType().Contains("MC")){
            referenceCentralityBin.first += 4;
            referenceCentralityBin.second += 4;
          }

          if(card[0]->GetDataType().Contains("pp") || card[0]->GetDataType().Contains("pPb")){
            iCentrality = 0;
            iCentralityReference = 0;
          } else {
            iCentrality = card[0]->FindBinIndexCentrality(centralityBin);
            iCentralityReference = individualCentrality ? iCentrality : card[0]->FindBinIndexCentrality(referenceCentralityBin);
          }

          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecRatio%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[0][iCentralityReference][iJetPtReference][iTrackPtReference]);

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // File loop
  
  // ==========================================================================
  //                Draw all the distribution in the same figure
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  drawer->SetLogX(true);

  TString compactCentralityString = "";
  TString compactJetPtString = "";
  TString compactTrackPtString = "";
  TString comparedVariableString = "";
  TString ratioName = "";
  TString energyWeightString = (card[0]->GetWeightExponent() == 1) ? "Nominal energy weight" : "Energy weight squared";
  TString legendString;
  int markerStyle[5] = {kFullCircle, kOpenSquare, kOpenCross, kFullStar, kFullCross};
  int color[] = {kBlack,kRed,kBlue,kGreen+3,kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};

  // Binning vectors
  std::vector<int> currentCentralityIndices;
  std::vector<int> currentJetPtIndices;
  std::vector<int> currentTrackPtIndices;
  bool colorWithCentrality = false;
  bool colorWithJetPt = false;
  bool colorWithTrackPt = false;
  int colorFinder = 0;
  int firstCentralityBin = 0;
  int firstTrackPtBin = 0;
  int firstJetPtBin = 0;
  int legendCentralityIndex = 0;
  int legendJetPtIndex = 0;
  int legendTrackPtIndex = 0;
  double minimumCandidate, maximumCandidate;
  TString individualLegend;
  TString yAxisName = "EEC (raw)";

  for(auto plottedBin : binningInformation){

    // Find the indices for each drawn bin
    currentCentralityIndices.clear();
    for(auto centralityBin : std::get<kCentrality>(plottedBin)){

      if(card[0]->GetDataType().Contains("MC")){
        centralityBin.first += 4;
        centralityBin.second += 4;
      }
      if(card[0]->GetDataType().Contains("pp") || card[0]->GetDataType().Contains("pPb")){
        currentCentralityIndices.push_back(0);
      } else {
        currentCentralityIndices.push_back(card[0]->FindBinIndexCentrality(centralityBin));
      }
    }

    currentJetPtIndices.clear();
    for(auto jetPtBin : std::get<kJetPt>(plottedBin)){
      currentJetPtIndices.push_back(card[0]->FindBinIndexJetPtEEC(jetPtBin));
    }

    currentTrackPtIndices.clear();
    for(auto trackPtBin : std::get<kTrackPt>(plottedBin)){
      currentTrackPtIndices.push_back(card[0]->GetBinIndexTrackPtEEC(trackPtBin));
    }

    // If the lenght of the indices array is larger than 1, we are plotting all these bins to the same plot
    colorWithCentrality = (currentCentralityIndices.size() > 1);
    colorWithTrackPt = (currentTrackPtIndices.size() > 1);
    colorWithJetPt = (currentJetPtIndices.size() > 1);

    // Remember the first drawn bin index
    firstCentralityBin = currentCentralityIndices.at(0);
    firstTrackPtBin = currentTrackPtIndices.at(0);
    firstJetPtBin = currentJetPtIndices.at(0);
          
    // Create a new canvas for the plot
    drawer->CreateSplitCanvas();

    // Logarithmic EEC axis
    drawer->SetLogY(true);

    TLegend* legend = new TLegend(0.18,0.04,0.45,0.58);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

    // Add the file description to the legend if coloring with any variable
    if(colorWithCentrality || colorWithJetPt || colorWithTrackPt){
      legend->AddEntry((TObject*) 0, fileDescription.at(0), "");
    }

    // Add common legend variables and define figure naming in case figures are saved
    if(!colorWithCentrality){ 
      if(card[0]->GetDataType().Contains("pp") || card[0]->GetDataType().Contains("pPb")){
        compactCentralityString = "";
      } else {
        legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second), "");
        compactCentralityString = Form("_C=%.0f-%.0f", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second);
      }
    } else {
      compactCentralityString = "";
      comparedVariableString = "_centralityComparison";
    }

    if(!colorWithJetPt) {
      legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", std::get<kJetPt>(plottedBin).at(0).first, std::get<kJetPt>(plottedBin).at(0).second), "");
      compactJetPtString = Form("_J=%.0f-%.0f", std::get<kJetPt>(plottedBin).at(0).first, std::get<kJetPt>(plottedBin).at(0).second);
    } else {
      compactJetPtString = "";
      comparedVariableString = "_jetPtComparison";
    }

    if(!colorWithTrackPt){ 
      legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", std::get<kTrackPt>(plottedBin).at(0)), "");
      compactTrackPtString = Form("_T>%.1f",std::get<kTrackPt>(plottedBin).at(0));
      compactTrackPtString.ReplaceAll(".","v");
    } else {
      compactTrackPtString = "";
      comparedVariableString = "_trackPtComparison";
    }

    legend->AddEntry((TObject*) 0, energyWeightString, "");

    // Set drawing style for all histograms
    colorFinder = 0;
    for(int iFile = 0; iFile < nComparisonFiles; iFile++){
      for(int iJetPt : currentJetPtIndices){
        for(int iTrackPt : currentTrackPtIndices){
          for(int iCentrality : currentCentralityIndices){
            if(colorWithCentrality) colorFinder++;
            hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFile+colorFinder]);
            hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile+colorFinder]);
            hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile+colorFinder]); 
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFile+colorFinder]);
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile+colorFinder]);
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile+colorFinder]);
          } // Centrality binning
          if(colorWithTrackPt) colorFinder++;
        } // Track pT binning 
        if(colorWithJetPt) colorFinder++;
      } // Jet pT binning
    } // File loop

    // Automatic zooming for the drawn histograms
    if(automaticZoom){
      eecZoom.first = 10000;
      eecZoom.second = 0;
      for(int iCentrality : currentCentralityIndices){
        for(int iJetPt : currentJetPtIndices){
          for(int iTrackPt : currentTrackPtIndices){
            for(int iFile = 0; iFile < nComparisonFiles; iFile++){
              hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              minimumCandidate = hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->GetMinimum();
              maximumCandidate = hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->GetMaximum();
              if(minimumCandidate < eecZoom.first) eecZoom.first = minimumCandidate;
              if(maximumCandidate > eecZoom.second) eecZoom.second = maximumCandidate;
            } // File loop
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
      eecZoom.first = eecZoom.first / 2.0;
      eecZoom.second = eecZoom.second * 2.0;
    }

    // Set the x- and y-axis drawing ranges
    hEnergyEnergyCorrelator[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
    hEnergyEnergyCorrelator[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);
          
    // Draw the histograms to the upper canvas
    if(drawnEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
      yAxisName = Form("EEC %s", histograms[0]->GetEnergyEnergyCorrelatorProcessSaveName(drawnEnergyEnergyCorrelator));
    }
    drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin], "#Deltar", yAxisName, " ");

    for(int iCentrality : currentCentralityIndices){
      for(int iJetPt : currentJetPtIndices){
        for(int iTrackPt : currentTrackPtIndices){
          for(int iFile = 0; iFile < nComparisonFiles; iFile++){
            if(iFile == 0 && iCentrality == firstCentralityBin && iJetPt == firstJetPtBin && iTrackPt == firstTrackPtBin) continue;
            hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
          } // File loop
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

    // Add legends for drawn histograms
    legendCentralityIndex = 0;
    legendJetPtIndex = 0;
    legendTrackPtIndex = 0;
    for(int iFile = 0; iFile < nComparisonFiles; iFile++){
      if(nComparisonFiles > 1){
        individualLegend = fileDescription.at(iFile);
      } else {
        individualLegend = "";
      }
      for(int iJetPt : currentJetPtIndices){
        if(colorWithJetPt) individualLegend = Form(" %.0f < jet p_{T} < %.0f GeV", std::get<kJetPt>(plottedBin).at(legendJetPtIndex).first, std::get<kJetPt>(plottedBin).at(legendJetPtIndex).second);
        legendJetPtIndex++;
        for(int iTrackPt : currentTrackPtIndices){
          if(colorWithTrackPt) individualLegend = Form(" track p_{T} > %.1f GeV", std::get<kTrackPt>(plottedBin).at(legendTrackPtIndex++));
          for(int iCentrality : currentCentralityIndices){
            if(colorWithCentrality) individualLegend = Form(" Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(legendCentralityIndex).first, std::get<kCentrality>(plottedBin).at(legendCentralityIndex).second);
            legendCentralityIndex++;
            legend->AddEntry(hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt], individualLegend.Data(), "p");
          } // Centrality loop 
        } // Track pT loop
      } // Jet pT loop
    } // File loop
  
    // Draw the legends to the upper pad
    legend->Draw();
          
    // Linear scale for the ratio
    drawer->SetLogY(false);
          
    // Set the drawing ranges
    if(includeRelativeUncertainty){
      hRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
      hRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
      hRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin]->SetFillColorAlpha(kBlack, 0.2);
      hRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin]->SetMarkerStyle(9);
      hRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin]->SetMarkerSize(0);
    }
    hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);          
    hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

    // Draw the histograms
    drawer->SetGridY(true);
    if(nComparisonFiles > 1){
      ratioName = fileDescription.at(0);
    } else if (colorWithCentrality){
      ratioName = Form("Cent: %.0f-%.0f%%", comparedCentralityBin.at(0).first, comparedCentralityBin.at(0).second);
    } else if (colorWithTrackPt){
      ratioName = Form("track p_{T} > %.1f GeV", comparedTrackPtBin.at(0));
    } else if (colorWithJetPt){
      ratioName = Form("%.0f < jet p_{T} < %.0f GeV", comparedJetPtBin.at(0).first, comparedJetPtBin.at(0).second);
    }

    if(includeRelativeUncertainty){
      drawer->DrawHistogramToLowerPad(hRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin], "#Deltar", Form("#frac{Color}{%s}", ratioName.Data()), " ", "e2");
    } else {
      //drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin], "#Deltar", Form("#frac{Color}{%s}", ratioName.Data()), " ");
      drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin], "#Deltar", "#frac{Color}{Truth signal}", " ");
    }
    for(int iCentrality : currentCentralityIndices){
      for(int iJetPt : currentJetPtIndices){
        for(int iTrackPt : currentTrackPtIndices){
          for(int iFile = 0; iFile < nComparisonFiles; iFile++){
            if(iFile == 0 && iCentrality == firstCentralityBin && iJetPt == firstJetPtBin && iTrackPt == firstTrackPtBin) continue;
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
          } // File loop
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
    drawer->SetGridY(false);
          
    // Save the figures to a file
    if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/eecComparison%s%s%s%s%s.%s", saveComment, comparedVariableString.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
    }
  }

}

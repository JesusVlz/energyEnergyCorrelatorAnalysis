#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "EECUnfoldConfiguration.h"
#include "JDrawer.h"

/*
 * Reset underflow and overflow bins for a one-dimensional histogram
 */
void removeOutOfRange(TH1D* histogramInNeedOfTrimming){
   histogramInNeedOfTrimming->SetBinContent(0, 0);
   histogramInNeedOfTrimming->SetBinContent(histogramInNeedOfTrimming->GetNbinsX()+1, 0);
}

/*
 * Reset underflow and overflow bins for a two-dimensional histogram
 */
void removeOutOfRange(TH2D* histogramInNeedOfTrimming)
{
   int NX = histogramInNeedOfTrimming->GetNbinsX();
   int NY = histogramInNeedOfTrimming->GetNbinsY();

   for(int iX = 0; iX <= NX + 1; iX++)
   {
      histogramInNeedOfTrimming->SetBinContent(iX, 0, 0);
      histogramInNeedOfTrimming->SetBinContent(iX, NY+1, 0);
   }
   for(int iY = 0; iY <= NY + 1; iY++)
   {
      histogramInNeedOfTrimming->SetBinContent(0, iY, 0);
      histogramInNeedOfTrimming->SetBinContent(NX+1, iY, 0);
   }
}

/*
 * Randomize a one-dimensional histogram within the statistical uncertainties of the histogram
 */
void randomizeWithinUncertainties(TH1D* randomizedHistogram){
  TRandom3* randomizer = new TRandom3();
  randomizer->SetSeed(0);
  double currentBinContent, currentBinError;
  for(int iBin = 1; iBin <= randomizedHistogram->GetNbinsX(); iBin++){
    currentBinContent = randomizedHistogram->GetBinContent(iBin);
    if(currentBinContent == 0) continue;
    currentBinError = randomizedHistogram->GetBinError(iBin);
    randomizedHistogram->SetBinContent(iBin, randomizer->Gaus(currentBinContent, currentBinError));
  }
  delete randomizer;
}

/*
 * Randomize a two-dimensional histogram within the statistical uncertainties of the histogram
 */
void randomizeWithinUncertainties(TH2D* randomizedHistogram){
  TRandom3* randomizer = new TRandom3();
  randomizer->SetSeed(0);
  double currentBinContent, currentBinError;
  for(int xBin = 1; xBin <= randomizedHistogram->GetNbinsX(); xBin++){
    for(int yBin = 1; yBin <= randomizedHistogram->GetNbinsY(); yBin++){
      currentBinContent = randomizedHistogram->GetBinContent(xBin, yBin);
      if(currentBinContent == 0) continue;
      currentBinError = randomizedHistogram->GetBinError(xBin, yBin);
      randomizedHistogram->SetBinContent(xBin, yBin, randomizer->Gaus(currentBinContent, currentBinError));
    }
  }
  delete randomizer;
}

/*
 * Macro for doing the jet pT unfolding for the 
 *
 *  Arguments:
 *   TString fileName = File from which the histograms are read and to which the processed histograms are written
 *   TString outputFileName = If given, the processed histograms are written to this file. Otherwise the fileName file is updated.
 *   const int iSplit = The Monte Carlo split used to derive the response matrix
 *                      0: Use the full Monte Carlo statistics to construct the response matrix
 *                      1: Use the first half of the Monte Carlo statistics to construct the response matrix
 *                      2: Use the second half of the Monte Carlo statistics to construct the response matrix
 *   const int iSystematic = Index for unfolding systematic uncertainty study
 *                           0: Default results, no systematic study
 *                           1: Evaluate down systematic uncertainties for jet pT resolution
 *                           2: Evaluate up systematic uncertainties for jet pT resolution
 *                           3: Evaluate down systematic uncertainties for jet energy scale
 *                           4: Evaluate up systematic uncertainties for jet energy scale
 *                           5: Evaluate systematic uncertainties for jet pT prior shape
 *                           6: Evaluate systematic uncertainty from 2% centrality shift
 *                           7: Evaluate systematic uncertainty from 6% centrality shift
 *                           8: Evaluate systematic uncertainty from fewer number of iterations in unfolding
 *                           9: Evaluate systematic uncertainty from larger number of iterations in unfolding
 *                          10: Evaluate systematic uncertainty from finite statistics of the MC sample
 *   const int iEnergyEnergyCorrelator = Energy-energy correlator index for the unfolded correlator. Indices are explained in EECHistogramManager.h
 *.  int iWeightExponent = Exponent given for the energy weights for energy-energy correaltors. 
 *                               0: Read the used weight exponent from the input file
 *                               >0: Use manually defined exponent. Currently 1 and 2 are implemented.
 *   TString responseFileName = Explicit processed response file; required for OO/ppRef.
 *   const int numberOfIterations = Explicit positive iteration count for that file, chosen with closure tests.
 *                                 Explicit file/count overrides automatic legacy split/systematic selection.
 */
void unfoldEEChistograms(TString dataFileName, TString outputFileName, const int iSplit = 0, const int iSystematic = 0, const int iEnergyEnergyCorrelator = EECHistogramManager::kEnergyEnergyCorrelator, int iWeightExponent = 0, TString responseFileName = "", const int numberOfIterations = 0){

  // **********************************
  //       Open the input files
  // **********************************

  // First, open the data iput file and obtain the EECCard from the file
  TFile* dataInputFile = TFile::Open(dataFileName);

  if(dataInputFile == NULL || dataInputFile->IsZombie()) {
    cout << "Error! The file " << dataFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* dataCard = new EECCard(dataInputFile);

  // OO, like PbPb, has distinct centrality bins.
  TString collisionSystem = dataCard->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  bool isOOData = collisionSystem.Contains("OO");
  bool hasCentrality = isPbPbData || isOOData;

  // Flag for generator level jets. If data file has generator level jets, do not unfold.
  // Running it still through this macro gives consistent histogram naming, which makes some things easier later.
  bool genJetsInData = (collisionSystem.Contains("GenGen") || collisionSystem.Contains("GenReco"));

  // If the weight exponent is 0, read the value that has been used in the input data file
  if(iWeightExponent == 0) iWeightExponent = dataCard->GetWeightExponent();

  // Use the information from the dataCard to create an unfolding configuration
  EECUnfoldConfiguration* unfoldConfigurationProvider = new EECUnfoldConfiguration(dataCard, iSplit, iSystematic, iWeightExponent, responseFileName, numberOfIterations);

  // Open the response matric file based on the information on the unfolding configuration
  TFile* responseInputFile = TFile::Open(unfoldConfigurationProvider->GetResponseFileName());
  
  if(responseInputFile == NULL || responseInputFile->IsZombie()) {
    cout << "Error! The file " << unfoldConfigurationProvider->GetResponseFileName().Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // Get the EECCard also from the response matrix file
  EECCard* responseCard = new EECCard(responseInputFile);
  
  
  // **********************************************************************************
  //       Check that the centrality, jet pT and track pT bins match with the files
  // **********************************************************************************

  if(isPbPbData != responseCard->GetDataType().Contains("PbPb") || isOOData != responseCard->GetDataType().Contains("OO") || collisionSystem.Contains("ppRef") != responseCard->GetDataType().Contains("ppRef")){
    cout << "Collision systems differ between the data and unfolding histograms!" << endl;
    cout << "Please ensure that the data and unfolding histograms are comparible!" << endl;
    return;
  }

  // PbPb legacy matrices can use shifted bins; explicit OO responses require matching bins.
  if(hasCentrality){
    const double centralityTolerance = isOOData ? 0.0001 : 6;
    if(responseCard->GetNCentralityBins() != dataCard->GetNCentralityBins()) {
      cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    for(int iCentrality = 0; iCentrality < responseCard->GetNCentralityBins(); iCentrality++){
      if(TMath::Abs(responseCard->GetLowBinBorderCentrality(iCentrality) - dataCard->GetLowBinBorderCentrality(iCentrality)) > centralityTolerance){
        cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      if(TMath::Abs(responseCard->GetHighBinBorderCentrality(iCentrality) - dataCard->GetHighBinBorderCentrality(iCentrality)) > centralityTolerance){
        cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
    }
  }  // Checking compatible centrality bins

  if(responseCard->GetNTrackPtBinsEEC() != dataCard->GetNTrackPtBinsEEC()){
    cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iTrackPt = 0; iTrackPt < responseCard->GetNTrackPtBinsEEC(); iTrackPt++){
    if(TMath::Abs(responseCard->GetLowBinBorderTrackPtEEC(iTrackPt) - dataCard->GetLowBinBorderTrackPtEEC(iTrackPt)) > 0.01){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(responseCard->GetHighBinBorderTrackPtEEC(iTrackPt) - dataCard->GetHighBinBorderTrackPtEEC(iTrackPt)) > 0.01){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }  // Checking the compatibility of track pT bins

  if(responseCard->GetNJetPtBinsEEC() != dataCard->GetNJetPtBinsEEC()){
    cout << "Error! Measured jet pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iJetPt = 0; iJetPt < responseCard->GetNJetPtBinsEEC(); iJetPt++){
    if(TMath::Abs(responseCard->GetLowBinBorderJetPtEEC(iJetPt) - dataCard->GetLowBinBorderJetPtEEC(iJetPt)) > 0.01){
      cout << "Error! Measured jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(responseCard->GetHighBinBorderJetPtEEC(iJetPt) - dataCard->GetHighBinBorderJetPtEEC(iJetPt)) > 0.01){
      cout << "Error! Measured jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }  // Checking the compatibility of reconstructed jet pT bins

  // ********************************************************
  //       Binning configuration for the unfolding study
  // ********************************************************
  
  // Find the number of bins from the card
  const int nCentralityBins = hasCentrality ? responseCard->GetNCentralityBins() : 1;
  const int nTrackPtBins = responseCard->GetNTrackPtBinsEEC();
  const int nJetPtBins = dataCard->GetNJetPtBinsEEC();

  // Bin range to be studied
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = nCentralityBins-1;
  
  int firstStudiedTrackPtBinEEC = 1;
  int lastStudiedTrackPtBinEEC = 5;
  if(firstStudiedTrackPtBinEEC < 0 || lastStudiedTrackPtBinEEC >= nTrackPtBins){
    cout << "Error! Configure the studied track pT bin range for the current input cards before unfolding." << endl;
    return;
  }

  // Run-3 cards define the jet pT binning; select complete bins inside the
  // existing 120-200 GeV study window. Preserve the legacy four-bin selection.
  std::vector<std::pair<double,double>> unfoldedJetPtBins = EECUnfoldConfiguration::SelectUnfoldedJetPtBins(*dataCard);

  const int nUnfoldedJetPtBins = unfoldedJetPtBins.size();
  for(const auto& bin : unfoldedJetPtBins){
    if(dataCard->FindBinIndexJetPtEEC(bin) < 0){
      cout << "Error! Requested unfolded jet pT bin " << bin.first << "-" << bin.second << " is absent from the input card. Configure unfoldedJetPtBins before unfolding." << endl;
      return;
    }
  }

  bool includeCovariance = false;

  bool saveFigures = false;
  TString saveComment = "_matrixInversion";
  TString figureFormat = "pdf";
    
  // ***************************************************************
  //    Create histogram managers and load the needed histograms
  // ***************************************************************

  // Define a histogram manager for unfolding response histograms
  EECHistogramManager* responseHistograms = new EECHistogramManager(responseInputFile, responseCard);

  // Define a histogram manager for the data histograms that are unfolded
  EECHistogramManager* dataHistograms = new EECHistogramManager(dataInputFile, dataCard);

  // Histograms that are needed to create the unfolding response
  TH1D* hUnfoldingMeasured[nCentralityBins][nTrackPtBins];
  TH1D* hUnfoldingTruth[nCentralityBins][nTrackPtBins];
  TH2D* hUnfoldingResponse[nCentralityBins][nTrackPtBins];
  TH2D* hUnfoldingCovariance[nCentralityBins][nTrackPtBins];
  TH1D* hUnfoldedDistribution[nCentralityBins][nTrackPtBins];
  TH2D* hCovarianceAfterUnfolding[nCentralityBins][nTrackPtBins];
  TH1D* energyEnergyCorrelatorForUnfolding[nCentralityBins][nTrackPtBins];
  TH1D* energyEnergyCorrelatorsFromData[nCentralityBins][nJetPtBins][nTrackPtBins];

  // RooUnfold response object
  RooUnfoldResponse* rooResponse[nCentralityBins][nTrackPtBins];

  // Initialize all the histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      hUnfoldingMeasured[iCentrality][iTrackPt] = NULL;
      hUnfoldingTruth[iCentrality][iTrackPt] = NULL;
      hUnfoldingResponse[iCentrality][iTrackPt] = NULL;
      hUnfoldingCovariance[iCentrality][iTrackPt] = NULL;
      hCovarianceAfterUnfolding[iCentrality][iTrackPt] = NULL;
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt] = NULL;
      rooResponse[iCentrality][iTrackPt] = NULL;
      hUnfoldedDistribution[iCentrality][iTrackPt] = NULL;
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt] = NULL;
      }
    } // Track pT loop
  }  // Centrality loop

  // Read the histograms needed for the unfolding response and create the RooUnfold response objects
  double normalizationFactor;
  std::pair<double,double> jetPtBinBorders;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      hUnfoldingMeasured[iCentrality][iTrackPt] = responseHistograms->GetHistogramJetPtUnfoldingMeasured(iCentrality, iTrackPt);
      hUnfoldingTruth[iCentrality][iTrackPt] = responseHistograms->GetHistogramJetPtUnfoldingTruth(iCentrality, iTrackPt);
      hUnfoldingResponse[iCentrality][iTrackPt] = responseHistograms->GetHistogramJetPtUnfoldingResponse(iCentrality, iTrackPt);

      if(!hUnfoldingMeasured[iCentrality][iTrackPt] || !hUnfoldingTruth[iCentrality][iTrackPt] || !hUnfoldingResponse[iCentrality][iTrackPt]){
        cout << "Error! Missing measured/truth/response histogram for centrality bin " << iCentrality << ", track pT bin " << iTrackPt << ". Cannot unfold." << endl;
        return;
      }
      if(hUnfoldingResponse[iCentrality][iTrackPt]->GetNbinsX() != hUnfoldingMeasured[iCentrality][iTrackPt]->GetNbinsX() || hUnfoldingResponse[iCentrality][iTrackPt]->GetNbinsY() != hUnfoldingTruth[iCentrality][iTrackPt]->GetNbinsX()){
        cout << "Error! Response dimensions differ from measured/truth distributions. Cannot unfold." << endl;
        return;
      }

      // Revome the overflow and underflow bins from the histograms
      removeOutOfRange(hUnfoldingMeasured[iCentrality][iTrackPt]);
      removeOutOfRange(hUnfoldingTruth[iCentrality][iTrackPt]);
      removeOutOfRange(hUnfoldingResponse[iCentrality][iTrackPt]);

      // If we are estimating systematic uncertainty due to limited statistics in MC randomize response matrix, fakes, and misses within their uncertainties
      if(iSystematic == 10){
        cout << "Going to randomize" << endl;
        randomizeWithinUncertainties(hUnfoldingMeasured[iCentrality][iTrackPt]);
        randomizeWithinUncertainties(hUnfoldingTruth[iCentrality][iTrackPt]);
        randomizeWithinUncertainties(hUnfoldingResponse[iCentrality][iTrackPt]);
        cout << "Done!" << endl;
      }

      // Create the response object with trimmed histograms
      rooResponse[iCentrality][iTrackPt] = new RooUnfoldResponse(hUnfoldingMeasured[iCentrality][iTrackPt], hUnfoldingTruth[iCentrality][iTrackPt], hUnfoldingResponse[iCentrality][iTrackPt]);

      // Read the covariance matrix from the data file
      if(includeCovariance) {
        hUnfoldingCovariance[iCentrality][iTrackPt] = dataHistograms->GetHistogramJetPtUnfoldingCovariance(EECHistogramManager::kCovarianceMatrixMeasured, iCentrality, iTrackPt);
        hCovarianceAfterUnfolding[iCentrality][iTrackPt] = (TH2D*) hUnfoldingCovariance[iCentrality][iTrackPt]->Clone(Form("covarianceAfterUnfolding%d%d", iCentrality, iTrackPt));
        hCovarianceAfterUnfolding[iCentrality][iTrackPt]->Reset();
      }

      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt] = dataHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);

      } // Measured jet pT loop
    } // Track pT loop
  } // Centrality loop


  // Next, we need to transform the data histograms into a format that can be read by RooUnfold
  // For this, we will need to combine the jet pT and deltaR axes
  int nDeltaRBinsData = energyEnergyCorrelatorsFromData[firstStudiedCentralityBin][0][firstStudiedTrackPtBinEEC]->GetNbinsX();
  double jetPtLowerBound;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    if(genJetsInData) break; // No unfolding for gen jets
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt] = (TH1D*) hUnfoldingMeasured[iCentrality][iTrackPt]->Clone(Form("dataCorrelatorForUnfolding%d%d", iCentrality, iTrackPt));
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->Reset();
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        jetPtLowerBound = dataCard->GetLowBinBorderJetPtEEC(iJetPt);
        for(int iBin = 1; iBin <= nDeltaRBinsData; iBin++){
          energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->SetBinContent(iBin + nDeltaRBinsData*iJetPt, energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin) * energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinWidth(iBin));
          energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->SetBinError(iBin + nDeltaRBinsData*iJetPt, energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinError(iBin) * energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinWidth(iBin));
        } // DeltaR bin loop
      } // Jet pT loop
    } // Track pT loop
  }  // Centrality loop

  // ********************************************************
  //                     The actual unfolding
  // ********************************************************

  // Error treatment options for unfolding
  // RooUnfolding::kNoErrors;
  // RooUnfolding::kErrors;
  // RooUnfolding::kCovariance;
  // RooUnfolding::kCovToys;

  // Unfolding using Bayesian unfolding
  RooUnfoldBayes* bayesUnfold[nCentralityBins][nTrackPtBins];
  int nMatrixBins = 0;
  //int nIterations = 8;
  if(includeCovariance) nMatrixBins = hUnfoldingCovariance[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetNbinsX();
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    if(genJetsInData) break; // No unfolding for gen jets
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      bayesUnfold[iCentrality][iTrackPt] = new RooUnfoldBayes(rooResponse[iCentrality][iTrackPt], energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt], unfoldConfigurationProvider->GetNumberOfIterations(dataCard->GetBinBordersCentrality(iCentrality), dataCard->GetLowBinBorderTrackPtEEC(iTrackPt)));
      //bayesUnfold[iCentrality][iTrackPt] = new RooUnfoldBayes(rooResponse[iCentrality][iTrackPt], energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt], nIterations);

      // Create a matrix and set it as covariance
      if(includeCovariance){
        TMatrixD covarianceMatrix(nMatrixBins, nMatrixBins);
        for(int iBin = 0; iBin < nMatrixBins; iBin++){
          for(int jBin = 0; jBin < nMatrixBins; jBin++){
            covarianceMatrix(iBin, jBin) = hUnfoldingCovariance[iCentrality][iTrackPt]->GetBinContent(iBin+1, jBin+1);
          }
        }
        bayesUnfold[iCentrality][iTrackPt]->SetMeasuredCov(covarianceMatrix);

        // Do the actual unfolding
        hUnfoldedDistribution[iCentrality][iTrackPt] = (TH1D*)bayesUnfold[iCentrality][iTrackPt]->Hunfold(RooUnfolding::kCovariance);

        // After unfolding, extract the unfolded coveriance matrix
        TMatrixD unfoldedCovarianceMatrix(nMatrixBins, nMatrixBins);
        unfoldedCovarianceMatrix = bayesUnfold[iCentrality][iTrackPt]->GetMeasuredCov();
        
        // Put the information from TMatrixD to TH2D:
        for(int iBin = 0; iBin < nMatrixBins; iBin++){
          for(int jBin = 0; jBin < nMatrixBins; jBin++){
            hCovarianceAfterUnfolding[iCentrality][iTrackPt]->SetBinContent(iBin+1, jBin+1, unfoldedCovarianceMatrix(iBin, jBin));
          }
        }


      } else {

        // Do the actual unfolding
        hUnfoldedDistribution[iCentrality][iTrackPt] = (TH1D*)bayesUnfold[iCentrality][iTrackPt]->Hunfold(RooUnfolding::kErrors);
      }

    } // Track pT loop
  }  // Centrality loop


  // *********************************************************************
  //     Transforming unfolded histograms back to one dimensional ones
  // *********************************************************************

  // Dissect the unfolded histograms in jet pT bins from the big histograms
  TH1D* hUnfolded[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBins];

  // Initialize the dissected histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        hUnfolded[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Find the DeltaR binning within a single jet pT bin from the big histogram
  double previousBinWidth = 0;
  double currentBinWidth = 0;
  int nDeltaRBins = 0;
  for(int iBin = 1; iBin < hUnfoldingTruth[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetNbinsX(); iBin++){
    currentBinWidth = hUnfoldingTruth[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetBinWidth(iBin);
    if(currentBinWidth < previousBinWidth) break;
    previousBinWidth = currentBinWidth;
    nDeltaRBins++;
  }

  double deltaRBinning[nDeltaRBins+1];
  for(int iBin = 1; iBin <= nDeltaRBins+1; iBin++){
    deltaRBinning[iBin-1] = hUnfoldingTruth[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetBinLowEdge(iBin);
  }

  // Create new histograms for the chosen bin range
  int jetPtUnfoldIndex = 0;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        if(genJetsInData){
          // If no unfolding is done, just copy the measured histograms directly for the unfolded collection
          jetPtUnfoldIndex = dataCard->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));
          hUnfolded[iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndex][iTrackPt]->Clone(Form("hUnfolded%d%d%d", iCentrality, iTrackPt, iJetPt));
        } else {
          hUnfolded[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("hUnfolded%d%d%d", iCentrality, iTrackPt, iJetPt), Form("hUnfolded%d%d%d", iCentrality, iTrackPt, iJetPt), nDeltaRBins, deltaRBinning);
        }
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  // Dissect the big histograms and fill small histograms based on the information
  // Also normalize the bins to bin width
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    if(genJetsInData) break; // No unfolding for gen jets
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        jetPtUnfoldIndex = dataCard->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));

        for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
          hUnfolded[iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinContent(iBin + nDeltaRBins*jetPtUnfoldIndex) / hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinWidth(iBin));
          hUnfolded[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinError(iBin + nDeltaRBins*jetPtUnfoldIndex) / hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinWidth(iBin));
        } // DeltaR bin loop

      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop
  
  // *********************************************************************
  //                Save the unfolded histograms to a file
  // *********************************************************************
  
  // First, add the git hash used to unfold the histograms to the card
  int firstUnfoldedJetPtBin = dataCard->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(0));
  int lastUnfoldedJetPtBin = dataCard->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(nUnfoldedJetPtBins-1));
  const char* gitHash = "GITHASHHERE";
  dataCard->AddUnfoldingGitHash(gitHash);
  dataCard->AddFileName(EECCard::kResponseMatrixFile, unfoldConfigurationProvider->GetResponseFileName());
  dataCard->AddOneDimensionalVector(EECCard::kFirstUnfoldedCentralityBin, firstStudiedCentralityBin);
  dataCard->AddOneDimensionalVector(EECCard::kLastUnfoldedCentralityBin, lastStudiedCentralityBin);
  dataCard->AddOneDimensionalVector(EECCard::kFirstUnfoldedTrackPtBin, firstStudiedTrackPtBinEEC);
  dataCard->AddOneDimensionalVector(EECCard::kLastUnfoldedTrackPtBin, lastStudiedTrackPtBinEEC);
  dataCard->AddOneDimensionalVector(EECCard::kFirstUnfoldedJetPtBin, firstUnfoldedJetPtBin);
  dataCard->AddOneDimensionalVector(EECCard::kLastUnfoldedJetPtBin, lastUnfoldedJetPtBin);


  EECHistogramManager* histogramsSavedToFile = new EECHistogramManager(dataCard);
  histogramsSavedToFile->SetLoadEnergyEnergyCorrelators(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelator);
  histogramsSavedToFile->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus);
  histogramsSavedToFile->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus);
  histogramsSavedToFile->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus);
  histogramsSavedToFile->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus);
  histogramsSavedToFile->SetLoadJetPtUnfoldingCovariance(includeCovariance);
  histogramsSavedToFile->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  histogramsSavedToFile->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  histogramsSavedToFile->SetJetPtBinRangeEEC(0, dataCard->GetNJetPtBinsEEC());

  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      if(includeCovariance) histogramsSavedToFile->SetUnfoldedCoverianceMatrix(hCovarianceAfterUnfolding[iCentrality][iTrackPt], iCentrality, iTrackPt);
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        jetPtUnfoldIndex = dataCard->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));
        histogramsSavedToFile->SetUnfoldedEnergyEnergyCorrelator(hUnfolded[iCentrality][iJetPt][iTrackPt], iEnergyEnergyCorrelator, iCentrality, jetPtUnfoldIndex, iTrackPt);
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  histogramsSavedToFile->WriteUnfoldedEnergyEnergyCorrelators(outputFileName, "UPDATE");
  if(includeCovariance) histogramsSavedToFile->WriteCovarianceMatrixAfterUnfolding(outputFileName, "UPDATE");
}

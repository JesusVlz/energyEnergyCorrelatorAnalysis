/*
 * Implementation of the EECUnfoldConfiguration class
 */

// Own includes
#include "EECUnfoldConfiguration.h"

// C++ includes
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

// ROOT includes
#include <TMath.h>

/*
 * Contructor
 */
EECUnfoldConfiguration::EECUnfoldConfiguration():
  fIsPbPbData(true),
  fSystematicIndex(0),
  fSplitIndex(0),
  fWeightExponent(1),
  fExplicitNumberOfIterations(0),
  fResponseMatrixFileName("")

{
  // Just add four to all best iteration numbers
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
      fBestNumberOfIterations[iCentrality][iTrackPt] = 4;
    }
  }

  // Find the actual numbers based on the study 
  InitializeArrays();
}

/*
 * Custom constructor
 */
EECUnfoldConfiguration::EECUnfoldConfiguration(EECCard* card, const int iSplit, const int iSystematic, const int iWeightExponent, const TString responseFileName, const int numberOfIterations):
  fSystematicIndex(iSystematic),
  fSplitIndex(iSplit),
  fWeightExponent(iWeightExponent),
  fExplicitNumberOfIterations(numberOfIterations),
  fResponseMatrixFileName(responseFileName)
{
  if(card == nullptr) throw std::invalid_argument("EECUnfoldConfiguration requires a data card.");
  TString collisionSystem = card->GetDataType();
  fIsPbPbData = collisionSystem.Contains("PbPb");
  
  // Just add four to all best iteration numbers
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
      fBestNumberOfIterations[iCentrality][iTrackPt] = 4;
    }
  }

  // Check that the split number is sensible
  if(fSplitIndex < 0 || fSplitIndex >= kNDatasetSplits){
    std::cout << "EECUnfoldConfiguration::ERROR! Selected split index for the response matrix is out or range!" << std::endl;
    std::cout << "You selected " << fSplitIndex << " while the allowed range is 0-" << kNDatasetSplits-1 << std::endl;
    std::cout << "Next the code will crash. Please fix your code." << std::endl;
    throw std::invalid_argument("iSplit in EECUnfoldConfiguration is out of allowed range!");
  }

  // Check that the weight exponent is sensible
  if(fWeightExponent < 1 || fWeightExponent > 2){
    std::cout << "EECUnfoldConfiguration::ERROR! Selected weight exponent for the response matrix is out or range!" << std::endl;
    std::cout << "You selected " << fWeightExponent << " while currently implemented range is 1-2" << std::endl;
    std::cout << "Next the code will crash. Please fix your code." << std::endl;
    throw std::invalid_argument("iWeightExponent in EECUnfoldConfiguration is out of allowed range!");
  }

  // Run-3 responses and iteration counts must come from a dedicated study.
  // In particular, neither OO nor ppRef may fall through to pp2017 defaults.
  const bool isRun3 = collisionSystem.Contains("OO") || collisionSystem.Contains("ppRef");
  if(isRun3 && fResponseMatrixFileName.IsNull()){
    throw std::invalid_argument("OO/ppRef unfolding requires an explicit Run-3 response ROOT file and a positive iteration count; legacy pp2017 responses are not valid.");
  }
  if(!fResponseMatrixFileName.IsNull()){
    if(fExplicitNumberOfIterations <= 0){
      throw std::invalid_argument("An explicit unfolding response requires a positive, explicitly selected iteration count. Determine it with an independent closure study.");
    }
    ValidateResponseFile(card);
    return;
  }
  if(fExplicitNumberOfIterations != 0){
    throw std::invalid_argument("An explicit unfolding iteration count must be accompanied by an explicit response ROOT file.");
  }

  // Find the actual numbers based on the legacy study
  InitializeArrays();
}

/*
 * Validate a user-supplied processed response before any unfolding is attempted.
 * The selected split/systematic file is supplied by the caller, rather than inferred
 * from legacy filenames. Its iteration count is likewise an explicit study input.
 */
void EECUnfoldConfiguration::ValidateResponseFile(EECCard* card) const{
  TDirectory::TContext directoryContext;
  std::unique_ptr<TFile> responseFile(TFile::Open(fResponseMatrixFileName, "READ"));
  if(!responseFile || responseFile->IsZombie()){
    throw std::invalid_argument("Cannot open explicit unfolding response: " + std::string(fResponseMatrixFileName.Data()));
  }

  // EECCard expects these vectors. Check them before constructing it to avoid
  // a null dereference when a file contains ordinary EECs but no response card.
  const char* requiredEntries[] = {"DataType", "McCorrelationType", "JetType", "JetAxis", "JetRadius", "WeightExponent", "CentralityBinEdges", "JetPtBinEdgesEEC", "TrackPtBinEdgesEEC"};
  for(const char* entry : requiredEntries){
    auto* values = dynamic_cast<TVectorT<float>*>(responseFile->Get(Form("JCard/%s", entry)));
    if(!values || values->GetNoElements() == 0 || values->GetLwb() != 1){
      throw std::invalid_argument("Explicit unfolding response is missing a valid JCard/" + std::string(entry));
    }
  }

  const TString collisionSystem = card->GetDataType();
  const int expectedMCType = collisionSystem.Contains("OO") ? 10 : (collisionSystem.Contains("ppRef") ? 12 : (fIsPbPbData ? 3 : 2));
  auto* responseType = static_cast<TVectorT<float>*>(responseFile->Get("JCard/DataType"));
  auto* correlationType = static_cast<TVectorT<float>*>(responseFile->Get("JCard/McCorrelationType"));
  if((*responseType)[1] != expectedMCType){
    throw std::invalid_argument("Explicit unfolding response must be MC of the same collision system as the data.");
  }
  // FillUnfoldingResponse uses generator particles independently of this mode.
  // Require valid MC metadata, but allow responses produced alongside any MC EEC mode.
  const float correlationMode = (*correlationType)[1];
  if(correlationMode < 0 || correlationMode > 3 || std::floor(correlationMode) != correlationMode){
    throw std::invalid_argument("Explicit unfolding response has an invalid McCorrelationType (expected 0-3).");
  }

  EECCard responseCard(responseFile.get());
  if(responseCard.FindWeightExponentIndex(fWeightExponent) < 1 || card->FindWeightExponentIndex(fWeightExponent) < 1){
    throw std::invalid_argument("Requested EEC weight exponent is absent from the data or explicit response card.");
  }
  const double epsilon = 0.0001;
  if(responseCard.GetJetType() != card->GetJetType() || card->GetJetAxis() < 0 || responseCard.GetJetAxis() < 0 || responseCard.GetJetAxis() != card->GetJetAxis() || std::abs(responseCard.GetJetRadius() - card->GetJetRadius()) > epsilon){
    throw std::invalid_argument("Jet type, axis or radius differs between data and the explicit unfolding response.");
  }
  if(responseCard.GetNJetPtBinsEEC() != card->GetNJetPtBinsEEC() || responseCard.GetNTrackPtBinsEEC() != card->GetNTrackPtBinsEEC()){
    throw std::invalid_argument("Jet/track pT binning differs between data and the explicit unfolding response.");
  }
  for(int iBin = 0; iBin < card->GetNJetPtBinsEEC(); ++iBin){
    if(std::abs(responseCard.GetLowBinBorderJetPtEEC(iBin) - card->GetLowBinBorderJetPtEEC(iBin)) > epsilon || std::abs(responseCard.GetHighBinBorderJetPtEEC(iBin) - card->GetHighBinBorderJetPtEEC(iBin)) > epsilon){
      throw std::invalid_argument("Jet pT bin edges differ between data and the explicit unfolding response.");
    }
  }
  for(int iBin = 0; iBin < card->GetNTrackPtBinsEEC(); ++iBin){
    if(std::abs(responseCard.GetLowBinBorderTrackPtEEC(iBin) - card->GetLowBinBorderTrackPtEEC(iBin)) > epsilon || std::abs(responseCard.GetHighBinBorderTrackPtEEC(iBin) - card->GetHighBinBorderTrackPtEEC(iBin)) > epsilon){
      throw std::invalid_argument("Track pT bin edges differ between data and the explicit unfolding response.");
    }
  }
  if(fIsPbPbData || collisionSystem.Contains("OO")){
    if(responseCard.GetNCentralityBins() != card->GetNCentralityBins()){
      throw std::invalid_argument("Centrality binning differs between data and the explicit unfolding response.");
    }
    for(int iBin = 0; iBin < card->GetNCentralityBins(); ++iBin){
      if(std::abs(responseCard.GetLowBinBorderCentrality(iBin) - card->GetLowBinBorderCentrality(iBin)) > epsilon || std::abs(responseCard.GetHighBinBorderCentrality(iBin) - card->GetHighBinBorderCentrality(iBin)) > epsilon){
        throw std::invalid_argument("Centrality bin edges differ between data and the explicit unfolding response.");
      }
    }
  }
  const char* responseDirectories[] = {"jetPtUnfoldingMeasured", "jetPtUnfoldingTruth", "jetPtUnfoldingResponse"};
  for(const char* name : responseDirectories){
    if(!responseFile->GetDirectory(name)){
      throw std::invalid_argument("Explicit response is missing processed unfolding histograms: " + std::string(name));
    }
  }
}

/*
 * Destructor
 */
EECUnfoldConfiguration::~EECUnfoldConfiguration(){

}

/*
 * Initialization for the number of iterations with Bayesian iteration
 */
void EECUnfoldConfiguration::InitializeArrays(){
  
  // ======================================================== //
  //                                                          //
  //  M     M   EEEEEE   TTTTTTT   H    H    OOOOO    DDDD    //
  //  MMM MMM   E           T      H    H   O     O   D   D   //
  //  M  M  M   EEEEEE      T      HHHHHH   O     O   D    D  //
  //  M     M   E           T      H    H   O     O   D   D   //
  //  M     M   EEEEEE      T      H    H    OOOOO    DDDD    //
  //                                                          //
  // ======================================================== //
  
  // The best number of iterations is determined from doing a chi2 fit to the unfolded to truth ratio.
  // The fit is done separately for each jet pT bin, and then the chi2 values are summed together for each response matrix.
  // For these summer chi2 values, we find a minimum of the distribution.
  // The arrays below show the minimum of these chi2 distributions.
  // There is also the response matrix names for which these numbers are determined included in the configurations.

  // Helper to select response matrix based on the selected Monte Carlo split
  const char* splitName[kNDatasetSplits] = {"", "_part1", "_part2"};

  // Select different response matrx configurations based on the exponent given for the energy weight
  switch(fWeightExponent){

    // ======================================================== //
    //                                                          //
    //    File definitions for regular energy weight pT1*pT2    //
    //                                                          //
    // ======================================================== //

    case 1:

      // The configuration for PbPb
      if(fIsPbPbData){

        if(fSystematicIndex == kNominalSmear || fSystematicIndex == kNumberOfIterationsDown || fSystematicIndex == kNumberOfIterationsUp || fSystematicIndex == kMCStatistics){ 
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for default response matrix
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_split1_nominalSmear_4pCentShift_2023-07-12.root
          //        chi2Histograms_PbPb_split2_nominalSmear_4pCentShift_2023-07-12.root 

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_nominalSmear_responseMatrix%s_processed_2024-01-16.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 4;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 4;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtResolutionUncertaintyDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT resolution uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_split1_uncertaintySmearDown_4pCentShift_2023-07-12.root
          //        chi2Histograms_PbPb_split2_uncertaintySmearDown_4pCentShift_2023-07-12.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_uncertaintySmearDown_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 4;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtResolutionUncertaintyUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT resolution uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_split1_uncertaintySmearUp_4pCentShift_2023-07-12.root
          //        chi2Histograms_PbPb_split2_uncertaintySmearUp_4pCentShift_2023-07-12.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_uncertaintySmearUp_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 3;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 4;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetEnergyScaleUncertaintyDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_split1_minusJECuncertainty_4pCentShift_2023-07-12.root
          //        chi2Histograms_PbPb_split2_minusJECuncertainty_4pCentShift_2023-07-12.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_minusJECuncertainty_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          fBestNumberOfIterations[0][3] = 4;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 4;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 4;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 4;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 4;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 4;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 5;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 4;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetEnergyScaleUncertaintyUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_split1_plusJECuncertainty_4pCentShift_2023-07-12.root
          //        chi2Histograms_PbPb_split2_plusJECuncertainty_4pCentShift_2023-07-12.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_plusJECuncertainty_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 3;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 4;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 4;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtPriorUncertainty){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT prior uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_split1_nominalSmear_jetPtWeight_4pCentShift_2023-07-12.root
          //        chi2Histograms_PbPb_split2_nominalSmear_jetPtWeight_4pCentShift_2023-07-12.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_nominalSmear_jetPtWeight_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 4;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 4;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kCentralityShiftDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for centrality shift uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_split1_nominalSmear_2pCentShift_2023-07-12.root
          //        chi2Histograms_PbPb_split2_nominalSmear_2pCentShift_2023-07-12.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_2pCentShift_cutBadPhi_optimizedUnfoldingBins_nominalSmear_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 3;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kCentralityShiftUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for centrality shift uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_split1_nominalSmear_6pCentShift_2023-07-12.root
          //        chi2Histograms_PbPb_split2_nominalSmear_6pCentShift_2023-07-12.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_6pCentShift_cutBadPhi_optimizedUnfoldingBins_nominalSmear_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 4;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 4;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 4;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 4;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else {
          std::cout << "EECUnfoldConfiguration::ERROR! " << fSystematicIndex  << " is undefined index for unfolding parameter set!" << std::endl;
          std::cout << "Yout code will crash soon if you are reading the response matrix file name from EECUnfoldConfiguration." << std::endl;
          std::cout << "If not, default number 4 is used in all cases. Please check your code and define a number in range 0-" << kNParameterSets-1 << std::endl;
        }

      // The configuration for pp
      } else {

        if(fSystematicIndex == kNominalSmear || fSystematicIndex == kNumberOfIterationsDown || fSystematicIndex == kNumberOfIterationsUp || fSystematicIndex == kMCStatistics){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for default response matrix
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_split1_nominalSmear_2023-06-23.root
          //        chi2Histograms_pp_split2_nominalSmear_2023-06-23.root

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_nominalSmear_responseMatrix%s_processed_2024-01-11.root", splitName[fSplitIndex]);
          //fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_shiftedPt_nominalSmear_responseMatrix%s_processed_2024-02-14.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][3] = 9;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 9;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 9;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtResolutionUncertaintyDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT resolution uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_split1_uncertaintySmearDown_2023-06-23.root
          //        chi2Histograms_pp_split2_uncertaintySmearDown_2023-06-23.root

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_uncertaintySmearDown_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][3] = 7;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 7;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 7;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtResolutionUncertaintyUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT resolution uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_split1_uncertaintySmearUp_2023-06-23.root
          //        chi2Histograms_pp_split2_uncertaintySmearUp_2023-06-23.root

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_uncertaintySmearUp_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][3] = 7;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 8;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 8;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetEnergyScaleUncertaintyDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_split1_minusJECuncertainty_2023-06-23.root
          //        chi2Histograms_pp_split2_minusJECuncertainty_2023-06-23.root

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_minusJECuncertainty_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][3] = 11;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 12;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 14;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetEnergyScaleUncertaintyUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_split1_plusJECuncertainty_2023-06-23.root
          //        chi2Histograms_pp_split2_plusJECuncertainty_2023-06-23.root

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_plusJECuncertainty_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][3] = 7;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 7;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 8;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtPriorUncertainty){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_split1_nominalSmear_jetPtWeight_2023-06-23.root
          //        chi2Histograms_pp_split2_nominalSmear_jetPtWeight_2023-06-23.root

          // TODO: Recheck the iterations for the updated response matrix
          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_nominalSmear_jetPtWeight_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          fBestNumberOfIterations[0][3] = 6;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 6;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 6;  // track pT > 3 GeV

        } else {
          std::cout << "EECUnfoldConfiguration::ERROR! " << fSystematicIndex << " is undefined index for unfolding parameter set!" << std::endl;
          std::cout << "Yout code will crash soon if you are reading the response matrix file name from EECUnfoldConfiguration." << std::endl;
          std::cout << "If not, default number 4 is used in all cases. Please check your code and define a number in range 0-" << kNParameterSets-1 << std::endl;
        }
      }

      break;

    // ============================================================ //
    //                                                              //
    //    File definitions for squared energy weight pT1^2*pT2^2    //
    //                                                              //
    // ============================================================ //

    case 2:

        // The configuration for PbPb
      if(fIsPbPbData){

        if(fSystematicIndex == kNominalSmear || fSystematicIndex == kNumberOfIterationsDown || fSystematicIndex == kNumberOfIterationsUp || fSystematicIndex == kMCStatistics){ 
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for default response matrix
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_energyWeightSquared_nominalReflectedCone_split1_nominalSmear_4pCentShift_2023-12-04.root
          //        chi2Histograms_PbPb_energyWeightSquared_nominalReflectedCone_split2_nominalSmear_4pCentShift_2023-12-04.root 

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_responseMatrix%s_processed_2024-01-18.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 4;  // Centrality = 0-10, track pT > 1 GeV
          fBestNumberOfIterations[1][1] = 4;  // Centrality = 10-30, track pT > 1 GeV
          fBestNumberOfIterations[2][1] = 4;  // Centrality = 30-50, track pT > 1 GeV
          fBestNumberOfIterations[3][1] = 3;  // Centrality = 50-90, track pT > 1 GeV

          fBestNumberOfIterations[0][2] = 4;  // Centrality = 0-10, track pT > 1.5 GeV
          fBestNumberOfIterations[1][2] = 4;  // Centrality = 10-30, track pT > 1.5 GeV
          fBestNumberOfIterations[2][2] = 4;  // Centrality = 30-50, track pT > 1.5 GeV
          fBestNumberOfIterations[3][2] = 3;  // Centrality = 50-90, track pT > 1.5 GeV

          fBestNumberOfIterations[0][3] = 4;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 4;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 4;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 4;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 3;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtResolutionUncertaintyDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT resolution uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_energyWeightSquared_split1_uncertaintySmearDown_4pCentShift_2023-11-16.root
          //        chi2Histograms_PbPb_energyWeightSquared_split2_uncertaintySmearDown_4pCentShift_2023-11-16.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_uncertaintySmearDown_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // Centrality = 0-10, track pT > 1 GeV
          fBestNumberOfIterations[1][1] = 3;  // Centrality = 10-30, track pT > 1 GeV
          fBestNumberOfIterations[2][1] = 3;  // Centrality = 30-50, track pT > 1 GeV
          fBestNumberOfIterations[3][1] = 3;  // Centrality = 50-90, track pT > 1 GeV

          fBestNumberOfIterations[0][2] = 3;  // Centrality = 0-10, track pT > 1.5 GeV
          fBestNumberOfIterations[1][2] = 3;  // Centrality = 10-30, track pT > 1.5 GeV
          fBestNumberOfIterations[2][2] = 3;  // Centrality = 30-50, track pT > 1.5 GeV
          fBestNumberOfIterations[3][2] = 3;  // Centrality = 50-90, track pT > 1.5 GeV

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 3;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 3;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtResolutionUncertaintyUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT resolution uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_energyWeightSquared_split1_uncertaintySmearUp_4pCentShift_2023-11-16.root
          //        chi2Histograms_PbPb_energyWeightSquared_split2_uncertaintySmearUp_4pCentShift_2023-11-16.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_uncertaintySmearUp_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // Centrality = 0-10, track pT > 1 GeV
          fBestNumberOfIterations[1][1] = 3;  // Centrality = 10-30, track pT > 1 GeV
          fBestNumberOfIterations[2][1] = 4;  // Centrality = 30-50, track pT > 1 GeV
          fBestNumberOfIterations[3][1] = 3;  // Centrality = 50-90, track pT > 1 GeV

          fBestNumberOfIterations[0][2] = 3;  // Centrality = 0-10, track pT > 1.5 GeV
          fBestNumberOfIterations[1][2] = 3;  // Centrality = 10-30, track pT > 1.5 GeV
          fBestNumberOfIterations[2][2] = 4;  // Centrality = 30-50, track pT > 1.5 GeV
          fBestNumberOfIterations[3][2] = 3;  // Centrality = 50-90, track pT > 1.5 GeV

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 3;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetEnergyScaleUncertaintyDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_energyWeightSquared_split1_minusJECuncertainty_4pCentShift_2023-11-16.root
          //        chi2Histograms_PbPb_energyWeightSquared_split2_minusJECuncertainty_4pCentShift_2023-11-16.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_minusJECuncertainty_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // Centrality = 0-10, track pT > 1 GeV
          fBestNumberOfIterations[1][1] = 4;  // Centrality = 10-30, track pT > 1 GeV
          fBestNumberOfIterations[2][1] = 4;  // Centrality = 30-50, track pT > 1 GeV
          fBestNumberOfIterations[3][1] = 4;  // Centrality = 50-90, track pT > 1 GeV

          fBestNumberOfIterations[0][2] = 3;  // Centrality = 0-10, track pT > 1.5 GeV
          fBestNumberOfIterations[1][2] = 4;  // Centrality = 10-30, track pT > 1.5 GeV
          fBestNumberOfIterations[2][2] = 4;  // Centrality = 30-50, track pT > 1.5 GeV
          fBestNumberOfIterations[3][2] = 4;  // Centrality = 50-90, track pT > 1.5 GeV

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 4;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 4;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 4;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 4;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 3;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 4;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetEnergyScaleUncertaintyUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_energyWeightSquared_split1_plusJECuncertainty_4pCentShift_2023-11-16.root
          //        chi2Histograms_PbPb_energyWeightSquared_split2_plusJECuncertainty_4pCentShift_2023-11-16.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_plusJECuncertainty_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // Centrality = 0-10, track pT > 1 GeV
          fBestNumberOfIterations[1][1] = 3;  // Centrality = 10-30, track pT > 1 GeV
          fBestNumberOfIterations[2][1] = 4;  // Centrality = 30-50, track pT > 1 GeV
          fBestNumberOfIterations[3][1] = 3;  // Centrality = 50-90, track pT > 1 GeV

          fBestNumberOfIterations[0][2] = 3;  // Centrality = 0-10, track pT > 1.5 GeV
          fBestNumberOfIterations[1][2] = 3;  // Centrality = 10-30, track pT > 1.5 GeV
          fBestNumberOfIterations[2][2] = 4;  // Centrality = 30-50, track pT > 1.5 GeV
          fBestNumberOfIterations[3][2] = 3;  // Centrality = 50-90, track pT > 1.5 GeV

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 3;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtPriorUncertainty){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT prior uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_energyWeightSquared_split1_nominalSmear_jetPtWeight_4pCentShift_2023-11-16.root
          //        chi2Histograms_PbPb_energyWeightSquared_split2_nominalSmear_jetPtWeight_4pCentShift_2023-11-16.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_jetPtWeight_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // Centrality = 0-10, track pT > 1 GeV
          fBestNumberOfIterations[1][1] = 3;  // Centrality = 10-30, track pT > 1 GeV
          fBestNumberOfIterations[2][1] = 4;  // Centrality = 30-50, track pT > 1 GeV
          fBestNumberOfIterations[3][1] = 3;  // Centrality = 50-90, track pT > 1 GeV

          fBestNumberOfIterations[0][2] = 3;  // Centrality = 0-10, track pT > 1.5 GeV
          fBestNumberOfIterations[1][2] = 4;  // Centrality = 10-30, track pT > 1.5 GeV
          fBestNumberOfIterations[2][2] = 4;  // Centrality = 30-50, track pT > 1.5 GeV
          fBestNumberOfIterations[3][2] = 3;  // Centrality = 50-90, track pT > 1.5 GeV

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 4;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 3;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 4;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 3;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kCentralityShiftDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for centrality shift uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_energyWeightSquared_split1_nominalSmear_2pCentShift_2023-11-16.root
          //        chi2Histograms_PbPb_energyWeightSquared_split2_nominalSmear_2pCentShift_2023-11-16.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_2pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // Centrality = 0-10, track pT > 1 GeV
          fBestNumberOfIterations[1][1] = 3;  // Centrality = 10-30, track pT > 1 GeV
          fBestNumberOfIterations[2][1] = 4;  // Centrality = 30-50, track pT > 1 GeV
          fBestNumberOfIterations[3][1] = 3;  // Centrality = 50-90, track pT > 1 GeV

          fBestNumberOfIterations[0][2] = 3;  // Centrality = 0-10, track pT > 1.5 GeV
          fBestNumberOfIterations[1][2] = 3;  // Centrality = 10-30, track pT > 1.5 GeV
          fBestNumberOfIterations[2][2] = 4;  // Centrality = 30-50, track pT > 1.5 GeV
          fBestNumberOfIterations[3][2] = 3;  // Centrality = 50-90, track pT > 1.5 GeV

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 4;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 3;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else if(fSystematicIndex == kCentralityShiftUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for centrality shift uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_PbPb_energyWeightSquared_split1_nominalSmear_6pCentShift_2023-11-16.root
          //        chi2Histograms_PbPb_energyWeightSquared_split2_nominalSmear_6pCentShift_2023-11-16.root

          fResponseMatrixFileName = Form("data/unfolding/PbPbMC2018_GenGen_akFlowJets_6pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // Centrality = 0-10, track pT > 1 GeV
          fBestNumberOfIterations[1][1] = 3;  // Centrality = 10-30, track pT > 1 GeV
          fBestNumberOfIterations[2][1] = 3;  // Centrality = 30-50, track pT > 1 GeV
          fBestNumberOfIterations[3][1] = 3;  // Centrality = 50-90, track pT > 1 GeV

          fBestNumberOfIterations[0][2] = 3;  // Centrality = 0-10, track pT > 1.5 GeV
          fBestNumberOfIterations[1][2] = 3;  // Centrality = 10-30, track pT > 1.5 GeV
          fBestNumberOfIterations[2][2] = 3;  // Centrality = 30-50, track pT > 1.5 GeV
          fBestNumberOfIterations[3][2] = 3;  // Centrality = 50-90, track pT > 1.5 GeV

          fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
          fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
          fBestNumberOfIterations[2][3] = 3;  // Centrality = 30-50, track pT > 2 GeV
          fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

          fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
          fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
          fBestNumberOfIterations[2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
          fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

          fBestNumberOfIterations[0][5] = 3;  // Centrality = 0-10, track pT > 3 GeV
          fBestNumberOfIterations[1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV
          fBestNumberOfIterations[2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV
          fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

        } else {
          std::cout << "EECUnfoldConfiguration::ERROR! " << fSystematicIndex  << " is undefined index for unfolding parameter set!" << std::endl;
          std::cout << "Yout code will crash soon if you are reading the response matrix file name from EECUnfoldConfiguration." << std::endl;
          std::cout << "If not, default number 4 is used in all cases. Please check your code and define a number in range 0-" << kNParameterSets-1 << std::endl;
        }

      // The configuration for pp
      } else {

        if(fSystematicIndex == kNominalSmear || fSystematicIndex == kNumberOfIterationsDown || fSystematicIndex == kNumberOfIterationsUp || fSystematicIndex == kMCStatistics){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for default response matrix
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_energyWeightSquared_split1_nominalSmear_2023-10-31.root
          //        chi2Histograms_pp_energyWeightSquared_split2_nominalSmear_2023-10-31.root
          // To determine these numbers, the 120 < jet pT < 140 GeV bin is ignored, because the results differ a lot from all other bins

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_responseMatrix%s_processed_2024-01-10.root", splitName[fSplitIndex]);
          //fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_shiftedPt_energyWeightSquared_nominalSmear_responseMatrix%s_processed_2024-02-14.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // track pT > 1 GeV
          fBestNumberOfIterations[0][2] = 3;  // track pT > 1.5 GeV
          fBestNumberOfIterations[0][3] = 3;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 3;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 3;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtResolutionUncertaintyDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT resolution uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_energyWeightSquared_split1_uncertaintySmearDown_2023-11-10.root
          //        chi2Histograms_pp_energyWeightSquared_split2_uncertaintySmearDown_2023-11-10.root
          // To determine these numbers, the 120 < jet pT < 140 GeV bin is ignored, because the results differ a lot from all other bins

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_uncertaintySmearDown_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // track pT > 1 GeV
          fBestNumberOfIterations[0][2] = 3;  // track pT > 1.5 GeV
          fBestNumberOfIterations[0][3] = 3;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 3;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 3;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtResolutionUncertaintyUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet pT resolution uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_energyWeightSquared_split1_uncertaintySmearUp_2023-11-10.root
          //        chi2Histograms_pp_energyWeightSquared_split2_uncertaintySmearUp_2023-11-10.root
          // To determine these numbers, the 120 < jet pT < 140 GeV bin is ignored, because the results differ a lot from all other bins

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_uncertaintySmearUp_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 4;  // track pT > 1 GeV
          fBestNumberOfIterations[0][2] = 4;  // track pT > 1.5 GeV
          fBestNumberOfIterations[0][3] = 4;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 4;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 4;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetEnergyScaleUncertaintyDown){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_energyWeightSquared_split1_minusJECuncertainty_2023-11-10.root
          //        chi2Histograms_pp_energyWeightSquared_split2_minusJECuncertainty_2023-11-10.root
          // To determine these numbers, the 120 < jet pT < 140 GeV bin is ignored, because the results differ a lot from all other bins

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_minusJECuncertainty_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // track pT > 1 GeV
          fBestNumberOfIterations[0][2] = 3;  // track pT > 1.5 GeV
          fBestNumberOfIterations[0][3] = 3;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 3;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 3;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetEnergyScaleUncertaintyUp){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_energyWeightSquared_split1_plusJECuncertainty_2023-11-10.root
          //        chi2Histograms_pp_energyWeightSquared_split2_plusJECuncertainty_2023-11-10.root
          // To determine these numbers, the 120 < jet pT < 140 GeV bin is ignored, because the results differ a lot from all other bins

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_plusJECuncertainty_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // track pT > 1 GeV
          fBestNumberOfIterations[0][2] = 3;  // track pT > 1.5 GeV
          fBestNumberOfIterations[0][3] = 3;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 3;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 3;  // track pT > 3 GeV

        } else if(fSystematicIndex == kJetPtPriorUncertainty){
          // TODO: Response matrix updated, still need to double check number of iterations
          // Configuration for jet energy scale uncertainty evaluation
          // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
          // Input: chi2Histograms_pp_energyWeightSquared_split1_nominalSmear_jetPtWeight_2023-11-10.root
          //        chi2Histograms_pp_energyWeightSquared_split2_nominalSmear_jetPtWeight_2023-11-10.root
          // To determine these numbers, the 120 < jet pT < 140 GeV bin is ignored, because the results differ a lot from all other bins

          fResponseMatrixFileName = Form("data/unfolding/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_jetPtWeight_responseMatrix%s_processed_2024-01-12.root", splitName[fSplitIndex]);

          // TODO: Recheck the iterations for the updated response matrix
          fBestNumberOfIterations[0][1] = 3;  // track pT > 1 GeV
          fBestNumberOfIterations[0][2] = 3;  // track pT > 1.5 GeV
          fBestNumberOfIterations[0][3] = 3;  // track pT > 2 GeV
          fBestNumberOfIterations[0][4] = 3;  // track pT > 2.5 GeV
          fBestNumberOfIterations[0][5] = 3;  // track pT > 3 GeV

        } else {
          std::cout << "EECUnfoldConfiguration::ERROR! " << fSystematicIndex << " is undefined index for unfolding parameter set!" << std::endl;
          std::cout << "Yout code will crash soon if you are reading the response matrix file name from EECUnfoldConfiguration." << std::endl;
          std::cout << "If not, default number 4 is used in all cases. Please check your code and define a number in range 0-" << kNParameterSets-1 << std::endl;
        }

      }

      break;

    default:

      std::cout << "EECUnfoldConfiguration::ERROR!" << std::endl;
      std::cout << "It should be impossible to reach this message. Something in the code is very wrong." << std::endl;
      std::cout << "Good luck debugging! Hopefully we will never meet again!" << std::endl;

  } // End of switch-case for different energy weight exponents
  
} 

/*
 * Get the optimal number of unfolding iterations for the given bin
 *
 *  const std::pair<double,double> centralityBinBorders = Centrality bin borders for the scaling factor
 *  double trackPtBorderLow = Lower track pT border for the studied track pT bin
 *
 *  return: Scaling factor corresponding to the bin with the given bin borders
 */
double EECUnfoldConfiguration::GetNumberOfIterations(const std::pair<double,double> centralityBinBorders, double trackPtBorderLow) const{

  if(fExplicitNumberOfIterations > 0) return fExplicitNumberOfIterations;

  // For now, just return 4. I will redo the logic later
  if(fSystematicIndex == kNumberOfIterationsDown) return 3;
  if(fSystematicIndex == kNumberOfIterationsUp) return 5;
  return 4;

  // ******************************************************************** //
  // First, find the bin indices that correspond to the input bin borders //
  // ******************************************************************** //

  // Small number
  const double epsilon = 0.0001;

  // Search if the given centrality bin borders are included in the background scale table
  int centralityIndex = -1;
  if(fIsPbPbData){
    for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
      if(TMath::Abs(fCentralityBinBorderLow[iCentrality] - centralityBinBorders.first) < epsilon){
        if(TMath::Abs(fCentralityBinBorderHigh[iCentrality] - centralityBinBorders.second) < epsilon){
        centralityIndex = iCentrality;
          break;
        }
      }

      if(TMath::Abs(fCentralityBinBorderLowShifted[iCentrality] - centralityBinBorders.first) < epsilon){
        if(TMath::Abs(fCentralityBinBorderHighShifted[iCentrality] - centralityBinBorders.second) < epsilon){
        centralityIndex = iCentrality;
          break;
        }
      }

      // We use 2% and 6% sentrality shifts for systematic estimation. Take those into account when finding centrality bin
      if(TMath::Abs(fCentralityBinBorderLowShifted[iCentrality] - centralityBinBorders.first - 2) < epsilon){
        if(TMath::Abs(fCentralityBinBorderHighShifted[iCentrality] - centralityBinBorders.second - 2) < epsilon){
        centralityIndex = iCentrality;
          break;
        }
      }

      if(TMath::Abs(fCentralityBinBorderLowShifted[iCentrality] - centralityBinBorders.first + 2) < epsilon){
        if(TMath::Abs(fCentralityBinBorderHighShifted[iCentrality] - centralityBinBorders.second + 2) < epsilon){
        centralityIndex = iCentrality;
          break;
        }
      }
    } // Loop over centrality bins included in the scaling table

    // If centrality bin is not found, print an error message and return -1 to show we did not find a proper scaling factor.
    if(centralityIndex == -1){
      std::cout << "EECUnfoldConfiguration::Error! Centrality bin " << centralityBinBorders.first << "-" << centralityBinBorders.second << " % was not found from the table." << std::endl;
      std::cout << "Cannot provide a background scaling factor. Please check your code for errors." << std::endl;
      return -1;
    }
  } else {
    centralityIndex = 0;
  }

  // Search if the given track pT bin borders are included in the background scale table
  int trackPtIndex = -1;
  for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
    if(TMath::Abs(fTrackPtBinBorderLow[iTrackPt] - trackPtBorderLow) < epsilon){
      trackPtIndex = iTrackPt;
      break;
    }
  } // Loop over track pT bins included in the scaling table

  // If track pT bin is not found, print an error message and return -1 to show we did not find a proper scaling factor.
  if(trackPtIndex == -1){
    std::cout << "EECUnfoldConfiguration::Error! Track pT lower border " << trackPtBorderLow << " GeV was not found from the table." << std::endl;
    std::cout << "Cannot provide a background scaling factor. Please check your code for errors." << std::endl;
    return -1;
  }

  // *********************************************************************************** //
  // After the correct bin has been found, return the corresponding number of iterations //
  // *********************************************************************************** //

  return fBestNumberOfIterations[centralityIndex][trackPtIndex];

}

// Getter for the name of the file containing the response matrices
TString EECUnfoldConfiguration::GetResponseFileName() const{
  return fResponseMatrixFileName;
}

// Select complete native Run-3 bins within the existing 120-200 GeV study window.
std::vector<std::pair<double,double>> EECUnfoldConfiguration::SelectUnfoldedJetPtBins(const EECCard& card){
  const TString collisionSystem = card.GetDataType();
  if(!collisionSystem.Contains("OO") && !collisionSystem.Contains("ppRef")){
    return {{120,140}, {140,160}, {160,180}, {180,200}};
  }

  std::vector<std::pair<double,double>> selectedBins;
  const double epsilon = 0.0001;
  for(int iBin = 0; iBin < card.GetNJetPtBinsEEC(); ++iBin){
    const auto bin = card.GetBinBordersJetPtEEC(iBin);
    if(bin.first >= 120 - epsilon && bin.second <= 200 + epsilon && bin.second > bin.first){
      selectedBins.push_back(bin);
    }
  }
  if(selectedBins.empty()){
    throw std::invalid_argument("No complete Run-3 jet pT bins lie inside the 120-200 GeV unfolding study window. Configure the input card or study window before unfolding.");
  }
  return selectedBins;
}

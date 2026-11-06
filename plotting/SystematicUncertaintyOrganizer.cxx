/*
 * Implementation of SystematicUncertaintyOrganizer
 */

#include "SystematicUncertaintyOrganizer.h"

#include <stdexcept>

/*
 * Default constructor
 */
SystematicUncertaintyOrganizer::SystematicUncertaintyOrganizer()
{
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] = NULL;
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Setup the flags for systematics grouping
  fSystematicsGroupFlag[kJetEnergyResolution] = kCorrelatedInDeltaR;
  fSystematicsGroupFlag[kJetEnergyScale] = kCorrelatedInDeltaR;
  fSystematicsGroupFlag[kUnfoldingTruth] = kCorrelatedInDeltaR;
  fSystematicsGroupFlag[kUnfoldingIterations] = kCorrelatedInDeltaR;
  fSystematicsGroupFlag[kTrackSelection] = kUncorrelatedInDeltaR;
  fSystematicsGroupFlag[kSingleTrackEfficiency] = kSkipped;
  fSystematicsGroupFlag[kTrackPairEfficiency] = kUncorrelatedInDeltaR;
  fSystematicsGroupFlag[kBackgroundSubtraction] = kUncorrelatedInDeltaR;
  fSystematicsGroupFlag[kSignalToBackgroundRatio] = kCorrelatedInDeltaR;
  fSystematicsGroupFlag[kCentralityShift] = kCorrelatedInDeltaR;
  fSystematicsGroupFlag[kMonteCarloStatistics] = kUncorrelatedInDeltaR;
  fSystematicsGroupFlag[kMonteCarloNonClosure] = kCorrelatedInDeltaR;
  fSystematicsGroupFlag[kAll] = kGroupForAll;

  // Setup the standardized color scheme for different uncertainty sources
  fUncertaintyColor[kJetEnergyResolution] = kBlue;
  fUncertaintyColor[kJetEnergyScale] = kRed;
  fUncertaintyColor[kUnfoldingTruth] = kGreen+3;
  fUncertaintyColor[kUnfoldingIterations] = kGray+1;
  fUncertaintyColor[kTrackSelection] = kMagenta;
  fUncertaintyColor[kSingleTrackEfficiency] = kWhite;
  fUncertaintyColor[kTrackPairEfficiency] = kCyan;
  fUncertaintyColor[kBackgroundSubtraction] = kViolet-6;
  fUncertaintyColor[kSignalToBackgroundRatio] = kYellow+1;
  fUncertaintyColor[kCentralityShift] = kSpring;
  fUncertaintyColor[kMonteCarloStatistics] = kAzure-7;
  fUncertaintyColor[kMonteCarloNonClosure] = kOrange+7;
  fUncertaintyColor[kAll] = kBlack;

  // Define which sources are relevant for pp
  fIsRelevant[0][kJetEnergyResolution] = true;
  fIsRelevant[0][kJetEnergyScale] = true;
  fIsRelevant[0][kUnfoldingTruth] = true;
  fIsRelevant[0][kUnfoldingIterations] = true;
  fIsRelevant[0][kTrackSelection] = true;
  fIsRelevant[0][kSingleTrackEfficiency] = false;
  fIsRelevant[0][kTrackPairEfficiency] = true;
  fIsRelevant[0][kBackgroundSubtraction] = true;
  fIsRelevant[0][kSignalToBackgroundRatio] = false;
  fIsRelevant[0][kCentralityShift] = false;
  fIsRelevant[0][kMonteCarloStatistics] = true;
  fIsRelevant[0][kMonteCarloNonClosure] = false;
  fIsRelevant[0][kAll] = true;

  // Define which sources are relevant for PbPb
  fIsRelevant[1][kJetEnergyResolution] = true;
  fIsRelevant[1][kJetEnergyScale] = true;
  fIsRelevant[1][kUnfoldingTruth] = true;
  fIsRelevant[1][kUnfoldingIterations] = true;
  fIsRelevant[1][kTrackSelection] = true;
  fIsRelevant[1][kSingleTrackEfficiency] = false;
  fIsRelevant[1][kTrackPairEfficiency] = true;
  fIsRelevant[1][kBackgroundSubtraction] = true;
  fIsRelevant[1][kSignalToBackgroundRatio] = true;
  fIsRelevant[1][kCentralityShift] = true;
  fIsRelevant[1][kMonteCarloStatistics] = true;
  fIsRelevant[1][kMonteCarloNonClosure] = false;
  fIsRelevant[1][kAll] = true;
}

/*
 * Constructor
 */
SystematicUncertaintyOrganizer::SystematicUncertaintyOrganizer(TFile* inputFile) :
  SystematicUncertaintyOrganizer()
{
  ReadInputFile(inputFile);
}

/*
 * Copy constructor
 */
SystematicUncertaintyOrganizer::SystematicUncertaintyOrganizer(const SystematicUncertaintyOrganizer& in)
{
  // Copy constructor
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] = in.fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty];
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
}

// Read the input file containing the uncertainty histograms
void SystematicUncertaintyOrganizer::ReadInputFile(TFile* inputFile){
  
  // The systematic uncertainty estimation is done for all unfolded bins. Read the indices for those from the card
  EECCard* card = new EECCard(inputFile);
  
  int firstStudiedCentralityBin = card->GetFirstUnfoldedCentralityBin();
  int lastStudiedCentralityBin = card->GetLastUnfoldedCentralityBin();
  
  int firstStudiedJetPtBinEEC = card->GetFirstUnfoldedJetPtBin();
  int lastStudiedJetPtBinEEC = card->GetLastUnfoldedJetPtBin();
  
  int firstStudiedTrackPtBinEEC = card->GetFirstUnfoldedTrackPtBin();
  int lastStudiedTrackPtBinEEC = card->GetLastUnfoldedTrackPtBin();

  // Helper variables
  TString saveName;
  
  // Read the uncertainties from the file
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          saveName = Form("systematicUncertainty_%s_C%dJ%dT%d", fSystematicUncertaintyName[iUncertainty].Data(), iCentrality, iJetPt, iTrackPt);
          fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] = (TH1D*) inputFile->Get(saveName);
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // ============================================================================================= //
  // For all the uncertainty sources that are not defined in the file, set the uncertainty to zero //
  // ============================================================================================= //

  // First, find at least one histogram that is not NULL
  TH1D* exampleHistogram = nullptr;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin && exampleHistogram == nullptr; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC && exampleHistogram == nullptr; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC && exampleHistogram == nullptr; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          if(fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] != NULL){
            exampleHistogram = fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty];
            break;
          }
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Without a template histogram there is no binning for the zero-error fallback.
  if(exampleHistogram == nullptr){
    throw std::runtime_error(Form("SystematicUncertaintyOrganizer::ReadInputFile: no systematic uncertainty histograms found in '%s' for C[%d,%d], J[%d,%d], T[%d,%d]; cannot infer histogram binning",
                                  inputFile->GetName(), firstStudiedCentralityBin, lastStudiedCentralityBin,
                                  firstStudiedJetPtBinEEC, lastStudiedJetPtBinEEC,
                                  firstStudiedTrackPtBinEEC, lastStudiedTrackPtBinEEC));
  }

  // Find all histograms that are NULL, and set them as error zero histograms
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          if(fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] == NULL){
            fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] = (TH1D*) exampleHistogram->Clone(Form("zeroErrorHistogram%d%d%d%d", iCentrality, iJetPt, iTrackPt, iUncertainty));
            for(int iBin = 1; iBin <= exampleHistogram->GetNbinsX(); iBin++){
              fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty]->SetBinError(iBin, 0);
            } 
          }
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
}

// Getter for absolute systematic uncertainty for energy-energy correlators
TH1D* SystematicUncertaintyOrganizer::GetSystematicUncertainty(const int iCentrality, const int iJetPt, const int iTrackPt, const int iUncertainty) const{
  
  // Return the uncertainty in the selected bin
  return fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty];
}

// Combine a predefined group of systematic uncertainty sources
TH1D* SystematicUncertaintyOrganizer::CombineUncertaintySources(const int iCentrality, const int iJetPt, const int iTrackPt, const int iGroup, const char* newName) const{
  TH1D* uncertaintyHistogram = (TH1D*) fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][kAll]->Clone(Form("%s%d%d%d", newName, iCentrality, iJetPt, iTrackPt));

  // Calculate a sum of squared from all evaluated uncertainties
  double sumOfSquares = 0;
  for(int iBin = 1; iBin <= uncertaintyHistogram->GetNbinsX(); iBin++){
    sumOfSquares = 0;
    for(int iUncertainty = 0; iUncertainty < kAll; iUncertainty++){
      if(fSystematicsGroupFlag[iUncertainty] == iGroup){
        sumOfSquares = sumOfSquares + TMath::Power(fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty]->GetBinError(iBin), 2);
      }
    }  // Uncertainty type loop
    uncertaintyHistogram->SetBinError(iBin, TMath::Sqrt(sumOfSquares));
  }  // Bin loop

  return uncertaintyHistogram;
}

// Getter for the absolute systematic uncertainty from all sources correlated in DeltaR
TH1D* SystematicUncertaintyOrganizer::GetCorrelatedSystematicUncertainty(const int iCentrality, const int iJetPt, const int iTrackPt) const{
  return CombineUncertaintySources(iCentrality, iJetPt, iTrackPt, kCorrelatedInDeltaR, "correlatedUncertainties");
}

// Getter for the absolute systematic uncertainty from all sources not correlated in DeltaR
TH1D* SystematicUncertaintyOrganizer::GetUncorrelatedSystematicUncertainty(const int iCentrality, const int iJetPt, const int iTrackPt) const{
  return CombineUncertaintySources(iCentrality, iJetPt, iTrackPt, kUncorrelatedInDeltaR, "uncorrelatedUncertainties");
}

// Getter for a name for the source of systematic uncertainty
TString SystematicUncertaintyOrganizer::GetSystematicUncertaintyName(const int iUncertainty) const{
  return fSystematicUncertaintyName[iUncertainty];
}

// Getter for a name suitable for a legend for the source of systematic uncertainty
TString SystematicUncertaintyOrganizer::GetSystematicUncertaintyLegendName(const int iUncertainty) const{
  return fSystematicUncertaintyLegendName[iUncertainty];
}


// Getter for an axis name for the source of systematic uncertainty
TString SystematicUncertaintyOrganizer::GetUncertaintyAxisName(const int iUncertainty) const{
  return fUncertaintyAxisName[iUncertainty];
}

// Getter for an axis name for the source of systematic uncertainty
int SystematicUncertaintyOrganizer::GetNUncertaintySources() const{
  return knUncertaintySources;
}

// Getter for standardized color for each systematic uncertainty source
int SystematicUncertaintyOrganizer::GetUncertaintyColor(const int iUncertainty) const{
  return fUncertaintyColor[iUncertainty];
}

// Getter for information if a systematic uncertainty is relevant for pp
bool SystematicUncertaintyOrganizer::GetSystematicUncertaintyRelevancyForPp(const int iUncertainty) const{
  return GetSystematicUncertaintyRelevancy(iUncertainty, false);
}

// Getter for information if a systematic uncertainty is relevant for PbPb
bool SystematicUncertaintyOrganizer::GetSystematicUncertaintyRelevancyForPbPb(const int iUncertainty) const{
  return GetSystematicUncertaintyRelevancy(iUncertainty, true);
}

// Getter for information if a systematic uncertainty is relevant for a system (pp or PbPb)
bool SystematicUncertaintyOrganizer::GetSystematicUncertaintyRelevancy(const int iUncertainty, const bool isPbPb) const{
  return fIsRelevant[isPbPb][iUncertainty];
}
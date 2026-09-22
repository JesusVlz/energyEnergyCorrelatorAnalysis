/*
 * Implementation of the TrackPairEfficiencyCorrector class
 */

// Own includes
#include "TrackPairEfficiencyCorrector.h"

/*
 * Contructor
 */
TrackPairEfficiencyCorrector::TrackPairEfficiencyCorrector() :
  fCard(NULL),
  fInputFile(NULL),
  fnCentralityBins(0),
  fnTrackPtBins(0),
  fnJetPtBins(0),
  fDisableCorrection(false),
  fUseSmoothedCorrection(false),
  fIsPbPbData(false)
{
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < kMaxTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt < kMaxTrackPtBins; iAssociatedPt++){
        fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt] = NULL;
        for(int iJetPt = 0; iJetPt < kMaxJetPtBins; iJetPt++){
          fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = NULL;
        } // Jet pT loop
      } // Associated pT loop
    } // Trigger pT loop
  } // Centrality loop
}

/*
 * Custom constructor
 */
TrackPairEfficiencyCorrector::TrackPairEfficiencyCorrector(TString inputFileName, bool useSmoothedCorrection, bool isPbPbData) :
  fDisableCorrection(false)
{
  fInputFile = TFile::Open(inputFileName);
  fUseSmoothedCorrection = useSmoothedCorrection;
  fIsPbPbData = isPbPbData;
  ReadCorrectionTables();
}

/*
 * Read the track pair efficiency correction tables from the file
 */
void TrackPairEfficiencyCorrector::ReadCorrectionTables(){
  
  // Read card from inputfile
  fCard = new BinningCard(fInputFile);
  
  // Read the dimensions of the arrays from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  fnTrackPtBins = fCard->GetNTrackPairPtBins();
  fnJetPtBins = fCard->GetNJetPtBins();
  
  // Read the centrality bin borders from the card
  for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fCentralityBinBorders[fnCentralityBins] = fCard->GetHighBinBorderCentrality(fnCentralityBins-1);
  
  // Read the track pT bin borders from the card
  for(int iTrackPt = 0; iTrackPt < fnTrackPtBins; iTrackPt++){
    fTrackPtBinBorders[iTrackPt] = fCard->GetLowBinBorderTrackPairPt(iTrackPt);
  }
  fTrackPtBinBorders[fnTrackPtBins] = fCard->GetHighBinBorderTrackPairPt(fnTrackPtBins-1);

  // Read the jet pT bin borders from the card
  for(int iJetPt = 0; iJetPt < fnJetPtBins; iJetPt++){
    fJetPtBinBorders[iJetPt] = fCard->GetLowBinBorderJetPt(iJetPt);
  }
  fJetPtBinBorders[fnJetPtBins] = fCard->GetHighBinBorderJetPt(fnJetPtBins-1);
  
  // Select the correct name for the corrections
  const char* smoothString = fUseSmoothedCorrection ? "smoothedT" : "t";

  // Read all the tables from the file
  for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < fnTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
        fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt] = (TH1D*) fInputFile->Get(Form("%srackPairEfficiencyCorrection_C%dT%dA%d", smoothString, iCentrality, iTriggerPt, iAssociatedPt));
        for(int iJetPt = 0; iJetPt < fnJetPtBins; iJetPt++){
          fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = (TH1D*) fInputFile->Get(Form("%srackPairEfficiencyCorrectionCloseToJet_C%dT%dA%dJ%d", smoothString, iCentrality, iTriggerPt, iAssociatedPt, iJetPt));
        }
      }
    }
  }
  
}

/*
 * Find a bin index from an array
 */
int TrackPairEfficiencyCorrector::FindBinIndex(const double* array, const int nBins, const double value) const{
  
  // If the value is below the lowest index, just use the lowest index
  if(value < array[0]) return 0;
  
  // Find the value from the arrey and return the bin index
  for(int iBin = 0; iBin < nBins; iBin++){
    if(value < array[iBin+1]) return iBin;
  }
  
  // If the value is higher than the upper bin border, just use the highest index
  return nBins-1;
}


/*
 * Find centrality bin index
 */
int TrackPairEfficiencyCorrector::FindCentralityBin(const double centrality) const{
  return FindBinIndex(fCentralityBinBorders, fnCentralityBins, centrality);
}

/*
 * Find track pT bin index
 */
int TrackPairEfficiencyCorrector::FindTrackPtBin(const double trackPt) const{
  return FindBinIndex(fTrackPtBinBorders, fnTrackPtBins, trackPt);
}

/*
 * Find jet pT bin index
 */
int TrackPairEfficiencyCorrector::FindJetPtBin(const double jetPt) const{
  return FindBinIndex(fJetPtBinBorders, fnJetPtBins, jetPt);
}

/*
 * Get the weight to be given to the reflected cone track based on a data study
 *
 *  const double deltaR = DeltaR between the two tracks
 *  const double centrality = Centrality of the event
 *  const double triggerPt = Trigger particle pT
 *  const double associatedPt = Associated particle pT
 *
 *  return: std::pair<double,double> where the first number is the correction, and the second is the error of the track pair efficiency
 */
std::pair<double,double> TrackPairEfficiencyCorrector::GetTrackPairEfficiencyCorrection(const double deltaR, const double centrality, const double triggerPt, const double associatedPt, const double jetPt) const{
  
  if(fDisableCorrection) return std::make_pair(1,0);

  // We do not need to correct anything above 0.4, since that is outside of the analysis range
  if(deltaR > 0.4) return std::make_pair(1,0);

  // Find the higher and lower pT from trigger and associated particles
  double higherPt = triggerPt;
  double lowerPt = associatedPt;
  if(associatedPt > triggerPt){
    higherPt = associatedPt;
    lowerPt = triggerPt;
  }

  // Determine the correct bin from the centrality and track pT values
  int iCentrality = FindCentralityBin(centrality);
  int iTriggerPt = FindTrackPtBin(higherPt);
  int iAssociatedPt = FindTrackPtBin(lowerPt);
  double correction = 1;
  double error = 0;
  bool skipJetCorrection = false;

  // If jet pT is given, use the correction with jet pT bins
  if(jetPt > 0) {

    // For pp, in some of the highest bins there is not enough statistics close to the jets. Use the correction without jets in these cases
    if(!fIsPbPbData){
      if(higherPt > 50 && lowerPt > 20){
        skipJetCorrection = true;
      } else if (higherPt > 100){
        skipJetCorrection = true;
      }
    }

    // For PbPb, in some of the highest bins there is not enough statistics close to the jets. Use the correction without jets in these cases
    if(fIsPbPbData){
      if(higherPt > 50 && lowerPt > 12){
        skipJetCorrection = true;
      } else if (higherPt > 20 && lowerPt > 20){
        skipJetCorrection = true;
      }
    }

    // If we are in bins where the statistics is not a problem, do the correction in jet pT bins
    if(!skipJetCorrection){

      // For small deltaR, proceed with the correction
      int iJetPt = FindJetPtBin(jetPt);

      correction = 1.0/fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->GetBinContent(fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->FindBin(deltaR));
      error = fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->GetBinError(fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->FindBin(deltaR));

      // Check that the correction is sane
      if(1.0/correction > 0.01 && 1.0/correction < 2){
        return std::make_pair(correction,error);
      }

      // If the correction does not make sense, use the correction without jet pT bins instead
      std::cout << "Insane track pair efficiency correction in bin: triggerPt: " << iTriggerPt << " assoc pT: " << iAssociatedPt << " jetPt: " << iJetPt  << " deltaR: " << deltaR << std::endl;
      std::cout << "Using the correction without jets instead" << std::endl;

    } // If for not skipping the corrections
    
  } // If for correction with jet pT bins

  // For pp, in some of the highest bins we skip the correction since we only see fluctuations there
  if(!fIsPbPbData){
    if(higherPt > 100 && lowerPt > 50){
      if(deltaR > 0.02) return std::make_pair(1,0);
    }
  }

  // Also for PbPb, we do not apply the correction in the highest bins due to not having enough statistics
  if(fIsPbPbData){
    if(higherPt > 50){
      if(lowerPt > 50 && deltaR > 0.028) return std::make_pair(1,0);
      if(lowerPt > 20 && deltaR > 0.05) return std::make_pair(1,0);
      if(lowerPt > 16 && deltaR > 0.1) return std::make_pair(1,0);
    }
  }
    
  // Ensure that the trigger pT is higher than associated pT
  correction = 1.0/fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt]->GetBinContent(fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt]->FindBin(deltaR));
  error = fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt]->GetBinError(fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt]->FindBin(deltaR));
 
  return std::make_pair(correction, error);
}

// Getter for flag to disable the track pair efficiency correction
bool TrackPairEfficiencyCorrector::GetDisableCorrection(){
  return fDisableCorrection;
}

// Setter for disabling the correction
void TrackPairEfficiencyCorrector::SetDisableCorrection(const bool disable){
  fDisableCorrection = disable;
}
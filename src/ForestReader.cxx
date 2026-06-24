// Implementation for ForestReader

// Own includes
#include "ForestReader.h"

/*
 * Default constructor
 */
ForestReader::ForestReader() :
  fDataType(0),
  fUseJetTrigger(0),
  fJetType(0),
  fJetAxis(0),
  fMatchJets(false),
  fReadTrackTree(true),
  fIsMiniAOD(false),
  fMixingMode(false),
  fMegaSkimMode(false),
  fIsMC(false),
  fIsPPb(false),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fHFPlusBranch(0),
  fHFMinusBranch(0),
  fPtHatBranch(0),
  fEventNumberBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetWTAPhiBranch(0),
  fJetEtaBranch(0),
  fJetWTAEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fJetRefPtBranch(0),
  fJetRefEtaBranch(0),
  fJetRefPhiBranch(0),
  fJetRefFlavorBranch(0),
  fJetMatchedPtBranch(0),
  fJetMatchedEtaBranch(0),
  fJetMatchedPhiBranch(0),
  fEventWeightBranch(0),
  fCaloJet60FilterBranch(0),
  fCaloJet80FilterBranch(0),
  fCaloJet100FilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fHBHENoiseBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fPileupFilterBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fHFPlus(0),
  fHFMinus(0),
  fPtHat(0),
  fEventNumber(0),
  fnJets(0),
  fnMatchedJets(0),
  fEventWeight(1),
  fCaloJet15FilterBit(0),
  fCaloJet30FilterBit(0),
  fCaloJet40FilterBit(0),
  fCaloJet60FilterBit(0),
  fCaloJet80FilterBit(0),
  fCaloJet100FilterBit(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fHBHENoiseFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0),
  fPileupFilterBit(0),
  fnTracks(0)

{
  // Default constructor
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC, 4 = pPb p -> -eta, 5 = pPb p -> +eta, 6 = pPb 5 TeV, 7 = pPb MC p -> -eta, 8 = pPb MC p -> +eta, 9 = OO, 10 = OO MC
 *   Int_t useJetTrigger: 0 = Do not use any triggers, > 0 = Require jet trigger
 *   Int_t jetType: 0 = Calo jets, 1 = PF jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = WTA axis
 *   Int_t matchJets: non-0 = Do matching for reco and gen jets. 0 = Do not require matching
 *   Bool_t readTrackTree: Read the track tree from the forest. Optimizes speed if tracks are not needed
 *   Bool_t mixingMode: Read only trees needed for mixed events
 *   Bool_t megaSkimMode: Assume that the file contains only the information strictly necessary for event mixing
 */
ForestReader::ForestReader(Int_t dataType, Int_t useJetTrigger, Int_t jetType, Int_t jetAxis, Int_t matchJets, Bool_t readTrackTree, Bool_t mixingMode, Bool_t megaSkimMode) :
  fUseJetTrigger(useJetTrigger),
  fJetType(jetType),
  fJetAxis(jetAxis),
  fMatchJets(matchJets),
  fReadTrackTree(readTrackTree),
  fIsMiniAOD(false),
  fMixingMode(mixingMode),
  fMegaSkimMode(megaSkimMode),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fHFPlusBranch(0),
  fHFMinusBranch(0),
  fPtHatBranch(0),
  fEventNumberBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetWTAPhiBranch(0),
  fJetEtaBranch(0),
  fJetWTAEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fJetRefPtBranch(0),
  fJetRefEtaBranch(0),
  fJetRefPhiBranch(0),
  fJetRefFlavorBranch(0),
  fJetMatchedPtBranch(0),
  fJetMatchedEtaBranch(0),
  fJetMatchedPhiBranch(0),
  fEventWeightBranch(0),
  fCaloJet60FilterBranch(0),
  fCaloJet80FilterBranch(0),
  fCaloJet100FilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fHBHENoiseBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fPileupFilterBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fHFPlus(0),
  fHFMinus(0),
  fPtHat(0),
  fEventNumber(0),
  fnJets(0),
  fnMatchedJets(0),
  fEventWeight(1),
  fCaloJet15FilterBit(0),
  fCaloJet30FilterBit(0),
  fCaloJet40FilterBit(0),
  fCaloJet60FilterBit(0),
  fCaloJet80FilterBit(0),
  fCaloJet100FilterBit(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fHBHENoiseFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0),
  fPileupFilterBit(0),
  fnTracks(0)
{
  // Custom constructor
  
  SetDataType(dataType);

  // We do not use trigger in mixed event mode, so disable jet trigger if mixing mode is used
  if(fMixingMode) fUseJetTrigger = false;
  
}

/*
 * Copy constructor
 */
ForestReader::ForestReader(const ForestReader& in) :
  fDataType(in.fDataType),
  fUseJetTrigger(in.fUseJetTrigger),
  fJetType(in.fJetType),
  fJetAxis(in.fJetAxis),
  fMatchJets(in.fMatchJets),
  fReadTrackTree(in.fReadTrackTree),
  fIsMiniAOD(in.fIsMiniAOD),
  fMixingMode(in.fMixingMode),
  fMegaSkimMode(in.fMegaSkimMode),
  fIsMC(in.fIsMC),
  fIsPPb(in.fIsPPb),
  fHiVzBranch(in.fHiVzBranch),
  fHiBinBranch(in.fHiBinBranch),
  fHFPlusBranch(in.fHFPlusBranch),
  fHFMinusBranch(in.fHFMinusBranch),
  fPtHatBranch(in.fPtHatBranch),
  fEventNumberBranch(in.fEventNumberBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetWTAPhiBranch(in.fJetWTAPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fJetWTAEtaBranch(in.fJetWTAEtaBranch),
  fJetRawPtBranch(in.fJetRawPtBranch),
  fJetMaxTrackPtBranch(in.fJetMaxTrackPtBranch),
  fJetRefPtBranch(in.fJetRefPtBranch),
  fJetRefEtaBranch(in.fJetRefEtaBranch),
  fJetRefPhiBranch(in.fJetRefPhiBranch),
  fJetRefFlavorBranch(in.fJetRefFlavorBranch),
  fJetMatchedPtBranch(in.fJetMatchedPtBranch),
  fJetMatchedEtaBranch(in.fJetMatchedEtaBranch),
  fJetMatchedPhiBranch(in.fJetMatchedPhiBranch),
  fEventWeightBranch(in.fEventWeightBranch),
  fCaloJet60FilterBranch(in.fCaloJet60FilterBranch),
  fCaloJet80FilterBranch(in.fCaloJet80FilterBranch),
  fCaloJet100FilterBranch(in.fCaloJet100FilterBranch),
  fPrimaryVertexBranch(in.fPrimaryVertexBranch),
  fBeamScrapingBranch(in.fBeamScrapingBranch),
  fHBHENoiseBranch(in.fHBHENoiseBranch),
  fHfCoincidenceBranch(in.fHfCoincidenceBranch),
  fClusterCompatibilityBranch(in.fClusterCompatibilityBranch),
  fPileupFilterBranch(in.fPileupFilterBranch),
  fTrackPtBranch(in.fTrackPtBranch),
  fTrackPtErrorBranch(in.fTrackPtErrorBranch),
  fTrackPhiBranch(in.fTrackPhiBranch),
  fTrackEtaBranch(in.fTrackEtaBranch),
  fHighPurityTrackBranch(in.fHighPurityTrackBranch),
  fTrackVertexDistanceZBranch(in.fTrackVertexDistanceZBranch),
  fTrackVertexDistanceZErrorBranch(in.fTrackVertexDistanceZErrorBranch),
  fTrackVertexDistanceXYBranch(in.fTrackVertexDistanceXYBranch),
  fTrackVertexDistanceXYErrorBranch(in.fTrackVertexDistanceXYErrorBranch),
  fTrackChi2Branch(in.fTrackChi2Branch),
  fnTrackDegreesOfFreedomBranch(in.fnTrackDegreesOfFreedomBranch),
  fnHitsTrackerLayerBranch(in.fnHitsTrackerLayerBranch),
  fnHitsTrackBranch(in.fnHitsTrackBranch),
  fTrackEnergyEcalBranch(in.fTrackEnergyEcalBranch),
  fTrackEnergyHcalBranch(in.fTrackEnergyHcalBranch),
  fVertexZ(in.fVertexZ),
  fHiBin(in.fHiBin),
  fHFPlus(in.fHFPlus),
  fHFMinus(in.fHFMinus),
  fPtHat(in.fPtHat),
  fEventNumber(in.fEventNumber),
  fnJets(in.fnJets),
  fnMatchedJets(in.fnMatchedJets),
  fEventWeight(in.fEventWeight),
  fCaloJet15FilterBit(in.fCaloJet15FilterBit),
  fCaloJet30FilterBit(in.fCaloJet30FilterBit),
  fCaloJet40FilterBit(in.fCaloJet40FilterBit),
  fCaloJet60FilterBit(in.fCaloJet60FilterBit),
  fCaloJet80FilterBit(in.fCaloJet80FilterBit),
  fCaloJet100FilterBit(in.fCaloJet100FilterBit),
  fPrimaryVertexFilterBit(in.fPrimaryVertexFilterBit),
  fBeamScrapingFilterBit(in.fBeamScrapingFilterBit),
  fHBHENoiseFilterBit(in.fHBHENoiseFilterBit),
  fHfCoincidenceFilterBit(in.fHfCoincidenceFilterBit),
  fClusterCompatibilityFilterBit(in.fClusterCompatibilityFilterBit),
  fPileupFilterBit(in.fPileupFilterBit),
  fnTracks(in.fnTracks)
{
  // Copy constructor
}

/*
 * Assignment operator
 */
ForestReader& ForestReader::operator=(const ForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  fDataType = in.fDataType;
  fUseJetTrigger = in.fUseJetTrigger;
  fJetType = in.fJetType;
  fJetAxis = in.fJetAxis;
  fMatchJets = in.fMatchJets;
  fReadTrackTree = in.fReadTrackTree;
  fIsMiniAOD = in.fIsMiniAOD;
  fMixingMode = in.fMixingMode;
  fMegaSkimMode = in.fMegaSkimMode;
  fIsMC = in.fIsMC;
  fIsPPb = in.fIsPPb;
  fHiVzBranch = in.fHiVzBranch;
  fHiBinBranch = in.fHiBinBranch;
  fHFPlusBranch = in.fHFPlusBranch;
  fHFMinusBranch = in.fHFMinusBranch;
  fPtHatBranch = in.fPtHatBranch;
  fEventNumberBranch = in.fEventNumberBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetWTAPhiBranch = in.fJetWTAPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fJetWTAEtaBranch = in.fJetWTAEtaBranch;
  fJetRawPtBranch = in.fJetRawPtBranch;
  fJetMaxTrackPtBranch = in.fJetMaxTrackPtBranch;
  fJetRefPtBranch = in.fJetRefPtBranch;
  fJetRefEtaBranch = in.fJetRefEtaBranch;
  fJetRefPhiBranch = in.fJetRefPhiBranch;
  fJetRefFlavorBranch = in.fJetRefFlavorBranch;
  fJetMatchedPtBranch = in.fJetMatchedPtBranch;
  fJetMatchedEtaBranch = in.fJetMatchedEtaBranch;
  fJetMatchedPhiBranch = in.fJetMatchedPhiBranch;
  fEventWeightBranch = in.fEventWeightBranch;
  fCaloJet60FilterBranch = in.fCaloJet60FilterBranch;
  fCaloJet80FilterBranch = in.fCaloJet80FilterBranch;
  fCaloJet100FilterBranch = in.fCaloJet100FilterBranch;
  fPrimaryVertexBranch = in.fPrimaryVertexBranch;
  fBeamScrapingBranch = in.fBeamScrapingBranch;
  fHBHENoiseBranch = in.fHBHENoiseBranch;
  fHfCoincidenceBranch = in.fHfCoincidenceBranch;
  fClusterCompatibilityBranch = in.fClusterCompatibilityBranch;
  fPileupFilterBranch = in.fPileupFilterBranch;
  fTrackPtBranch = in.fTrackPtBranch;
  fTrackPtErrorBranch = in.fTrackPtErrorBranch;
  fTrackPhiBranch = in.fTrackPhiBranch;
  fTrackEtaBranch = in.fTrackEtaBranch;
  fHighPurityTrackBranch = in.fHighPurityTrackBranch;
  fTrackVertexDistanceZBranch = in.fTrackVertexDistanceZBranch;
  fTrackVertexDistanceZErrorBranch = in.fTrackVertexDistanceZErrorBranch;
  fTrackVertexDistanceXYBranch = in.fTrackVertexDistanceXYBranch;
  fTrackVertexDistanceXYErrorBranch = in.fTrackVertexDistanceXYErrorBranch;
  fTrackChi2Branch = in.fTrackChi2Branch;
  fnTrackDegreesOfFreedomBranch = in.fnTrackDegreesOfFreedomBranch;
  fnHitsTrackerLayerBranch = in.fnHitsTrackerLayerBranch;
  fnHitsTrackBranch = in.fnHitsTrackBranch;
  fTrackEnergyEcalBranch = in.fTrackEnergyEcalBranch;
  fTrackEnergyHcalBranch = in.fTrackEnergyHcalBranch;
  fVertexZ = in.fVertexZ;
  fHiBin = in.fHiBin;
  fHFPlus = in.fHFPlus;
  fHFMinus = in.fHFMinus;
  fPtHat = in.fPtHat;
  fEventNumber = in.fEventNumber;
  fnJets = in.fnJets;
  fnMatchedJets = in.fnMatchedJets;
  fEventWeight = in.fEventWeight;
  fCaloJet15FilterBit = in.fCaloJet15FilterBit;
  fCaloJet30FilterBit = in.fCaloJet30FilterBit;
  fCaloJet40FilterBit = in.fCaloJet40FilterBit;
  fCaloJet60FilterBit = in.fCaloJet60FilterBit;
  fCaloJet80FilterBit = in.fCaloJet80FilterBit;
  fCaloJet100FilterBit = in.fCaloJet100FilterBit;
  fPrimaryVertexFilterBit = in.fPrimaryVertexFilterBit;
  fBeamScrapingFilterBit = in.fBeamScrapingFilterBit;
  fHBHENoiseFilterBit = in.fHBHENoiseFilterBit;
  fHfCoincidenceFilterBit = in.fHfCoincidenceFilterBit;
  fClusterCompatibilityFilterBit = in.fClusterCompatibilityFilterBit;
  fPileupFilterBit = in.fPileupFilterBit;
  fnTracks = in.fnTracks;
  
  return *this;
}

/*
 * Destructor
 */
ForestReader::~ForestReader(){
  // destructor
}

/*
 * Setter for fDataType
 */
void ForestReader::SetDataType(Int_t dataType){
  
  //Sanity check for given data type
  if(dataType < 0 || dataType > knDataTypes-1){
    cout << "ERROR: Data type input " << dataType << " is invalid in ForestReader.cxx!" << endl;
    cout << "Please give integer between 0 and " << knDataTypes-1 << "." << endl;
    cout << "Setting data type to 0 (pp)." << endl;
    fDataType = 0;
  } else {
    
    // If the sanity check passes, set the given data type
    fDataType = dataType;

    // Also set the relevant flags for this data type
    fIsMC = (fDataType == kPpMC) || (fDataType == kPbPbMC) || (fDataType == kPPbMC_pToMinusEta) || (fDataType == kPPbMC_pToPlusEta);
    fIsPPb = (fDataType == kPPb_pToMinusEta) || (fDataType == kPPb_pToPlusEta) || (fDataType == kPPb_pToMinusEta_5TeV) || (fDataType == kPPbMC_pToMinusEta) || (fDataType == kPPbMC_pToPlusEta);
  }
}

// Getter for number of events in the tree
Int_t ForestReader::GetNEvents() const{
  return fHiVzBranch->GetEntries();
}

// Getter for number of jets in an event
Int_t ForestReader::GetNJets() const{
  return fnJets;
}

// Getter for vertex z position
Float_t ForestReader::GetVz() const{
  return fVertexZ;
}

// Getter for centrality. CMS has integer centrality bins from 0 to 200, thus division by 2.
Float_t ForestReader::GetCentrality() const{
  return fHiBin/2.0;
}

// Getter for hiBin. Return 1 for negative values (for easier handling of tracking efficiency correction)
Int_t ForestReader::GetHiBin() const{
  if(fHiBin < 0) return 1;
  return fHiBin;
}

 // Getter for energy in HF plus calorimeters
Float_t ForestReader::GetHFPlus() const{
  return fHFPlus;
}

// Getter for energy in HF minus calorimeters   
Float_t ForestReader::GetHFMinus() const{
  return fHFMinus;
}

// Getter for energy in HF calorimeters        
Float_t ForestReader::GetHFSum() const{
  return fHFPlus + fHFMinus;
}        

// Getter for pT hat
Float_t ForestReader::GetPtHat() const{
  return fPtHat;
}

// Getter for event number
ULong64_t ForestReader::GetEventNumber() const{
  return fEventNumber;
}

// Getter for pT hat
Float_t ForestReader::GetEventWeight() const{
  return fEventWeight;
}

// Getter for calorimeter jet filter bit with threshold 15 GeV.
Int_t ForestReader::GetCaloJet15FilterBit() const{
  return fCaloJet15FilterBit;
}

// Getter for calorimeter jet filter bit with threshold 30 GeV.
Int_t ForestReader::GetCaloJet30FilterBit() const{
  return fCaloJet30FilterBit;
}

// Getter for calorimeter jet filter bit with threshold 40 GeV.
Int_t ForestReader::GetCaloJet40FilterBit() const{
  return fCaloJet40FilterBit;
}

// Getter for calorimeter jet filter bit with threshold 60 GeV. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetCaloJet60FilterBit() const{
  return fCaloJet60FilterBit;
}

// Getter for calorimeter jet filter bit with threshold 80 GeV. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetCaloJet80FilterBit() const{
  return fCaloJet80FilterBit;
}

// Getter for calorimeter jet filter bit with threshold 100 GeV. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetCaloJet100FilterBit() const{
  return fCaloJet100FilterBit;
}

// Getter for primary vertex filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetPrimaryVertexFilterBit() const{
  return fPrimaryVertexFilterBit;
}

// Getter for beam scraping filter bit. Always 1 for MC and PbPb (set in the initializer).
Int_t ForestReader::GetBeamScrapingFilterBit() const{
  return fBeamScrapingFilterBit;
}

// Getter for HB/HE noise filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetHBHENoiseFilterBit() const{
  return fHBHENoiseFilterBit;
}

// Getter for HF energy coincidence filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetHfCoincidenceFilterBit() const{
  return fHfCoincidenceFilterBit;
}

// Getter for cluster compatibility filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetClusterCompatibilityFilterBit() const{
  return fClusterCompatibilityFilterBit;
}

// Getter for cluster pile-up filter bit.
Int_t ForestReader::GetPileupFilterBit() const{
  return fPileupFilterBit;
}

// Getter for number of tracks in an event
Int_t ForestReader::GetNTracks() const{
  return fnTracks;
}

// Getter for track algorithm
Int_t ForestReader::GetTrackAlgorithm(Int_t iTrack) const{
  return 0;
}

// Getter for the original track algorithm
Int_t ForestReader::GetTrackOriginalAlgorithm(Int_t iTrack) const{
  return 0;
}

// Getter for track MVA
Float_t ForestReader::GetTrackMVA(Int_t iTrack) const{
  return 0;
}

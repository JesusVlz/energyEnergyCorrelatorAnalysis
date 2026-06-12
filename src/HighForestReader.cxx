// Implementation for HighForestReader

//
// Modified: Jesus Corral, 09-May-2026
// ===========================================================
// Added support for OO data and MC
// including the new jet trigger, branches and filters for OO and the centrality definition for OO. 
// 
//===========================================================
// Own includes
#include "HighForestReader.h"

/*
 * Default constructor
 */
HighForestReader::HighForestReader() :
  ForestReader(),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fParticleTree(0),
  fnJetsBranch(0),
  fnTracksBranch(0),
  fMatchedJetWTAEtaBranch(0),
  fMatchedJetWTAPhiBranch(0),
  fTrackAlgorithmBranch(0),
  fTrackOriginalAlgorithmBranch(0),
  fTrackMVABranch(0),
  fParticlePtBranch(0),
  fParticleEtaBranch(0),
  fParticlePhiBranch(0),
  fParticleChargeBranch(0),
  fParticleSubeventBranch(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fJetRawPtArray(),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fHighPurityTrackArray(),
  fTrackVertexDistanceZArray(),
  fTrackVertexDistanceZErrorArray(),
  fTrackVertexDistanceXYArray(),
  fTrackVertexDistanceXYErrorArray(),
  fTrackChi2Array(),
  fnTrackDegreesOfFreedomArray(),
  fnHitsTrackerLayerArray(),
  fnHitsTrackArray(),
  fTrackEnergyEcalArray(),
  fTrackEnergyHcalArray(),
  fTrackAlgorithmArray(),
  fTrackOriginalAlgorithmArray(),
  fTrackMVAArray(),
  fTrackPtVector(0),
  fTrackPtErrorVector(0),
  fTrackPhiVector(0),
  fTrackEtaVector(0),
  fHighPurityTrackVector(0),
  fTrackVertexDistanceZVector(0),
  fTrackVertexDistanceZErrorVector(0),
  fTrackVertexDistanceXYVector(0),
  fTrackVertexDistanceXYErrorVector(0),
  fTrackNormalizedChi2Vector(0),
  fnHitsTrackerLayerVector(0),
  fnHitsTrackVector(0),
  fTrackEnergyEcalVector(0),
  fTrackEnergyHcalVector(0),
  fParticlePtVector(0),
  fParticlePhiVector(0),
  fParticleEtaVector(0),
  fParticleChargeVector(0),
  fParticleSubeventVector(0)
{
  // Default constructor
  
  // Initialize fJetMaxTrackPtArray to -1
  for(int i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC, 4 = pPb p -> -eta, 5 = pPb p -> +eta 6 = pPb 5 TeV, 7 = pPb MC p -> -eta, 8 = pPb MC p -> +eta
 *   Int_t useJetTrigger: 0 = Do not use any triggers, > 0 = Require jet triggers
 *   Int_t jetType: 0 = Calo jets, 1 = PF jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = WTA axis
 *   Bool_t matchJets: non-0 = Do matching for reco and gen jets. 0 = Do not require matching
 *   Bool_t readTrackTree: Read the track trees from the forest. Optimizes speed if tracks are not needed
 *   Bool_t mixingMode: Flag for mixed events more (false = regular mode, true = mixed event mode)
 *   Bool_t megaSkimMode: Assume that the file contains only the information strictly necessary for event mixing
 */
HighForestReader::HighForestReader(Int_t dataType, Int_t useJetTrigger, Int_t jetType, Int_t jetAxis, Int_t matchJets, Bool_t readTrackTree, Bool_t mixingMode, Bool_t megaSkimMode) :
  ForestReader(dataType,useJetTrigger,jetType,jetAxis,matchJets,readTrackTree,mixingMode,megaSkimMode),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fParticleTree(0),
  fnJetsBranch(0),
  fnTracksBranch(0),
  fMatchedJetWTAEtaBranch(0),
  fMatchedJetWTAPhiBranch(0),
  fTrackAlgorithmBranch(0),
  fTrackOriginalAlgorithmBranch(0),
  fTrackMVABranch(0),
  fParticlePtBranch(0),
  fParticleEtaBranch(0),
  fParticlePhiBranch(0),
  fParticleChargeBranch(0),
  fParticleSubeventBranch(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fHighPurityTrackArray(),
  fTrackVertexDistanceZArray(),
  fTrackVertexDistanceZErrorArray(),
  fTrackVertexDistanceXYArray(),
  fTrackVertexDistanceXYErrorArray(),
  fTrackChi2Array(),
  fnTrackDegreesOfFreedomArray(),
  fnHitsTrackerLayerArray(),
  fnHitsTrackArray(),
  fTrackEnergyEcalArray(),
  fTrackEnergyHcalArray(),
  fTrackAlgorithmArray(),
  fTrackOriginalAlgorithmArray(),
  fTrackMVAArray(),
  fTrackPtVector(0),
  fTrackPtErrorVector(0),
  fTrackPhiVector(0),
  fTrackEtaVector(0),
  fHighPurityTrackVector(0),
  fTrackVertexDistanceZVector(0),
  fTrackVertexDistanceZErrorVector(0),
  fTrackVertexDistanceXYVector(0),
  fTrackVertexDistanceXYErrorVector(0),
  fTrackNormalizedChi2Vector(0),
  fnHitsTrackerLayerVector(0),
  fnHitsTrackVector(0),
  fTrackEnergyEcalVector(0),
  fTrackEnergyHcalVector(0),
  fParticlePtVector(0),
  fParticlePhiVector(0),
  fParticleEtaVector(0),
  fParticleChargeVector(0),
  fParticleSubeventVector(0)
{
  // Custom constructor
  
  // Initialize fJetMaxTrackPtArray to -1
  for(int i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
  
}

/*
 * Copy constructor
 */
HighForestReader::HighForestReader(const HighForestReader& in) :
  ForestReader(in),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fHltTree(in.fHltTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fParticleTree(in.fParticleTree),
  fnJetsBranch(in.fnJetsBranch),
  fnTracksBranch(in.fnTracksBranch),
  fMatchedJetWTAEtaBranch(in.fMatchedJetWTAEtaBranch),
  fMatchedJetWTAPhiBranch(in.fMatchedJetWTAPhiBranch),
  fTrackAlgorithmBranch(in.fTrackAlgorithmBranch),
  fTrackOriginalAlgorithmBranch(in.fTrackOriginalAlgorithmBranch),
  fTrackMVABranch(in.fTrackMVABranch),
  fParticlePtBranch(in.fParticlePtBranch),
  fParticleEtaBranch(in.fParticleEtaBranch),
  fParticlePhiBranch(in.fParticlePhiBranch),
  fParticleChargeBranch(in.fParticleChargeBranch),
  fParticleSubeventBranch(in.fParticleSubeventBranch),
  fTrackPtVector(in.fTrackPtVector),
  fTrackPtErrorVector(in.fTrackPtErrorVector),
  fTrackPhiVector(in.fTrackPhiVector),
  fTrackEtaVector(in.fTrackEtaVector),
  fHighPurityTrackVector(in.fHighPurityTrackVector),
  fTrackVertexDistanceZVector(in.fTrackVertexDistanceZVector),
  fTrackVertexDistanceZErrorVector(in.fTrackVertexDistanceZErrorVector),
  fTrackVertexDistanceXYVector(in.fTrackVertexDistanceXYVector),
  fTrackVertexDistanceXYErrorVector(in.fTrackVertexDistanceXYErrorVector),
  fTrackNormalizedChi2Vector(in.fTrackNormalizedChi2Vector),
  fnHitsTrackerLayerVector(in.fnHitsTrackerLayerVector),
  fnHitsTrackVector(in.fnHitsTrackVector),
  fTrackEnergyEcalVector(in.fTrackEnergyEcalVector),
  fTrackEnergyHcalVector(in.fTrackEnergyHcalVector),
  fParticlePtVector(in.fParticlePtVector),
  fParticlePhiVector(in.fParticlePhiVector),
  fParticleEtaVector(in.fParticleEtaVector),
  fParticleChargeVector(in.fParticleChargeVector),
  fParticleSubeventVector(in.fParticleSubeventVector)
{
  // Copy constructor
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
  }
  
  for(Int_t i = 0; i < fnMaxTrack; i++){
    fTrackPtArray[i] = in.fTrackPtArray[i];
    fTrackPtErrorArray[i] = in.fTrackPtErrorArray[i];
    fTrackPhiArray[i] = in.fTrackPhiArray[i];
    fTrackEtaArray[i] = in.fTrackEtaArray[i];
    fHighPurityTrackArray[i] = in.fHighPurityTrackArray[i];
    fTrackVertexDistanceZArray[i] = in.fTrackVertexDistanceZArray[i];
    fTrackVertexDistanceZErrorArray[i] = in.fTrackVertexDistanceZErrorArray[i];
    fTrackVertexDistanceXYArray[i] = in.fTrackVertexDistanceXYArray[i];
    fTrackVertexDistanceXYErrorArray[i] = in.fTrackVertexDistanceXYErrorArray[i];
    fTrackChi2Array[i] = in.fTrackChi2Array[i];
    fnTrackDegreesOfFreedomArray[i] = in.fnTrackDegreesOfFreedomArray[i];
    fnHitsTrackerLayerArray[i] = in.fnHitsTrackerLayerArray[i];
    fnHitsTrackArray[i] = in.fnHitsTrackArray[i];
    fTrackEnergyEcalArray[i] = in.fTrackEnergyEcalArray[i];
    fTrackEnergyHcalArray[i] = in.fTrackEnergyHcalArray[i];
    fTrackAlgorithmArray[i] = in.fTrackAlgorithmArray[i];
    fTrackOriginalAlgorithmArray[i] = in.fTrackOriginalAlgorithmArray[i];
    fTrackMVAArray[i] = in.fTrackMVAArray[i];
  }
}

/*
 * Assignment operator
 */
HighForestReader& HighForestReader::operator=(const HighForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  ForestReader::operator=(in);
  
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fHltTree = in.fHltTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  fParticleTree = in.fParticleTree;
  fnJetsBranch = in.fnJetsBranch;
  fnTracksBranch = in.fnTracksBranch;
  fMatchedJetWTAEtaBranch = in.fMatchedJetWTAEtaBranch;
  fMatchedJetWTAPhiBranch = in.fMatchedJetWTAPhiBranch;
  fTrackAlgorithmBranch = in.fTrackAlgorithmBranch;
  fTrackOriginalAlgorithmBranch = in.fTrackOriginalAlgorithmBranch;
  fTrackMVABranch = in.fTrackMVABranch;
  fParticlePtBranch = in.fParticlePtBranch;
  fParticleEtaBranch = in.fParticleEtaBranch;
  fParticlePhiBranch = in.fParticlePhiBranch;
  fParticleChargeBranch = in.fParticleChargeBranch;
  fParticleSubeventBranch = in.fParticleSubeventBranch;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
  }
  
  for(Int_t i = 0; i < fnMaxTrack; i++){
    fTrackPtArray[i] = in.fTrackPtArray[i];
    fTrackPtErrorArray[i] = in.fTrackPtErrorArray[i];
    fTrackPhiArray[i] = in.fTrackPhiArray[i];
    fTrackEtaArray[i] = in.fTrackEtaArray[i];
    fHighPurityTrackArray[i] = in.fHighPurityTrackArray[i];
    fTrackVertexDistanceZArray[i] = in.fTrackVertexDistanceZArray[i];
    fTrackVertexDistanceZErrorArray[i] = in.fTrackVertexDistanceZErrorArray[i];
    fTrackVertexDistanceXYArray[i] = in.fTrackVertexDistanceXYArray[i];
    fTrackVertexDistanceXYErrorArray[i] = in.fTrackVertexDistanceXYErrorArray[i];
    fTrackChi2Array[i] = in.fTrackChi2Array[i];
    fnTrackDegreesOfFreedomArray[i] = in.fnTrackDegreesOfFreedomArray[i];
    fnHitsTrackerLayerArray[i] = in.fnHitsTrackerLayerArray[i];
    fnHitsTrackArray[i] = in.fnHitsTrackArray[i];
    fTrackEnergyEcalArray[i] = in.fTrackEnergyEcalArray[i];
    fTrackEnergyHcalArray[i] = in.fTrackEnergyHcalArray[i];
    fTrackAlgorithmArray[i] = in.fTrackAlgorithmArray[i];
    fTrackOriginalAlgorithmArray[i] = in.fTrackOriginalAlgorithmArray[i];
    fTrackMVAArray[i] = in.fTrackMVAArray[i];
  }
  
  fTrackPtVector = in.fTrackPtVector;
  fTrackPtErrorVector = in.fTrackPtErrorVector;
  fTrackPhiVector = in.fTrackPhiVector;
  fTrackEtaVector = in.fTrackEtaVector;
  fHighPurityTrackVector = in.fHighPurityTrackVector;
  fTrackVertexDistanceZVector = in.fTrackVertexDistanceZVector;
  fTrackVertexDistanceZErrorVector = in.fTrackVertexDistanceZErrorVector;
  fTrackVertexDistanceXYVector = in.fTrackVertexDistanceXYVector;
  fTrackVertexDistanceXYErrorVector = in.fTrackVertexDistanceXYErrorVector;
  fTrackNormalizedChi2Vector = in.fTrackNormalizedChi2Vector;
  fnHitsTrackerLayerVector = in.fnHitsTrackerLayerVector;
  fnHitsTrackVector = in.fnHitsTrackVector;
  fTrackEnergyEcalVector = in.fTrackEnergyEcalVector;
  fTrackEnergyHcalVector = in.fTrackEnergyHcalVector;

  fParticlePtVector = in.fParticlePtVector;
  fParticlePhiVector = in.fParticlePhiVector;
  fParticleEtaVector = in.fParticleEtaVector;
  fParticleChargeVector = in.fParticleChargeVector;
  fParticleSubeventVector = in.fParticleSubeventVector;
  
  return *this;
}

/*
 * Destructor
 */
HighForestReader::~HighForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void HighForestReader::Initialize(){

  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*", 0);
  fHeavyIonTree->SetBranchStatus("vz", 1);
  fHeavyIonTree->SetBranchAddress("vz", &fVertexZ, &fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("evt", 1);
  fHeavyIonTree->SetBranchAddress("evt", &fEventNumber, &fEventNumberBranch);

  if(fDataType == kPp || fDataType == kPpMC){
    // We do not have HF tower information for pp. In this case find HF like energy from particle flow candidates
    fHeavyIonTree->SetBranchStatus("hiHFPlus_pf", 1);
    fHeavyIonTree->SetBranchAddress("hiHFPlus_pf", &fHFPlus, &fHFPlusBranch);
    fHeavyIonTree->SetBranchStatus("hiHFMinus_pf", 1);
    fHeavyIonTree->SetBranchAddress("hiHFMinus_pf", &fHFMinus, &fHFMinusBranch);
  } else {
    fHeavyIonTree->SetBranchStatus("hiHFplus", 1);
    fHeavyIonTree->SetBranchAddress("hiHFplus", &fHFPlus, &fHFPlusBranch);
    fHeavyIonTree->SetBranchStatus("hiHFminus", 1);
    fHeavyIonTree->SetBranchAddress("hiHFminus", &fHFMinus, &fHFMinusBranch);
  }

  if(fDataType == kPbPb || fDataType == kPbPbMC || fDataType == kOO || fDataType == kOOMC){
    fHeavyIonTree->SetBranchStatus("hiBin", 1);
    fHeavyIonTree->SetBranchAddress("hiBin", &fHiBin, &fHiBinBranch);
  } else {
    fHiBin = -1;  // No centrality definition for pp or pPb
  }
  if(fIsMC){
    fHeavyIonTree->SetBranchStatus("pthat", 1);
    fHeavyIonTree->SetBranchAddress("pthat", &fPtHat, &fPtHatBranch); // pT hat only for MC
    fHeavyIonTree->SetBranchStatus("weight", 1);
    fHeavyIonTree->SetBranchAddress("weight", &fEventWeight, &fEventWeightBranch);
  } else {
    fPtHat = 0;        // We do not have pT hat information for real data
    fEventWeight = 1;  // No event weighting for real data
  }
  
  // Connect the branches to the jet tree
  const char* jetAxis[2] = {"jt","WTA"};
  const char* branchName;

  // Jet trees are not needed for mixed events
  if(fMixingMode){

    fnJets = 0;
    fnMatchedJets = 0;

  } else {
  
    fJetTree->SetBranchStatus("*", 0);

    // The corrected jet pT does not exist in the pPb forest
    if(!fIsPPb){
      fJetTree->SetBranchStatus("jtpt", 1);
      fJetTree->SetBranchAddress("jtpt", &fJetPtArray, &fJetPtBranch);
    }
  
    // If specified, select WTA axis for jet phi
    branchName = Form("%sphi",jetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName, 1);
    fJetTree->SetBranchAddress(branchName, &fJetPhiArray, &fJetPhiBranch);
  
    // If specified, select WTA axis for jet eta
    branchName = Form("%seta",jetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName, 1);
    fJetTree->SetBranchAddress(branchName, &fJetEtaArray, &fJetEtaBranch);
  
    fJetTree->SetBranchStatus("nref", 1);
    fJetTree->SetBranchAddress("nref", &fnJets, &fnJetsBranch);

    // OO forests use different branch names for raw pT and maximum track pT inside a jet
    // In OO forests, the raw pT is stored in the "jtptUncorrected" branch, while in other forests it is stored in the "rawpt" branch. 
    // Similarly, the maximum track pT inside a jet is stored in the "chargedMax" branch in OO forests and in the "trackMax" branch in other forests.
    if(fDataType == kOO || fDataType == kOOMC){
      fJetTree->SetBranchStatus("jtptUncorrected", 1);
      fJetTree->SetBranchAddress("jtptUncorrected", &fJetRawPtArray, &fJetRawPtBranch);
      fJetTree->SetBranchStatus("chargedMax", 1);
      fJetTree->SetBranchAddress("chargedMax", &fJetMaxTrackPtArray, &fJetMaxTrackPtBranch);
    } else {
      fJetTree->SetBranchStatus("rawpt", 1);
      fJetTree->SetBranchAddress("rawpt", &fJetRawPtArray, &fJetRawPtBranch);
      fJetTree->SetBranchStatus("trackMax", 1);
      fJetTree->SetBranchAddress("trackMax", &fJetMaxTrackPtArray, &fJetMaxTrackPtBranch);
    }
    
  
    // If we are looking at Monte Carlo, connect the reference pT and parton arrays
    if(fIsMC){
      fJetTree->SetBranchStatus("refpt", 1);
      fJetTree->SetBranchAddress("refpt", &fJetRefPtArray, &fJetRefPtBranch);
      fJetTree->SetBranchStatus("refeta", 1);
      fJetTree->SetBranchAddress("refeta", &fJetRefEtaArray, &fJetRefEtaBranch);
      fJetTree->SetBranchStatus("refphi", 1);
      fJetTree->SetBranchAddress("refphi", &fJetRefPhiArray, &fJetRefPhiBranch);
      fJetTree->SetBranchStatus("refparton_flavorForB", 1);
      fJetTree->SetBranchAddress("refparton_flavorForB", &fJetRefFlavorArray, &fJetRefFlavorBranch);
      fJetTree->SetBranchStatus("genpt", 1);
      fJetTree->SetBranchAddress("genpt", &fMatchedJetPtArray, &fJetMatchedPtBranch);
    
      // Read jet phi and eta
      fJetTree->SetBranchStatus("genphi", 1);
      fJetTree->SetBranchAddress("genphi", &fMatchedJetPhiArray, &fJetMatchedPhiBranch);
      fJetTree->SetBranchStatus("geneta", 1);
      fJetTree->SetBranchAddress("geneta", &fMatchedJetEtaArray, &fJetMatchedEtaBranch);
    
      // WTA axis is not available in pPb files
      if(!fIsPPb){
        fJetTree->SetBranchStatus("WTAgenphi", 1);
        fJetTree->SetBranchAddress("WTAgenphi", &fMatchedJetWTAPhiArray, &fMatchedJetWTAPhiBranch);
        fJetTree->SetBranchStatus("WTAgeneta", 1);
        fJetTree->SetBranchAddress("WTAgeneta", &fMatchedJetWTAEtaArray, &fMatchedJetWTAEtaBranch);
      }
    
      fJetTree->SetBranchStatus("ngen", 1);
      fJetTree->SetBranchAddress("ngen", &fnMatchedJets, &fnMatchedJetsBranch);
    }
  }
  
  // Event selection summary
  //
  //         tree                      branch                         What it is
  //  hltanalysis/HltTree   HLT_HIPuAK4CaloJet100_Eta5p1_v1      Event selection for PbPb
  //  hltanalysis/HltTree   HLT_OxyL1SingleJet20_v1              Event selection for OO
  //  hltanalysis/HltTree      HLT_AK4CaloJet80_Eta5p1_v1         Event selection for pp
  // skimanalysis/HltTree         pprimaryVertexFilter           Event selection for PbPb
  // skimanalysis/HltTree    HBHENoiseFilterResultRun2Loose   Event selection for pp and PbPb
  // skimanalysis/HltTree         pPAprimaryVertexFilter          Event selection for pp
  // skimanalysis/HltTree           pBeamScrapingFilter           Event selection for pp
  
  // Connect the branches to the HLT tree
  if(fUseJetTrigger){
    
    fHltTree->SetBranchStatus("*", 0);
    
    if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC

      // Calo jet 15 trigger
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet15_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet15_v1", &fCaloJet15FilterBit, &fCaloJet15FilterBranch);

      // Calo jet 30 trigger
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet30_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet30_v1", &fCaloJet30FilterBit, &fCaloJet30FilterBranch);

      // Calo jet 40 trigger
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet40_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet40_v1", &fCaloJet40FilterBit, &fCaloJet40FilterBranch);

      // Calo jet 60 trigger
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet60_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet60_v1", &fCaloJet60FilterBit, &fCaloJet60FilterBranch);

      // Calo jet 80 trigger
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet80_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet80_v1", &fCaloJet80FilterBit, &fCaloJet80FilterBranch);
      
      // Calo jet 100 trigger
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet100_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet100_v1", &fCaloJet100FilterBit, &fCaloJet100FilterBranch);
      
    } else if(fIsPPb) { // pPb data or MC

      // No low jet pT triggers in pPb forests
      fCaloJet15FilterBit = 1;
      fCaloJet30FilterBit = 1;
      fCaloJet40FilterBit = 1;

      // Calo jet 60 trigger
      fHltTree->SetBranchStatus("HLT_PAAK4CaloJet60_Eta5p1_v3", 1);
      fHltTree->SetBranchAddress("HLT_PAAK4CaloJet60_Eta5p1_v3", &fCaloJet60FilterBit, &fCaloJet60FilterBranch);

      // Calo jet 80 trigger
      fHltTree->SetBranchStatus("HLT_PAAK4CaloJet80_Eta5p1_v3", 1);
      fHltTree->SetBranchAddress("HLT_PAAK4CaloJet80_Eta5p1_v3", &fCaloJet80FilterBit, &fCaloJet80FilterBranch);
      
      // Calo jet 100 trigger
      fHltTree->SetBranchStatus("HLT_PAAK4CaloJet100_Eta5p1_v3", 1);
      fHltTree->SetBranchAddress("HLT_PAAK4CaloJet100_Eta5p1_v3", &fCaloJet100FilterBit, &fCaloJet100FilterBranch);

    } else if(fDataType == kOO || fDataType == kOOMC) { //  OO data or MC

      // Calo jet 20 trigger
      fHltTree->SetBranchStatus("HLT_OxyL1SingleJet20_v1", 1);
      fHltTree->SetBranchAddress("HLT_OxyL1SingleJet20_v1", &fCaloJet15FilterBit, &fCaloJet15FilterBranch);

      // Calo jet 35 trigger
      fHltTree->SetBranchStatus("HLT_OxyL1SingleJet35_v1", 1);
      fHltTree->SetBranchAddress("HLT_OxyL1SingleJet35_v1", &fCaloJet30FilterBit, &fCaloJet30FilterBranch);

      // Calo jet 44 trigger
      fHltTree->SetBranchStatus("HLT_OxyL1SingleJet44_v1", 1);
      fHltTree->SetBranchAddress("HLT_OxyL1SingleJet44_v1", &fCaloJet40FilterBit, &fCaloJet40FilterBranch);

       // Calo jet 60 trigger
      fHltTree->SetBranchStatus("HLT_OxyL1SingleJet60_v1", 1);
      fHltTree->SetBranchAddress("HLT_OxyL1SingleJet60_v1", &fCaloJet60FilterBit, &fCaloJet60FilterBranch);
      
      // No high jet pT triggers in OO forests
      fCaloJet80FilterBit = 1;
      fCaloJet100FilterBit = 1;

    } else { // PbPb data or MC

      // No low jet pT triggers in PbPb forests
      fCaloJet15FilterBit = 1;
      fCaloJet30FilterBit = 1;

      // Calo jet 60 trigger
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet40Eta5p1_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet40Eta5p1_v1", &fCaloJet40FilterBit, &fCaloJet40FilterBranch);

      // Calo jet 60 trigger
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet60Eta5p1_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet60Eta5p1_v1", &fCaloJet60FilterBit, &fCaloJet60FilterBranch);

      // Calo jet 80 trigger
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet80Eta5p1_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80Eta5p1_v1", &fCaloJet80FilterBit, &fCaloJet80FilterBranch);
      
      // Calo jet 100 trigger
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1", &fCaloJet100FilterBit, &fCaloJet100FilterBranch);
      
    }
  } else {
    fCaloJet60FilterBit = 1;
    fCaloJet80FilterBit = 1;
    fCaloJet100FilterBit = 1;
  }
  
  // Connect the branches to the skim tree (different for pp and PbPb data and Monte Carlo)
  if(fMegaSkimMode){

    // In mega skim mode all the event selections are already done. Set all flags to 1.
    fPrimaryVertexFilterBit = 1;
    fBeamScrapingFilterBit = 1;
    fHBHENoiseFilterBit = 1;
    fHfCoincidenceFilterBit = 1;
    fClusterCompatibilityFilterBit = 1;

  } else {

    fSkimTree->SetBranchStatus("*", 0);
    // pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter
    if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC
      
      // HBHE noise of beam scraping filter bits are not available in the MiniAOD forests.
      // Also the name of the vertex filter is changed for MiniAOD forests
      if(fIsMiniAOD){
        fHBHENoiseFilterBit = 1;
        fBeamScrapingFilterBit = 1;
        fSkimTree->SetBranchStatus("pprimaryVertexFilter", 1);
        fSkimTree->SetBranchAddress("pprimaryVertexFilter", &fPrimaryVertexFilterBit, &fPrimaryVertexBranch);
      } else {
        fSkimTree->SetBranchStatus("pPAprimaryVertexFilter", 1);
        fSkimTree->SetBranchAddress("pPAprimaryVertexFilter", &fPrimaryVertexFilterBit, &fPrimaryVertexBranch);
        fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
        fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &fHBHENoiseFilterBit, &fHBHENoiseBranch);
        fSkimTree->SetBranchStatus("pBeamScrapingFilter", 1);
        fSkimTree->SetBranchAddress("pBeamScrapingFilter", &fBeamScrapingFilterBit, &fBeamScrapingBranch);
      }
      fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
      fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
      fPileupFilterBit = 1; // No pileup filter for pp

    } else if(fIsPPb){ // pPb data or MC

      fSkimTree->SetBranchStatus("pPAprimaryVertexFilter", 1);
      fSkimTree->SetBranchAddress("pPAprimaryVertexFilter", &fPrimaryVertexFilterBit, &fPrimaryVertexBranch);
      fSkimTree->SetBranchStatus("pBeamScrapingFilter", 1);
      fSkimTree->SetBranchAddress("pBeamScrapingFilter", &fBeamScrapingFilterBit, &fBeamScrapingBranch);
      fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
      fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &fHBHENoiseFilterBit, &fHBHENoiseBranch);
      fSkimTree->SetBranchStatus("phfCoincFilter", 1);
      fSkimTree->SetBranchAddress("phfCoincFilter", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
      fSkimTree->SetBranchStatus("pVertexFilterCutdz1p0", 1);
      fSkimTree->SetBranchAddress("pVertexFilterCutdz1p0", &fPileupFilterBit, &fPileupFilterBranch);

      fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pPb

     } else if(fDataType == kOO || fDataType == kOOMC){ // OO data or MC

      // Primary vertex has at least two tracks, is within 25 cm in z-direction and within 2 cm in xy-direction
      fSkimTree->SetBranchStatus("pprimaryVertexFilter", 1);
      fSkimTree->SetBranchAddress("pprimaryVertexFilter", &fPrimaryVertexFilterBit, &fPrimaryVertexBranch);
     
      // Have at least two towers on both of the HF calorimerter to have energies above 4 GeV
      // accumulated by the energies of PF (particle-flow)candidates
      fSkimTree->SetBranchStatus("pphfCoincFilterPF2Th4", 1);
      fSkimTree->SetBranchAddress("pphfCoincFilterPF2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);

      fHBHENoiseFilterBit = 1; // HBHE noise filter bit is not available in the OO MiniAOD forests.
      fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for OO
      fBeamScrapingFilterBit = 1;  // No beam scraping filter for OO
      fPileupFilterBit = 1;        // No pile-up filter for OO

    } else { // PbPb data or MC
    
      // Primary vertex has at least two tracks, is within 25 cm in z-direction and within 2 cm in xy-direction
      fSkimTree->SetBranchStatus("pprimaryVertexFilter", 1);
      fSkimTree->SetBranchAddress("pprimaryVertexFilter", &fPrimaryVertexFilterBit, &fPrimaryVertexBranch);
    
      // Cut on noise on HCAL
      if(fIsMiniAOD){
        fHBHENoiseFilterBit = 1; // HBHE noise filter bit is not available in the MiniAOD forests.
      
        // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
        fSkimTree->SetBranchStatus("pphfCoincFilter2Th4", 1);
        fSkimTree->SetBranchAddress("pphfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
      
      } else {
        //fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
        //fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &fHBHENoiseFilterBit, &fHBHENoiseBranch);
        fHBHENoiseFilterBit = 1; // For testing purpose, set this to 1
      
        // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
        fSkimTree->SetBranchStatus("phfCoincFilter2Th4", 1);
        fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
      }
    
      // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible
      fSkimTree->SetBranchStatus("pclusterCompatibilityFilter", 1);
      fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fClusterCompatibilityFilterBit, &fClusterCompatibilityBranch);
    
      fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
      fPileupFilterBit = 1;        // No pile-up filter for PbPb
    }
  }
  
  // Connect the branches to the track tree
  if(fReadTrackTree){
    
    fTrackTree->SetBranchStatus("*", 0);
    
    // We need to read the forest to vectors for MiniAODs and to arrays for AODs
    if(fIsMiniAOD){
      
      fTrackTree->SetBranchStatus("trkPt", 1);
      fTrackTree->SetBranchAddress("trkPt", &fTrackPtVector, &fTrackPtBranch);
      fTrackTree->SetBranchStatus("trkPhi", 1);
      fTrackTree->SetBranchAddress("trkPhi", &fTrackPhiVector, &fTrackPhiBranch);
      fTrackTree->SetBranchStatus("trkEta", 1);
      fTrackTree->SetBranchAddress("trkEta", &fTrackEtaVector, &fTrackEtaBranch);
      fTrackTree->SetBranchStatus("nTrk", 1);
      fTrackTree->SetBranchAddress("nTrk", &fnTracks, &fnTracksBranch);

      // In mega skim mode, only basic track kinematics are available
      if(!fMegaSkimMode){
        fTrackTree->SetBranchStatus("trkPtError", 1);
        fTrackTree->SetBranchAddress("trkPtError", &fTrackPtErrorVector, &fTrackPtErrorBranch);
        fTrackTree->SetBranchStatus("highPurity", 1);
        fTrackTree->SetBranchAddress("highPurity", &fHighPurityTrackVector, &fHighPurityTrackBranch);
        fTrackTree->SetBranchStatus("trkDzFirstVtx", 1);
        fTrackTree->SetBranchAddress("trkDzFirstVtx", &fTrackVertexDistanceZVector, &fTrackVertexDistanceZBranch);
        fTrackTree->SetBranchStatus("trkDzErrFirstVtx", 1);
        fTrackTree->SetBranchAddress("trkDzErrFirstVtx", &fTrackVertexDistanceZErrorVector, &fTrackVertexDistanceZErrorBranch);
        fTrackTree->SetBranchStatus("trkDxyFirstVtx", 1);
        fTrackTree->SetBranchAddress("trkDxyFirstVtx", &fTrackVertexDistanceXYVector, &fTrackVertexDistanceXYBranch);
        fTrackTree->SetBranchStatus("trkDxyErrFirstVtx", 1);
        fTrackTree->SetBranchAddress("trkDxyErrFirstVtx", &fTrackVertexDistanceXYErrorVector, &fTrackVertexDistanceXYErrorBranch);
        fTrackTree->SetBranchStatus("trkNormChi2", 1);
        fTrackTree->SetBranchAddress("trkNormChi2", &fTrackNormalizedChi2Vector, &fTrackChi2Branch);
        fTrackTree->SetBranchStatus("trkNLayers", 1);
        fTrackTree->SetBranchAddress("trkNLayers", &fnHitsTrackerLayerVector, &fnHitsTrackerLayerBranch);
        fTrackTree->SetBranchStatus("trkNHits", 1);
        fTrackTree->SetBranchAddress("trkNHits", &fnHitsTrackVector, &fnHitsTrackBranch);
        fTrackTree->SetBranchStatus("pfEcal", 1);
        fTrackTree->SetBranchAddress("pfEcal", &fTrackEnergyEcalVector, &fTrackEnergyEcalBranch);
        fTrackTree->SetBranchStatus("pfHcal", 1);
        fTrackTree->SetBranchAddress("pfHcal", &fTrackEnergyHcalVector, &fTrackEnergyHcalBranch);
      }
      
    } else { // Read the tree from AOD files
      
      fTrackTree->SetBranchStatus("trkPt", 1);
      fTrackTree->SetBranchAddress("trkPt", &fTrackPtArray, &fTrackPtBranch);
      fTrackTree->SetBranchStatus("trkPhi", 1);
      fTrackTree->SetBranchAddress("trkPhi", &fTrackPhiArray, &fTrackPhiBranch);
      fTrackTree->SetBranchStatus("trkEta", 1);
      fTrackTree->SetBranchAddress("trkEta", &fTrackEtaArray, &fTrackEtaBranch);
      fTrackTree->SetBranchStatus("nTrk", 1);
      fTrackTree->SetBranchAddress("nTrk", &fnTracks, &fnTracksBranch);

      // In mega skim mode, only basic track kinematics are available
      if(!fMegaSkimMode){
        fTrackTree->SetBranchStatus("trkPtError", 1);
        fTrackTree->SetBranchAddress("trkPtError", &fTrackPtErrorArray, &fTrackPtErrorBranch);
        fTrackTree->SetBranchStatus("highPurity", 1);
        fTrackTree->SetBranchAddress("highPurity", &fHighPurityTrackArray, &fHighPurityTrackBranch);
        fTrackTree->SetBranchStatus("trkDz1", 1);
        fTrackTree->SetBranchAddress("trkDz1", &fTrackVertexDistanceZArray, &fTrackVertexDistanceZBranch);
        fTrackTree->SetBranchStatus("trkDzError1", 1);
        fTrackTree->SetBranchAddress("trkDzError1", &fTrackVertexDistanceZErrorArray, &fTrackVertexDistanceZErrorBranch);
        fTrackTree->SetBranchStatus("trkDxy1", 1);
        fTrackTree->SetBranchAddress("trkDxy1", &fTrackVertexDistanceXYArray, &fTrackVertexDistanceXYBranch);
        fTrackTree->SetBranchStatus("trkDxyError1", 1);
        fTrackTree->SetBranchAddress("trkDxyError1", &fTrackVertexDistanceXYErrorArray, &fTrackVertexDistanceXYErrorBranch);
        fTrackTree->SetBranchStatus("pfEcal", 1);
        fTrackTree->SetBranchAddress("pfEcal", &fTrackEnergyEcalArray, &fTrackEnergyEcalBranch);
        fTrackTree->SetBranchStatus("pfHcal", 1);
        fTrackTree->SetBranchAddress("pfHcal", &fTrackEnergyHcalArray, &fTrackEnergyHcalBranch);

        // Branches that do not exist for the pPb forests
        if(!fIsPPb){
          fTrackTree->SetBranchStatus("trkChi2", 1);
          fTrackTree->SetBranchAddress("trkChi2", &fTrackChi2Array, &fTrackChi2Branch);
          fTrackTree->SetBranchStatus("trkNdof", 1);
          fTrackTree->SetBranchAddress("trkNdof", &fnTrackDegreesOfFreedomArray, &fnTrackDegreesOfFreedomBranch);
          fTrackTree->SetBranchStatus("trkNlayer", 1);
          fTrackTree->SetBranchAddress("trkNlayer", &fnHitsTrackerLayerArray, &fnHitsTrackerLayerBranch);
          fTrackTree->SetBranchStatus("trkNHit", 1);
          fTrackTree->SetBranchAddress("trkNHit", &fnHitsTrackArray, &fnHitsTrackBranch);
      
          // Additional information for track cuts
          fTrackTree->SetBranchStatus("trkAlgo", 1);
          fTrackTree->SetBranchAddress("trkAlgo", &fTrackAlgorithmArray, &fTrackAlgorithmBranch);
          fTrackTree->SetBranchStatus("trkOriginalAlgo", 1);
          fTrackTree->SetBranchAddress("trkOriginalAlgo", &fTrackOriginalAlgorithmArray, &fTrackOriginalAlgorithmBranch);
        }
      
        // Track MVA only in PbPb trees
        if(fDataType == kPbPb || fDataType == kPbPbMC){
          fTrackTree->SetBranchStatus("trkMVA", 1);
          fTrackTree->SetBranchAddress("trkMVA", &fTrackMVAArray, &fTrackMVABranch);
        }
      }
    }
  } // Reading track trees

  // Generator level particles needed only is we do mixing, and we are working on Monte Carlo simulations
  if(fIsMC){

    fParticleTree->SetBranchStatus("*", 0);
    fParticleTree->SetBranchStatus("pt", 1);
    fParticleTree->SetBranchAddress("pt", &fParticlePtVector, &fParticlePtBranch);
    fParticleTree->SetBranchStatus("phi", 1);
    fParticleTree->SetBranchAddress("phi", &fParticlePhiVector, &fParticlePhiBranch);
    fParticleTree->SetBranchStatus("eta", 1);
    fParticleTree->SetBranchAddress("eta", &fParticleEtaVector, &fParticleEtaBranch);
    fParticleTree->SetBranchStatus("sube", 1);
    fParticleTree->SetBranchAddress("sube", &fParticleSubeventVector, &fParticleSubeventBranch);

    // The track charge is not available in mega skims
    if(!fMegaSkimMode){
      fParticleTree->SetBranchStatus("chg", 1);
      fParticleTree->SetBranchAddress("chg", &fParticleChargeVector, &fParticleChargeBranch);
    }

  }

  // We need to load one event to initialize TChains properly. Not sure why, but this is how things seem to work
  GetEvent(0);
  
}

/*
 * Connect a new tree to the reader
 */
void HighForestReader::ReadForestFromFile(TFile* inputFile){
  std::vector<TString> fileList;
  fileList.push_back(inputFile->GetName());
  ReadForestFromFileList(fileList);
}

/*
 * Connect a new tree to the reader
 */
void HighForestReader::ReadForestFromFileList(std::vector<TString> fileList){

  // Open one file to determine if we are dealing with miniAOD or regular AOD
  TFile* inputFile = TFile::Open(fileList.at(0));
  TTree* miniAODcheck;
  // The track tree has different name in miniAOD and AOD
  if(fDataType == kPbPb || fDataType == kPbPbMC || fDataType == kOO || fDataType == kOOMC){
    miniAODcheck = (TTree*)inputFile->Get("PbPbTracks/trackTree");
  } else {
    miniAODcheck = (TTree*)inputFile->Get("ppTracks/trackTree");
  }
  fIsMiniAOD = !(miniAODcheck == NULL);
  inputFile->Close();

  // Helper variable for finding the correct tree
  const char *treeName[4] = {"none","none","none","none"};

  // Connect a trees from the file to the reader
  fHeavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
  if(fUseJetTrigger) fHltTree = new TChain("hltanalysis/HltTree");
  if(!fMegaSkimMode) fSkimTree = new TChain("skimanalysis/HltTree"); // Mega skims have all event selection built-in

  // The jet tree has different name for PbPb data and other datasets
  if(fDataType == kPbPb || fDataType == kPbPbMC){
    treeName[0] = "akPu4CaloJetAnalyzer/t";     // Tree for calo jets
    treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for csPF jets
    treeName[2] = "akPu4PFJetAnalyzer/t";       // Tree for puPF jets
    treeName[3] = "akFlowPuCs4PFJetAnalyzer/t"; // Tree for flow subtracted csPF jets
  } else if(fDataType == kOO || fDataType == kOOMC){
    treeName[0] = "akCs2PFJetAnalyzer/t";       // Tree for cs R=0.2 PF jets
    treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for cs R=0.4 PF jets
    treeName[2] = "akCs6PFJetAnalyzer/t";       // Tree for cs R=0.6 PF jets
    treeName[3] = "akCs8PFJetAnalyzer/t";       // Tree for cs R=0.8 PF jets
  } else { // pp or pPb data or MC
    treeName[0] = "ak4CaloJetAnalyzer/t"; // Tree for calo jets
    treeName[1] = "ak4PFJetAnalyzer/t";   // Tree for PF jets
    treeName[2] = "akCs4PFJetAnalyzer/t"; // Tree for constituent subtracted PF jets
  }

  if(!fMixingMode) fJetTree = new TChain(treeName[fJetType]);

  if(fReadTrackTree){
    if(fIsMiniAOD && (fDataType == kPbPb || fDataType == kPbPbMC || fDataType == kOO || fDataType == kOOMC)){
      fTrackTree = new TChain("PbPbTracks/trackTree");
    } else {
      if(fIsMiniAOD){
        fTrackTree = new TChain("ppTracks/trackTree");
      } else {
        fTrackTree = new TChain("ppTrack/trackTree");
      }
    }
  }

  if(fIsMC){
    fParticleTree = new TChain("HiGenParticleAna/hi");
  }
  
  for(std::vector<TString>::iterator listIterator = fileList.begin(); listIterator != fileList.end(); listIterator++){
    fHeavyIonTree->Add(*listIterator);
    if(!fMegaSkimMode) fSkimTree->Add(*listIterator);
    if(fUseJetTrigger) fHltTree->Add(*listIterator);
    if(!fMixingMode) fJetTree->Add(*listIterator);
    if(fReadTrackTree) fTrackTree->Add(*listIterator);
    if(fIsMC){
      fParticleTree->Add(*listIterator);
    }
  }
  
  // Connect branches to trees
  Initialize();
}

/*
 * Burn the current forest.
 */
void HighForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fHltTree->Delete();
  fSkimTree->Delete();
  fJetTree->Delete();
  fTrackTree->Delete();
  fParticleTree->Delete();
}

/*
 * Check if there are problems in the file list
 */
bool HighForestReader::CheckFileProblems(){

  // Loop over all the files that are in TChains and check they are doing fine
  TObjArray *fileElements = fHeavyIonTree->GetListOfFiles();
  for (TObject *op: *fileElements) {
    auto chainElement = static_cast<TChainElement *>(op);
    TFile *testFile = TFile::Open(chainElement->GetTitle());

    // Check that the file is still open
    if(!testFile->IsOpen()) {
      std::cout << "The file " << testFile->GetName() << " is not open anymore!" << std::endl;
      return true;
    }

    // Check that the file did not turn into a zombie
    if(testFile->IsZombie()){
      std::cout << "The file " << testFile->GetName() << " turned into a zombie!" << std::endl;
      return true;
    } 

    // If everything is fine with the file, we can close the test file
    testFile->Close();
  }

  // If no problems are found from any of the files, return false
  return false;

}

/*
 * Load an event to memory
 */
void HighForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  if(!fMixingMode) fJetTree->GetEntry(nEvent);
  if(fUseJetTrigger) fHltTree->GetEntry(nEvent);
  if(!fMegaSkimMode) fSkimTree->GetEntry(nEvent);
  if(fReadTrackTree) {
    fTrackTree->GetEntry(nEvent);
    
    // For miniAOD, the number of tracks can be deduced from the vectors
    if(fIsMiniAOD){
      fnTracks = fTrackPtVector->size();
    }
  }
  if(fIsMC){
    fParticleTree->GetEntry(nEvent);
  }
}

// Getter for jet pT
Float_t HighForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t HighForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray[iJet];
}

// Getter for jet eta
Float_t HighForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray[iJet];
}

// Getter for jet raw pT
Float_t HighForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetRawPtArray[iJet];
}

// Getter for maximum track pT inside a jet
Float_t HighForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return fJetMaxTrackPtArray[iJet];
}

// Getter for track pT
Float_t HighForestReader::GetTrackPt(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackPtVector->at(iTrack);
  return fTrackPtArray[iTrack];
}

// Getter for track pT error
Float_t HighForestReader::GetTrackPtError(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackPtErrorVector->at(iTrack);
  return fTrackPtErrorArray[iTrack];
}

// Getter for track phi
Float_t HighForestReader::GetTrackPhi(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackPhiVector->at(iTrack);
  return fTrackPhiArray[iTrack];
}

// Getter for track eta
Float_t HighForestReader::GetTrackEta(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackEtaVector->at(iTrack);
  return fTrackEtaArray[iTrack];
}

// Getter for high purity of the track
Bool_t HighForestReader::GetTrackHighPurity(Int_t iTrack) const{
  if(fIsMiniAOD) return fHighPurityTrackVector->at(iTrack);
  return fHighPurityTrackArray[iTrack];
}

// Getter for track distance from primary vertex in z-direction
Float_t HighForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceZVector->at(iTrack);
  return fTrackVertexDistanceZArray[iTrack];
}

// Getter for error of track distance from primary vertex in z-direction
Float_t HighForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceZErrorVector->at(iTrack);
  return fTrackVertexDistanceZErrorArray[iTrack];
}

// Getter for track distance from primary vertex in xy-direction
Float_t HighForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceXYVector->at(iTrack);
  return fTrackVertexDistanceXYArray[iTrack];
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t HighForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceXYErrorVector->at(iTrack);
  return fTrackVertexDistanceXYErrorArray[iTrack];
}

// Getter for normalized track chi2 value from reconstruction fit
Float_t HighForestReader::GetTrackNormalizedChi2(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackNormalizedChi2Vector->at(iTrack);
  return GetTrackChi2(iTrack) / (1.0*GetNTrackDegreesOfFreedom(iTrack));
}

// Getter for track chi2 value from reconstruction fit
Float_t HighForestReader::GetTrackChi2(Int_t iTrack) const{
  if(fIsMiniAOD) return -1; // Does not exist in MiniAOD forest
  return fTrackChi2Array[iTrack];
}

// Getter for number of degrees of freedom in reconstruction fit
Int_t HighForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  if(fIsMiniAOD) return -1; // Does not exist in MiniAOD forest
  return fnTrackDegreesOfFreedomArray[iTrack];
}

// Getter for number of hits in tracker layers
Int_t HighForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  if(fIsMiniAOD) return fnHitsTrackerLayerVector->at(iTrack);
  return fnHitsTrackerLayerArray[iTrack];
}

// Getter for number of hits for the track
Int_t HighForestReader::GetNHitsTrack(Int_t iTrack) const{
  if(fIsMiniAOD) return fnHitsTrackVector->at(iTrack);
  return fnHitsTrackArray[iTrack];
}

// Getter for track energy in ECal
Float_t HighForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackEnergyEcalVector->at(iTrack);
  return fTrackEnergyEcalArray[iTrack];
}

// Getter for track energy in HCal
Float_t HighForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackEnergyHcalVector->at(iTrack);
  return fTrackEnergyHcalArray[iTrack];
}

// Getter for track charge. Charge is only relevant for generator level tracks.
Int_t HighForestReader::GetTrackCharge(Int_t iTrack) const{
  return 1;
}

// Getter for track subevent index. Relevant only for generator level tracks.
Int_t HighForestReader::GetTrackSubevent(Int_t iTrack) const{
  return -1;
}

// Getter for track MC status. Relevant only for generator level tracks.
Int_t HighForestReader::GetTrackMCStatus(Int_t iTrack) const{
  return 1;
}

// Getter for track algorithm
Int_t HighForestReader::GetTrackAlgorithm(Int_t iTrack) const{
  if(fIsMiniAOD) return 0; // Does not exist in MiniAOD forest
  return fTrackAlgorithmArray[iTrack];
}

// Getter for track original algorithm
Int_t HighForestReader::GetTrackOriginalAlgorithm(Int_t iTrack) const{
  if(fIsMiniAOD) return 0; // Does not exist in MiniAOD forest
  return fTrackOriginalAlgorithmArray[iTrack];
}

// Getter for track MVA
Float_t HighForestReader::GetTrackMVA(Int_t iTrack) const{
  if(fIsMiniAOD) return 1; // Does not exist in MiniAOD forest
  return fTrackMVAArray[iTrack];
}

// Getter for the number of generator level particles
Int_t HighForestReader::GetNParticles() const{
  return fParticlePtVector->size();
}

// Getter for generator level particle pT
Float_t HighForestReader::GetParticlePt(Int_t iParticle) const{
  return fParticlePtVector->at(iParticle);
} 

// Getter for generator level particle phi
Float_t HighForestReader::GetParticlePhi(Int_t iParticle) const{
  return fParticlePhiVector->at(iParticle);
} 

// Getter for generator level particle eta
Float_t HighForestReader::GetParticleEta(Int_t iParticle) const{
  return fParticleEtaVector->at(iParticle);
} 

// Getter for generator level particle charge       
Int_t HighForestReader::GetParticleCharge(Int_t iParticle) const{
  return fParticleChargeVector->at(iParticle);
}

// Getter for generator level particle subevent index
Int_t HighForestReader::GetParticleSubevent(Int_t iParticle) const{
  return fParticleSubeventVector->at(iParticle);
}

// Check if generator level jet has a matching reconstructed jet
Bool_t HighForestReader::HasMatchingJet(Int_t iJet) const{
  
  // If we are not matching jets or are considering real data, everything is fine
  if(!fMatchJets || !fIsMC) return true;
  
  // For each reconstructed jet there is a reference pT, which tells the the pT of a matched generator level jet
  // If this number is -999, it means that there are no generator level jets matching the reconstructed jet
  if(fJetRefPtArray[iJet] < 0) return false;
  return true;
}

// Getter for matched generator level jet pT
Float_t HighForestReader::GetMatchedPt(Int_t iJet) const{
  
  // If we are not matching jets or are considering real data, this value has no meaning
  if(!fMatchJets || !fIsMC) return 0;
  
  // Find the index of the matching generator level jet
  Int_t matchedIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;

  // Return matched gen pT
  return fMatchedJetPtArray[matchedIndex];
}

// Parton flavor for the parton initiating the jet
Int_t HighForestReader::GetPartonFlavor(Int_t iJet) const{
  
  // If we are considering real data, this value has no meaning
  if(!fIsMC) return 0;
  
  // Return the matching parton flavor
  return fJetRefFlavorArray[iJet];
}

// Get the matching generator level jet index for the given reconstructed jet
Int_t HighForestReader::GetMatchingIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, it has a matching reconstructed jet
  Int_t matchedIndex = -1;
  for(Int_t iMatched = 0; iMatched < fnMatchedJets; iMatched++){
    if(TMath::Abs(fJetRefEtaArray[iJet] - fMatchedJetEtaArray[iMatched]) < 0.015){
      if(TMath::Abs(fJetRefPhiArray[iJet] - fMatchedJetPhiArray[iMatched]) < 0.015){
        if(TMath::Abs(fJetRefPtArray[iJet] - fMatchedJetPtArray[iMatched]) < 0.03*fMatchedJetPtArray[iMatched]){
          matchedIndex = iMatched;
          break;
        }
      }
    }
  }
  
  return matchedIndex;
  
}

// Get the eta of the matched generator level jet
Float_t HighForestReader::GetMatchedEta(Int_t iJet) const{
  
  // If not matching jets, just return something because this has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching generator level jet
  Int_t matchedIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  if(fJetAxis == 0) return fMatchedJetEtaArray[matchedIndex];
  return fMatchedJetWTAEtaArray[matchedIndex];
  
}

// Get the pT of the matched reconstructed jet
Float_t HighForestReader::GetMatchedPhi(Int_t iJet) const{
  
  // If not matching jets, just return something because this has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching generator level jet
  Int_t matchedIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  if(fJetAxis == 0) return fMatchedJetPhiArray[matchedIndex];
  return fMatchedJetWTAPhiArray[matchedIndex];
  
}

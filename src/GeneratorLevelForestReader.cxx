// Implementation for GeneratorLevelForestReader

// Own includes
#include "GeneratorLevelForestReader.h"

/*
 * Default constructor
 */
GeneratorLevelForestReader::GeneratorLevelForestReader() :
  ForestReader(),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fTrackPtArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fTrackChargeArray(0),
  fTrackSubeventArray(0),
  fTrackStatusArray(0)
{
  // Default constructor
}

/*
 * Custom constructor
 *
 *  Arguments:
 // 0 = pp Ref 5 TeV,
 // 1 = PbPb 5 TeV, 
 // 2 = pp MC, 
 // 3 = PbPb MC, 
 // 4 = pPb p -> -eta, 
 // 5 = pPb p -> +eta, 
 // 6 = pPb 5 TeV, 
 // 7 = pPb MC p -> -eta, 
 // 8 = pPb MC p -> +eta, 
 // 9 = OO, 
 // 10 = OO MC, 
 // 11 = pp Ref 5.36 TeV, 
 // 12 = pp Ref 5.36 TeV MC
 // 
 //  *   Int_t useJetTrigger: 0 = Do not use any triggers, > 0 = Require jet trigger
 *   Int_t jetType: 0 = Calo jets, 1 = PF jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = WTA axis
 *   Int_t matchJets: non-0 = Do matching for reco and gen jets. 0 = Do not require matching
 *   Bool_t readTrackTree: Read the track trees from the forest. Optimizes speed if tracks are not needed
 *   Bool_t mixingMode: Flag for mixed events more (false = regular mode, true = mixed event mode)
 *   Bool_t megaSkimMode: Assume that the file contains only the information strictly necessary for event mixing
 */
GeneratorLevelForestReader::GeneratorLevelForestReader(Int_t dataType, Int_t useJetTrigger, Int_t jetType, Int_t jetAxis, Int_t matchJets, Bool_t readTrackTree, Bool_t mixingMode, Bool_t megaSkimMode) :
  ForestReader(dataType,useJetTrigger,jetType,jetAxis,matchJets,readTrackTree,mixingMode,megaSkimMode),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fJetRefFlavorArray(),
  fRecoJetPtArray(),
  fRecoJetEtaArray(),
  fRecoJetPhiArray(),
  fRecoJetRawPtArray(),
  fTrackPtArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fTrackChargeArray(0),
  fTrackSubeventArray(0),
  fTrackStatusArray(0)
{
  // Custom constructor
  
}

/*
 * Copy constructor
 */
GeneratorLevelForestReader::GeneratorLevelForestReader(const GeneratorLevelForestReader& in) :
  ForestReader(in),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fHltTree(in.fHltTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fTrackPtArray(in.fTrackPtArray),
  fTrackPhiArray(in.fTrackPhiArray),
  fTrackEtaArray(in.fTrackEtaArray),
  fTrackChargeArray(in.fTrackChargeArray),
  fTrackSubeventArray(in.fTrackSubeventArray),
  fTrackStatusArray(in.fTrackStatusArray)
{
  // Copy constructor
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRefFlavorArray[i] = in.fJetRefFlavorArray[i];
    fRecoJetPtArray[i] = in.fRecoJetPtArray[i];
    fRecoJetEtaArray[i] = in.fRecoJetEtaArray[i];
    fRecoJetPhiArray[i] = in.fRecoJetPhiArray[i];
    fRecoJetRawPtArray[i] = in.fRecoJetRawPtArray[i];
  }
}

/*
 * Assignment operator
 */
GeneratorLevelForestReader& GeneratorLevelForestReader::operator=(const GeneratorLevelForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  ForestReader::operator=(in);
  
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fHltTree = in.fHltTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  fTrackPtArray = in.fTrackPtArray;
  fTrackPhiArray = in.fTrackPhiArray;
  fTrackEtaArray = in.fTrackEtaArray;
  fTrackChargeArray = in.fTrackChargeArray;
  fTrackSubeventArray = in.fTrackSubeventArray;
  fTrackStatusArray = in.fTrackStatusArray;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRefFlavorArray[i] = in.fJetRefFlavorArray[i];
    fRecoJetPtArray[i] = in.fRecoJetPtArray[i];
    fRecoJetEtaArray[i] = in.fRecoJetEtaArray[i];
    fRecoJetPhiArray[i] = in.fRecoJetPhiArray[i];
    fRecoJetRawPtArray[i] = in.fRecoJetRawPtArray[i];
  }
  
  return *this;
}

/*
 * Destructor
 */
GeneratorLevelForestReader::~GeneratorLevelForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void GeneratorLevelForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("evt",1);
  fHeavyIonTree->SetBranchAddress("evt", &fEventNumber, &fEventNumberBranch);

  if(fDataType == kPpMC || fDataType == kOOMC || fDataType == kPpRef5p36TeVMC){
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

  if(fDataType == kPbPbMC){
    fHeavyIonTree->SetBranchStatus("hiBin", 1);
    fHeavyIonTree->SetBranchAddress("hiBin", &fHiBin, &fHiBinBranch);
  } else {
    fHiBin = -1;  // No centrality definition for pp or pPb
  }

  fHeavyIonTree->SetBranchStatus("pthat",1);
  fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  
  // Read event weight for simulation
  fHeavyIonTree->SetBranchStatus("weight",1);
  fHeavyIonTree->SetBranchAddress("weight",&fEventWeight,&fEventWeightBranch);
  
  // Connect the branches to the jet tree
  const char* jetAxis[2] = {"jt","WTA"};
  const char* branchName;

  // Jet trees are not needed for mixed events
  if(fMixingMode){

    fnJets = 0;
    fnMatchedJets = 0;

  } else {
  
    fJetTree->SetBranchStatus("*",0);
    fJetTree->SetBranchStatus("genpt",1);
    fJetTree->SetBranchAddress("genpt",&fJetPtArray,&fJetPtBranch);
  
    // Read jet eta and phi
    fJetTree->SetBranchStatus("genphi",1);
    fJetTree->SetBranchAddress("genphi",&fJetPhiArray,&fJetPhiBranch);
    fJetTree->SetBranchStatus("geneta",1);
    fJetTree->SetBranchAddress("geneta",&fJetEtaArray,&fJetEtaBranch);

    // WTA axis is not available in pPb files
    if(!fIsPPb){
      fJetTree->SetBranchStatus("WTAgenphi",1);
      fJetTree->SetBranchAddress("WTAgenphi",&fJetWTAPhiArray,&fJetWTAPhiBranch);
      fJetTree->SetBranchStatus("WTAgeneta",1);
      fJetTree->SetBranchAddress("WTAgeneta",&fJetWTAEtaArray,&fJetWTAEtaBranch);
    }
  
    fJetTree->SetBranchStatus("ngen",1);
    fJetTree->SetBranchAddress("ngen",&fnJets,&fJetRawPtBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  
    // If we are matching jets, connect the ref pT array
    if(fMatchJets){
      fJetTree->SetBranchStatus("refpt",1);
      fJetTree->SetBranchAddress("refpt",&fJetRefPtArray,&fJetRefPtBranch);
      fJetTree->SetBranchStatus("refeta",1);
      fJetTree->SetBranchAddress("refeta",&fJetRefEtaArray,&fJetRefEtaBranch);
      fJetTree->SetBranchStatus("refphi",1);
      fJetTree->SetBranchAddress("refphi",&fJetRefPhiArray,&fJetRefPhiBranch);
      fJetTree->SetBranchStatus("refparton_flavorForB",1);
      fJetTree->SetBranchAddress("refparton_flavorForB",&fJetRefFlavorArray,&fJetRefFlavorBranch);
      fJetTree->SetBranchStatus("jtpt",1);
      fJetTree->SetBranchAddress("jtpt",&fRecoJetPtArray,&fJetMatchedPtBranch);
      fJetTree->SetBranchStatus("rawpt",1);
      fJetTree->SetBranchAddress("rawpt",&fRecoJetRawPtArray,&fHighPurityTrackBranch); // Reuse a branch from ForestReader that is not otherwise needed here
    
      // If specified, select WTA axis for jet phi
      branchName = Form("%sphi",jetAxis[fJetAxis]);
      fJetTree->SetBranchStatus(branchName,1);
      fJetTree->SetBranchAddress(branchName,&fRecoJetPhiArray,&fJetMatchedPhiBranch);
    
      // If specified, select WTA axis for jet eta
      branchName = Form("%seta",jetAxis[fJetAxis]);
      fJetTree->SetBranchStatus(branchName,1);
      fJetTree->SetBranchAddress(branchName,&fRecoJetEtaArray,&fJetMatchedEtaBranch);
    
      fJetTree->SetBranchStatus("nref",1);
      fJetTree->SetBranchAddress("nref",&fnMatchedJets,&fnTrackDegreesOfFreedomBranch); // Reuse a branch from ForestReader that is not otherwise needed here
    }
  }
  
  // Connect the branches to the HLT tree
  if(fUseJetTrigger){
    
    fHltTree->SetBranchStatus("*",0);
    
    if(fDataType == kPpMC){ // pp MC

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
      
    } else if(fIsPPb) { // pPb MC

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

    } else if(fDataType == kOOMC) { //  OO data or MC

      // Run-3 OO minimum-bias trigger, required for both data and MC.
      //
      // The first legacy trigger slot is reused as the generic Run-3
      // minimum-bias trigger bit. TriggerSelection 8 reads this slot.
      fHltTree->SetBranchStatus("HLT_MinimumBiasHF_OR_BptxAND_v1", 1);
      fHltTree->SetBranchAddress("HLT_MinimumBiasHF_OR_BptxAND_v1", &fCaloJet15FilterBit, &fCaloJet15FilterBranch);

      // The Run-3 minimum-bias analysis does not use jet-trigger paths.
      fCaloJet30FilterBit = 0;
      fCaloJet40FilterBit = 0;
      fCaloJet60FilterBit = 0;
      fCaloJet80FilterBit = 0;
      fCaloJet100FilterBit = 0;

    } else if(fDataType == kPpRef5p36TeVMC){

      // Run-3 pp-reference ZeroBias trigger, required for both data and MC.
      //
      // The first legacy trigger slot is reused as the generic Run-3
      // minimum-bias/zero-bias trigger bit. TriggerSelection 8 reads it.
      fHltTree->SetBranchStatus("HLT_PPRefZeroBias_v6", 1);
      fHltTree->SetBranchAddress("HLT_PPRefZeroBias_v6", &fCaloJet15FilterBit, &fCaloJet15FilterBranch);

      // The Run-3 ZeroBias analysis does not use jet-trigger paths.
      fCaloJet30FilterBit = 0;
      fCaloJet40FilterBit = 0;
      fCaloJet60FilterBit = 0;
      fCaloJet80FilterBit = 0;
      fCaloJet100FilterBit = 0;

    } else { // PbPb MC

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
  
  // Connect the branches to the skim tree (different for pp and PbPb Monte Carlo)
  if(fMegaSkimMode){

    // In mega skim mode all the event selections are already done. Set all flags to 1.
    fPrimaryVertexFilterBit = 1;
    fBeamScrapingFilterBit = 1;
    fHBHENoiseFilterBit = 1;
    fHfCoincidenceFilterBit = 1;
    fClusterCompatibilityFilterBit = 1;
    fPileupFilterBit = 1;

  } else {

    fSkimTree->SetBranchStatus("*",0);
    if(fDataType == kPpMC){ // pp MC
      fSkimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
      fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
      fSkimTree->SetBranchStatus("pBeamScrapingFilter",1);
      fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    
      if(fIsMiniAOD){
        fHBHENoiseFilterBit = 1; // This filter bit is not available in MiniAODs
      } else {
        fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
        fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
      }
    
      fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
      fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
      fPileupFilterBit = 1;        // No pile-up filter for pp
    } else if(fIsPPb){ // pPb  MC

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

    } else if(fDataType == kOOMC){ // OO data or MC

      // Primary vertex has at least two tracks, is within 25 cm in z-direction and within 2 cm in xy-direction
      fSkimTree->SetBranchStatus("pprimaryVertexFilter", 1);
      fSkimTree->SetBranchAddress("pprimaryVertexFilter", &fPrimaryVertexFilterBit, &fPrimaryVertexBranch);
     
      // Have at least two towers on both of the HF calorimerter to have energies above 4 GeV
      // accumulated by the energies of PF (particle-flow)candidates
      fSkimTree->SetBranchStatus("OOphfCoincFilterPF2Th4", 1);
      fSkimTree->SetBranchAddress("OOphfCoincFilterPF2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);

      fSkimTree->SetBranchStatus("pileupVertexFilter", 1);
      fSkimTree->SetBranchAddress("pileupVertexFilter", &fPileupFilterBit, &fPileupFilterBranch);

      // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible
      //fSkimTree->SetBranchStatus("pclusterCompatibilityFilter", 1);
      //fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fClusterCompatibilityFilterBit, &fClusterCompatibilityBranch);

      //fPrimaryVertexFilterBit = 1; // No primary vertex filter for OO
      //fHfCoincidenceFilterBit = 1; // No HF coincidence filter for OO
      //fPileupFilterBit = 1; // No pile-up filter for OO
      fHBHENoiseFilterBit = 1; // HBHE noise filter bit is not available in the OO MiniAOD forests.
      fBeamScrapingFilterBit = 1;  // No beam scraping filter for OO
      fClusterCompatibilityFilterBit = 1; // No Cluster compatibility, no good for 2025 data

     } else if(fDataType == kPpRef5p36TeV || fDataType == kPpRef5p36TeVMC){ // ppRef data or MC 5.36TeV

      // Primary vertex has at least two tracks, is within 25 cm in z-direction and within 2 cm in xy-direction
      fSkimTree->SetBranchStatus("pprimaryVertexFilter", 1);
      fSkimTree->SetBranchAddress("pprimaryVertexFilter", &fPrimaryVertexFilterBit, &fPrimaryVertexBranch);
     
      // Have at least two towers on both of the HF calorimerter to have energies above 4 GeV
      // accumulated by the energies of PF (particle-flow)candidates
      //fSkimTree->SetBranchStatus("OOpfCoincFilterPF2Th4", 1);
      //fSkimTree->SetBranchAddress("OOpfCoincFilterPF2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);

      //fSkimTree->SetBranchStatus("pileupVertexFilter", 1);
      //fSkimTree->SetBranchAddress("pileupVertexFilter", &fPileupFilterBit, &fPileupFilterBranch);

      // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible
      //fSkimTree->SetBranchStatus("pclusterCompatibilityFilter", 1);
      //fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fClusterCompatibilityFilterBit, &fClusterCompatibilityBranch);

      //fPrimaryVertexFilterBit = 1; // No primary vertex filter for OO
      fHfCoincidenceFilterBit = 1; // No HF coincidence filter for OO
      fPileupFilterBit = 1; // No pile-up filter for OO
      fHBHENoiseFilterBit = 1; // HBHE noise filter bit is not available in the OO MiniAOD forests.
      fBeamScrapingFilterBit = 1;  // No beam scraping filter for OO
      fClusterCompatibilityFilterBit = 1; // No Cluster compatibility, no good for 2025 data

    } else { // PbPb MC
    
      // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction
      fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
      fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    
      if(fIsMiniAOD){
      
        fHBHENoiseFilterBit = 1; // This filter bit is not available in MiniAODs
      
        // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
        fSkimTree->SetBranchStatus("pphfCoincFilter2Th4",1);
        fSkimTree->SetBranchAddress("pphfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
      
      } else {
      
        // Cut on noise on HCAL
        fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
        fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
      
        // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
        fSkimTree->SetBranchStatus("phfCoincFilter2Th4",1);
        fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
      
      }

      // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible
      fSkimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
      fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fClusterCompatibilityFilterBit, &fClusterCompatibilityBranch);
    
      fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
      fPileupFilterBit = 1;        // No pile-up filter for PbPb
    }
  }
  
  // Connect the branches to the track tree
  if(fReadTrackTree){
    fTrackTree->SetBranchStatus("*",0);
    fTrackTree->SetBranchStatus("pt",1);
    fTrackTree->SetBranchAddress("pt",&fTrackPtArray,&fTrackPtBranch);
    fTrackTree->SetBranchStatus("phi",1);
    fTrackTree->SetBranchAddress("phi",&fTrackPhiArray,&fTrackPhiBranch);
    fTrackTree->SetBranchStatus("eta",1);
    fTrackTree->SetBranchAddress("eta",&fTrackEtaArray,&fTrackEtaBranch);
    fTrackTree->SetBranchStatus("sube",1);
    fTrackTree->SetBranchAddress("sube",&fTrackSubeventArray,&fTrackChi2Branch);  // Reuse a branch from ForestReader that is not otherwise needed here

    // The track charge is not available in mega skims
    if(!fMegaSkimMode){
      fTrackTree->SetBranchStatus("chg",1);
      fTrackTree->SetBranchAddress("chg",&fTrackChargeArray,&fTrackPtErrorBranch);  // Reuse a branch from ForestReader that is not otherwise needed here
    }
  } // Reading track trees

  // We need to load one event to initialize TChains properly. Not sure why, but this is how things seem to work
  GetEvent(0);

}

/*
 * Connect a new tree to the reader
 */
void GeneratorLevelForestReader::ReadForestFromFile(TFile *inputFile){
  std::vector<TString> fileList;
  fileList.push_back(inputFile->GetName());
  ReadForestFromFileList(fileList);
}

/*
 * Connect a new tree to the reader
 */
void GeneratorLevelForestReader::ReadForestFromFileList(std::vector<TString> fileList){

  // Open one file to determine if we are dealing with miniAOD or regular AOD
  TFile *inputFile = TFile::Open(fileList.at(0));
  TTree* miniAODcheck = (TTree*)inputFile->Get("HiForestInfo/HiForest");
  fIsMiniAOD = !(miniAODcheck == NULL);
  inputFile->Close();

  // Helper variable for finding the correct tree
  const char *treeName[4] = {"none","none","none","none"};

  // Connect a trees from the file to the reader
  fHeavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
  if(fUseJetTrigger) fHltTree = new TChain("hltanalysis/HltTree");
  if(!fMegaSkimMode) fSkimTree = new TChain("skimanalysis/HltTree");

  // The jet tree has different name in different datasets
  if(fDataType == kPp || fDataType == kPpMC || fDataType == kPpRef5p36TeV || fDataType == kPpRef5p36TeVMC){
    treeName[0] = "ak4CaloJetAnalyzer/t"; // Tree for calo jets
    treeName[1] = "ak4PFJetAnalyzer/t";   // Tree for PF jets
  } else { // PbPb data or MC
    treeName[0] = "akPu4CaloJetAnalyzer/t";     // Tree for calo jets
    treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for csPF jets
    treeName[2] = "akPu4PFJetAnalyzer/t";       // Tree for puPF jets
    treeName[3] = "akFlowPuCs4PFJetAnalyzer/t"; // Tree for flow subtracted csPF jets
  }

  if(!fMixingMode) fJetTree = new TChain(treeName[fJetType]);

  // The track tree is different than in other types of forests
  if(fReadTrackTree) fTrackTree = new TChain("HiGenParticleAna/hi");

  for(std::vector<TString>::iterator listIterator = fileList.begin(); listIterator != fileList.end(); listIterator++){
    fHeavyIonTree->Add(*listIterator);
    if(!fMegaSkimMode) fSkimTree->Add(*listIterator);
    if(fUseJetTrigger) fHltTree->Add(*listIterator);
    if(!fMixingMode) fJetTree->Add(*listIterator);
    if(fReadTrackTree) fTrackTree->Add(*listIterator);
  }
  
  // Connect branches to trees
  Initialize();
}

/*
 * Burn the current forest.
 */
void GeneratorLevelForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fJetTree->Delete();
  fHltTree->Delete();
  fSkimTree->Delete();
  fTrackTree->Delete();
}

/*
 * Check if there are problems in the file list
 */
bool GeneratorLevelForestReader::CheckFileProblems(){

  // Loop over all the files that are in TChains and check they are doing fine
  TObjArray *fileElements = fHeavyIonTree->GetListOfFiles();
  for (TObject *op: *fileElements) {
    auto chainElement = static_cast<TChainElement *>(op);
    TFile *testFile = TFile::Open(chainElement->GetTitle());

    // Check that the file is still open
    if(!testFile->IsOpen()) return true;

    // Check that the file did not turn into a zombie
    if(testFile->IsZombie()) return true;

    // If everything is fine with the file, we can close the test file
    testFile->Close();
  }

  // If no problems are found from any of the files, return false
  return false;

}

/*
 * Load an event to memory
 */
void GeneratorLevelForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  if(!fMixingMode) fJetTree->GetEntry(nEvent);
  if(fUseJetTrigger) fHltTree->GetEntry(nEvent);
  if(!fMegaSkimMode) fSkimTree->GetEntry(nEvent);
  if(fReadTrackTree) {
    fTrackTree->GetEntry(nEvent);
    
    // Read the numbers of tracks for this event
    fnTracks = fTrackPtArray->size();
  }
}

// Getter for jet pT
Float_t GeneratorLevelForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t GeneratorLevelForestReader::GetJetPhi(Int_t iJet) const{
  if(fJetAxis == 0) return fJetPhiArray[iJet];
  return fJetWTAPhiArray[iJet];
}

// Getter for jet eta
Float_t GeneratorLevelForestReader::GetJetEta(Int_t iJet) const{
  if(fJetAxis == 0) return fJetEtaArray[iJet];
  return fJetWTAEtaArray[iJet];
}

// Getter for jet raw pT (not relevant for generator jets, just return regular jet pT)
Float_t GeneratorLevelForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetPtArray[iJet];;
}

// Getter for maximum track pT inside a jet (not relevant for generator jets, return negative value to show this)
Float_t GeneratorLevelForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return -1;
}

// Getter for track pT
Float_t GeneratorLevelForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtArray->at(iTrack);
}

// Getter for track phi
Float_t GeneratorLevelForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiArray->at(iTrack);
}

// Getter for track eta
Float_t GeneratorLevelForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaArray->at(iTrack);
}

// Getter for track charge
Int_t GeneratorLevelForestReader::GetTrackCharge(Int_t iTrack) const{
  return fTrackChargeArray->at(iTrack);
}

// Getter for track subevent index
Int_t GeneratorLevelForestReader::GetTrackSubevent(Int_t iTrack) const{
  return fTrackSubeventArray->at(iTrack);
}

// Getter for track MC status
Int_t GeneratorLevelForestReader::GetTrackMCStatus(Int_t iTrack) const{
  return 1;
}

// Getter for the number of generator level particles. Same as nTracks for generator level forest reader
Int_t GeneratorLevelForestReader::GetNParticles() const{
  return fnTracks;
}

// Getter for generator level particle pT. Same as regular track pT for generator level forest
Float_t GeneratorLevelForestReader::GetParticlePt(Int_t iParticle) const{
  return GetTrackPt(iParticle);
} 

// Getter for generator level particle phi. Same as regular track phi for generator level forest
Float_t GeneratorLevelForestReader::GetParticlePhi(Int_t iParticle) const{
  return GetTrackPhi(iParticle);
} 

// Getter for generator level particle eta. Same as regular track eta for generator level forest
Float_t GeneratorLevelForestReader::GetParticleEta(Int_t iParticle) const{
  return GetTrackEta(iParticle);
} 

// Getter for generator level particle charge       
Int_t GeneratorLevelForestReader::GetParticleCharge(Int_t iParticle) const{
  return GetTrackCharge(iParticle);
}

// Getter for generator level particle subevent index. Same as regular subevent index for generator level forest
Int_t GeneratorLevelForestReader::GetParticleSubevent(Int_t iParticle) const{
  return GetTrackSubevent(iParticle);
}

// Getter for track pT error (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackPtError(Int_t iTrack) const{
  return 0; // Setting all errors to 0 always passes the track quality cut
}

// Getter for high purity of the track (not relevant for generator tracks)
Bool_t GeneratorLevelForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return true; // All the generator tracks are of high purity
}

// Getter for track distance from primary vertex in z-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in z-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for track distance from primary vertex in xy-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in xy-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for normalized track chi2 value from reconstruction fit (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackNormalizedChi2(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for track chi2 value from reconstruction fit (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackChi2(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of degrees of freedom in reconstruction fit (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of hits in tracker layers (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of hits for the track (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNHitsTrack(Int_t iTrack) const{
  return 1; // Note: The cut on NHits is disabled for generator level tracks in the main analysis code
}

// Getter for track energy in ECal (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et is disabled for generator level tracks in the main analysis code
}

// Getter for track energy in HCal (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et is disabled for generator level tracks in the main analysis code
}

// Check if generator level jet has a matching reconstructed jet
Bool_t GeneratorLevelForestReader::HasMatchingJet(Int_t iJet) const{
  
  // If not matching jets, just tell that everything is fine
  if(!fMatchJets) return true;
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, check also eta and phi
  // If all values are close by, we must have a matching jet
  Double_t jetPt = GetJetPt(iJet);
  for(Int_t iRef = 0; iRef < fnMatchedJets; iRef++){
    if(TMath::Abs(fJetEtaArray[iJet] - fJetRefEtaArray[iRef]) < 0.015){
      if(TMath::Abs(fJetPhiArray[iJet] - fJetRefPhiArray[iRef]) < 0.015){
        if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.03*jetPt) {
          return true;
        }
      }
    }
  }
  
  return false;
}

// Get the pT of the matched reconstructed jet
Int_t GeneratorLevelForestReader::GetMatchingIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, check also eta and phi
  // If all values are close by, we must have a matching jet
  Double_t jetPt = GetJetPt(iJet);
  Int_t matchingIndex = -1;
  for(Int_t iRef = 0; iRef < fnMatchedJets; iRef++){
    if(TMath::Abs(fJetEtaArray[iJet] - fJetRefEtaArray[iRef]) < 0.015){
      if(TMath::Abs(fJetPhiArray[iJet] - fJetRefPhiArray[iRef]) < 0.015){
        if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.03*jetPt) {
          matchingIndex = iRef;
          break;
        }
      }
    }
  }
  
  // Return the matching index
  return matchingIndex;
}

// Get the pT of the matched reconstructed jet
Float_t GeneratorLevelForestReader::GetMatchedPt(Int_t iJet) const{
  
  // If not matching jets, just return 0 as this value has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet pT. Need raw pT for reco jets before jet corrections are in the forest
  return fRecoJetRawPtArray[matchingIndex]; // Use this if doing jet pT correction manually
  //return fRecoJetPtArray[matchingIndex]; // Use this if jet pT corrected in the forest
}

// Get the phi of the matched reconstructed jet
Float_t GeneratorLevelForestReader::GetMatchedPhi(Int_t iJet) const{
  
  // If not matching jets, just return 0 as this value has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet phi
  return fRecoJetPhiArray[matchingIndex];
}

// Get the eta of the matched reconstructed jet
Float_t GeneratorLevelForestReader::GetMatchedEta(Int_t iJet) const{
  
  // If not matching jets, just return 0 as this value has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet eta
  return fRecoJetEtaArray[matchingIndex];
}

// Check if generator level jet has a matching reconstructed jet
Int_t GeneratorLevelForestReader::GetPartonFlavor(Int_t iJet) const{
  
  // If not matching jets, just return 0 as this value has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching parton flavor
  return fJetRefFlavorArray[matchingIndex];
}

// Class for the main analysis algorithms for the energy-energy correlator analysis

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "EECAnalyzer.h"

using namespace std;

/*
 * Total weight function to match multiplicity in MC to that in data.
 *
 * Derived using the multiplicityPlotter macro. First the shape of the multiplicity distribution is determined
 * from the non-efficiency corrected distribution. Then the multiplicity boundaries are obtained from the
 * efficiency corrected distribution. This gives the base value below. The weight comes from the difference
 * between the multiplicity distribution in MC from a flat line.
 *
 * The used data files:
 *  dijetPbPb2018_akCaloJet_onlyJets_jet100TriggerEta1v3_withTrackEff_processed_2022-01-19.root
 *  PbPbMC2018_RecoGen_akCaloJet_onlyJets_noCentShift_noCentWeight_jetEta1v3_processed_2022-01-21.root
 */
double totalMultiplicityWeight(double *x, double *par){
  
  double weight = 0;
  double base = 0;
  
  if(x[0] < 600){
    weight = (0.0309424+0.110670*TMath::Exp(3000*(-0.00138397))) / (0.0685758+0.374178*TMath::Exp(x[0]*(-0.00550382)));
  } else {
    weight = (0.0309424+0.110670*TMath::Exp(3000*(-0.00138397))) / (0.0309424+0.110670*TMath::Exp(x[0]*(-0.00138397)));
  }
  
  // Gauss function, if multiplicity > 3000
  if(x[0] > 3000) {
    base = TMath::Exp(-0.5*TMath::Power((3000-x[0])/215,2));
  } else if(x[0] < 300) {
    // Second order polynomial, if multiplicity < 300
    base = (-10.5233 + 1.47193 * x[0] -0.0024 * x[0] * x[0]) / 503.0;
  } else {
    base = (0.10665541 * x[0] + 183.00338) / 503.0;
  }
  
  return base * weight;
  
}

/*
 * Default constructor
 */
EECAnalyzer::EECAnalyzer() :
  fFileNames(0),
  fCard(0),
  fHistograms(0),
  fVzWeightFunction(0),
  fCentralityWeightFunctionCentral(0),
  fCentralityWeightFunctionPeripheral(0),
  fMultiplicityWeightFunction(0),
  fPtWeightFunction(0),
  fSmearingFunction(0),
  fTrackEfficiencyCorrector(),
  fJetCorrector(),
  fJetUncertainty(),
  fRng(0),
  fDataType(-1),
  fTriggerSelection(0),
  fJetType(0),
  fMatchJets(0),
  fDebugLevel(0),
  fIsRealData(true),
  fIsPPbData(false),
  fIsPpData(false),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1),
  fJetAxis(0),
  fVzCut(0),
  fMinimumPtHat(0),
  fMaximumPtHat(0),
  fJetEtaCut(0),
  fJetMinimumPtCut(0),
  fJetMaximumPtCut(0),
  fCutBadPhiRegion(false),
  fCutWithJetEtaCM(false),
  fMinimumMaxTrackPtFraction(0),
  fMaximumMaxTrackPtFraction(0),
  fJetUncertaintyMode(0),
  fTrackEtaCut(0),
  fTrackMinPtCut(0),
  fTrackMaxPtCut(0),
  fMaxTrackPtRelativeError(0),
  fMaxTrackDistanceToVertex(0),
  fCalorimeterSignalLimitPt(0),
  fHighPtEtFraction(0),
  fChi2QualityCut(0),
  fMinimumTrackHits(0),
  fSubeventCut(0),
  fTrackEfficiencyVariation(0.024),
  fTrackSelectionVariation(""),
  fJetPtWeightConfiguration(0),
  fReconstructedJetMinimumPtCut(0),
  fGeneratorJetMinimumPtCut(0),
  fSkipCovarianceMatrix(false),
  fMcCorrelationType(0),
  fJetRadius(0.4),
  fWeightExponent(0),
  fSmearDeltaR(false),
  fSmearEnergyWeight(false),
  fDoReflectedCone(false),
  fDoMixedCone(false),
  fMegaSkimMode(false),
  fDoPerpendicularCone(false),
  fLocalRun(0),
  fLocalMixing(0),
  fMixingListIndex(1),
  fMixingStartIndex(0),
  fRunningMixingIndex(0),
  fnEventsInMixingFile(0),
  fMixedEventVz(0),
  fMixedEventHiBin(0),
  fMixedEventMultiplicity(0),
  fMixedEventHFEnergy(0),
  fMixedEventEventNumber(0),
  fMixedEventVzBin(0),
  fMixedEventHiBinBin(0),
  fMixedEventMultiplicityBin(0),
  fMixedEventHFEnergyBin(0),
  fEventMatchForMixing(0),
  fBinnedMixing(false),
  fVzMixingBinBorders(0),
  fHiBinMixingBinBorders(0),
  fHFEnergyMixingBinBorders(0),
  fMultiplicityMixingBinBorders(0),
  fHFEnergyMatchingShift(0),
  fFillEventInformation(false),
  fFillJetHistograms(false),
  fFillTrackHistograms(false),
  fFillJetConeHistograms(false),
  fFillEnergyEnergyCorrelators(false),
  fFillEnergyEnergyCorrelatorsSystematics(false),
  fFillJetPtClosure(false),
  fFillJetPtUnfoldingResponse(false),
  fFillTrackParticleMatchingHistograms(false),
  fMultiplicityMode(false)
{
  // Default constructor
  fHistograms = new EECHistograms();
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fJetReader = NULL;
  fRecoJetReader = NULL;
  fTrackReader = NULL;
  fMixedEventReader = NULL;
  
  // Create a corrector for track pair efficiency
  fTrackPairEfficiencyCorrector = new TrackPairEfficiencyCorrector();

  // Create a manager for jet energy resolution smearing in MC
  fEnergyResolutionSmearingFinder = new JetMetScalingFactorManager();

  // Create smearing providers for DeltaR and energy weights
  fDeltaRSmearer = new SmearingProvider();
  fEnergyWeightSmearer = new SmearingProvider();

  for(int oBin = 0; oBin < 5; oBin++){
    for(int iBin = 0; iBin < 10; iBin++){
      fThisEventCorrelator[oBin][iBin] = NULL;
    }
  }
}

/*
 * Custom constructor
 */
EECAnalyzer::EECAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard, Int_t runLocal, Int_t mixingListIndex) :
  fFileNames(fileNameVector),
  fCard(newCard),
  fHistograms(0),
  fJetCorrector(),
  fJetUncertainty(),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1),
  fMixingListIndex(mixingListIndex)
{
  // Custom constructor
  fHistograms = new EECHistograms(fCard);
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fJetReader = NULL;
  fRecoJetReader = NULL;
  fTrackReader = NULL;
  fMixedEventReader = NULL;

  // Decode indices for local running to read mixing files from proper location
  if(runLocal > 1){
    fLocalRun = 0;
    fLocalMixing = 1;
  } else {
    fLocalRun = runLocal;
    fLocalMixing = 0;
  }
  
  // Configurure the analyzer from input card
  ReadConfigurationFromCard();
  
  // pT weight function for Pythia to match 2017 MC and data pT spectra. Derived from all jets above 120 GeV
  fPtWeightFunction = new TF1("fPtWeightFunction","pol3",0,500);
  //fPtWeightFunction->SetParameters(0.699073,0.00287672,-6.35568e-06,5.02292e-09);
  //fPtWeightFunction->SetParameters(0.708008,0.0032891,-1.05716e-05,1.16656e-08); // From JECv4
  fPtWeightFunction->SetParameters(0.79572,0.0021861,-6.35407e-06,6.66435e-09); // From JECv6
  
  // Function for smearing the jet pT for systemtic uncertainties
  fSmearingFunction = new TF1("fSmearingFunction","pol4",0,500);
  
  // Find the correct folder for track correction tables based on data type
  if(fIsPpData){
    
    // Track correction for 2017 pp data
    fTrackEfficiencyCorrector = new TrkEff2017pp(false, "trackCorrectionTables/pp2017/");
    
    // Weight function for 2017 MC
    fVzWeightFunction = new TF1("fvz","pol6",-15,15);  // Weight function for 2017 MC
    fVzWeightFunction->SetParameters(0.973805, 0.00339418, 0.000757544, -1.37331e-06, -2.82953e-07, -3.06778e-10, 3.48615e-09);
    fCentralityWeightFunctionCentral = NULL;
    fCentralityWeightFunctionPeripheral = NULL;
    fMultiplicityWeightFunction = NULL;
    
    // Track pair efficiency corrector for pp
    fTrackPairEfficiencyCorrector = new TrackPairEfficiencyCorrector("trackCorrectionTables/trackPairEfficiencyCorrection_pp2017_32DeltaRBins_2023-07-11.root", false, false);

    // Jet energy resolution smearing scale factor manager
    if(fJetUncertaintyMode > 0){
      fEnergyResolutionSmearingFinder = new JetMetScalingFactorManager(false, fJetUncertaintyMode-3);
    } else {
      fEnergyResolutionSmearingFinder = new JetMetScalingFactorManager();
    }

    // Smearing configuration for DeltaR and energy weight in energy-energy correlators
    if(fSmearDeltaR){
      const char* deltaRFileName = "ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_trackParticleResponseMatrixWithTrackPairEfficiency_processed_2023-10-30.root";
      fDeltaRSmearer = new SmearingProvider(Form("smearingFiles/%s", deltaRFileName), "particleDeltaRResponseMatrix", false, false);
      cout << "Smearing DeltaR in energy-energy correlators using the file " << deltaRFileName << endl;
    } else {
      fDeltaRSmearer = NULL;
    }

    if(fSmearEnergyWeight){
      const char* energyWeightFileName = "ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_trackParticleResponseMatrixWithTrackPairEfficiency_processed_2023-10-30.root";
      fEnergyWeightSmearer = new SmearingProvider(Form("smearingFiles/%s", energyWeightFileName), "particlePairPtClosure", true, false);
      cout << "Smearing energy weights for energy-energy correlators using the file " << energyWeightFileName << endl;
    } else {
      fEnergyWeightSmearer = NULL;
    } 
  
  } else if(fIsOOData) {
    
   // 1. Skip Track Corrections: Return a dummy corrector or NULL
    fTrackEfficiencyCorrector = NULL; 

    // 2. Apply PbPb Weights: Use the parameters you provided for Vz and Centrality
    fVzWeightFunction = new TF1("fvz","pol6",-15,15);
    fVzWeightFunction->SetParameters(1.00591, -0.0193751, 0.000961142, -2.44303e-05, -8.24443e-06, 1.66679e-07, 1.11028e-08);
    
    fCentralityWeightFunctionCentral = new TF1("fCentralWeight","pol6",0,30);
    fCentralityWeightFunctionCentral->SetParameters(4.73421, -0.0477343, -0.0332804, 0.00355699, -0.00017427, 4.18398e-06, -3.94746e-08);
    
    fCentralityWeightFunctionPeripheral = new TF1("fPeripheralWeight","pol6",30,90);
    fCentralityWeightFunctionPeripheral->SetParameters(3.38091, -0.0609601, -0.00228529, 9.43076e-05, -1.39593e-06, 9.85435e-09, -2.77153e-11);
    
    fMultiplicityWeightFunction = new TF1("fMultiWeight", totalMultiplicityWeight, 0, 5000, 0);

    // 3. Skip JEC: Do not initialize fJetCorrector here
    // In the JEC section, ensure OO does not call AddLevel()

    // 4. Skip Mixing: Disable flags
    fDoMixedCone = false;
    fDoPerpendicularCone = false;

    // Set remaining helpers to NULL to ensure they are skipped
    fTrackPairEfficiencyCorrector = NULL;
    fEnergyResolutionSmearingFinder = NULL;
    fDeltaRSmearer = NULL;
    fEnergyWeightSmearer = NULL;

  } 

  else if(fIsPPbData){
    
    // Track correction for 2016 pPb data
    fTrackEfficiencyCorrector = new TrkEff2016pPb(false, "trackCorrectionTables/pPb2016/");
    
    // Weight function for 2016 pPb MC. Copied from Dener's analysis
    fVzWeightFunction = new TF1("VzWeightFunction", "pol8", -15.1, 15.1, TF1::EAddToList::kNo);
    fVzWeightFunction->SetParameters(0.856516, -0.0159813, 0.00436628, -0.00012862, 2.61129e-05, -4.16965e-07, 1.73711e-08, -3.11953e-09, 6.24993e-10);
    fCentralityWeightFunctionCentral = NULL;
    fCentralityWeightFunctionPeripheral = NULL;
    fMultiplicityWeightFunction = NULL;
    
    // Track pair efficiency corrector for pp TODO: Pair efficiency correction copied from pp!
    fTrackPairEfficiencyCorrector = new TrackPairEfficiencyCorrector("trackCorrectionTables/trackPairEfficiencyCorrection_pp2017_32DeltaRBins_2023-07-11.root", false, false);

    // Jet energy resolution smearing scale factor manager
    if(fJetUncertaintyMode > 0){
      fEnergyResolutionSmearingFinder = new JetMetScalingFactorManager(false, fJetUncertaintyMode-3);
    } else {
      fEnergyResolutionSmearingFinder = new JetMetScalingFactorManager();
    }

    // Smearing configuration for DeltaR and energy weight in energy-energy correlators TODO: Smearing here copied from pp
    if(fSmearDeltaR){
      const char* deltaRFileName = "ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_trackParticleResponseMatrixWithTrackPairEfficiency_processed_2023-10-30.root";
      fDeltaRSmearer = new SmearingProvider(Form("smearingFiles/%s", deltaRFileName), "particleDeltaRResponseMatrix", false, false);
      cout << "Smearing DeltaR in energy-energy correlators using the file " << deltaRFileName << endl;
    } else {
      fDeltaRSmearer = NULL;
    }

    if(fSmearEnergyWeight){
      const char* energyWeightFileName = "ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_trackParticleResponseMatrixWithTrackPairEfficiency_processed_2023-10-30.root";
      fEnergyWeightSmearer = new SmearingProvider(Form("smearingFiles/%s", energyWeightFileName), "particlePairPtClosure", true, false);
      cout << "Smearing energy weights for energy-energy correlators using the file " << energyWeightFileName << endl;
    } else {
      fEnergyWeightSmearer = NULL;
    } 
    
  } else if (fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPbPbMC){
    
    // Track correction for 2018 PbPb data
    fTrackEfficiencyCorrector = new TrkEff2018PbPb("general", fTrackSelectionVariation, false, "trackCorrectionTables/PbPb2018/");
    
    // The vz weight function is rederived from the miniAOD dataset.
    // Macro used for derivation: deriveMonteCarloWeights.C, Git hash: d4eab1cd188da72f5a81b8902cb6cc55ea1baf23
    // Input files: eecAnalysis_akFlowJet_onlyJets_weightEventInfo_combinedTriggers_processed_2023-03-06.root
    //              PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_wtaAxis_onlyJets_noTrigger_ptHatWeight_processed_2023-03-06.root
    fVzWeightFunction = new TF1("fvz","pol6",-15,15);
    fVzWeightFunction->SetParameters(1.00591, -0.0193751, 0.000961142, -2.44303e-05, -8.24443e-06, 1.66679e-07, 1.11028e-08);
    
    // The centrality weight function is rederived for the miniAOD dataset.
    // Macro used for derivation: deriveMonteCarloWeights.C, Git hash: d4eab1cd188da72f5a81b8902cb6cc55ea1baf23
    // Input files: eecAnalysis_akFlowJet_onlyJets_weightEventInfo_combinedTriggers_processed_2023-03-06.root
    //              PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_wtaAxis_onlyJets_noTrigger_ptHatWeight_processed_2023-03-06.root
    fCentralityWeightFunctionCentral = new TF1("fCentralWeight","pol6",0,30);
    fCentralityWeightFunctionCentral->SetParameters(4.73421, -0.0477343, -0.0332804, 0.00355699, -0.00017427, 4.18398e-06, -3.94746e-08);
    fCentralityWeightFunctionPeripheral = new TF1("fPeripheralWeight","pol6",30,90);
    fCentralityWeightFunctionPeripheral->SetParameters(3.38091, -0.0609601, -0.00228529, 9.43076e-05, -1.39593e-06, 9.85435e-09, -2.77153e-11);
    
    // Multiplicity based weight function
    fMultiplicityWeightFunction = new TF1("fMultiWeight", totalMultiplicityWeight, 0, 5000, 0);
    
    // Track pair efficiency corrector for PbPb
    fTrackPairEfficiencyCorrector = new TrackPairEfficiencyCorrector("trackCorrectionTables/trackPairEfficiencyCorrection_PbPb2018_32DeltaRBins_2023-07-11.root", true, true);

    // Jet energy resolution smearing scale factor manager
    if(fJetUncertaintyMode > 0){
      fEnergyResolutionSmearingFinder = new JetMetScalingFactorManager(true, fJetUncertaintyMode-3);
    } else {
      fEnergyResolutionSmearingFinder = new JetMetScalingFactorManager();
    }

    // Smearing configuration for DeltaR and energy weight in energy-energy correlators
    if(fSmearDeltaR){
      const char* deltaRFileName = "PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_trackParticleResponseMatrixWithTrackPairEfficiency_processed_2023-10-30.root";
      fDeltaRSmearer = new SmearingProvider(Form("smearingFiles/%s", deltaRFileName), "particleDeltaRResponseMatrix", false, true);
      cout << "Smearing DeltaR in energy-energy correlators using the file " << deltaRFileName << endl;
    } else {
      fDeltaRSmearer = NULL;
    }

    if(fSmearEnergyWeight){
      const char* energyWeightFileName = "PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_trackParticleResponseMatrixWithTrackPairEfficiency_processed_2023-10-30.root";
      fEnergyWeightSmearer = new SmearingProvider(Form("smearingFiles/%s", energyWeightFileName), "particlePairPtClosure", true, true);
      cout << "Smearing energy weights for energy-energy correlators using the file " << energyWeightFileName << endl;
    } else {
      fEnergyWeightSmearer = NULL;
    }
    
  } else {
    fVzWeightFunction = NULL;
    fCentralityWeightFunctionCentral = NULL;
    fCentralityWeightFunctionPeripheral = NULL;
    fMultiplicityWeightFunction = NULL;
    fTrackPairEfficiencyCorrector = NULL;
    fEnergyResolutionSmearingFinder = NULL;
    fDeltaRSmearer = NULL;
    fEnergyWeightSmearer = NULL;
  }
  
  // Initialize the random number generator with a random seed
  fRng = new TRandom3();
  fRng->SetSeed(0);

  // Energy-energy correlator constructed only from this event
  const Int_t nDeltaRBins = fHistograms->GetNDeltaRBinsEEC();
  Double_t* deltaRBins = new Double_t[nDeltaRBins+1];
  for(Int_t iBin = 1; iBin <= nDeltaRBins; iBin++){
    deltaRBins[iBin-1] = fHistograms->GetDeltaRBinBorderLowEEC(iBin);
  }
  deltaRBins[nDeltaRBins] = fHistograms->GetDeltaRBinBorderHighEEC(nDeltaRBins);

  const Int_t nTrackPtBinsEEC = fCard->GetNBin("TrackPtBinEdgesEEC");
  for(auto iWeightExponent = 0; iWeightExponent < fWeightExponent.size(); iWeightExponent++){
    for(Int_t iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      fThisEventCorrelator[iWeightExponent][iTrackPt] = new TH1D(Form("thisEventCorrelator%d%d", iWeightExponent, iTrackPt), Form("thisEventCorrelator%d%d", iWeightExponent, iTrackPt), nDeltaRBins, deltaRBins);
    }
  }

  // Release the memory reserved with new
  delete[] deltaRBins;
  
  //**********************************************************
  //    Disable/enable the track pair efficiency correction
  //**********************************************************

  bool disableTrackPairEfficiencyCorrection = (fCard->Get("DisableTrackPairEfficiencyCorrection") == 1);
  if((fMcCorrelationType == kGenGen) || (fMcCorrelationType == kRecoGen)) disableTrackPairEfficiencyCorrection = true; // Disable the track pair efficiency correction for generator level particles
  if((fDataType == ForestReader::kOO || fDataType == ForestReader::kOOMC)) disableTrackPairEfficiencyCorrection = true; // Disable the track pair efficiency correction for OO data sets

  if(fTrackPairEfficiencyCorrector != NULL) {
    fTrackPairEfficiencyCorrector->SetDisableCorrection(fCard->Get("DisableTrackPairEfficiencyCorrection"));
  }
}

/*
 * Copy constructor
 */
EECAnalyzer::EECAnalyzer(const EECAnalyzer& in) :
  fJetReader(in.fJetReader),
  fRecoJetReader(in.fRecoJetReader),
  fTrackReader(in.fTrackReader),
  fMixedEventReader(in.fMixedEventReader),
  fFileNames(in.fFileNames),
  fCard(in.fCard),
  fHistograms(in.fHistograms),
  fVzWeightFunction(in.fVzWeightFunction),
  fCentralityWeightFunctionCentral(in.fCentralityWeightFunctionCentral),
  fCentralityWeightFunctionPeripheral(in.fCentralityWeightFunctionPeripheral),
  fMultiplicityWeightFunction(in.fMultiplicityWeightFunction),
  fPtWeightFunction(in.fPtWeightFunction),
  fSmearingFunction(in.fSmearingFunction),
  fRng(in.fRng),
  fDataType(in.fDataType),
  fTriggerSelection(in.fTriggerSelection),
  fJetType(in.fJetType),
  fMatchJets(in.fMatchJets),
  fDebugLevel(in.fDebugLevel),
  fIsRealData(in.fIsRealData),
  fIsPPbData(in.fIsPPbData),
  fIsPpData(in.fIsPpData),
  fVzWeight(in.fVzWeight),
  fCentralityWeight(in.fCentralityWeight),
  fPtHatWeight(in.fPtHatWeight),
  fTotalEventWeight(in.fTotalEventWeight),
  fJetAxis(in.fJetAxis),
  fVzCut(in.fVzCut),
  fMinimumPtHat(in.fMinimumPtHat),
  fMaximumPtHat(in.fMaximumPtHat),
  fJetEtaCut(in.fJetEtaCut),
  fJetMinimumPtCut(in.fJetMinimumPtCut),
  fJetMaximumPtCut(in.fJetMaximumPtCut),
  fCutBadPhiRegion(in.fCutBadPhiRegion),
  fCutWithJetEtaCM(in.fCutWithJetEtaCM),
  fMinimumMaxTrackPtFraction(in.fMinimumMaxTrackPtFraction),
  fMaximumMaxTrackPtFraction(in.fMaximumMaxTrackPtFraction),
  fJetUncertaintyMode(in.fJetUncertaintyMode),
  fTrackEtaCut(in.fTrackEtaCut),
  fTrackMinPtCut(in.fTrackMinPtCut),
  fTrackMaxPtCut(in.fTrackMaxPtCut),
  fMaxTrackPtRelativeError(in.fMaxTrackPtRelativeError),
  fMaxTrackDistanceToVertex(in.fMaxTrackDistanceToVertex),
  fCalorimeterSignalLimitPt(in.fCalorimeterSignalLimitPt),
  fHighPtEtFraction(in.fHighPtEtFraction),
  fChi2QualityCut(in.fChi2QualityCut),
  fMinimumTrackHits(in.fMinimumTrackHits),
  fSubeventCut(in.fSubeventCut),
  fTrackEfficiencyVariation(in.fTrackEfficiencyVariation),
  fTrackSelectionVariation(in.fTrackSelectionVariation),
  fJetPtWeightConfiguration(in.fJetPtWeightConfiguration),
  fReconstructedJetMinimumPtCut(in.fReconstructedJetMinimumPtCut),
  fGeneratorJetMinimumPtCut(in.fGeneratorJetMinimumPtCut),
  fSkipCovarianceMatrix(in.fSkipCovarianceMatrix),
  fMcCorrelationType(in.fMcCorrelationType),
  fJetRadius(in.fJetRadius),
  fWeightExponent(in.fWeightExponent),
  fSmearDeltaR(in.fSmearDeltaR),
  fSmearEnergyWeight(in.fSmearEnergyWeight),
  fDoReflectedCone(in.fDoReflectedCone),
  fDoMixedCone(in.fDoMixedCone),
  fMegaSkimMode(in.fMegaSkimMode),
  fDoPerpendicularCone(in.fDoPerpendicularCone),
  fLocalRun(in.fLocalRun),
  fLocalMixing(in.fLocalMixing),
  fMixingListIndex(in.fMixingListIndex),
  fMixingStartIndex(in.fMixingStartIndex),
  fRunningMixingIndex(in.fRunningMixingIndex),
  fnEventsInMixingFile(in.fnEventsInMixingFile),
  fMixedEventVz(in.fMixedEventVz),
  fMixedEventHiBin(in.fMixedEventHiBin),
  fMixedEventMultiplicity(in.fMixedEventMultiplicity),
  fMixedEventHFEnergy(in.fMixedEventHFEnergy),
  fMixedEventEventNumber(in.fMixedEventEventNumber),
  fMixedEventVzBin(in.fMixedEventVzBin),
  fMixedEventHiBinBin(in.fMixedEventHiBinBin),
  fMixedEventMultiplicityBin(in.fMixedEventMultiplicityBin),
  fMixedEventHFEnergyBin(in.fMixedEventHFEnergyBin),
  fEventMatchForMixing(in.fEventMatchForMixing),
  fBinnedMixing(in.fBinnedMixing),
  fVzMixingBinBorders(in.fVzMixingBinBorders),
  fHiBinMixingBinBorders(in.fHiBinMixingBinBorders),
  fHFEnergyMixingBinBorders(in.fHFEnergyMixingBinBorders),
  fMultiplicityMixingBinBorders(in.fMultiplicityMixingBinBorders),
  fHFEnergyMatchingShift(in.fHFEnergyMatchingShift),
  fFillEventInformation(in.fFillEventInformation),
  fFillJetHistograms(in.fFillJetHistograms),
  fFillTrackHistograms(in.fFillTrackHistograms),
  fFillJetConeHistograms(in.fFillJetConeHistograms),
  fFillEnergyEnergyCorrelators(in.fFillEnergyEnergyCorrelators),
  fFillEnergyEnergyCorrelatorsSystematics(in.fFillEnergyEnergyCorrelatorsSystematics),
  fFillJetPtClosure(in.fFillJetPtClosure),
  fFillJetPtUnfoldingResponse(in.fFillJetPtUnfoldingResponse),
  fFillTrackParticleMatchingHistograms(in.fFillTrackParticleMatchingHistograms),
  fMultiplicityMode(in.fMultiplicityMode)
{
  // Copy constructor

  for(int oBin = 0; oBin < 5; oBin++){
    for(int iBin = 0; iBin < 10; iBin++){
      fThisEventCorrelator[oBin][iBin] = in.fThisEventCorrelator[oBin][iBin];
    }
  }
  
}

/*
 * Assingment operator
 */
EECAnalyzer& EECAnalyzer::operator=(const EECAnalyzer& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fJetReader = in.fJetReader;
  fRecoJetReader = in.fRecoJetReader;
  fTrackReader = in.fTrackReader;
  fMixedEventReader = in.fMixedEventReader;
  fFileNames = in.fFileNames;
  fCard = in.fCard;
  fHistograms = in.fHistograms;
  fVzWeightFunction = in.fVzWeightFunction;
  fCentralityWeightFunctionCentral = in.fCentralityWeightFunctionCentral;
  fCentralityWeightFunctionPeripheral = in.fCentralityWeightFunctionPeripheral;
  fMultiplicityWeightFunction = in.fMultiplicityWeightFunction;
  fPtWeightFunction = in.fPtWeightFunction;
  fSmearingFunction = in.fSmearingFunction;
  fRng = in.fRng;
  fDataType = in.fDataType;
  fTriggerSelection = in.fTriggerSelection;
  fJetType = in.fJetType;
  fMatchJets = in.fMatchJets;
  fDebugLevel = in.fDebugLevel;
  fIsRealData = in.fIsRealData;
  fIsPPbData = in.fIsPPbData;
  fIsPpData = in.fIsPpData;
  fVzWeight = in.fVzWeight;
  fCentralityWeight = in.fCentralityWeight;
  fPtHatWeight = in.fPtHatWeight;
  fTotalEventWeight = in.fTotalEventWeight;
  fJetAxis = in.fJetAxis;
  fVzCut = in.fVzCut;
  fMinimumPtHat = in.fMinimumPtHat;
  fMaximumPtHat = in.fMaximumPtHat;
  fJetEtaCut = in.fJetEtaCut;
  fJetMinimumPtCut = in.fJetMinimumPtCut;
  fJetMaximumPtCut = in.fJetMaximumPtCut;
  fCutBadPhiRegion = in.fCutBadPhiRegion;
  fCutWithJetEtaCM = in.fCutWithJetEtaCM;
  fMinimumMaxTrackPtFraction = in.fMinimumMaxTrackPtFraction;
  fMaximumMaxTrackPtFraction = in.fMaximumMaxTrackPtFraction;
  fJetUncertaintyMode = in.fJetUncertaintyMode;
  fTrackEtaCut = in.fTrackEtaCut;
  fTrackMinPtCut = in.fTrackMinPtCut;
  fTrackMaxPtCut = in.fTrackMaxPtCut;
  fMaxTrackPtRelativeError = in.fMaxTrackPtRelativeError;
  fMaxTrackDistanceToVertex = in.fMaxTrackDistanceToVertex;
  fCalorimeterSignalLimitPt = in.fCalorimeterSignalLimitPt;
  fHighPtEtFraction = in.fHighPtEtFraction;
  fChi2QualityCut = in.fChi2QualityCut;
  fMinimumTrackHits = in.fMinimumTrackHits;
  fSubeventCut = in.fSubeventCut;
  fTrackEfficiencyVariation = in.fTrackEfficiencyVariation;
  fTrackSelectionVariation = in.fTrackSelectionVariation;
  fJetPtWeightConfiguration = in.fJetPtWeightConfiguration;
  fReconstructedJetMinimumPtCut = in.fReconstructedJetMinimumPtCut;
  fGeneratorJetMinimumPtCut = in.fGeneratorJetMinimumPtCut;
  fSkipCovarianceMatrix = in.fSkipCovarianceMatrix;
  fMcCorrelationType = in.fMcCorrelationType;
  fJetRadius = in.fJetRadius;
  fWeightExponent = in.fWeightExponent;
  fSmearDeltaR = in.fSmearDeltaR;
  fSmearEnergyWeight = in.fSmearEnergyWeight;
  fDoReflectedCone = in.fDoReflectedCone;
  fDoMixedCone = in.fDoMixedCone;
  fMegaSkimMode = in.fMegaSkimMode;
  fDoPerpendicularCone = in.fDoPerpendicularCone;
  fLocalRun = in.fLocalRun;
  fLocalMixing = in.fLocalMixing;
  fMixingListIndex = in.fMixingListIndex;
  fMixingStartIndex = in.fMixingStartIndex;
  fRunningMixingIndex = in.fRunningMixingIndex;
  fnEventsInMixingFile = in.fnEventsInMixingFile;
  fMixedEventVz = in.fMixedEventVz;
  fMixedEventHiBin = in.fMixedEventHiBin;
  fMixedEventMultiplicity = in.fMixedEventMultiplicity;
  fMixedEventHFEnergy = in.fMixedEventHFEnergy;
  fMixedEventEventNumber = in.fMixedEventEventNumber;
  fMixedEventVzBin = in.fMixedEventVzBin;
  fMixedEventHiBinBin = in.fMixedEventHiBinBin;
  fMixedEventMultiplicityBin = in.fMixedEventMultiplicityBin;
  fMixedEventHFEnergyBin = in.fMixedEventHFEnergyBin;
  fEventMatchForMixing = in.fEventMatchForMixing;
  fBinnedMixing = in.fBinnedMixing;
  fVzMixingBinBorders = in.fVzMixingBinBorders;
  fHiBinMixingBinBorders = in.fHiBinMixingBinBorders;
  fHFEnergyMixingBinBorders = in.fHFEnergyMixingBinBorders;
  fMultiplicityMixingBinBorders = in.fMultiplicityMixingBinBorders;
  fHFEnergyMatchingShift = in.fHFEnergyMatchingShift;
  fFillEventInformation = in.fFillEventInformation;
  fFillJetHistograms = in.fFillJetHistograms;
  fFillTrackHistograms = in.fFillTrackHistograms;
  fFillJetConeHistograms = in.fFillJetConeHistograms;
  fFillEnergyEnergyCorrelators = in.fFillEnergyEnergyCorrelators;
  fFillEnergyEnergyCorrelatorsSystematics = in.fFillEnergyEnergyCorrelatorsSystematics;
  fFillJetPtClosure = in.fFillJetPtClosure;
  fFillJetPtUnfoldingResponse = in.fFillJetPtUnfoldingResponse;
  fFillTrackParticleMatchingHistograms = in.fFillTrackParticleMatchingHistograms;
  fMultiplicityMode = in.fMultiplicityMode;

  for(int oBin = 0; oBin < 5; oBin++){
    for(int iBin = 0; iBin < 10; iBin++){
      fThisEventCorrelator[oBin][iBin] = in.fThisEventCorrelator[oBin][iBin];
    }
  }
  
  return *this;
}

/*
 * Destructor
 */
EECAnalyzer::~EECAnalyzer(){
  // destructor
  delete fHistograms;
  if(fVzWeightFunction) delete fVzWeightFunction;
  if(fTrackEfficiencyCorrector) delete fTrackEfficiencyCorrector;
  if(fJetCorrector) delete fJetCorrector;
  if(fJetUncertainty) delete fJetUncertainty;
  if(fTrackPairEfficiencyCorrector) delete fTrackPairEfficiencyCorrector;
  if(fEnergyResolutionSmearingFinder) delete fEnergyResolutionSmearingFinder;
  if(fCentralityWeightFunctionCentral) delete fCentralityWeightFunctionCentral;
  if(fCentralityWeightFunctionPeripheral) delete fCentralityWeightFunctionPeripheral;
  if(fMultiplicityWeightFunction) delete fMultiplicityWeightFunction;
  if(fPtWeightFunction) delete fPtWeightFunction;
  if(fSmearingFunction) delete fSmearingFunction;
  if(fRng) delete fRng;
  if(fJetReader) delete fJetReader;
  if(fRecoJetReader) delete fRecoJetReader;
  if(fTrackReader && (fMcCorrelationType == kGenReco || fMcCorrelationType == kRecoGen)) delete fTrackReader;
}

/*
 * Read all the configuration from the input card
 */
void EECAnalyzer::ReadConfigurationFromCard(){
  
  //****************************************
  //     Analyzed data type and trigger
  //****************************************
  fDataType = fCard->Get("DataType");
  fTriggerSelection = fCard->Get("TriggerSelection");

  // Determine the helper data types
  fIsRealData = (fDataType != ForestReader::kPpMC && fDataType != ForestReader::kPbPbMC && fDataType != ForestReader::kPPbMC_pToMinusEta && fDataType != ForestReader::kPPbMC_pToPlusEta && fDataType != ForestReader::kOOMC);
  fIsPPbData = (fDataType == ForestReader::kPPb_pToMinusEta || fDataType == ForestReader::kPPb_pToPlusEta || fDataType == ForestReader::kPPb_pToMinusEta_5TeV || fDataType == ForestReader::kPPbMC_pToMinusEta || fDataType == ForestReader::kPPbMC_pToPlusEta);
  fIsPpData = (fDataType == ForestReader::kPp || fDataType == ForestReader::kPpMC);
  fIsOOData = (fDataType == ForestReader::kOO || fDataType == ForestReader::kOOMC);

  
  //****************************************
  //         Event selection cuts
  //****************************************
  
  fVzCut = fCard->Get("ZVertexCut");          // Event cut vor the z-position of the primary vertex
  fMinimumPtHat = fCard->Get("LowPtHatCut");  // Minimum accepted pT hat value
  fMaximumPtHat = fCard->Get("HighPtHatCut"); // Maximum accepted pT hat value
  
  //****************************************
  //          Jet selection cuts
  //****************************************
  
  fJetEtaCut = fCard->Get("JetEtaCut");           // Eta cut around midrapidity
  fJetMinimumPtCut = fCard->Get("MinJetPtCut");   // Minimum pT cut for jets
  fJetMaximumPtCut = fCard->Get("MaxJetPtCut");   // Maximum pT accepted for jets (and tracks)
  fMinimumMaxTrackPtFraction = fCard->Get("MinMaxTrackPtFraction");  // Cut for jets consisting only from soft particles
  fMaximumMaxTrackPtFraction = fCard->Get("MaxMaxTrackPtFraction");  // Cut for jets consisting only from one high pT particle
  fCutBadPhiRegion = (fCard->Get("CutBadPhi") == 1);   // Flag for cutting the phi region with bad tracking efficiency from the analysis
  fCutWithJetEtaCM = (fCard->Get("CutWithEtaCM") == 1);   // Flag for cutting symmetric region in center-of-mass jet eta that is within defined jet eta acceptence
  fJetUncertaintyMode = fCard->Get("JetUncertainty");  // Select whether to use nominal jet pT or vary it within uncertainties
  
  //****************************************
  //   Extra jet cuts for unfolding study
  //****************************************

  fReconstructedJetMinimumPtCut = fCard->Get("MinJetPtUnfoldingReco");
  fGeneratorJetMinimumPtCut = fCard->Get("MinJetPtUnfoldingTruth");

  //****************************************
  //        Track selection cuts
  //****************************************
  
  fTrackEtaCut = fCard->Get("TrackEtaCut");     // Eta cut around midrapidity
  fTrackMinPtCut = fCard->Get("MinTrackPtCut"); // Minimum track pT cut
  fTrackMaxPtCut = fCard->Get("MaxTrackPtCut"); // Maximum track pT cut
  fMaxTrackPtRelativeError = fCard->Get("MaxTrackPtRelativeError");   // Maximum relative error for pT
  fMaxTrackDistanceToVertex = fCard->Get("VertexMaxDistance");        // Maximum distance to primary vetrex
  fCalorimeterSignalLimitPt = fCard->Get("CalorimeterSignalLimitPt"); // Require signal in calorimeters for track above this pT
  fHighPtEtFraction = fCard->Get("HighPtEtFraction"); // For high pT tracks, minimum required Et as a fraction of track pT
  fChi2QualityCut = fCard->Get("Chi2QualityCut");     // Quality cut for track reconstruction
  fMinimumTrackHits = fCard->Get("MinimumTrackHits"); // Quality cut for track hits
  
  fSubeventCut = fCard->Get("SubeventCut");     // Required subevent type

  //****************************************
  //        Systematic variations
  //****************************************

  fTrackEfficiencyVariation = fCard->Get("TrackEfficiencyVariation"); // Relative amount with which the tracking efficiency corrections are varied to estimate systematic uncertainties

  // Determine the correct efficiency table to use based on track selection cuts
  fTrackSelectionVariation = "";
  if(TMath::Abs(fMaxTrackPtRelativeError-0.05) < 0.001 && TMath::Abs(fMaxTrackDistanceToVertex-2) < 0.001 && TMath::Abs(fChi2QualityCut-0.15) < 0.001) fTrackSelectionVariation = "Tight";
  if(TMath::Abs(fMaxTrackPtRelativeError-0.15) < 0.001 && TMath::Abs(fMaxTrackDistanceToVertex-5) < 0.001) fTrackSelectionVariation = "Loose";

  // Select the jet pT weighting configuration
  fJetPtWeightConfiguration = fCard->Get("JetPtWeight");
  

  //****************************************
  //    Correlation type for Monte Carlo
  //****************************************
  if(fIsRealData) {
    fMcCorrelationType = -100;
  } else {
    fMcCorrelationType = fCard->Get("McCorrelationType");         // Correlation type for Monte Carlo
  }
  fMultiplicityMode = (fCard->Get("MultiplicityMode") == 1);

  //****************************************
  //       Jet selection and matching
  //****************************************
  fJetType = fCard->Get("JetType");              // Select the type of analyzed jets (Calo, CSPF, PuPF, FlowPF)
  fJetAxis = fCard->Get("JetAxis");              // Select between escheme and WTA axes
  fMatchJets = fCard->Get("MatchJets");          // Match flag between reconstructed and generator level jets
  
  //********************************************************************
  //    Turn off certain track cuts for generated tracks, pp, and pPb 
  //********************************************************************
  
  if(fMcCorrelationType == kGenGen || fMcCorrelationType == kRecoGen || fIsPpData || fIsPPbData){
    fCalorimeterSignalLimitPt = 10000;
    fChi2QualityCut = 10000;
    fMinimumTrackHits = 0;
  }
  
  //************************************************
  //   Configuration for energy-energy correlators
  //************************************************
  fJetRadius = fCard->Get("JetRadius");
  fSmearDeltaR = (fCard->Get("SmearDeltaR") == 1);
  fSmearEnergyWeight = (fCard->Get("SmearEnergyWeight") == 1);
  fSkipCovarianceMatrix = (fCard->Get("SkipCovarianceMatrix") == 1);

  const Int_t nWeightExponents = fCard->GetN("WeightExponent");
  fWeightExponent.clear();
  for(Int_t iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    fWeightExponent.push_back(fCard->Get("WeightExponent", iWeightExponent));
  }

  //******************************************************************
  //   Variables for energy-energy correlator background evaluation
  //******************************************************************
  int backgroundMethods = fCard->Get("BackgroundMethods");
  std::bitset<3> backgroundBitChecker(backgroundMethods);
  fDoReflectedCone = backgroundBitChecker.test(0);
  fDoPerpendicularCone = backgroundBitChecker.test(1);
  fDoMixedCone = backgroundBitChecker.test(2);
  fMegaSkimMode = (fCard->Get("MegaSkimMixing") == 1);

  // Configuration for event mixing
  fEventMatchForMixing = fCard->Get("EventMatchForMixing");
  fBinnedMixing = (fCard->Get("BinnedMixing") == 1);

  // If exactly two numbers are defined for mixing bins, those are interpreted as example bin border and bin width
  // Otherwise these are interpreted directly as bin borders
  const Int_t nVzMixingBins = fCard->GetN("VzMixingBins");
  fVzMixingBinBorders.clear();
  for(Int_t iBin = 0; iBin < nVzMixingBins; iBin++){
    fVzMixingBinBorders.push_back(fCard->Get("VzMixingBins", iBin));
  }

  const Int_t nHiBinMixingbins = fCard->GetN("HiBinMixingBins");
  fHiBinMixingBinBorders.clear();
  for(Int_t iBin = 0; iBin < nHiBinMixingbins; iBin++){
    fHiBinMixingBinBorders.push_back(fCard->Get("HiBinMixingBins", iBin));
  }

  const Int_t nMixingHFEnergyBins = fCard->GetN("HFEnergyMixingBins");
  fHFEnergyMixingBinBorders.clear();
  for(Int_t iBin = 0; iBin < nMixingHFEnergyBins; iBin++){
    fHFEnergyMixingBinBorders.push_back(fCard->Get("HFEnergyMixingBins", iBin));
  }

  const Int_t nMultiplicityMixingBins = fCard->GetN("MultiplicityMixingBins");
  fMultiplicityMixingBinBorders.clear();
  for(Int_t iBin = 0; iBin < nMultiplicityMixingBins; iBin++){
    fMultiplicityMixingBinBorders.push_back(fCard->Get("MultiplicityMixingBins", iBin));
  }

  // Shift in HF energy when finding events to mix with
  fHFEnergyMatchingShift = fCard->Get("HFEnergyShift");
  
  //************************************************
  //         Read which histograms to fill
  //************************************************
  int filledHistograms = fCard->Get("FilledHistograms");
  std::bitset<knFillTypes> bitChecker(filledHistograms);
  fFillEventInformation = bitChecker.test(kFillEventInformation);
  fFillJetHistograms = bitChecker.test(kFillJets);
  fFillTrackHistograms = bitChecker.test(kFillTracks);
  fFillJetConeHistograms = bitChecker.test(kFillJetConeHistograms);
  fFillEnergyEnergyCorrelators = bitChecker.test(kFillEnergyEnergyCorrelators);
  fFillEnergyEnergyCorrelatorsSystematics = bitChecker.test(kFillEnergyEnergyCorrelatorsSystematics);
  fFillJetPtClosure = bitChecker.test(kFillJetPtClosure);
  fFillJetPtUnfoldingResponse = bitChecker.test(kFillJetPtUnfoldingResponse);
  fFillTrackParticleMatchingHistograms = bitChecker.test(kFillTrackParticleMatchingHistograms);
  
  //************************************************
  //              Debug messages
  //************************************************
  fDebugLevel = fCard->Get("DebugLevel");
}

/*
 * Main analysis loop
 */
void EECAnalyzer::RunAnalysis(){
  
  //************************************************
  //  Define variables needed in the analysis loop
  //************************************************
  
  // Input files and forest readers for analysis
  TFile* inputFile;
  TFile* copyInputFile;      // If we read forest for tracks and jets with different readers, we need to give different file pointers to them
  TFile* unfoldingInputFile; // If using several readers, avoid problems by using different file pointers for different readers
  
  // Event variables
  Int_t nEvents = 0;                // Number of events
  Double_t vz = 0;                  // Vertex z-position
  Double_t centrality = 0;          // Event centrality
  Int_t hiBin = 0;                  // CMS hiBin (centrality * 2)
  Double_t ptHat = 0;               // pT hat for MC events
  
  // Combining bools to make the code more readable
  Bool_t useDifferentReaderForJetsAndTracks = (fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenReco); // Use different forest reader for jets and tracks
  
  // Variables for trigger flags
  Bool_t caloJet15Trigger;
  Bool_t caloJet30Trigger;
  Bool_t caloJet40Trigger;
  Bool_t caloJet60Trigger;
  Bool_t caloJet80Trigger;
  Bool_t caloJet100Trigger;
  
  // Variables for jets
  Double_t jetPt = 0;               // pT of the i:th jet in the event
  Double_t jetPtCorrected = 0;      // Jet pT corrected with the JFF correction
  Double_t jetPhi = 0;              // phi of the i:th jet in the event
  Double_t jetEta = 0;              // eta of the i:th jet in the event
  Double_t jetEtaCM = 0;            // eta of the i:th jet in the center-of-mass frame
  Int_t jetFlavor = 0;              // Flavor of the jet. 0 = Quark jet. 1 = Gluon jet.
  Double_t jetPtWeight = 1;         // Weighting for jet pT
  Bool_t jetOver80Found = false;    // Flag for finding a jet above 80 GeV from the event
  Bool_t jetOver120Found = false;   // Flag for finding a jet above 120 GeV from the event
  Double_t leadingJetPt = 0;        // Leading jet pT in an event
  Int_t nGoodJets = 0;              // Nuber of good jets in the event

  // Variables for reflected cone QA study
  Double_t reflectedConeJetPt = 0;   // pT of a jet that might be in the reflected cone
  Double_t reflectedConeJetPhi = 0;  // phi of a jet that might be in the reflected cone
  Double_t reflectedConeJetEta = 0;  // eta of a jet that might be in the reflected cone
  Int_t nJetsInReflectedCone = 0;    // number of jets found inside the reflected cone
  
  // Variables for smearing study
  Double_t smearingFactor = 0;       // Larger of the JEC uncertainties
  
  // Variables for jet matching and closure
  Int_t unmatchedCounter = 0;       // Number of jets that fail the matching
  Int_t matchedCounter = 0;         // Number of jets that are matched
  Bool_t antiMatchPtVeto = false;   // Reject matching if the pT is too far off
  Int_t partonFlavor = -999;        // Code for parton flavor in Monte Carlo
  Int_t matchRejectedDueToPtCounter = 0; // Counter for too large pT difference between matched jets
  
  // Variables for tracks
  Double_t fillerTrack[4];                // Track histogram filler
  Double_t trackEfficiencyCorrection;     // Track efficiency correction
  Int_t nTracks;                          // Number of tracks in an event
  Double_t trackPt = 0;                   // Track pT
  Double_t trackEta = 0;                  // Track eta
  Double_t trackPhi = 0;                  // Track phi
  Double_t trackMultiplicity = 0;         // Multiplicity
  Double_t trackMultiplicityWeighted = 0; // Weighted multiplicity
  Int_t trackSubevent = 0;                // Subevent index in Pythia+Hydjet simulation
  Int_t trackSubeventIndex = 0;           // Simplified subevent index
  Double_t maxTrackPtInJetSignal = 0;     // Maximum signal track pT in the jet
  Double_t maxTrackPtInJetBackground = 0; // Maximum background track pT in the jet

  // Variables for covariance study
  const Int_t nDeltaRBins = fHistograms->GetNDeltaRBinsEEC(); // Number of DeltaR bins
  Double_t eecCovariance = 0;             // Covariance between two deltaR bins in energy-energy correlators
  Double_t deltaRBinContent1 = 0;         // Bin content of the i:th DeltaR bin
  Double_t deltaRBinContent2 = 0;         // Bin content of the j:th DeltaR bin
  Double_t transformedDeltaR1 = 0;        // Transformed bin contents of the i:th DeltaR bin
  Double_t transformedDeltaR2 = 0;        // Transformed bin contents of the j:th DeltaR bin
  
  // Study for track multiplicity inside the jet cone
  const Int_t nTrackPtBinsEEC = fCard->GetNBin("TrackPtBinEdgesEEC");
  Double_t* trackPtBinsEEC = new double[nTrackPtBinsEEC];
  for(Int_t iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
    trackPtBinsEEC[iTrackPt] = fCard->Get("TrackPtBinEdgesEEC",iTrackPt);
  }

  // Mapping a bin index to multiplicity number
  std::map<std::tuple<Int_t,Int_t,Int_t>, Int_t> uncorrectedMultiplicityInJetCone;
  std::map<std::tuple<Int_t,Int_t,Int_t>, Int_t> multiplicityInJetCone;
  std::tuple<Int_t,Int_t,Int_t> multiplicityBin;

  //Int_t uncorrectedMultiplicityInJetCone[nTrackPtBinsEEC][EECHistograms::knJetConeTypes][EECHistograms::knSubeventTypes+1]; // Particle multiplicity within a jet cone, no tracking efficiency correction
  //Double_t multiplicityInJetCone[nTrackPtBinsEEC][EECHistograms::knJetConeTypes][EECHistograms::knSubeventTypes+1];;        // Efficiency corrected particle multiplicity within a jet cone

  // Event mixing variables
  std::vector<TString> mixingFiles;      // List of mixing files
  Int_t currentMultiplicity;             // Multiplicity in the current event
  Int_t mixedEventFillIndex;             // Index used to fill mixed event track information into vectors
  std::vector<Int_t> mixedEventIndices;  // Vector for holding the events we have already mixed with
  Int_t eventNumber;                     // Event number
  Double_t hfSum;                        // Energy sum in HF calorimeters

  // HF calorimeters
  Double_t hfPlus;                       // Energy sum in HF plus calorimeters
  Double_t hfMinus;                      // Energy sum in HF minus calorimeters
  
  // Variables for energy-energy correlators
  vector<double> selectedTrackPt[EECHistograms::knJetConeTypes];     // Track pT for tracks selected for energy-energy correlator analysis
  vector<double> relativeTrackEta[EECHistograms::knJetConeTypes];    // Track eta relative to the jet axis
  vector<double> relativeTrackPhi[EECHistograms::knJetConeTypes];    // Track phi relative to the jet axis
  vector<int> selectedTrackSubevent[EECHistograms::knJetConeTypes];  // Track subevent for tracks selected for energy-energy correlator analysis
  Double_t jetReflectedEta = 0;          // Reflected jet eta to be used for background estimation
  Double_t jetPerpendicularPhi = 0;      // Perpendicular jet phi to be used for background estimation
  Double_t deltaRTrackJet = 0;           // DeltaR between tracks and jet axis
  
  // File name helper variables
  TString currentFile;
  TFile *testFile;
  
  // Histograms that should be filled
  Int_t filledHistograms = fCard->Get("FilledHistograms");
  Bool_t readTrackTree = true;
  if(filledHistograms < 4) {
    useDifferentReaderForJetsAndTracks = false;
    readTrackTree = false;  // Do not read track trees if only jets are used. Makes things much faster in this case.
  }
  
  // Fillers for THnSparses
  const Int_t nFillJet = 6;
  const Int_t nFillMultiplicity = 3;
  const Int_t nFillMultiplicityInJetCone = 5;
  const Int_t nFillParticleDensityInJetCone = 6;
  const Int_t nFillMaxParticlePtInJetCone = 4;
  const Int_t nFillReflectedConeQA = 2;
  const Int_t nFillUnfoldingCovariance = 5;
  const Int_t nFillHFEnergy = 5;
  Double_t fillerJet[nFillJet];
  Double_t fillerMultiplicity[nFillMultiplicity];
  Double_t fillerMultiplicityInJetCone[nFillMultiplicityInJetCone];
  Double_t fillerParticleDensityInJetCone[nFillParticleDensityInJetCone];
  Double_t fillerMaxParticlePtInJetCone[nFillMaxParticlePtInJetCone];
  Double_t fillerReflectedConeQA[nFillReflectedConeQA];
  Double_t fillerUnfoldingCovariance[nFillUnfoldingCovariance];
  Double_t fillerHFEnergy[nFillHFEnergy];
  
  // Setup the jet energy correction files based on the analyzed data type
  std::string correctionFileRelative[ForestReader::knDataTypes];
  correctionFileRelative[ForestReader::kPp] = "jetEnergyCorrections/Spring18_ppRef5TeV_V6_DATA_L2Relative_AK4PF.txt";
  correctionFileRelative[ForestReader::kPbPb] = "jetEnergyCorrections/Autumn18_HI_V8_DATA_L2Relative_AK4PF.txt";
  correctionFileRelative[ForestReader::kPpMC] = "jetEnergyCorrections/Spring18_ppRef5TeV_V6_MC_L2Relative_AK4PF.txt";
  correctionFileRelative[ForestReader::kPbPbMC] = "jetEnergyCorrections/Autumn18_HI_V8_MC_L2Relative_AK4PF.txt";
  correctionFileRelative[ForestReader::kPPb_pToMinusEta] = "jetEnergyCorrections/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt";
  correctionFileRelative[ForestReader::kPPb_pToPlusEta] = "jetEnergyCorrections/Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt";
  correctionFileRelative[ForestReader::kPPb_pToMinusEta_5TeV] = "jetEnergyCorrections/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt";
  correctionFileRelative[ForestReader::kPPbMC_pToMinusEta] = "jetEnergyCorrections/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt";
  correctionFileRelative[ForestReader::kPPbMC_pToPlusEta] = "jetEnergyCorrections/Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt";
  correctionFileRelative[ForestReader::kOO] = "";
  correctionFileRelative[ForestReader::kOOMC] = "";

  std::string correctionFileResidual[ForestReader::knDataTypes];
  correctionFileResidual[ForestReader::kPp] = "jetEnergyCorrections/Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt";
  correctionFileResidual[ForestReader::kPbPb] = "jetEnergyCorrections/Autumn18_HI_V8_DATA_L2L3Residual_AK4PF.txt";
  correctionFileResidual[ForestReader::kPpMC] = "CorrectionNotAppliedPF.txt";
  correctionFileResidual[ForestReader::kPbPbMC] = "CorrectionNotAppliedPF.txt";
  correctionFileResidual[ForestReader::kPPb_pToMinusEta] = "jetEnergyCorrections/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt";
  correctionFileResidual[ForestReader::kPPb_pToPlusEta] = "jetEnergyCorrections/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt";
  correctionFileResidual[ForestReader::kPPb_pToMinusEta_5TeV] = "jetEnergyCorrections/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt";
  correctionFileResidual[ForestReader::kPPbMC_pToMinusEta] = "CorrectionNotAppliedPF.txt";
  correctionFileResidual[ForestReader::kPPbMC_pToPlusEta] = "CorrectionNotAppliedPF.txt";
  correctionFileResidual[ForestReader::kOO] = "CorrectionNotAppliedPF.txt"; // Use the dummy string
  correctionFileResidual[ForestReader::kOOMC] = "CorrectionNotAppliedPF.txt";

  std::string uncertaintyFile[ForestReader::knDataTypes];
  uncertaintyFile[ForestReader::kPp] = "jetEnergyCorrections/Spring18_ppRef5TeV_V6_DATA_Uncertainty_AK4PF.txt";
  uncertaintyFile[ForestReader::kPbPb] = "jetEnergyCorrections/Autumn18_HI_V8_DATA_Uncertainty_AK4PF.txt";
  uncertaintyFile[ForestReader::kPpMC] = "jetEnergyCorrections/Spring18_ppRef5TeV_V6_MC_Uncertainty_AK4PF.txt";
  uncertaintyFile[ForestReader::kPbPbMC] = "jetEnergyCorrections/Autumn18_HI_V8_MC_Uncertainty_AK4PF.txt";
  uncertaintyFile[ForestReader::kPPb_pToMinusEta] = "jetEnergyCorrections/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF_modifiedtopPb.txt";
  uncertaintyFile[ForestReader::kPPb_pToPlusEta] = "jetEnergyCorrections/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF_modifiedtopPb.txt";
  uncertaintyFile[ForestReader::kPPb_pToMinusEta_5TeV] = "jetEnergyCorrections/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF_modifiedtopPb.txt";
  uncertaintyFile[ForestReader::kPPbMC_pToMinusEta] = "jetEnergyCorrections/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF_modifiedtopPb.txt";
  uncertaintyFile[ForestReader::kPPbMC_pToPlusEta] = "jetEnergyCorrections/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF_modifiedtopPb.txt";
  uncertaintyFile[ForestReader::kOO] = "";
  uncertaintyFile[ForestReader::kOOMC] = "";
  
  // For calo jets, use the correction files for calo jets (otherwise same name, but replace PF with Calo)
  if(fJetType == 0 && !fIsOOData){
    size_t pfIndex = 0;
    pfIndex = correctionFileRelative[fDataType].find("PF", pfIndex);
    correctionFileRelative[fDataType].replace(pfIndex, 2, "Calo");
    pfIndex = 0;
    pfIndex = correctionFileResidual[fDataType].find("PF", pfIndex);
    correctionFileResidual[fDataType].replace(pfIndex, 2, "Calo");
    pfIndex = 0;
    pfIndex = uncertaintyFile[fDataType].find("PF", pfIndex);
    uncertaintyFile[fDataType].replace(pfIndex, 2, "Calo");
    
  }
  
  vector<string> correctionFiles;
  correctionFiles.push_back(correctionFileRelative[fDataType]);
  if(fIsRealData)  correctionFiles.push_back(correctionFileResidual[fDataType]);
  
  fJetCorrector = new JetCorrector(correctionFiles);
  fJetUncertainty = new JetUncertainty(uncertaintyFile[fDataType]);
  
  //************************************************
  //      Find forest readers for data files
  //************************************************
  
  if(fMcCorrelationType == kGenReco || fMcCorrelationType == kGenGen){
    fJetReader = new GeneratorLevelForestReader(fDataType,fTriggerSelection,fJetType,fJetAxis,fMatchJets,readTrackTree);
  } else {
    fJetReader = new HighForestReader(fDataType,fTriggerSelection,fJetType,fJetAxis,fMatchJets,readTrackTree);
  }
  
  // Select the reader for tracks based on forest and MC correlation type
  if(fMcCorrelationType == kRecoGen){
    fTrackReader = new GeneratorLevelForestReader(fDataType,fTriggerSelection,fJetType,fJetAxis,fMatchJets);
  } else if (fMcCorrelationType == kGenReco){
    fTrackReader = new HighForestReader(fDataType,fTriggerSelection,fJetType,fJetAxis,fMatchJets);
  } else {
    fTrackReader = fJetReader;
  }

  // Select the reader for mixed events
  if(fDoMixedCone){
    if(fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenGen){
      fMixedEventReader = new GeneratorLevelForestReader(fDataType,fTriggerSelection,fJetType,fJetAxis,false,true,true,fMegaSkimMode);
    } else {
      fMixedEventReader = new HighForestReader(fDataType,fTriggerSelection,fJetType,fJetAxis,false,true,true,fMegaSkimMode);
    }
  }

  if(fFillJetPtUnfoldingResponse || fFillTrackParticleMatchingHistograms){
    fUnfoldingForestReader = new UnfoldingForestReader(fDataType,fJetType,fJetAxis);
  }

  // Prepare mixed event files and a forest reader for reflected cone mixing
  const char* fileListName[2][ForestReader::knDataTypes][2];

  // Transform the CRAB job index into a accepted range of mixing job indices
  if(fMixingListIndex < 0){
    fMixingListIndex = fMixingListIndex * -1;
    if((fDataType == ForestReader::kPbPbMC && !fMegaSkimMode) || fDataType == ForestReader::kPPb_pToMinusEta){
      // There are four copies of non-mega skimmed mixing file list for PbPb MC
      // There is also four parts for the megaskimmed files for pPb p -> -eta at 8 TeV
      fMixingListIndex = (fMixingListIndex % 4) + 1;
    } else if (ForestReader::kPPb_pToMinusEta) {
      // There is seven parts for the megaskimmed files for pPb p -> +eta at 8 TeV
      fMixingListIndex = (fMixingListIndex % 7) + 1;
    } else {
      // There are eight copies all other mixing file lists
      fMixingListIndex = (fMixingListIndex % 8) + 1;
    }
  }

  // CRAB running, regular mixing forest
  fileListName[0][ForestReader::kPp][0] = "none";  // only mega skimmed mixing available for pp
  fileListName[0][ForestReader::kPbPb][0] = Form("mixingFileList/PbPbData2018_minBiasFiles_copy%d.txt", fMixingListIndex); // PbPb data for CRAB
  fileListName[0][ForestReader::kPpMC][0] = "none";  // only mega skimmed mixing available for pp MC
  fileListName[0][ForestReader::kPbPbMC][0] = Form("mixingFileList/PbPbMC2018_minBiasHydjetFiles_copy%d.txt", fMixingListIndex); // PbPb MC for CRAB
  fileListName[0][ForestReader::kPPb_pToMinusEta][0] = "none"; // only mega skimmed mixing available for pPb p -> -eta
  fileListName[0][ForestReader::kPPb_pToPlusEta][0] = "none"; // only mega skimmed mixing available for pPb p -> +eta
  fileListName[0][ForestReader::kPPb_pToMinusEta_5TeV][0] = "none"; // only mega skimmed mixing available for pPb p -> -eta 5 TeV
  fileListName[0][ForestReader::kPPbMC_pToMinusEta][0] = "none"; // currently no mixing implemented for pPb MC
  fileListName[0][ForestReader::kPPbMC_pToPlusEta][0] = "none"; // currently no mixing implemented for pPb MC
  fileListName[0][ForestReader::kOO][0] = "none"; // currently no mixing implemented for OO
  fileListName[0][ForestReader::kOOMC][0] = "none"; // currently no mixing implemented for OOMC


  // Local test, regular mixing forest
  fileListName[1][ForestReader::kPp][0] = "none";  // only mega skimmed mixing available for pp
  fileListName[1][ForestReader::kPbPb][0] = "mixingFileList/mixingFilesPbPb.txt"; // PbPb data for local test
  fileListName[1][ForestReader::kPpMC][0] = "none";  // only mega skimmed mixing available for pp MC
  fileListName[1][ForestReader::kPbPbMC][0] = "mixingFileList/mixingFilesPbPbMC.txt";  // PbPb MC for local test
  fileListName[1][ForestReader::kPPb_pToMinusEta][0] = "none"; // only mega skimmed mixing available for pPb p -> -eta
  fileListName[1][ForestReader::kPPb_pToPlusEta][0] = "none"; // only mega skimmed mixing available for pPb p -> +eta
  fileListName[1][ForestReader::kPPb_pToMinusEta_5TeV][0] = "none"; // only mega skimmed mixing available for pPb p -> -eta 5 TeV
  fileListName[1][ForestReader::kPPbMC_pToMinusEta][0] = "none"; // currently no mixing implemented for pPb MC
  fileListName[1][ForestReader::kPPbMC_pToPlusEta][0] = "mixingFileList/mixingFilesPPb_pToPlusEta.txt"; // low statistics test mixing file for pPb MC p -> + eta
  fileListName[1][ForestReader::kOO][0] = "none"; // currently no mixing implemented for OO
  fileListName[1][ForestReader::kOOMC][0] = "none"; // currently no mixing implemented for OOMC

  // CRAB running, mega skimmed mixing forest
  fileListName[0][ForestReader::kPp][1] = "mixingFileList/zeroBiasPp2017_5TeV_megaSkim_2025-06-24.txt"; // pp data for CRAB
  fileListName[0][ForestReader::kPbPb][1] = Form("mixingFileList/PbPbData2018_minBiasMegaSkims_copy%d.txt", fMixingListIndex); // PbPb data for CRAB
  fileListName[0][ForestReader::kPpMC][1] = "none";  // no mixing currently implemented for pp MC
  fileListName[0][ForestReader::kPbPbMC][1] = Form("mixingFileList/PbPbMC2018_minBiasHydjetMegaSkims_copy%d.txt", fMixingListIndex); // PbPb MC for CRAB
  fileListName[0][ForestReader::kPPb_pToMinusEta][1] = Form("mixingFileList/minimumBiasPPb2016_8TeV_megaSkim_pToMinusEta_part%d.txt", fMixingListIndex); // pPb p -> -eta for CRAB
  fileListName[0][ForestReader::kPPb_pToPlusEta][1] = Form("mixingFileList/minimumBiasPPb2016_8TeV_megaSkim_pToPlusEta_part%d.txt", fMixingListIndex); // pPb p -> +eta for CRAB
  fileListName[0][ForestReader::kPPb_pToMinusEta_5TeV][1] = "mixingFileList/minimumBiasPPb2016_5TeV_megaSkim_pToMinusEta.txt"; // pPb pToMinusEta 5 TeV for CRAB
  fileListName[0][ForestReader::kPPbMC_pToMinusEta][1] = "none"; // currently no mixing implemented for pPb MC
  fileListName[0][ForestReader::kPPbMC_pToPlusEta][1] = "mixingFileList/EPOS_pToPlusEta_pPb816Summer16DR_megaSkim_2026-02-13.txt"; // currently no mixing implemented for pPb MC
  fileListName[0][ForestReader::kOO][0] = "none"; // currently no mixing implemented for OO
  fileListName[0][ForestReader::kOOMC][0] = "none"; // currently no mixing implemented for OOMC


  // Local test, mega skimmed mixing forest
  fileListName[1][ForestReader::kPp][1] = "mixingFileList/mixingFilesPp_zeroBias_megaSkim.txt";  // pp data for local test
  fileListName[1][ForestReader::kPbPb][1] = "mixingFileList/mixingFilesPbPb_megaSkim.txt"; // PbPb data for local test
  fileListName[1][ForestReader::kPpMC][1] = "none";  // no mixing currently implemented for pp MC
  fileListName[1][ForestReader::kPbPbMC][1] = "mixingFileList/mixingFilesPbPbMC_megaSkim.txt";  // PbPb MC for local test
  fileListName[1][ForestReader::kPPb_pToMinusEta][1] = "mixingFileList/mixingFilesPPb_pToMinusEta_megaSkim.txt"; // pPb p -> -eta for local test
  fileListName[1][ForestReader::kPPb_pToPlusEta][1] = "mixingFileList/mixingFilesPPb_pToPlusEta_megaSkim.txt";  // pPb p -> +eta for local test
  fileListName[1][ForestReader::kPPb_pToMinusEta_5TeV][1] = "mixingFileList/mixingFilesPPb_pToMinusEta_megaSkim.txt"; // pPb p -> -eta 5 TeV for local test
  fileListName[1][ForestReader::kPPbMC_pToMinusEta][1] = "none"; // currently no mixing implemented for pPb MC
  fileListName[1][ForestReader::kPPbMC_pToPlusEta][1] = "mixingFileList/mixingFilesEpos_pToPlusEta_megaSkim.txt"; // pPb MC p -> +eta
  fileListName[1][ForestReader::kOO][0] = "none"; // currently no mixing implemented for OO
  fileListName[1][ForestReader::kOOMC][0] = "none"; // currently no mixing implemented for OOMC

        
  // Read the mixing files if defined and a file list exists
  if(fDoMixedCone){
    if(std::string(fileListName[fLocalRun][fDataType][fMegaSkimMode]) == "none"){

      cout << "EECAnalyzer::Note: You are trying to do mixed cone background estimation for a dataset for which no mixed event file list is defined in the code! The input file list is used also for mixing." << endl;
      mixingFiles = fFileNames;
    
    } else {

      std::string lineInFile;
      std::ifstream mixingFileStream(fileListName[fLocalRun][fDataType][fMegaSkimMode]);
      while (std::getline(mixingFileStream,lineInFile)) {

        // For slurm running, remove xrootd pointer from mixing lists
        if(fLocalMixing){
          // Find the first occurrence of "//"
          size_t firstDoubleSlash = lineInFile.find("//");
          if(firstDoubleSlash != std::string::npos){
          // Find the second occurrence of "//"
            size_t secondDoubleSlash = lineInFile.find("//", firstDoubleSlash + 2);
            if(secondDoubleSlash != std::string::npos){
              // Remove the part with xrootd pointer from mixing files
              mixingFiles.push_back(lineInFile.substr(secondDoubleSlash + 1));
            }
          }
        } else{
          // If not running on slurm, we can directly use whatever is in the file lists
          mixingFiles.push_back(lineInFile);
        }
      }
    } // Mixing file list exists

    // Check that all the mixing files we are supposed to read exist and can be opened
    for(auto currentFile : mixingFiles){
      testFile = TFile::Open(currentFile);

      // Check that the file exists
      if(!testFile){
        cout << "Error! Could not find the file: " << currentFile.Data() << endl;
        assert(0);
      }

      // Check that the file is open
      if(!testFile->IsOpen()){
        cout << "Error! Could not open the file: " << currentFile.Data() << endl;
        assert(0);
      }
    
      // Check that the file is not zombie
      if(testFile->IsZombie()){
        cout << "Error! The following file is a zombie: " << currentFile.Data() << endl;
        assert(0);
      }

      // Close the test file after the test 
      testFile->Close();
    }

    // If everything is fine, we can read the mixing information from the minimum bias files
    fMixedEventReader->ReadForestFromFileList(mixingFiles);

    // Print the used mixing files
    if(fDebugLevel > 0){
      cout << "Mixing files: " << endl;
      for(std::vector<TString>::iterator mixIterator = mixingFiles.begin(); mixIterator != mixingFiles.end(); mixIterator++){
        cout << *mixIterator << endl;
      }
      cout << endl;
    }

    // Scourge the mixed events to find centrality and vz values that are matched for the mixed events
    fnEventsInMixingFile = fMixedEventReader->GetNEvents();  

    if(fBinnedMixing){

      // Find the event matching variables from mixed event list
      PrepareBinnedMixingVectors();

      // Make a two dimensional illustration of the mixed event statistics
      std::vector<Double_t> vzMixingCountVector = GetBinnedMixingBinBorders(fVzMixingBinBorders, fMixedEventVzBin);
      std::vector<Double_t> eventActivityCountVector;
      std::vector<Double_t> eventActivityBinBorders;

      if(fEventMatchForMixing == kMatchHiBin){
        eventActivityCountVector = GetBinnedMixingBinBorders(fHiBinMixingBinBorders, fMixedEventHiBinBin);
        eventActivityBinBorders = fHiBinMixingBinBorders;
      } else if (fEventMatchForMixing == kMatchHFEnergy){
        eventActivityCountVector = GetBinnedMixingBinBorders(fHFEnergyMixingBinBorders, fMixedEventHFEnergyBin);
        eventActivityBinBorders = fHFEnergyMixingBinBorders;
      } else {
        eventActivityCountVector = GetBinnedMixingBinBorders(fMultiplicityMixingBinBorders, fMixedEventMultiplicityBin);
        eventActivityBinBorders = fMultiplicityMixingBinBorders;
      }
      

      fHistograms->fhMixedEventCounts->SetBins(vzMixingCountVector.size() - 1, vzMixingCountVector.data(), eventActivityCountVector.size() - 1, eventActivityCountVector.data());

      for(int iMixedEvent = 0; iMixedEvent < fMixedEventVzBin.size(); iMixedEvent++){
        fHistograms->fhMixedEventCounts->Fill(FindExampleValue(fMixedEventVzBin.at(iMixedEvent), fVzMixingBinBorders), FindExampleValue(eventActivityCountVector.at(iMixedEvent), eventActivityBinBorders));
      }

    } else {

      // Find the event matching variables from mixed event list
      PrepareMixingVectors();

      // Make a two dimensional illustration of the mixed event statistics
      std::vector<Double_t> vzMixingCountVector = GetMixingBinBorders(fVzMixingBinBorders, fMixedEventVz);
      std::vector<Double_t> eventActivityCountVector;
      std::vector<Double_t> eventActivityVector;

      if(fEventMatchForMixing == kMatchHiBin){
        for(auto number : fMixedEventHiBin) eventActivityVector.push_back(number);
        eventActivityCountVector = GetMixingBinBorders(fHiBinMixingBinBorders, eventActivityVector);
      } else if (fEventMatchForMixing == kMatchHFEnergy){
        eventActivityCountVector = GetMixingBinBorders(fHFEnergyMixingBinBorders, fMixedEventHFEnergy);
        eventActivityVector = fMixedEventHFEnergy;
      } else {
        for(auto number : fMixedEventMultiplicity) eventActivityVector.push_back(number);
        eventActivityCountVector = GetMixingBinBorders(fMultiplicityMixingBinBorders, eventActivityVector);
      }

      fHistograms->fhMixedEventCounts->SetBins(vzMixingCountVector.size() - 1, vzMixingCountVector.data(), eventActivityCountVector.size() - 1, eventActivityCountVector.data());

      for(int iMixedEvent = 0; iMixedEvent < fMixedEventVz.size(); iMixedEvent++){
        fHistograms->fhMixedEventCounts->Fill(fMixedEventVz.at(iMixedEvent), eventActivityVector.at(iMixedEvent));
      }
    }

  } // Prepare event mixing for PbPb data and MC 
  
  //************************************************
  //       Main analysis loop over all files
  //************************************************
  
  // Loop over files
  Int_t nFiles = fFileNames.size();
  for(Int_t iFile = 0; iFile < nFiles; iFile++) {
    
    //************************************************
    //              Find and open files
    //************************************************
    
    // Find the filename and open the input file
    currentFile = fFileNames.at(iFile);
    inputFile = TFile::Open(currentFile);
    if(useDifferentReaderForJetsAndTracks) copyInputFile = TFile::Open(currentFile);
    if(fFillJetPtUnfoldingResponse || fFillTrackParticleMatchingHistograms) unfoldingInputFile = TFile::Open(currentFile);
    
    // Check that the file exists
    if(!inputFile){
      cout << "Error! Could not find the file: " << currentFile.Data() << endl;
      assert(0);
    }

    // Check that the file is open
    if(!inputFile->IsOpen()){
      cout << "Error! Could not open the file: " << currentFile.Data() << endl;
      assert(0);
    }
    
    // Check that the file is not zombie
    if(inputFile->IsZombie()){
      cout << "Error! The following file is a zombie: " << currentFile.Data() << endl;
      assert(0);
    }
    

    // Print the used files
    if(fDebugLevel > 0) cout << "Reading from file: " << currentFile.Data() << endl;
    
    //************************************************
    //            Read forest from file
    //************************************************
    
    // If file is good, read the forest from the file
    fJetReader->ReadForestFromFile(inputFile);  // There might be a memory leak in handling the forest...
    if(useDifferentReaderForJetsAndTracks) fTrackReader->ReadForestFromFile(copyInputFile); // If we mix reco and gen, the reader for jets and tracks is different
    if(fFillJetPtUnfoldingResponse || fFillTrackParticleMatchingHistograms) fUnfoldingForestReader->ReadForestFromFile(unfoldingInputFile);
    nEvents = fJetReader->GetNEvents();


    //**********************************************************************************************
    //    Special case for pPb MC. We need to limit ptHat range to those in the generated files
    //**********************************************************************************************

    // TODO: Figure out proper event weighting once the new MC is ready
    /*if(fDataType == ForestReader::kPPbMC_pToMinusEta || fDataType == ForestReader::kPPbMC_pToPlusEta){

      // The bin borders for pThat for which pPb MC is generated
      std::vector<std::pair<int,double>> ptHatBinBordersAndWeights;
      ptHatBinBordersAndWeights.push_back(std::make_pair(15, 1.0404701e-06 * 961104));
      ptHatBinBordersAndWeights.push_back(std::make_pair(30, 7.7966624e-08 * 952110));
      ptHatBinBordersAndWeights.push_back(std::make_pair(50, 1.0016052e-08 * 952554));
      ptHatBinBordersAndWeights.push_back(std::make_pair(80, 1.3018269e-09 * 996844));
      ptHatBinBordersAndWeights.push_back(std::make_pair(120, 2.2648493e-10 * 964681));
      ptHatBinBordersAndWeights.push_back(std::make_pair(170, 4.0879112e-11 * 999260));
      ptHatBinBordersAndWeights.push_back(std::make_pair(220, 1.1898939e-11 * 964336));
      ptHatBinBordersAndWeights.push_back(std::make_pair(280, 3.3364433e-12 * 995036));
      ptHatBinBordersAndWeights.push_back(std::make_pair(370, 7.6612402e-13 * 958160));
      ptHatBinBordersAndWeights.push_back(std::make_pair(460, 2.1341026e-13 * 981427));
      ptHatBinBordersAndWeights.push_back(std::make_pair(540, 7.9191586e-14 * 1000000));
      ptHatBinBordersAndWeights.push_back(std::make_pair(8160, 0));
      
      // The expected pThat range is defined in the file name
      for(int iPtHat = 0; iPtHat < ptHatBinBordersAndWeights.size()-1; iPtHat++){
        if(fFileNames.at(iFile).Contains(Form("pthat%d", ptHatBinBordersAndWeights.at(iPtHat).first))){
          fMinimumPtHat = ptHatBinBordersAndWeights.at(iPtHat).first;
          fMaximumPtHat = ptHatBinBordersAndWeights.at(iPtHat+1).first;
          fPtHatWeight = ptHatBinBordersAndWeights.at(iPtHat).second / nEvents;
          break;
        }
      }

    }*/

    //************************************************
    //         Main event loop for each file
    //************************************************
    
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){ // nEvents

      // For each event, chack that the file stays open:
      // This is to try to combat file read errors occasionally happening during CRAB running.
      // Will need to monitor the situation and see if this really works.
      if(!inputFile->IsOpen() || inputFile->IsZombie()){
        cout << "Error! Lost access to the file: " << currentFile.Data() << endl;
        assert(0);
      }

      // Do the same check for all mixing files if we are mixing
      //if(fDoMixedCone){
      //  if(fMixedEventReader->CheckFileProblems()){
      //    cout << "Error! Lost access to mixing files!" << endl;
      //    assert(0);
      //  }
      //}

      // Reset the flags that show if high pT jets are found in the events
      jetOver80Found = false;
      jetOver120Found = false;
      
      //************************************************
      //         Read basic event information
      //************************************************
      
      // Print to console how the analysis is progressing
      if(fDebugLevel > 1 && iEvent % 1000 == 0) cout << "Analyzing event " << iEvent << endl;

      // Read the event to memory
      fJetReader->GetEvent(iEvent);
      
      // If track reader is not the same as jet reader, read the event to memory in trackReader
      if(useDifferentReaderForJetsAndTracks) fTrackReader->GetEvent(iEvent);
      if(fFillJetPtUnfoldingResponse || fFillTrackParticleMatchingHistograms) fUnfoldingForestReader->GetEvent(iEvent);

      // Get vz, centrality and pT hat information
      vz = fJetReader->GetVz();
      centrality = fJetReader->GetCentrality();
      hiBin = fJetReader->GetHiBin();
      ptHat = fJetReader->GetPtHat();
      eventNumber = fJetReader->GetEventNumber();
      hfSum = fJetReader->GetHFSum();
      hfPlus = fJetReader->GetHFPlus();
      hfMinus = fJetReader->GetHFMinus();
      
      // We need to apply pT hat cuts before getting pT hat weight. There might be rare events above the upper
      // limit from which the weights are calculated, which could cause the code to crash.
      if(ptHat < fMinimumPtHat || ptHat >= fMaximumPtHat) continue;
      
      // Get the weighting for the event
      fVzWeight = GetVzWeight(vz);
      if(fMultiplicityMode){
        // Multiplicity based weight
        trackMultiplicity = GetMultiplicity();
        fCentralityWeight = GetMultiplicityWeight(trackMultiplicity);
        centrality = GetCentralityFromMultiplicity(trackMultiplicity);
      } else {
        // Regular centrality based weight
        fCentralityWeight = GetCentralityWeight(hiBin);
      }

      // Event weight for 2018 MC
      if(!(fDataType == ForestReader::kPPbMC_pToMinusEta) && !(fDataType == ForestReader::kPPbMC_pToPlusEta)){

        // Event-by-event weight for pp and PbPb MC
        fPtHatWeight = fJetReader->GetEventWeight();
      } 
      fTotalEventWeight = fVzWeight*fCentralityWeight*fPtHatWeight;
      
      // Fill event counter histogram
      fHistograms->fhEvents->Fill(EECHistograms::kAll);          // All the events looped over
      
      //  ============================================
      //  ===== Apply all the event quality cuts =====
      //  ============================================
      
      // Jet trigger combinations
      caloJet15Trigger = (fJetReader->GetCaloJet15FilterBit() == 1);
      caloJet30Trigger = (fJetReader->GetCaloJet30FilterBit() == 1);
      caloJet40Trigger = (fJetReader->GetCaloJet40FilterBit() == 1);
      caloJet60Trigger = (fJetReader->GetCaloJet60FilterBit() == 1);
      caloJet80Trigger = (fJetReader->GetCaloJet80FilterBit() == 1);
      caloJet100Trigger = (fJetReader->GetCaloJet100FilterBit() == 1);
      
      // Fill the trigger histograms before the trigger selection
      if(fFillEventInformation){
        if(!caloJet60Trigger && !caloJet80Trigger && !caloJet100Trigger) fHistograms->fhTriggers->Fill(EECHistograms::kNoTrigger);
        if(caloJet60Trigger && !caloJet80Trigger && !caloJet100Trigger) fHistograms->fhTriggers->Fill(EECHistograms::kOnlyCaloJet60);
        if(!caloJet60Trigger && caloJet80Trigger && !caloJet100Trigger) fHistograms->fhTriggers->Fill(EECHistograms::kOnlyCaloJet80);
        if(!caloJet60Trigger && !caloJet80Trigger && caloJet100Trigger) fHistograms->fhTriggers->Fill(EECHistograms::kOnlyCaloJet100);
        if(caloJet60Trigger && caloJet80Trigger && !caloJet100Trigger) fHistograms->fhTriggers->Fill(EECHistograms::kCaloJet60And80);
        if(caloJet60Trigger && !caloJet80Trigger && caloJet100Trigger) fHistograms->fhTriggers->Fill(EECHistograms::kCaloJet60And100);
        if(!caloJet60Trigger && caloJet80Trigger && caloJet100Trigger) fHistograms->fhTriggers->Fill(EECHistograms::kCaloJet80And100);
        if(caloJet60Trigger && caloJet80Trigger && caloJet100Trigger) fHistograms->fhTriggers->Fill(EECHistograms::kCaloJet60And80And100);
      }
      
      // Make jet trigger selection
      if(fTriggerSelection == 1 && !caloJet80Trigger) continue;   // Select events with CaloJet80 trigger
      if(fTriggerSelection == 2 && !caloJet100Trigger) continue;  // Select events with CaloJet100 trigger
      if(fTriggerSelection == 3 && (!caloJet80Trigger && !caloJet100Trigger)) continue; // Select events with CaloJet80 OR CaloJet100 triggers
      if(fTriggerSelection == 4 && !caloJet60Trigger) continue;   // Select events with CaloJet60 trigger
      if(fTriggerSelection == 5 && (!caloJet80Trigger && !caloJet100Trigger)) continue; // Select events with CaloJet60 OR CaloJet80 OR CaloJet100 triggers. This selection is used with the sample forested filtering with CaloJet80 and CaloJet100 trigger.
      if(fTriggerSelection == 6 && (!caloJet60Trigger || caloJet80Trigger || caloJet100Trigger)) continue; // Select events with CaloJet60 OR CaloJet80 OR CaloJet100 triggers. This selection is used with the sample forested filtering with CaloJet60 trigger. Any events containing CaloJet80 or CaloJet100 triggers must be vetoed to avoid double counting when combining the two samples.
      if(fTriggerSelection == 7 && (!caloJet60Trigger && !caloJet80Trigger)) continue; // Select events with CaloJet60 OR CaloJet80 triggers
      if(fTriggerSelection == 15 && !caloJet15Trigger) continue; // Select events with CaloJet15 trigger
      if(fTriggerSelection == 30 && !caloJet30Trigger) continue; // Select events with CaloJet30 trigger
      if(fTriggerSelection == 40 && !caloJet40Trigger) continue; // Select events with CaloJet40 trigger

      // Fill the histogram for triggered events
      fHistograms->fhEvents->Fill(EECHistograms::kTriggered);

      // Check event cuts
      if(!PassEventCuts(fJetReader,fFillEventInformation)) continue;
      
      // If combining triggers, need to include event weight for events that only fire the lower trigger
      // The weight used here is the inverse of the average effective prescale in the whole sample
      // The effective prescale number is determined using the macro findEffectivePrescale
      // Git hash for the version used to get the number below is 2db6ffe6f5433b6b9f56a30fa2cd16b21d76561d
      // The input file used is eecAnalysis_akFlowJet_findEffectivePrescale_processed_2023-02-28.root
      if(fTriggerSelection == 3 && caloJet80Trigger && !caloJet100Trigger){
        fTotalEventWeight = fTotalEventWeight * 2.56248;
      }

      // The weight for the CaloJet60 trigger is much bigger, since the prescale is also much bigger
      // Trigger selection bit 6 already strips away events where also CaloJet80 or CaloJet100 fires, se we do not need to repeat that condition here
      // The number below is extracted from the following files:
      // eecAnalysis_akFlowJet_onlyJets_findTrigger60Weight_60sample_2023-04-26.root
      // eecAnalysis_akFlowJet_onlyJets_findTrigger60Weight_80and100sample_2023-04-26.root
      if(fTriggerSelection == 6){
        fTotalEventWeight = fTotalEventWeight * 35.143313;
      }
      
      // If combining triggers, need to include event weight for events that only fire the lower trigger
      // The weight used here is the inverse of the average effective prescale in the whole sample
      // The effective prescale number is calculated directly from the entries of the trigger histogram
      // The input file used is ppData_pfJets_wtaAxis_triggerCountsForCaloJet60Prescale_2023-08-07.root
      if(fTriggerSelection == 7 && caloJet60Trigger && !caloJet80Trigger){
        fTotalEventWeight = fTotalEventWeight * 6.3567036;
      }
      
      // Fill the event information histograms for the events that pass the event cuts
      if(fFillEventInformation){
        fHistograms->fhVertexZ->Fill(vz, fPtHatWeight);                         // z vertex distribution from all events
        fHistograms->fhVertexZWeighted->Fill(vz, fTotalEventWeight);            // z-vertex distribution weighted with the weight function
        fHistograms->fhCentrality->Fill(centrality, fPtHatWeight);              // Centrality filled from all events
        fHistograms->fhCentralityWeighted->Fill(centrality, fTotalEventWeight); // Centrality weighted with the centrality weighting function
        fHistograms->fhPtHat->Fill(ptHat);                                      // pT hat histogram
        fHistograms->fhPtHatWeighted->Fill(ptHat, fTotalEventWeight);           // pT het histogram weighted with corresponding cross section and event number

        // Fill the trigger histograms after the trigger selection
        if(!caloJet60Trigger && !caloJet80Trigger && !caloJet100Trigger) fHistograms->fhTriggersAfterSelection->Fill(EECHistograms::kNoTrigger);
        if(caloJet60Trigger && !caloJet80Trigger && !caloJet100Trigger) fHistograms->fhTriggersAfterSelection->Fill(EECHistograms::kOnlyCaloJet60);
        if(!caloJet60Trigger && caloJet80Trigger && !caloJet100Trigger) fHistograms->fhTriggersAfterSelection->Fill(EECHistograms::kOnlyCaloJet80);
        if(!caloJet60Trigger && !caloJet80Trigger && caloJet100Trigger) fHistograms->fhTriggersAfterSelection->Fill(EECHistograms::kOnlyCaloJet100);
        if(caloJet60Trigger && caloJet80Trigger && !caloJet100Trigger) fHistograms->fhTriggersAfterSelection->Fill(EECHistograms::kCaloJet60And80);
        if(caloJet60Trigger && !caloJet80Trigger && caloJet100Trigger) fHistograms->fhTriggersAfterSelection->Fill(EECHistograms::kCaloJet60And100);
        if(!caloJet60Trigger && caloJet80Trigger && caloJet100Trigger) fHistograms->fhTriggersAfterSelection->Fill(EECHistograms::kCaloJet80And100);
        if(caloJet60Trigger && caloJet80Trigger && caloJet100Trigger) fHistograms->fhTriggersAfterSelection->Fill(EECHistograms::kCaloJet60And80And100);
      }
      
      // ======================================
      // ===== Event quality cuts applied =====
      // ======================================
      
      // The unfolding study is completely separate from regular analysis to make code easier to follow
      if(fFillJetPtUnfoldingResponse) FillUnfoldingResponse();
      if(fFillTrackParticleMatchingHistograms) ConstructParticleResponses(); // Construct DeltaR and pT1*pT2 response matrices

      //****************************************************
      //    Loop over all jets and find tracks within jets
      //****************************************************
      
      const Bool_t circleJet = false; // Instead of actual jets, just draw a random circle somewhere to the event to mimic a jet
      Int_t nJets = circleJet ? 1 : fJetReader->GetNJets();

      // Jet loop
      leadingJetPt = 0;
      nGoodJets = 0;
      for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++){
        
        // Only do actual jet stuff with actual jets
        if(circleJet){
          jetPt = 125;
          jetPhi = fRng->Uniform(-TMath::Pi(), TMath::Pi());
          jetEta = fRng->Uniform(-fJetEtaCut, fJetEtaCut);
          jetReflectedEta = GetReflectedEta(jetEta);
          jetFlavor = 0;
        } else {
          jetPt = fJetReader->GetJetRawPt(jetIndex);  // Get the raw pT and do manual correction later
          jetPhi = fJetReader->GetJetPhi(jetIndex);
          jetEta = fJetReader->GetJetEta(jetIndex);
          jetReflectedEta = GetReflectedEta(jetEta);
          jetFlavor = 0;

          // For symmetric systems, eta in center-of-mass frame is the same as eta in the lab frame
          jetEtaCM = jetEta;

          // For asymmetric systems, we will need to boost to get to center-of-mass frame
          if(fDataType == ForestReader::kPPb_pToMinusEta) jetEtaCM += 0.465; 
          if(fDataType == ForestReader::kPPb_pToPlusEta) jetEtaCM -= 0.465;
          if(fDataType == ForestReader::kPPb_pToMinusEta_5TeV) jetEtaCM += 0.465; 
          if(fDataType == ForestReader::kPPbMC_pToMinusEta) jetEtaCM += 0.465; 
          if(fDataType == ForestReader::kPPbMC_pToPlusEta) jetEtaCM -= 0.465;
          
          // For data, instead of jet flavor, mark positive vz with 1 and negative with 0
          // This is used in one of the systematic checks for long range correlations
          if((fIsRealData) && vz > 0) jetFlavor = 1;
          
          //  ========================================
          //  ======== Apply jet quality cuts ========
          //  ========================================
          
          if(TMath::Abs(jetEta) >= fJetEtaCut) continue; // Cut for jet eta
          if(fCutWithJetEtaCM && TMath::Abs(jetEtaCM) >= fJetEtaCut - 0.465) continue;
          if(fCutBadPhiRegion && (jetPhi > -0.1 && jetPhi < 1.2)) continue; // Cut the area of large inefficiency in tracker
          
          // No jet quality cuts for generator level jets
          if(!(fMcCorrelationType == kGenGen || fMcCorrelationType == kGenReco)){
            
            // I forgot to include the trackMax branch in the "NewRelease" MC production
            // This branch in correctly included in default MC and in the data production
            // But in case we want to use the NewRelease-files, we need to calculate trackMax manually
            // The trackMax array is initialized to -1, so if trackMax is exactly -1, we know this branch is not included in the forest.
            // This piece of code can be removed if NewRelease part of the MC is reforested or is not used anymore.
            bool isNewReleaseMC = (fJetReader->GetJetMaxTrackPt(jetIndex) == -1);
            
            if(isNewReleaseMC){
              
              double manualMaxTrackPt = 0;
              
              nTracks = fJetReader->GetNTracks();
              for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
                
                // Check that all the track cuts are passed
                if(!PassTrackCuts(fJetReader,iTrack,fHistograms->fhTrackCuts,true)) continue;
                
                // Find the track pT, eta and phi
                trackPt = fJetReader->GetTrackPt(iTrack);
                trackEta = fJetReader->GetTrackEta(iTrack);
                trackPhi = fJetReader->GetTrackPhi(iTrack);
                
                // If the track is close to a jet, change the track eta-phi coordinates to a system where the jet axis is at origin
                deltaRTrackJet = GetDeltaR(jetEta, jetPhi, trackEta, trackPhi);
                
                if(deltaRTrackJet < 0.4){
                  if(trackPt > manualMaxTrackPt) manualMaxTrackPt = trackPt;
                }
                
              } // Track pT loop for manual calculation
              
              if(fMinimumMaxTrackPtFraction >= manualMaxTrackPt/fJetReader->GetJetRawPt(jetIndex)){
                continue; // Cut for jets with only very low pT particles
              }
              if(fMaximumMaxTrackPtFraction <= manualMaxTrackPt/fJetReader->GetJetRawPt(jetIndex)){
                continue; // Cut for jets where all the pT is taken by one track
              }
              
            } else {
              // Default jet quality cut, used in all other cases expect for NewRelease MC test
              
              if(fMinimumMaxTrackPtFraction >= fJetReader->GetJetMaxTrackPt(jetIndex)/fJetReader->GetJetRawPt(jetIndex)){
                continue; // Cut for jets with only very low pT particles
              }
              if(fMaximumMaxTrackPtFraction <= fJetReader->GetJetMaxTrackPt(jetIndex)/fJetReader->GetJetRawPt(jetIndex)){
                continue; // Cut for jets where all the pT is taken by one track
              }
            }
          }
          
          // Jet matching between reconstructed and generator level jets
          if((fMatchJets == 1) && !fJetReader->HasMatchingJet(jetIndex)){
            unmatchedCounter++;
            continue;
          }      
          
          // Require also reference parton flavor to be quark [-6,-1] U [1,6] or gluon (21)
          // We need to match gen jets to reco to get the parton flavor, but for reco jets it is always available in the forest
          // Here should implement an option if only quark and gluon tagged jets should be allowed in final results!
          if(fMatchJets == 1){
            
            matchedCounter++; // For debugging purposes, count the number of matched jets
            jetFlavor = 0;    // Jet flavor. 0 = Quark jet.
            
            partonFlavor = fJetReader->GetPartonFlavor(jetIndex);
            if(TMath::Abs(partonFlavor) == 21) jetFlavor = 1; // 1 = Gluon jet
            
          }
          
          //  ========================================
          //  ======= Jet quality cuts applied =======
          //  ========================================
          
          // Do a correction for the jet pT
          fJetCorrector->SetJetPT(jetPt);
          fJetCorrector->SetJetEta(jetEta);
          fJetCorrector->SetJetPhi(jetPhi);
          
          fJetUncertainty->SetJetPT(jetPt);
          fJetUncertainty->SetJetEta(jetEta);
          fJetUncertainty->SetJetPhi(jetPhi);
          
          jetPtCorrected = fJetCorrector->GetCorrectedPT();
          
          // No jet pT correction for generator level jets
          if(!(fMcCorrelationType == kGenGen || fMcCorrelationType == kGenReco)) {
            jetPt = jetPtCorrected;
            
            // If we are using smearing scenario, modify the jet pT using gaussian smearing
            // Note that for MC, we need to smear the jet pT to match the resolution in data
            if(fJetUncertaintyMode > 0){
              smearingFactor = GetSmearingFactor(jetPt, jetEta, centrality);
              jetPt = jetPt * fRng->Gaus(1,smearingFactor);
            }

            // If we are making runs using variation of jet pT within uncertainties, modify the jet pT here
            // Notice that we still need to apply the nominal jet pT smearing in MC before varying the scale
            if(fJetUncertaintyMode == 1) jetPt = jetPt * (1 - fJetUncertainty->GetUncertainty().first);
            if(fJetUncertaintyMode == 2) jetPt = jetPt * (1 + fJetUncertainty->GetUncertainty().second);
            
          }
          
          // After the jet pT can been corrected, apply analysis jet pT cuts
          if(jetPt < fJetMinimumPtCut) continue;
          if(jetPt > fJetMaximumPtCut) continue;

          // Check if a high pT jet is found from the event
          if(jetPt > 80 && !jetOver80Found){
            jetOver80Found = true;
            fHistograms->fhEvents->Fill(EECHistograms::kJetOver80);
          }

          if(jetPt > 120 && !jetOver120Found){
            jetOver120Found = true;
            fHistograms->fhEvents->Fill(EECHistograms::kJetOver120);
          }

          
          // If we are matching jets, require that the matched jet has at least half of the reconstructed pT
          if(fMatchJets > 0){
            antiMatchPtVeto = false;
            if(jetPt*0.5 > fJetReader->GetMatchedPt(jetIndex) || fJetReader->GetMatchedPt(jetIndex) * 0.5 > jetPt){
              matchRejectedDueToPtCounter++;
              antiMatchPtVeto = true;
              if(fMatchJets == 1) continue;
            }
          }

          // Check if an anti-matching jet exists after jet pT cuts
          // Leave these debugging messages here for now
          /*if(fMatchJets == 2){
            if(!fJetReader->HasMatchingJet(jetIndex)){
              cout << "There is an anti-mathcing jet in event: " << iEvent << endl;
            }
            if(antiMatchPtVeto){
              cout << "There is a large pT gap in events: " << iEvent << endl;
            }
          }*/

          //cout << "Jet pT: " << jetPt << endl;
          //cout << "Matched pT: " << fJetReader->GetMatchedPt(jetIndex) << endl;

          // Anti-matching between reconstructed and generator level jets
          if((fMatchJets == 2) && fJetReader->HasMatchingJet(jetIndex) && !antiMatchPtVeto){
            continue;
          }    
          
          //************************************************
          //       Fill histograms for jet pT closure
          //************************************************
          
          // Only fill is matching is enabled and histograms are selected for filling
          if(fFillJetPtClosure && (fMatchJets == 1)) FillJetPtClosureHistograms(jetIndex);
          
        } // Circle jet if

        // Once the jet quality cuts are applied, determine the leading jet pT
        if(jetPt > leadingJetPt) leadingJetPt = jetPt;

        // We have a good jet!
        nGoodJets++;
        
        //************************************************
        //         Fill histograms for all jets
        //************************************************
        
        // Only fill the any jet histogram if selected
        if(fFillJetHistograms){
          
          jetPtCorrected = jetPt;
          
          // Find the pT weight for the jet
          jetPtWeight = GetJetPtWeight(jetPtCorrected);
          
          // Fill the axes in correct order
          fillerJet[0] = jetPtCorrected;          // Axis 0 = inclusive jet pT
          fillerJet[1] = jetPhi;                  // Axis 1 = inclusive jet phi
          fillerJet[2] = jetEta;                  // Axis 2 = inclusive jet eta
          fillerJet[3] = jetEtaCM;                // Axis 3 = inclusive jet eta_CM
          fillerJet[4] = centrality;              // Axis 4 = centrality
          fillerJet[5] = jetFlavor;               // Axis 5 = flavor of the jet
          
          fHistograms->fhInclusiveJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight); // Fill the data point to histogram
          
          
        } // Check if we want to fill any jet histograms
        
        //************************************************
        //   Do energy-energy correlation within jets
        //************************************************
        
        if(fFillEnergyEnergyCorrelators || fFillEnergyEnergyCorrelatorsSystematics || fFillJetConeHistograms){
          
          // Clear the vectors of track kinematics for tracks selected for energy-energy correlators
          for(int iJetConeType = 0; iJetConeType < EECHistograms::knJetConeTypes; iJetConeType++){
            selectedTrackPt[iJetConeType].clear();
            relativeTrackEta[iJetConeType].clear();
            relativeTrackPhi[iJetConeType].clear();
            selectedTrackSubevent[iJetConeType].clear();
          }
          
          // Clear the multiplicity arrays
          for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
            for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
              for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
                multiplicityBin = std::make_tuple(iTrackPt,iJetCone,iSubevent);
                uncorrectedMultiplicityInJetCone[multiplicityBin] = 0;
                multiplicityInJetCone[multiplicityBin] = 0;
              }
            }
          }
          
          // Loop over tracks and check which are within the jet radius
          nTracks = fTrackReader->GetNTracks();
          trackMultiplicity = 0;
          trackMultiplicityWeighted = 0;
          for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
            
            // Check that all the track cuts are passed
            if(!PassTrackCuts(fTrackReader,iTrack,fHistograms->fhTrackCuts,true)) continue;
            
            // Find the track eta and phi
            trackPt = fTrackReader->GetTrackPt(iTrack);
            trackEta = fTrackReader->GetTrackEta(iTrack);
            trackPhi = fTrackReader->GetTrackPhi(iTrack);
            trackEfficiencyCorrection = GetTrackEfficiencyCorrection(trackPt, trackEta, hiBin);
            trackSubevent = fTrackReader->GetTrackSubevent(iTrack);
            trackSubeventIndex = GetSubeventIndex(trackSubevent);

            trackMultiplicity += 1;
            trackMultiplicityWeighted += trackEfficiencyCorrection;
            
            // If the track is close to a jet, change the track eta-phi coordinates to a system where the jet axis is at origin
            deltaRTrackJet = GetDeltaR(jetEta, jetPhi, trackEta, trackPhi);
            if(deltaRTrackJet < fJetRadius){
              relativeTrackPhi[EECHistograms::kSignalCone].push_back(trackPhi-jetPhi);
              relativeTrackEta[EECHistograms::kSignalCone].push_back(trackEta-jetEta);
              
              // Also remember track pT and subevent
              selectedTrackPt[EECHistograms::kSignalCone].push_back(trackPt);
              selectedTrackSubevent[EECHistograms::kSignalCone].push_back(trackSubevent);
              
              // If we are calculating multiplicity within the jet cone, update the multiplicity arrays
              if(fFillJetConeHistograms){
                for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
                  if(trackPt < trackPtBinsEEC[iTrackPt]) break;
                  multiplicityBin = std::make_tuple(iTrackPt,EECHistograms::kSignalCone,EECHistograms::knSubeventTypes);
                  uncorrectedMultiplicityInJetCone[multiplicityBin]++;
                  multiplicityInJetCone[multiplicityBin] += trackEfficiencyCorrection;
                  if(trackSubeventIndex >= 0){
                    multiplicityBin = std::make_tuple(iTrackPt,EECHistograms::kSignalCone,trackSubeventIndex);
                    uncorrectedMultiplicityInJetCone[multiplicityBin]++;
                    multiplicityInJetCone[multiplicityBin] += trackEfficiencyCorrection;
                  }
                } // Track pT loop
                
                // Also if the jet cone histograms are filled, find the maximum pT within the jet
                if(trackSubeventIndex > 0){
                  if(trackPt > maxTrackPtInJetBackground) maxTrackPtInJetBackground = trackPt;
                } else {
                  if(trackPt > maxTrackPtInJetSignal) maxTrackPtInJetSignal = trackPt;
                }
                
              } // Fill jet cone histograms
              
            } // Track close to jet
            
            // For the track density, use a fixed cone size around the jet axis. TODO: Synchronize the cone size with EECHistograms
            if(fFillJetConeHistograms){
              
              if(deltaRTrackJet < 0.8){
                fillerParticleDensityInJetCone[0] = deltaRTrackJet;     // Axis 0: DeltaR between the track and the jet
                fillerParticleDensityInJetCone[1] = jetPt;              // Axis 1: jet pT
                fillerParticleDensityInJetCone[2] = trackPt;            // Axis 2: track pT
                fillerParticleDensityInJetCone[3] = centrality;         // Axis 3: centrality
                fillerParticleDensityInJetCone[4] = 0;                  // Axis 4: 0 is the index for signal cone
                fillerParticleDensityInJetCone[5] = trackSubeventIndex; // Axis 5: Subevent index for the track
                fHistograms->fhParticleDensityAroundJet->Fill(fillerParticleDensityInJetCone, fTotalEventWeight * trackEfficiencyCorrection);
                fHistograms->fhParticlePtDensityAroundJet->Fill(fillerParticleDensityInJetCone, fTotalEventWeight * trackEfficiencyCorrection * trackPt);
              }
            }
            
            // If the track is close to a reflected jet axis, change the track eta-phi coordinates to a system where the reflected jet axis is at origin
            if(fDoReflectedCone){
              deltaRTrackJet = GetDeltaR(jetReflectedEta, jetPhi, trackEta, trackPhi);
              if(deltaRTrackJet < fJetRadius){
                relativeTrackPhi[EECHistograms::kReflectedCone].push_back(trackPhi-jetPhi);
                relativeTrackEta[EECHistograms::kReflectedCone].push_back(trackEta-jetReflectedEta);
                
                // Also remember track pT and subevent
                selectedTrackPt[EECHistograms::kReflectedCone].push_back(trackPt);
                selectedTrackSubevent[EECHistograms::kReflectedCone].push_back(trackSubevent);
                
                // If we are calculating multiplicity within the reflected cone, update the multiplicity arrays
                if(fFillJetConeHistograms){
                  for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
                    if(trackPt < trackPtBinsEEC[iTrackPt]) break;
                    multiplicityBin = std::make_tuple(iTrackPt,EECHistograms::kReflectedCone,EECHistograms::knSubeventTypes);
                    uncorrectedMultiplicityInJetCone[multiplicityBin]++;
                    multiplicityInJetCone[multiplicityBin] += trackEfficiencyCorrection;
                    if(trackSubeventIndex >= 0){
                      multiplicityBin = std::make_tuple(iTrackPt,EECHistograms::kReflectedCone,trackSubeventIndex);
                      uncorrectedMultiplicityInJetCone[multiplicityBin]++;
                      multiplicityInJetCone[multiplicityBin] += trackEfficiencyCorrection;
                    }
                  } // Track pT loop
                } // Fill jet cone histograms
              }
              
              // For the track density, use a fixed cone size around the jet axis. TODO: Synchronize the cone size with EECHistograms
              if(fFillJetConeHistograms){
                if(deltaRTrackJet < 0.8){
                  fillerParticleDensityInJetCone[0] = deltaRTrackJet;     // Axis 0: DeltaR between the track and the jet
                  fillerParticleDensityInJetCone[1] = jetPt;              // Axis 1: jet pT
                  fillerParticleDensityInJetCone[2] = trackPt;            // Axis 2: track pT
                  fillerParticleDensityInJetCone[3] = centrality;         // Axis 3: centrality
                  fillerParticleDensityInJetCone[4] = 1;                  // Axis 4: 1 is the index for reflected cone
                  fillerParticleDensityInJetCone[5] = trackSubeventIndex; // Axis 5: Subevent index for the track
                  fHistograms->fhParticleDensityAroundJet->Fill(fillerParticleDensityInJetCone, fTotalEventWeight * trackEfficiencyCorrection);
                  fHistograms->fhParticlePtDensityAroundJet->Fill(fillerParticleDensityInJetCone, fTotalEventWeight * trackEfficiencyCorrection * trackPt);
                }
              }
            } // Search for track from reflected cone

            // Collect the selected track information for ALICE-style perpendicular cones
            if(fDoPerpendicularCone){

              for(int iPerpendicularCone = 0; iPerpendicularCone <= 1; iPerpendicularCone++){

                jetPerpendicularPhi = GetPerpendicularPhi(jetPhi, iPerpendicularCone);

                deltaRTrackJet = GetDeltaR(jetEta, jetPerpendicularPhi, trackEta, trackPhi);
                if(deltaRTrackJet < fJetRadius){
                  relativeTrackPhi[EECHistograms::kPerpendicularCone + iPerpendicularCone].push_back(trackPhi-jetPerpendicularPhi);
                  relativeTrackEta[EECHistograms::kPerpendicularCone + iPerpendicularCone].push_back(trackEta-jetEta);
                
                  // Also remember track pT and subevent
                  selectedTrackPt[EECHistograms::kPerpendicularCone + iPerpendicularCone].push_back(trackPt);
                  selectedTrackSubevent[EECHistograms::kPerpendicularCone + iPerpendicularCone].push_back(trackSubevent);
                
                  // If we are calculating multiplicity within the perpendicular cone, update the multiplicity arrays
                  if(fFillJetConeHistograms){
                    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
                      if(trackPt < trackPtBinsEEC[iTrackPt]) break;
                      multiplicityBin = std::make_tuple(iTrackPt,EECHistograms::kPerpendicularCone + iPerpendicularCone,EECHistograms::knSubeventTypes);
                      uncorrectedMultiplicityInJetCone[multiplicityBin]++;
                      multiplicityInJetCone[multiplicityBin] += trackEfficiencyCorrection;
                      if(trackSubeventIndex >= 0){
                        multiplicityBin = std::make_tuple(iTrackPt,EECHistograms::kPerpendicularCone + iPerpendicularCone,trackEfficiencyCorrection);
                        uncorrectedMultiplicityInJetCone[multiplicityBin]++;
                        multiplicityInJetCone[multiplicityBin] += trackEfficiencyCorrection;
                      }
                    } // Track pT loop
                  } // Fill jet cone histograms
                } // Found track to add to vectors
              } // Loop over two perpendicular cones
            } // Search for track from perpendicular cone
            
          } // Track loop

          //******************************************
          //           Reflected cone QA       
          //******************************************

          // Fill the quality assuranve histograms for reflected cone
          if(fDoReflectedCone){

            // Loop over the jets again
            nJetsInReflectedCone = 0;

            for(Int_t reflectedConeJetIndex = 0; reflectedConeJetIndex < nJets; reflectedConeJetIndex++){

              reflectedConeJetPt = fJetReader->GetJetRawPt(reflectedConeJetIndex);  // Get the raw pT and do manual correction later
              reflectedConeJetPhi = fJetReader->GetJetPhi(reflectedConeJetIndex);
              reflectedConeJetEta = fJetReader->GetJetEta(reflectedConeJetIndex);

              // First check if the jet is close to the reflected cone
              if(GetDeltaR(jetReflectedEta, jetPhi, reflectedConeJetEta, reflectedConeJetPhi) > fJetRadius) continue;

              // Jet quality cuts
              if(TMath::Abs(reflectedConeJetEta) >= fJetEtaCut) continue; // Cut for jet eta
              if(fCutBadPhiRegion && (reflectedConeJetPhi > -0.1 && reflectedConeJetPhi < 1.2)) continue; // Cut the area of large inefficiency in tracker

              if(fMinimumMaxTrackPtFraction >= fJetReader->GetJetMaxTrackPt(reflectedConeJetIndex)/fJetReader->GetJetRawPt(reflectedConeJetIndex)){
                continue; // Cut for jets with only very low pT particles
              }
              if(fMaximumMaxTrackPtFraction <= fJetReader->GetJetMaxTrackPt(reflectedConeJetIndex)/fJetReader->GetJetRawPt(reflectedConeJetIndex)){
                continue; // Cut for jets where all the pT is taken by one track
              }

              // For non-generator level jet pT: do a correction for the jet pT
              if(!(fMcCorrelationType == kGenGen || fMcCorrelationType == kGenReco)){              
                fJetCorrector->SetJetPT(reflectedConeJetPt);
                fJetCorrector->SetJetEta(reflectedConeJetEta);
                fJetCorrector->SetJetPhi(reflectedConeJetPhi);

                reflectedConeJetPt = fJetCorrector->GetCorrectedPT();
              }

              // After the jet pT can been corrected, apply jet pT cuts for the jets in reflected cone
              if(reflectedConeJetPt < 25) continue;
              if(reflectedConeJetPt > fJetMaximumPtCut) continue;

              // If all the cuts are passed, we have found a jet within the reflected cone
              nJetsInReflectedCone++;

              // Fill the histograms containing all the jet pT:s found within the reflected cone
              fillerReflectedConeQA[0] = reflectedConeJetPt;
              fillerReflectedConeQA[1] = centrality;
              fHistograms->fhJetPtInReflectedCone->Fill(fillerReflectedConeQA, fTotalEventWeight);
            }

            // Fill the histograms for the total number of jets found in reflected cone
            fillerReflectedConeQA[0] = nJetsInReflectedCone;
            fillerReflectedConeQA[1] = centrality;
            fHistograms->fhJetNumberInReflectedCone->Fill(fillerReflectedConeQA, fTotalEventWeight);

          } // Reflected cone QA if

          //***********************************************
          //            Mixed cone background      
          //***********************************************

          // Estimate background from cones dropped to mixed events.
          if(fDoMixedCone){

            // Event matching is different for different systems
            if(fEventMatchForMixing == kMatchHiBin){
              // Match HiBin
              FindHiBinMatchedEvents(mixedEventIndices, 2, vz, hiBin, eventNumber, iEvent, fBinnedMixing);
            } else if (fEventMatchForMixing == kMatchHFEnergy){
              // Match HF energy
              FindHFEnergyMatchedEvents(mixedEventIndices, 2, vz, hfSum, eventNumber, iEvent, fBinnedMixing);
            } else {
              // Match background multiplicity
              currentMultiplicity = GetGenMultiplicity(fTrackReader, kSubeventNonZero, false);
              FindMultiplicityMatchedEvents(mixedEventIndices, 2, vz, currentMultiplicity, eventNumber, iEvent, fBinnedMixing); 
            }

            // Make the pPb MC "minimum bias" by only accepting particles from the EPOS simulation
            //if(fDataType == ForestReader::kPPbMC_pToPlusEta || fDataType == ForestReader::kPPbMC_pToMinusEta){
            //  if(fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenGen){
            //    fSubeventCut = kSubeventNonZero;
            //  }
            //}   

            // Fill the pairing arrays with the found mixed event indices
            mixedEventFillIndex = EECHistograms::kMixedCone;

            for(int mixedEventIndex : mixedEventIndices){

              // If match vz and hiBin, then load the event from the mixed event tree and find particles from jet cone region
              fMixedEventReader->GetEvent(mixedEventIndex);

              nTracks = fMixedEventReader->GetNTracks();
              for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){

                // Check that all the track cuts are passed. In mega skim mode, only tracks passing the cuts are saved, so no need to check the track cuts again here.
                if(!fMegaSkimMode){
                  if(!PassTrackCuts(fMixedEventReader,iTrack,fHistograms->fhTrackCuts,true)) continue;
                }

                // Find the track kinematics
                trackPt = fMixedEventReader->GetTrackPt(iTrack);
                trackEta = fMixedEventReader->GetTrackEta(iTrack);
                trackPhi = fMixedEventReader->GetTrackPhi(iTrack);
                trackSubevent = fMixedEventReader->GetTrackSubevent(iTrack);

                // If the track is close to a jet, change the track eta-phi coordinates to a system where the jet axis is at origin
                deltaRTrackJet = GetDeltaR(jetEta, jetPhi, trackEta, trackPhi);
                if(deltaRTrackJet < fJetRadius){
                  relativeTrackPhi[mixedEventFillIndex].push_back(trackPhi-jetPhi);
                  relativeTrackEta[mixedEventFillIndex].push_back(trackEta-jetEta);

                  // Also remember track pT and subevent
                  selectedTrackPt[mixedEventFillIndex].push_back(trackPt);
                  selectedTrackSubevent[mixedEventFillIndex].push_back(trackSubevent);
                } // Track is close to a jet
              } // Mixed event track loop

              // Currently there are two events we mix with, fill the second mixed cone to the second index
              mixedEventFillIndex = EECHistograms::kSecondMixedCone;
      
            } // Loop over good events to mix with

            // Return to original cut after analyzing all particles for the mixing check
            //fSubeventCut = fCard->Get("SubeventCut");

          } // Event mixing for PbPb data and MC

          //************************
          //      Jet cone Q/A      
          //************************
          
          // Fill the multiplicity histograms within the jet and maximum track pT within the jet cone histograms
          if(fFillJetConeHistograms){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
              
              // Fill the subevent hisotgrams only for Pythia+Hydjet simulation
              if(fDataType != ForestReader::kPbPbMC && iSubevent < EECHistograms::knSubeventTypes) continue;
              
              fillerMultiplicityInJetCone[1] = jetPt;      // Axis 1: Jet pT
              fillerMultiplicityInJetCone[3] = centrality; // Axis 3: centrality
              fillerMultiplicityInJetCone[4] = iSubevent;  // Axis 4: Subevent index
              for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){

                // Multiplicity in jet cone
                multiplicityBin = std::make_tuple(iTrackPt,EECHistograms::kSignalCone,iSubevent);
                fillerMultiplicityInJetCone[0] = multiplicityInJetCone[multiplicityBin];                         // Axis 0: Multiplicity
                fillerMultiplicityInJetCone[2] = trackPtBinsEEC[iTrackPt]+0.01;                                  // Axis 2: Track pT
                fHistograms->fhParticleMultiplicityInJet->Fill(fillerMultiplicityInJetCone,fTotalEventWeight);   // Fill the multiplicity within the jet cone histogram

                fillerMultiplicityInJetCone[0] = uncorrectedMultiplicityInJetCone[multiplicityBin];                        // Axis 0: Multiplicity
                fHistograms->fhParticleMultiplicityInJetUncorrected->Fill(fillerMultiplicityInJetCone,fTotalEventWeight);  // Fill the uncorrected multiplicity within the jet cone histogram

                // Multiplicity in reflected cone
                multiplicityBin = std::make_tuple(iTrackPt,EECHistograms::kReflectedCone,iSubevent);
                fillerMultiplicityInJetCone[0] = multiplicityInJetCone[multiplicityBin];                                    // Axis 0: Multiplicity
                fHistograms->fhParticleMultiplicityInReflectedCone->Fill(fillerMultiplicityInJetCone,fTotalEventWeight);    // Fill the multiplicity within the reflected cone histogram

                fillerMultiplicityInJetCone[0] = uncorrectedMultiplicityInJetCone[multiplicityBin];                                  // Axis 0: Multiplicity
                fHistograms->fhParticleMultiplicityInReflectedConeUncorrected->Fill(fillerMultiplicityInJetCone,fTotalEventWeight);  // Fill the uncorrected multiplicity within the reflected cone histogram
              }
              
            } // Subevent loop
            
            // Fill the max particle pT within the jet cone histograms
            fillerMaxParticlePtInJetCone[0] = jetPt;                     // Axis 0: Jet pT
            fillerMaxParticlePtInJetCone[1] = maxTrackPtInJetSignal;     // Axis 1: Maximum particle pT for signal particles
            fillerMaxParticlePtInJetCone[2] = maxTrackPtInJetBackground; // Axis 2: Maximum particle pT for background particles
            fillerMaxParticlePtInJetCone[3] = centrality;                // Axis 3: centrality
            fHistograms->fhMaxPtParticleInJet->Fill(fillerMaxParticlePtInJetCone,fTotalEventWeight);
            
          } // Fill multiplicity in jets and maximum track pT within the jet cone histograms

          //************************************************
          //        Fill energy-energy correlators    
          //************************************************
          
          // Calculate the energy-energy correlator within this jet
          if(fFillEnergyEnergyCorrelators || fFillEnergyEnergyCorrelatorsSystematics){

            // Before calculating the energy-energy correlators, reset histograms for covariance matrix calculation
            // Covariance matrix calculation only for data, not needed in MC
            if(!fSkipCovarianceMatrix){
              for(auto iWeightExponent = 0; iWeightExponent < fWeightExponent.size(); iWeightExponent++){
                for(Int_t iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
                  fThisEventCorrelator[iWeightExponent][iTrackPt]->Reset();
                } // Track pT loop for covariance matrices
              }
            }

            // Then calculate the energy-energy correlator for this jet
            CalculateEnergyEnergyCorrelator(selectedTrackPt, relativeTrackEta, relativeTrackPhi, selectedTrackSubevent, jetPt);

            // After that is done, add the contents from this event to the total coveriance matrix
            // Covariance matrix calculation only for data, not needed in MC
            if(!fSkipCovarianceMatrix){
              fillerUnfoldingCovariance[3] = centrality; // Axis 3: Centrality
              for(auto iWeightExponent = 0; iWeightExponent < fWeightExponent.size(); iWeightExponent++){
                fillerUnfoldingCovariance[4] = iWeightExponent; // Axis 4: Weight exponent index
                for(Int_t iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
                  trackPt = fCard->Get("TrackPtBinEdgesEEC", iTrackPt) + 0.1;
                  fillerUnfoldingCovariance[2] = trackPt;  // Axis 2: Track pT corresponding to bin
                  for(Int_t iDeltaRBin = 1; iDeltaRBin <= nDeltaRBins; iDeltaRBin++){
                    deltaRBinContent1 = fThisEventCorrelator[iWeightExponent][iTrackPt]->GetBinContent(iDeltaRBin);
                    transformedDeltaR1 = TransformToUnfoldingAxis(fHistograms->GetDeltaRBinBorderLowEEC(iDeltaRBin)+0.0001, jetPt);

                    for(Int_t jDeltaRBin = 1; jDeltaRBin <= nDeltaRBins; jDeltaRBin++){
                      deltaRBinContent2 = fThisEventCorrelator[iWeightExponent][iTrackPt]->GetBinContent(jDeltaRBin);
                      transformedDeltaR2 = TransformToUnfoldingAxis(fHistograms->GetDeltaRBinBorderLowEEC(jDeltaRBin)+0.0001, jetPt);

                      // Calculate the covariance between the two bins
                      eecCovariance = deltaRBinContent1 * deltaRBinContent2;

                      // Fill the calculation result to correct bin
                      fillerUnfoldingCovariance[0] = transformedDeltaR1; // Axis 0: Transformed deltaR axis
                      fillerUnfoldingCovariance[1] = transformedDeltaR2; // Axis 1: Transformed deltaR axis

                      fHistograms->fhJetPtUnfoldingCovariance->Fill(fillerUnfoldingCovariance, eecCovariance);
                    
                    } // Inner DeltaR loop
                  } // Outer DeltaR loop
                } // Track pT loop for covariance matrices
              } // Weight exponent loop
            } // Fill covariance matrix histograms
          } // Fill energy-energy correlator histograms
          
        } // Fill anything in-jet particle pairing related
        
      } // End of jet loop

      // After the jet loop, we can fill the HF energy histograms for with leading jet pT
      if(fFillEventInformation){

        // Fill the axes for HF energy histogram
        fillerHFEnergy[0] = hfPlus;        // Axis 0 = HF plus
        fillerHFEnergy[1] = hfMinus;       // Axis 1 = HF minus
        fillerHFEnergy[2] = hfSum;         // Axis 2 = HF sum
        fillerHFEnergy[3] = leadingJetPt;  // Axis 3 = leading jet pT
        fillerHFEnergy[4] = nGoodJets;     // Axis 4 = Number of jets passing selection

        fHistograms->fhHFEnergy->Fill(fillerHFEnergy, fTotalEventWeight); 

      }
      
      //************************************************
      //       Fill histograms for inclusive tracks
      //************************************************
      
      // Inclusive track histograms
      if((fFillTrackHistograms || fFillEventInformation)){
        
        // Loop over all track in the event
        nTracks = fTrackReader->GetNTracks();
        trackMultiplicity = 0;
        trackMultiplicityWeighted = 0;
        for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
          
          // Check that all the track cuts are passed
          if(!PassTrackCuts(fTrackReader,iTrack,fHistograms->fhTrackCuts)) continue;
          
          // Get the efficiency correction
          trackPt = fTrackReader->GetTrackPt(iTrack);
          trackEta = fTrackReader->GetTrackEta(iTrack);
          trackPhi = fTrackReader->GetTrackPhi(iTrack);
          trackEfficiencyCorrection = GetTrackEfficiencyCorrection(iTrack);
          
          trackMultiplicity += 1;
          trackMultiplicityWeighted += trackEfficiencyCorrection;
          
          // Fill track histograms
          if(fFillTrackHistograms){
            fillerTrack[0] = trackPt;      // Axis 0: Track pT
            fillerTrack[1] = trackPhi;     // Axis 1: Track phi
            fillerTrack[2] = trackEta;     // Axis 2: Track eta
            fillerTrack[3] = centrality;   // Axis 3: Centrality
            fHistograms->fhTrack->Fill(fillerTrack,trackEfficiencyCorrection*fTotalEventWeight);  // Fill the track histogram
            fHistograms->fhTrackUncorrected->Fill(fillerTrack,fTotalEventWeight);                 // Fill the uncorrected track histogram
          }
          
        } // Track loop
        
        // Fill multiplicity histogram from all events
        if(fFillEventInformation){
          fillerMultiplicity[0] = trackMultiplicity;
          fillerMultiplicity[1] = trackMultiplicityWeighted;
          fillerMultiplicity[2] = centrality;
          fHistograms->fhMultiplicity->Fill(fillerMultiplicity, fTotalEventWeight);
        }
      }
      
    } // Event loop
    
    //************************************************
    //      Cleanup at the end of the file loop
    //************************************************
    
    // Close the input files after the event has been read
    inputFile->Close();
    if(useDifferentReaderForJetsAndTracks) copyInputFile->Close();
    if(fFillJetPtUnfoldingResponse || fFillTrackParticleMatchingHistograms) unfoldingInputFile->Close();
    
    // Write some debug messages if prompted
    if(fDebugLevel > 1){
      if(fMatchJets == 1){
        cout << "There were " << matchedCounter << " matched jets." << endl;
        cout << "There were " << unmatchedCounter << " unmatched jets." << endl;
        cout << "Too large pT gap between matched jets in " << matchRejectedDueToPtCounter << " pairs." << endl;
      }
    }
    
  } // File loop

  // Release the memory for arrays reserved with new operator
  delete[] trackPtBinsEEC;
  
}

/*
 * Method for calculating the energy-energy correlators for tracks within the jets
 *
 *  const vector<double> selectedTrackPt[EECHistograms::knJetConeTypes] = pT array for tracks selected for the analysis
 *  const vector<double> relativeTrackEta[EECHistograms::knJetConeTypes] = relative eta array for tracks selected for the analysis
 *  const vector<double> relativeTrackPhi[EECHistograms::knJetConeTypes] = relative phi array for tracks selected for the analysis
 *  const vector<int> selectedTrackSubevent[EECHistograms::knJetConeTypes] = subevent array for tracks selected for the analysis
 *  const double jetPt = pT of the jet the tracks are close to
 */
void EECAnalyzer::CalculateEnergyEnergyCorrelator(const vector<double> selectedTrackPt[EECHistograms::knJetConeTypes], const vector<double> relativeTrackEta[EECHistograms::knJetConeTypes], const vector<double> relativeTrackPhi[EECHistograms::knJetConeTypes], const vector<int> selectedTrackSubevent[EECHistograms::knJetConeTypes], const double jetPt){
  
  // Define a filler for THnSparse
  Double_t fillerEnergyEnergyCorrelator[7];       // Axes: deltaR, Jet pT, lower track pT, centrality, pairing type, subevent
  Double_t fillerParticleDeltaRResponseMatrix[5]; // Filler for validating deltaR smearing
  Double_t fillerParticlePtResponseMatrix[7];     // Filler for validating particle pT smearing
  
  // Indices for different pairing types
  Int_t firstParticleType[EECHistograms::knPairingTypes] =  {EECHistograms::kSignalCone, EECHistograms::kSignalCone,    EECHistograms::kReflectedCone, EECHistograms::kSignalCone, EECHistograms::kMixedCone, EECHistograms::kSignalCone, EECHistograms::kMixedCone, EECHistograms::kSecondMixedCone, EECHistograms::kSignalCone, EECHistograms::kPerpendicularCone, EECHistograms::kSignalCone, EECHistograms::kPerpendicularCone, EECHistograms::kSecondPerpendicularCone};
  Int_t secondParticleType[EECHistograms::knPairingTypes] = {EECHistograms::kSignalCone, EECHistograms::kReflectedCone, EECHistograms::kReflectedCone, EECHistograms::kMixedCone, EECHistograms::kMixedCone, EECHistograms::kSecondMixedCone, EECHistograms::kSecondMixedCone, EECHistograms::kSecondMixedCone, EECHistograms::kPerpendicularCone, EECHistograms::kPerpendicularCone, EECHistograms::kSecondPerpendicularCone, EECHistograms::kSecondPerpendicularCone, EECHistograms::kSecondPerpendicularCone};
  
  // Event information
  Double_t centrality = fTrackReader->GetCentrality();
  Int_t hiBin = fTrackReader->GetHiBin();
  if(fMultiplicityMode) centrality = GetCentralityFromMultiplicity(GetMultiplicity());
  
  // Variables for tracks
  Double_t trackPt1;       // Track pT for the first track
  Double_t trackEta1;      // Track eta for the first track
  Double_t trackPhi1;      // Track phi for the first track
  Double_t trackPt2;       // Track pT for the second track
  Double_t trackEta2;      // Track eta for the second track
  Double_t trackPhi2;      // Track phi for the second track
  Double_t trackDeltaR;    // DeltaR between the two tracks
  Double_t smearedDeltaR;  // Smeared DeltaR between two tracks
  Double_t lowerTrackPt;   // Lower of the two track pT:s
  Double_t trackEfficiencyCorrection1;    // Efficiency correction for the first track
  Double_t trackEfficiencyCorrection2;    // Efficiency correction for the first track
  Double_t trackPairEfficiencyCorrection; // Efficiency correction for track pair efficiency
  Double_t correlatorWeight;            // Weight given to the energy-energy correlator pT1*pT2
  Double_t smearedCorrelatorWeight;     // Smeared correlator weight
  Int_t subeventTrack1;    // Subevent index for the first track (0 = pythia, > 0 = hydjet)
  Int_t subeventTrack2;    // Subevent index for the second track (0 = pythia, > 0 = hydjet)
  Int_t subeventCombination;      // Subevent combination type (0 = pythia-pythia, 1 = pythia-hydjet, 2 = hydjet-pythia, 3 = hydjet-hydjet)
  Int_t startIndex;        // First index when looping over the second tracks
  Int_t weightIndex;       // Inxed in histogram axis where the weight exponent is saved

  // Variables for covariance matrix calculation
  const Int_t nTrackPtBinsEEC = fCard->GetNBin("TrackPtBinEdgesEEC"); // Number of track pT bins

  // Variables for systematic uncertainty study
  Double_t trackPairEfficiencyError;            // Error of the evaluated track pair efficiency
  Double_t trackPairInefficiency;               // Value for track pair inefficiency
  Double_t trackPairInefficiencyVariation;      // Systematic variation applied for the track pair inefficiency
  Double_t variedTrackPairEfficiencyCorrection; // Varied correction for track pair efficiency
  
  // Variables for the number of tracks
  Int_t nTracks1;
  Int_t nTracks2;
  
  // For MC: weight according to jet pT
  Double_t jetPtWeight = GetJetPtWeight(jetPt);
  
  // Loop over the pairing types (same jet/reflected cone jet)
  for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
    
    // Skip all the pairing types that are not specified in the configuration
    if(!fDoReflectedCone){
      if(iPairingType == EECHistograms::kSignalReflectedConePair) continue;
      if(iPairingType == EECHistograms::kReflectedConePair) continue;
    }

    if(!fDoMixedCone){
      if(iPairingType == EECHistograms::kSignalMixedConePair) continue;
      if(iPairingType == EECHistograms::kMixedConePair) continue;
      if(iPairingType == EECHistograms::kSignalSecondMixedConePair) continue;
      if(iPairingType == EECHistograms::kMixedMixedConePair) continue;
      if(iPairingType == EECHistograms::kSecondMixedConePair) continue;
    }

    if(!fDoPerpendicularCone){
      if(iPairingType == EECHistograms::kSignalPerpendicularConePair) continue;
      if(iPairingType == EECHistograms::kPerpendicularConePair) continue;
      if(iPairingType == EECHistograms::kSignalSecondPerpendicularConePair) continue;
      if(iPairingType == EECHistograms::kPerpendicularPerpendicularConePair) continue;
      if(iPairingType == EECHistograms::kSecondPerpendicularConePair) continue;
    }
    
    // Find the numbers of tracks for the track loops
    nTracks1 = selectedTrackPt[firstParticleType[iPairingType]].size();
    nTracks2 = selectedTrackPt[secondParticleType[iPairingType]].size();
    
    // Loop over all the tracks in the jet cone
    for(Int_t iFirstTrack = 0; iFirstTrack < nTracks1; iFirstTrack++){
      
      // Get the kinematics for the first track
      trackPt1 = selectedTrackPt[firstParticleType[iPairingType]].at(iFirstTrack);
      trackPhi1 = relativeTrackPhi[firstParticleType[iPairingType]].at(iFirstTrack);
      trackEta1 = relativeTrackEta[firstParticleType[iPairingType]].at(iFirstTrack);
      
      // Get the efficiency correction for the first track
      trackEfficiencyCorrection1 = GetTrackEfficiencyCorrection(trackPt1, trackEta1, hiBin);
      
      // Find the track subevent (only relevant for simulation)
      subeventTrack1 = selectedTrackSubevent[firstParticleType[iPairingType]].at(iFirstTrack);
      
      // Find the starting index for the second track. Avoid double counting of pairs if we pair tracks within the same jet cone.
      // To optimize statistics, loop over all tracks when the second track is taken from reflected jet cone.
      if(firstParticleType[iPairingType] == secondParticleType[iPairingType]){
        startIndex = iFirstTrack+1;
      } else {
        startIndex = 0;
      }
      
      // Loop over the tracks in either the same jet cone or reflected jet cone to create all possible pairs of tracks.
      for(Int_t iSecondTrack = startIndex; iSecondTrack < nTracks2; iSecondTrack++){
        
        // Get the kinematics for the second track
        trackPt2 = selectedTrackPt[secondParticleType[iPairingType]].at(iSecondTrack);
        trackPhi2 = relativeTrackPhi[secondParticleType[iPairingType]].at(iSecondTrack);
        trackEta2 = relativeTrackEta[secondParticleType[iPairingType]].at(iSecondTrack);
        
        // Get the efficiency correction for the second track
        trackEfficiencyCorrection2 = GetTrackEfficiencyCorrection(trackPt2, trackEta2, hiBin);
        
        // Find the track subevent (only relevant for simulation)
        subeventTrack2 = selectedTrackSubevent[secondParticleType[iPairingType]].at(iSecondTrack);
        
        // Find the subevent type based on the subevents of the tracks (only relevant for simulation)
        subeventCombination = GetSubeventCombination(subeventTrack1, subeventTrack2);
        
        // Find the deltaR between the tracks
        trackDeltaR = GetDeltaR(trackEta1, trackPhi1, trackEta2, trackPhi2);
        //trackDeltaR = SimpleSmearDeltaR(trackDeltaR); // Simple smearing for track deltaR
        
        // Find the lower of the two track pT:s
        lowerTrackPt = trackPt1;
        if(trackPt2 < trackPt1) lowerTrackPt = trackPt2;

        // Realistic smearing for track DeltaR if selected
        if(fSmearDeltaR){
          smearedDeltaR = fDeltaRSmearer->GetSmearedValue(trackDeltaR, centrality, jetPt, lowerTrackPt);
          fillerParticleDeltaRResponseMatrix[0] = smearedDeltaR;    // Axis 0: Smeared DeltaR between the two tracks
          fillerParticleDeltaRResponseMatrix[1] = trackDeltaR;      // Axis 1: DeltaR between the two tracks
          fillerParticleDeltaRResponseMatrix[2] = jetPt;            // Axis 2: pT of the jet the tracks are near of
          fillerParticleDeltaRResponseMatrix[3] = lowerTrackPt;     // Axis 3: Lower of the two track pT:s
          fillerParticleDeltaRResponseMatrix[4] = centrality;       // Axis 4: Event centrality
          fHistograms->fhParticleDeltaRResponse->Fill(fillerParticleDeltaRResponseMatrix, fTotalEventWeight);
          trackDeltaR = smearedDeltaR;
        }
        
        // Calculate the weights given to the energy-energy correlators for all defined weight values
        // Factor for smearing study: fRng->Gaus(1,0.0237) for PbPb fRng->Gaus(1,0.0244) for pp (weight = 1)

        weightIndex = 0;
        for(Double_t currentExponent : fWeightExponent){

          correlatorWeight = TMath::Power(trackPt1*trackPt2, currentExponent);

          // Realistic smearing for energy weight if selected
          if(fSmearEnergyWeight){
            smearedCorrelatorWeight = fEnergyWeightSmearer->GetSmearedValue(correlatorWeight, centrality, jetPt, lowerTrackPt);
            fillerParticlePtResponseMatrix[0] = smearedCorrelatorWeight;
            fillerParticlePtResponseMatrix[1] = correlatorWeight;
            fillerParticlePtResponseMatrix[2] = jetPt;
            fillerParticlePtResponseMatrix[3] = lowerTrackPt;
            fillerParticlePtResponseMatrix[4] = centrality;
            fillerParticlePtResponseMatrix[5] = smearedCorrelatorWeight / correlatorWeight;
            fillerParticlePtResponseMatrix[6] = weightIndex;   // Axis 6: Inxed for the current energy-energy correlator weight
            fHistograms->fhParticlePtResponse->Fill(fillerParticlePtResponseMatrix, fTotalEventWeight);
            correlatorWeight = smearedCorrelatorWeight;
          }

          // Find the pair efficiency correction for the track pair
          std::tie(trackPairEfficiencyCorrection, trackPairEfficiencyError) = fTrackPairEfficiencyCorrector->GetTrackPairEfficiencyCorrection(trackDeltaR, centrality, trackPt1, trackPt2, jetPt);

          // Fill the energy-energy correlator histograms
          fillerEnergyEnergyCorrelator[0] = trackDeltaR;               // Axis 0: DeltaR between the two tracks
          fillerEnergyEnergyCorrelator[1] = jetPt;                     // Axis 1: pT of the jet the tracks are near of
          fillerEnergyEnergyCorrelator[2] = lowerTrackPt;              // Axis 2: Lower of the two track pT:s
          fillerEnergyEnergyCorrelator[3] = centrality;                // Axis 3: Event centrality
          fillerEnergyEnergyCorrelator[4] = iPairingType;              // Axis 4: Track pairing type (signal-signal, signal-reflected cone, reflected cone-reflected cone)
          fillerEnergyEnergyCorrelator[5] = subeventCombination;       // Axis 5: Subevent combination type
          fillerEnergyEnergyCorrelator[6] = weightIndex;               // Axis 6: Inxed for the current energy-energy correlator weight

          if(fFillEnergyEnergyCorrelators){
            fHistograms->fhEnergyEnergyCorrelator->Fill(fillerEnergyEnergyCorrelator, trackEfficiencyCorrection1 * trackEfficiencyCorrection2 * fTotalEventWeight * correlatorWeight * trackPairEfficiencyCorrection * jetPtWeight);  // Fill the energy-energy correlator histogram

            // Fill also the histograms just from this jet to calculate covariances
            // Covariance matrix calculation only for data, not needed in MC
            if(!fSkipCovarianceMatrix){
              for(Int_t iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
                if(lowerTrackPt >= fCard->Get("TrackPtBinEdgesEEC", iTrackPt)){
                  fThisEventCorrelator[weightIndex][iTrackPt]->Fill(trackDeltaR, trackEfficiencyCorrection1 * trackEfficiencyCorrection2 * fTotalEventWeight * correlatorWeight * trackPairEfficiencyCorrection * jetPtWeight);
                }
              } // Track pT loop for covariance matrices
            } // Fill coveriance matrices

          } // Filling energy-energy correlators

          if(fFillEnergyEnergyCorrelatorsSystematics){
            // For first systematic variation, increase the track efficiency corrections by the defined amount
            fHistograms->fhEnergyEnergyCorrelatorEfficiencyVariationPlus->Fill(fillerEnergyEnergyCorrelator, (trackEfficiencyCorrection1 + trackEfficiencyCorrection1*fTrackEfficiencyVariation) * (trackEfficiencyCorrection2 + trackEfficiencyCorrection2*fTrackEfficiencyVariation) * fTotalEventWeight * correlatorWeight * trackPairEfficiencyCorrection * jetPtWeight);

            // For second systematic variation, decrease the track efficiency corrections by the defined amount
            fHistograms->fhEnergyEnergyCorrelatorEfficiencyVariationMinus->Fill(fillerEnergyEnergyCorrelator, (trackEfficiencyCorrection1 - trackEfficiencyCorrection1*fTrackEfficiencyVariation) * (trackEfficiencyCorrection2 - trackEfficiencyCorrection2*fTrackEfficiencyVariation) * fTotalEventWeight *  correlatorWeight * trackPairEfficiencyCorrection * jetPtWeight);

            // For the systematic uncertainties vary the track pair inefficiency either by twice the track efficiency uncertainty or by it's error, which one is bigger
            trackPairInefficiency = 1.0 - (1.0 / trackPairEfficiencyCorrection);
            trackPairInefficiencyVariation = (2*fTrackEfficiencyVariation*trackPairInefficiency > trackPairEfficiencyError) ? 2*fTrackEfficiencyVariation*trackPairInefficiency : trackPairEfficiencyError;
            variedTrackPairEfficiencyCorrection = 1.0 / (1.0 - (trackPairInefficiency + trackPairInefficiencyVariation));

            // Fill the histogram with decresed track pair efficiency
            fHistograms->fhEnergyEnergyCorrelatorPairEfficiencyVariationMinus->Fill(fillerEnergyEnergyCorrelator, trackEfficiencyCorrection1 * trackEfficiencyCorrection2 * fTotalEventWeight * correlatorWeight * variedTrackPairEfficiencyCorrection * jetPtWeight); 

            // Next, vary the track pair inefficiency to the other direction
            variedTrackPairEfficiencyCorrection = 1.0 / (1.0 - (trackPairInefficiency - trackPairInefficiencyVariation));

            // Fill the histogarm with increased track pair efficiency
            fHistograms->fhEnergyEnergyCorrelatorPairEfficiencyVariationPlus->Fill(fillerEnergyEnergyCorrelator, trackEfficiencyCorrection1 * trackEfficiencyCorrection2 * fTotalEventWeight * correlatorWeight * variedTrackPairEfficiencyCorrection * jetPtWeight);
          } // Fill systematic uncertainty evaluation correlators

          weightIndex++;
        } // Loop over all defined weight exponents
      } // Inner track loop
    } // Outer track loop
  } // Loop over pairing types (same jet/reflected cone jet)
  
}

/*
 * Fill the jet pT closure histograms
 *
 *  const Int_t jetIndex = Index of a jet for which the closure is filled
 */
void EECAnalyzer::FillJetPtClosureHistograms(const Int_t jetIndex){
  
  // Define a filler for the closure histogram
  const Int_t nAxesClosure = 7;
  Double_t fillerClosure[nAxesClosure];
  
  // Find the pT of the matched gen jet and flavor of reference parton
  Float_t matchedGenPt = fJetReader->GetMatchedPt(jetIndex);  // If not doing jet pT correction manually, change the implementation in GeneratorLevelForestReader!
  Float_t matchedGenEta = fJetReader->GetMatchedEta(jetIndex);
  Float_t matchedGenPhi = fJetReader->GetMatchedPhi(jetIndex);
  Int_t referencePartonFlavor = fJetReader->GetPartonFlavor(jetIndex);
  
  // Find the centrality of the event and the pT of the reconstructed jet
  Double_t recoPt = fJetReader->GetJetRawPt(jetIndex); // Read raw pT and do correction manually
  Double_t centrality = fJetReader->GetCentrality();
  Double_t jetEta = fJetReader->GetJetEta(jetIndex);
  Double_t jetPhi = fJetReader->GetJetPhi(jetIndex);
  
  // If we are using generator level jets, swap reco and gen variables
  if(fMcCorrelationType == kGenReco || fMcCorrelationType == kGenGen){
    Double_t swapper = matchedGenPt;
    matchedGenPt = recoPt;
    recoPt = swapper;
    swapper = matchedGenEta;
    matchedGenEta = jetEta;
    jetEta = swapper;
    swapper = matchedGenPhi;
    matchedGenPhi = jetPhi;
    jetPhi = swapper;
  }
  
  // Helper variable for smearing study
  Double_t smearingFactor;
  
  // For 2018 data, we need to correct the reconstructed pT with jet energy correction
  fJetCorrector->SetJetPT(recoPt);
  fJetCorrector->SetJetEta(jetEta);
  fJetCorrector->SetJetPhi(jetPhi);
  recoPt = fJetCorrector->GetCorrectedPT();
  
  // If we are using smearing scenario, modify the reconstructed jet pT using gaussian smearing
  if(fJetUncertaintyMode > 0){
    smearingFactor = GetSmearingFactor(recoPt, jetEta, centrality);
    recoPt = recoPt * fRng->Gaus(1,smearingFactor);
  }
  
  
  // Define index for parton flavor using algoritm: [-6,-1] U [1,6] -> kQuark, 21 -> kGluon, anything else -> -1
  Int_t referencePartonIndex = -1;
  if(referencePartonFlavor >= -6 && referencePartonFlavor <= 6 && referencePartonFlavor != 0) referencePartonIndex = EECHistograms::kQuark;
  if(referencePartonFlavor == 21) referencePartonIndex = EECHistograms::kGluon;
  
  // Fill the different axes for the filler
  fillerClosure[0] = matchedGenPt;         // Axis 0: pT of the matched generator level jet
  fillerClosure[1] = recoPt;               // Axis 1: pT of the matched reconstructed jet
  fillerClosure[2] = jetEta;               // Axis 2: eta of the jet under consideration
  fillerClosure[3] = centrality;           // Axis 3: Centrality of the event
  fillerClosure[4] = referencePartonIndex; // Axis 4: Reference parton type (quark/gluon)
  fillerClosure[5] = recoPt/matchedGenPt;  // Axis 5: Reconstructed level jet to generator level jet pT ratio
  fillerClosure[6] = jetPhi;               // Axis 6: phi of the jet under consideration
  
  // Fill the closure histogram
  fHistograms->fhJetPtClosure->Fill(fillerClosure,fTotalEventWeight);
  
}


/*
 * Fill all the histograms needed in the unfolding study
 *
 */
void EECAnalyzer::FillUnfoldingResponse(){

  // Variables for jets
  Double_t jetPt;             // Reconstructed jet pT
  Double_t jetPhi;            // Reconstructed jet phi
  Double_t jetEta;            // Reconstructed jet eta
  Double_t genPt;             // Matched generator level jet pT
  Int_t nJets;                // Number of jets
  Double_t smearingFactor;    // Smearing factor for systematic uncertainty study

  // Variables for tracks
  Double_t trackPt;           // Generator level particle pT
  Double_t trackEta;          // Generator level particle eta
  Double_t trackPhi;          // Generator level particle phi
  Int_t nTracks;              // Number of generator level particles
  
  // Variables for energy-energy correlators
  vector<double> selectedTrackPt;     // Track pT for tracks selected for energy-energy correlator analysis
  vector<double> relativeTrackEta;    // Track eta relative to the jet axis
  vector<double> relativeTrackPhi;    // Track phi relative to the jet axis
  vector<int> mathedGenIndices;       // Indices of the gen jets that have been matched with reconstructed jets
  Double_t deltaRTrackJet;            // Delta R between reconstructed jet, and generator level particles

  // Reconstructed jet loop
  nJets = fUnfoldingForestReader->GetNJets();
  mathedGenIndices.clear();
  Double_t centrality = fUnfoldingForestReader->GetCentrality();
  for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++){

    jetPt = fUnfoldingForestReader->GetJetRawPt(jetIndex);  // Get the raw pT and do manual correction later
    jetPhi = fUnfoldingForestReader->GetJetPhi(jetIndex);
    jetEta = fUnfoldingForestReader->GetJetEta(jetIndex);

    //  ========================================
    //  ======== Apply jet quality cuts ========
    //  ========================================

    if(TMath::Abs(jetEta) >= fJetEtaCut) continue;                     // Cut for jet eta
    if(fCutBadPhiRegion && (jetPhi > -0.1 && jetPhi < 1.2)) continue;  // Cut the area of large inefficiency in tracker

    if(fMinimumMaxTrackPtFraction >= fUnfoldingForestReader->GetJetMaxTrackPt(jetIndex) / fUnfoldingForestReader->GetJetRawPt(jetIndex)){
      continue;  // Cut for jets with only very low pT particles
    }
    if(fMaximumMaxTrackPtFraction <= fUnfoldingForestReader->GetJetMaxTrackPt(jetIndex) / fUnfoldingForestReader->GetJetRawPt(jetIndex)){
      continue;  // Cut for jets where all the pT is taken by one track
    }
    

    //  ========================================
    //  ======= Jet quality cuts applied =======
    //  ========================================

    // For 2018 data: do a correction for the jet pT
    fJetCorrector->SetJetPT(jetPt);
    fJetCorrector->SetJetEta(jetEta);
    fJetCorrector->SetJetPhi(jetPhi);

    fJetUncertainty->SetJetPT(jetPt);
    fJetUncertainty->SetJetEta(jetEta);
    fJetUncertainty->SetJetPhi(jetPhi);

    jetPt = fJetCorrector->GetCorrectedPT();

    // If we are using smearing scenario, modify the jet pT using gaussian smearing
    // Note that for MC, we need to smear the jet pT to match the resolution in data
    if(fJetUncertaintyMode > 0) {
      smearingFactor = GetSmearingFactor(jetPt, jetEta, centrality);
      jetPt = jetPt * fRng->Gaus(1, smearingFactor);
    }

    // If we are making runs using variation of jet pT within uncertainties, modify the jet pT here
    // Notice that we still need to use nominal jet pT smearing in MC before applying the shift
    if(fJetUncertaintyMode == 1) jetPt = jetPt * (1 - fJetUncertainty->GetUncertainty().first);
    if(fJetUncertaintyMode == 2) jetPt = jetPt * (1 + fJetUncertainty->GetUncertainty().second);

    // After the jet pT can been corrected, apply analysis jet pT cuts
    if(jetPt < fReconstructedJetMinimumPtCut) continue;
    if(jetPt > fJetMaximumPtCut) continue;

    //************************************************
    //   Do energy-energy correlation within jets
    //************************************************

    // Clear the vectors of track kinematics for tracks selected for energy-energy correlators
    selectedTrackPt.clear();
    relativeTrackEta.clear();
    relativeTrackPhi.clear();

    // For now, only use generator level particles for the unfolding study
    nTracks = fUnfoldingForestReader->GetNGenParticles();
    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

      // Check that all the generator level particle selections are passed
      if(!PassGenParticleSelection(fUnfoldingForestReader, iTrack)) continue;

      // Find the track kinematics
      trackPt = fUnfoldingForestReader->GetGenParticlePt(iTrack);
      trackEta = fUnfoldingForestReader->GetGenParticleEta(iTrack);
      trackPhi = fUnfoldingForestReader->GetGenParticlePhi(iTrack);

      // If the track is close to a jet, change the track eta-phi coordinates to a system where the jet axis is at origin
      deltaRTrackJet = GetDeltaR(jetEta, jetPhi, trackEta, trackPhi);
      if(deltaRTrackJet < fJetRadius) {
        relativeTrackPhi.push_back(trackPhi - jetPhi);
        relativeTrackEta.push_back(trackEta - jetEta);
        selectedTrackPt.push_back(trackPt);
      } // Track close to jet
    } // Track loop

    // Find a matching generator level jet for the reconstructed jet
    if(fUnfoldingForestReader->HasMatchingGenJet(jetIndex)){
      genPt = fUnfoldingForestReader->GetMatchedGenPt(jetIndex);
      mathedGenIndices.push_back(fUnfoldingForestReader->GetMatchingGenIndex(jetIndex));
    } else {
      genPt = -1;
    }

    // Calculate the energy-energy correlator within this jet
    CalculateEnergyEnergyCorrelatorForUnfolding(selectedTrackPt, relativeTrackEta, relativeTrackPhi, jetPt, genPt);

    // Fill also histograms to do jet pT unfolding in one dimension
    FillOneDimensionalJetPtUnfoldingHistograms(jetPt, genPt);

  }  // Reconstructed jet loop

  // Generator level jet loop
  nJets = fUnfoldingForestReader->GetNGeneratorJets();
  for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++){

    // If this index is already filled, do not fill it again
    auto searchResult = std::find(begin(mathedGenIndices), end(mathedGenIndices), jetIndex);
    if(searchResult != std::end(mathedGenIndices)) continue;

    // If the index is not jet filled, fill it to the true distribution
    genPt = fUnfoldingForestReader->GetGeneratorJetPt(jetIndex); 
    jetPhi = fUnfoldingForestReader->GetGeneratorJetPhi(jetIndex);
    jetEta = fUnfoldingForestReader->GetGeneratorJetEta(jetIndex);

    //  ==========================================
    //  ======== Apply jet kinematic cuts ========
    //  ==========================================

    if(TMath::Abs(jetEta) >= fJetEtaCut) continue;                     // Cut for jet eta
    if(fCutBadPhiRegion && (jetPhi > -0.1 && jetPhi < 1.2)) continue;  // Cut the area of large inefficiency in tracker
    if(genPt < fGeneratorJetMinimumPtCut) continue;                    // Cut for minimum jet pT
    if(genPt > fJetMaximumPtCut) continue;                             // Cut for maximum jet pT

    //************************************************
    //   Do energy-energy correlation within jets
    //************************************************

    // Clear the vectors of track kinematics for tracks selected for energy-energy correlators
    selectedTrackPt.clear();
    relativeTrackEta.clear();
    relativeTrackPhi.clear();

    // For now, only use generator level particles for the unfolding study
    nTracks = fUnfoldingForestReader->GetNGenParticles();
    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){

      // Check that all the generator level particle selections are passed
      if(!PassGenParticleSelection(fUnfoldingForestReader, iTrack)) continue;

      // Find the track kinematics
      trackPt = fUnfoldingForestReader->GetGenParticlePt(iTrack);
      trackEta = fUnfoldingForestReader->GetGenParticleEta(iTrack);
      trackPhi = fUnfoldingForestReader->GetGenParticlePhi(iTrack);

      // If the track is close to a jet, change the track eta-phi coordinates to a system where the jet axis is at origin
      deltaRTrackJet = GetDeltaR(jetEta, jetPhi, trackEta, trackPhi);
      if(deltaRTrackJet < fJetRadius) {
        relativeTrackPhi.push_back(trackPhi - jetPhi);
        relativeTrackEta.push_back(trackEta - jetEta);
        selectedTrackPt.push_back(trackPt);
      } // Track close to jet
    } // Track loop

    // Calculate the energy-energy correlator within this jet
    CalculateEnergyEnergyCorrelatorForUnfolding(selectedTrackPt, relativeTrackEta, relativeTrackPhi, -1, genPt);

    // Fill also histograms to do jet pT unfolding in one dimension
    FillOneDimensionalJetPtUnfoldingHistograms(-1, genPt);

  } // Generator level jet loop

}

/*
 * Method for calculating the energy-energy correlators for tracks within the jets for unfolding
 *
 *  const vector<double> selectedTrackPt = pT array for tracks selected for the analysis
 *  const vector<double> relativeTrackEta = relative eta array for tracks selected for the analysis
 *  const vector<double> relativeTrackPhi = relative phi array for tracks selected for the analysis
 *  const double jetPt = pT of the reconstructed jet the tracks are close to
 *  const double genPt = pT of the generator level jet the tracks are close to
 */
void EECAnalyzer::CalculateEnergyEnergyCorrelatorForUnfolding(const vector<double> selectedTrackPt, const vector<double> relativeTrackEta, const vector<double> relativeTrackPhi, const double jetPt, const double genPt){
  
  // Define fillers for THnSparses
  Double_t fillerDistribution[4]; // Axes: deltaR as a function of reco/gen jet pT, track pT cut, centrality, energy weight
  Double_t fillerResponse[5]; // Axes: deltaR as a function of reco jet pT, deltaR as a function of gen jet pT, track pT cut, centrality, energy weight
  Int_t weightIndex;       // Inxed in histogram axis where the weight exponent is saved
  
  // Event information
  Double_t centrality = fUnfoldingForestReader->GetCentrality();
  if(fMultiplicityMode) centrality = GetCentralityFromMultiplicity(GetMultiplicity());
  
  // Variables for generator level particles
  Double_t trackPt1;         // Track pT for the first track
  Double_t trackEta1;        // Track eta for the first track
  Double_t trackPhi1;        // Track phi for the first track
  Double_t trackPt2;         // Track pT for the second track
  Double_t trackEta2;        // Track eta for the second track
  Double_t trackPhi2;        // Track phi for the second track
  Double_t trackDeltaR;      // DeltaR between the two tracks
  Double_t lowerTrackPt;     // Lower of the two track pT:s
  Double_t correlatorWeight; // Weight given to the energy-energy correlator pT1*pT2
  Int_t nTracks;             // Number of tracks

  // Variables for energy-energy correlator response
  Double_t unfoldingDeltaRGeneratorLevel; // DeltaR transformed into generator level unfulding axis
  Double_t unfoldingDeltaRReconstructed;  // DeltaR transformed into reconstructed unfolding axis

  // For MC: reweight only true jet pT
  Double_t jetPtWeightReco = (fJetPtWeightConfiguration == 2) ? 1 : GetJetPtWeight(jetPt);
  Double_t jetPtWeightGen = GetJetPtWeight(genPt);

  // Find the numbers of tracks for the track loops
  nTracks = selectedTrackPt.size();

  // Loop over all the tracks in the jet cone
  for(Int_t iFirstTrack = 0; iFirstTrack < nTracks; iFirstTrack++){

    // Get the kinematics for the first track
    trackPt1 = selectedTrackPt.at(iFirstTrack);
    trackPhi1 = relativeTrackPhi.at(iFirstTrack);
    trackEta1 = relativeTrackEta.at(iFirstTrack);

    // Loop over the all track paris not already studied
    for(Int_t iSecondTrack = iFirstTrack + 1; iSecondTrack < nTracks; iSecondTrack++) {

      // Get the kinematics for the second track
      trackPt2 = selectedTrackPt.at(iSecondTrack);
      trackPhi2 = relativeTrackPhi.at(iSecondTrack);
      trackEta2 = relativeTrackEta.at(iSecondTrack);

      // Find the deltaR between the tracks
      trackDeltaR = GetDeltaR(trackEta1, trackPhi1, trackEta2, trackPhi2); // The actual thingy

      /*trackDeltaR = trackDeltaR * fRng->Gaus(1, 0.2);
      if(trackDeltaR < 0 || trackDeltaR > 0.8){
        trackDeltaR = fRng->Rndm()*0.8; // Random number between 0 and 0.8.
      }*/

      //trackDeltaR = fRng->Rndm()*0.8; // Random number between 0 and 0.8.

      // Find the lower of the two track pT:s
      lowerTrackPt = trackPt1;
      if(trackPt2 < trackPt1) lowerTrackPt = trackPt2;

      // Realistic smearing for DeltaR between two particles if selected
      if(fSmearEnergyWeight){
        trackDeltaR = fDeltaRSmearer->GetSmearedValue(trackDeltaR, centrality, jetPt, lowerTrackPt);
      }

      // Calculate the all the defined weights given to the energy-energy correlators
      weightIndex = 0;
      for(Double_t currentExponent : fWeightExponent){

        correlatorWeight = TMath::Power(trackPt1*trackPt2, currentExponent);

        // Realistic smearing for energy weight is selected
        if(fSmearEnergyWeight){
          correlatorWeight = fEnergyWeightSmearer->GetSmearedValue(correlatorWeight, centrality, jetPt, lowerTrackPt);
        }

        // Transform deltaR into deltaR as a function of reconstructed jet pT
        unfoldingDeltaRReconstructed = TransformToUnfoldingAxis(trackDeltaR, jetPt);

        // Transform deltaR into deltaR as a function of generator level jet pT
        unfoldingDeltaRGeneratorLevel = TransformToUnfoldingAxis(trackDeltaR, genPt);

        // If both reconstructed and generator level jet pT are given, fill the response matrix
        if(jetPt >= fReconstructedJetMinimumPtCut && genPt >= fGeneratorJetMinimumPtCut){
          fillerResponse[0] = unfoldingDeltaRReconstructed;
          fillerResponse[1] = unfoldingDeltaRGeneratorLevel;
          fillerResponse[2] = lowerTrackPt;
          fillerResponse[3] = centrality;
          fillerResponse[4] = weightIndex;
          fHistograms->fhUnfoldingResponse->Fill(fillerResponse, fTotalEventWeight * correlatorWeight * jetPtWeightReco * jetPtWeightGen);
        }

        // If the reconstructed jet pT is given, fill the reconstructed distribution for unfolding
        if(jetPt >= fReconstructedJetMinimumPtCut){
          fillerDistribution[0] = unfoldingDeltaRReconstructed;
          fillerDistribution[1] = lowerTrackPt;
          fillerDistribution[2] = centrality;
          fillerDistribution[3] = weightIndex;
          fHistograms->fhUnfoldingMeasured->Fill(fillerDistribution, fTotalEventWeight * correlatorWeight * jetPtWeightReco);
        }

        // If the generator level jet pT is given, fill the generator level distribution for unfolding
        if(genPt >= fGeneratorJetMinimumPtCut){
          fillerDistribution[0] = unfoldingDeltaRGeneratorLevel;
          fillerDistribution[1] = lowerTrackPt;
          fillerDistribution[2] = centrality;
          fillerDistribution[3] = weightIndex;
          fHistograms->fhUnfoldingTruth->Fill(fillerDistribution, fTotalEventWeight * correlatorWeight * jetPtWeightGen);
        }

        weightIndex++;
      } // Energy weight loop
    } // Inner track loop
  } // Outer track loop
}

/*
 * Method for filling histograms for one dimensional jet pT unfolding
 *
 *  const double jetPt = pT of the reconstructed jet 
 *  const double genPt = pT of the generator level jet
 */
void EECAnalyzer::FillOneDimensionalJetPtUnfoldingHistograms(const double jetPt, const double genPt){

  // Define fillers for THnSparses
  Double_t fillerDistribution[2]; // Axes: deltaR as a function of reco/gen jet pT, centrality
  Double_t fillerResponse[3]; // Axes: reco jet pT, gen jet pT, centrality
  
  // Event information
  Double_t centrality = fUnfoldingForestReader->GetCentrality();
  if(fMultiplicityMode) centrality = GetCentralityFromMultiplicity(GetMultiplicity());

  // For MC: reweight only true jet pT
  Double_t jetPtWeightReco = (fJetPtWeightConfiguration == 2) ? 1 : GetJetPtWeight(jetPt);
  Double_t jetPtWeightGen = GetJetPtWeight(genPt);

  // If both reconstructed and generator level jet pT are given, fill the response matrix
  if(jetPt >= fReconstructedJetMinimumPtCut && genPt >= fGeneratorJetMinimumPtCut){
    fillerResponse[0] = jetPt;
    fillerResponse[1] = genPt;
    fillerResponse[2] = centrality;
    fHistograms->fhJetPtUnfoldingResponse->Fill(fillerResponse, fTotalEventWeight * jetPtWeightReco * jetPtWeightGen);
  }

  // If the reconstructed jet pT is given, fill the reconstructed distribution for unfolding
  if(jetPt >= fReconstructedJetMinimumPtCut){
    fillerDistribution[0] = jetPt;
    fillerDistribution[1] = centrality;
    fHistograms->fhJetPtUnfoldingMeasured->Fill(fillerDistribution, fTotalEventWeight * jetPtWeightReco);
  }

  // If the generator level jet pT is given, fill the generator level distribution for unfolding
  if(genPt >= fGeneratorJetMinimumPtCut){
    fillerDistribution[0] = genPt;
    fillerDistribution[1] = centrality;
    fHistograms->fhJetPtUnfoldingTruth->Fill(fillerDistribution, fTotalEventWeight * jetPtWeightGen);
  }
}

/*
 * Construct responses for DeltaR and (pT1*pT2)^n
 *
 */
void EECAnalyzer::ConstructParticleResponses(){

  // Variables for jets
  Double_t jetPt;             // Reconstructed jet pT
  Double_t jetPhi;            // Reconstructed jet phi
  Double_t jetEta;            // Reconstructed jet eta
  Int_t nJets;                // Number of jets
  Double_t smearingFactor;    // Smearing factor for systematic uncertainty study

  // Variables for tracks
  Double_t trackPt;                       // Track or generator level particle pT
  Double_t trackEta;                      // Track or generator level particle eta
  Double_t trackPhi;                      // Track or generator level particle phi
  Double_t trackEfficiencyCorrection;     // Track efficiency correction
  Double_t trackPairEfficiencyCorrection; // Track pair efficiency correction
  Double_t trackPairEfficiencyError;      // Track pair efficiency error
  Int_t trackCharge;                      // Charge of the particle
  Int_t nTracks;                          // Number of generator level particles
  Bool_t trackPairCorrectionFlag;         // Flag about enabling track pair efficiency correction

  // Event variables
  Int_t hiBin = fUnfoldingForestReader->GetHiBin();
  Double_t centrality = fUnfoldingForestReader->GetCentrality();
  Int_t weightIndex;
  
  // Variables for tracks that are matched
  vector<std::tuple<double,double,double,double,int,int>> selectedTrackInformation; 
  vector<std::tuple<double,double,double,double,int,int>> selectedParticleInformation; 
  vector<std::tuple<double,int>> possibleMatches;
  Double_t deltaTrackPt;     // Difference between two track pT:s
  Double_t deltaRTrackJet;   // DeltaR between reconstructed jet, and generator level particles
  Double_t deltaRTracks;     // DeltaR between two tracks
  Double_t deltaRParticles;  // DeltaR between two particles
  Double_t trackMomentumProduct;    // pT1*pT2 for tracks
  Double_t particleMomentumProduct; // pT1*pT2 for particles

  // Filler for test
  double fillerParticleMatching[4];
  double fillerParticleDeltaRResponseMatrix[5];
  double fillerParticlePtResponseMatrix[7];

  // Reconstructed jet loop
  nJets = fUnfoldingForestReader->GetNJets();
  for(Int_t jetIndex = 0; jetIndex < nJets; jetIndex++){

    jetPt = fUnfoldingForestReader->GetJetRawPt(jetIndex);  // Get the raw pT and do manual correction later
    jetPhi = fUnfoldingForestReader->GetJetPhi(jetIndex);
    jetEta = fUnfoldingForestReader->GetJetEta(jetIndex);

    //  ========================================
    //  ======== Apply jet quality cuts ========
    //  ========================================

    if(TMath::Abs(jetEta) >= fJetEtaCut) continue;                     // Cut for jet eta
    if(fCutBadPhiRegion && (jetPhi > -0.1 && jetPhi < 1.2)) continue;  // Cut the area of large inefficiency in tracker

    if(fMinimumMaxTrackPtFraction >= fUnfoldingForestReader->GetJetMaxTrackPt(jetIndex) / fUnfoldingForestReader->GetJetRawPt(jetIndex)){
      continue;  // Cut for jets with only very low pT particles
    }
    if(fMaximumMaxTrackPtFraction <= fUnfoldingForestReader->GetJetMaxTrackPt(jetIndex) / fUnfoldingForestReader->GetJetRawPt(jetIndex)){
      continue;  // Cut for jets where all the pT is taken by one track
    }
    

    //  ========================================
    //  ======= Jet quality cuts applied =======
    //  ========================================

    // For 2018 data: do a correction for the jet pT
    fJetCorrector->SetJetPT(jetPt);
    fJetCorrector->SetJetEta(jetEta);
    fJetCorrector->SetJetPhi(jetPhi);

    fJetUncertainty->SetJetPT(jetPt);
    fJetUncertainty->SetJetEta(jetEta);
    fJetUncertainty->SetJetPhi(jetPhi);

    jetPt = fJetCorrector->GetCorrectedPT();

    // If we are using smearing scenario, modify the jet pT using gaussian smearing
    if(fJetUncertaintyMode > 0) {
      smearingFactor = GetSmearingFactor(jetPt, jetEta, centrality);
      jetPt = jetPt * fRng->Gaus(1, smearingFactor);
    }

    // If we are making runs using variation of jet pT within uncertainties, modify the jet pT here
    // Notice that we still need to use nominal jet pT smearing in MC before applying the shift
    if(fJetUncertaintyMode == 1) jetPt = jetPt * (1 - fJetUncertainty->GetUncertainty().first);
    if(fJetUncertaintyMode == 2) jetPt = jetPt * (1 + fJetUncertainty->GetUncertainty().second);

    // After the jet pT can been corrected, apply analysis jet pT cuts
    if(jetPt < fReconstructedJetMinimumPtCut) continue;
    if(jetPt > fJetMaximumPtCut) continue;

    //************************************************
    //   Do energy-energy correlation within jets
    //************************************************

    // Clear the vectors of track kinematics for tracks selected for energy-energy correlators
    selectedTrackInformation.clear();
    selectedParticleInformation.clear();

    // Loop over generator level particles and find all particles close to the jets
    nTracks = fUnfoldingForestReader->GetNTracks();
    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

      // Check that all the track cuts are passed
      if(!PassTrackCuts(fUnfoldingForestReader, iTrack)) continue;

      // Find the track kinematics
      trackPt = fUnfoldingForestReader->GetTrackPt(iTrack);
      trackEta = fUnfoldingForestReader->GetTrackEta(iTrack);
      trackPhi = fUnfoldingForestReader->GetTrackPhi(iTrack);
      trackEfficiencyCorrection = GetTrackEfficiencyCorrection(trackPt, trackEta, hiBin);
      trackCharge = fUnfoldingForestReader->GetTrackCharge(iTrack);

      // If the track is close to a jet, change the track eta-phi coordinates to a system where the jet axis is at origin
      deltaRTrackJet = GetDeltaR(jetEta, jetPhi, trackEta, trackPhi);
      if(deltaRTrackJet < fJetRadius && trackPt > fTrackMinPtCut){
        selectedTrackInformation.push_back(std::make_tuple(trackPt, trackEta - jetEta, trackPhi - jetPhi, trackEfficiencyCorrection, trackCharge, -1));
      } // Track close to jet
    } // Track loop

    // Loop over generator level particles and find all particles close to the jets
    nTracks = fUnfoldingForestReader->GetNGenParticles();
    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

      // Check that all the generator level particle selections are passed
      if(!PassGenParticleSelection(fUnfoldingForestReader, iTrack)) continue;

      // Find the track kinematics
      trackPt = fUnfoldingForestReader->GetGenParticlePt(iTrack);
      trackEta = fUnfoldingForestReader->GetGenParticleEta(iTrack);
      trackPhi = fUnfoldingForestReader->GetGenParticlePhi(iTrack);
      trackCharge = fUnfoldingForestReader->GetGenParticleCharge(iTrack);

      // If the track is close to a jet, change the track eta-phi coordinates to a system where the jet axis is at origin
      deltaRTrackJet = GetDeltaR(jetEta, jetPhi, trackEta, trackPhi);
      if(deltaRTrackJet < fJetRadius && trackPt > fTrackMinPtCut) {
        selectedParticleInformation.push_back(std::make_tuple(trackPt, trackEta - jetEta, trackPhi - jetPhi, 1, trackCharge, -1));
      } // Track close to jet
    } // Track loop

    // For each reconstructed track, calculate the number of generator level particles within Delta R<0.05
    // Use a 2 GeV pT cut for reconstructed tracks and 1 GeV cut for generator level particles
    int counter;
    fillerParticleMatching[1] = jetPt;
    fillerParticleMatching[3] = centrality;
    for(Int_t iTrack = 0; iTrack < selectedTrackInformation.size(); iTrack++){
      counter = 0;
      for(Int_t iParticle = 0; iParticle < selectedParticleInformation.size(); iParticle++){
        if(GetDeltaR(std::get<kTrackEta>(selectedTrackInformation.at(iTrack)), std::get<kTrackPhi>(selectedTrackInformation.at(iTrack)), std::get<kTrackEta>(selectedParticleInformation.at(iParticle)), std::get<kTrackPhi>(selectedParticleInformation.at(iParticle))) < 0.05){
          if(std::get<kTrackCharge>(selectedTrackInformation.at(iTrack)) == std::get<kTrackCharge>(selectedParticleInformation.at(iParticle))){
            if(TMath::Abs(std::get<kTrackPt>(selectedParticleInformation.at(iParticle))-std::get<kTrackPt>(selectedTrackInformation.at(iTrack))) < 0.5*std::get<kTrackPt>(selectedTrackInformation.at(iTrack))){
              counter++;
            }
          } 
        } 
      }
      fillerParticleMatching[0] = counter;
      fillerParticleMatching[2] = std::get<kTrackPt>(selectedTrackInformation.at(iTrack));
      fHistograms->fhParticlesCloseToTracks->Fill(fillerParticleMatching);
    }

    // Implement potential matching algorithm and see how many matches are found
    

    // Loop over the reconstructed track starting from the track with highest pT
    std::sort(selectedTrackInformation.begin(), selectedTrackInformation.end(), std::greater<std::tuple<double,double,double,double,int,int>>());
    for(Int_t iTrack = 0; iTrack < selectedTrackInformation.size(); iTrack++){
      possibleMatches.clear();
      for(Int_t iParticle = 0; iParticle < selectedParticleInformation.size(); iParticle++){

        // First, check that the generator level particle is not already matched to any reconstructed track
        if(std::get<kMatchIndex>(selectedParticleInformation.at(iParticle)) < 0) {

          // Then check that the generator level particle is within 0.05 from the reconstructed particle
          if(GetDeltaR(std::get<kTrackEta>(selectedTrackInformation.at(iTrack)), std::get<kTrackPhi>(selectedTrackInformation.at(iTrack)), std::get<kTrackEta>(selectedParticleInformation.at(iParticle)), std::get<kTrackPhi>(selectedParticleInformation.at(iParticle))) < 0.05) {

            // Then, check the charges match between reconstructed and generator level particles
            if(std::get<kTrackCharge>(selectedTrackInformation.at(iTrack)) == std::get<kTrackCharge>(selectedParticleInformation.at(iParticle))) {

              // Then, require the the generator level particle pT is within 50% of the reconstructed particle pT
              deltaTrackPt = TMath::Abs(std::get<kTrackPt>(selectedParticleInformation.at(iParticle)) - std::get<kTrackPt>(selectedTrackInformation.at(iTrack)));
              if(deltaTrackPt < 0.5 * std::get<kTrackPt>(selectedTrackInformation.at(iTrack))) {

                // If all the conditions are fulfilled, we have a potential match
                possibleMatches.push_back(std::make_tuple(deltaTrackPt, iParticle));
              }
            }
          } 
        } // Check that generator level particle has not already been matched
      } // Generator level particle loop

      // If there is more than one match, sort the possible matches based on the difference in pT
      if(possibleMatches.size() > 1) std::sort(possibleMatches.begin(), possibleMatches.end());

      // If there is at least one match, set the matching indices for the track and the particle
      if(possibleMatches.size() > 0){
        std::get<kMatchIndex>(selectedTrackInformation.at(iTrack)) = std::get<kPossibleMatchIndex>(possibleMatches.at(0));
        std::get<kMatchIndex>(selectedParticleInformation.at(std::get<kPossibleMatchIndex>(possibleMatches.at(0)))) = iTrack;
      }
    }

    // For each reconstructed track, fill a histogram telling if a matching generator level particle was found or not
    fillerParticleMatching[1] = jetPt;
    fillerParticleMatching[3] = centrality;
    for(Int_t iTrack = 0; iTrack < selectedTrackInformation.size(); iTrack++){
      fillerParticleMatching[0] = std::get<kMatchIndex>(selectedTrackInformation.at(iTrack)) < 0 ? 0 : 1;
      fillerParticleMatching[2] = std::get<kTrackPt>(selectedTrackInformation.at(iTrack));
      fHistograms->fhTracksWithMatchedParticle->Fill(fillerParticleMatching);
    }

    // After the matching has been done, we can construct the response matrices
    for(Int_t iTrack = 0; iTrack < selectedTrackInformation.size(); iTrack++){

      // Only loop over tracks that have matched generator level particle
      if(std::get<kMatchIndex>(selectedTrackInformation.at(iTrack)) < 0) continue;

      // Second track loop to make the track pairs
      for(Int_t jTrack = iTrack+1; jTrack < selectedTrackInformation.size(); jTrack++){

        // Only loop over tracks that have matched generator level particle
        if(std::get<kMatchIndex>(selectedTrackInformation.at(jTrack)) < 0) continue;

        // Calculate the deltaR from tracks and from the matched particles
        deltaRTracks = GetDeltaR(std::get<kTrackEta>(selectedTrackInformation.at(iTrack)), std::get<kTrackPhi>(selectedTrackInformation.at(iTrack)), std::get<kTrackEta>(selectedTrackInformation.at(jTrack)), std::get<kTrackPhi>(selectedTrackInformation.at(jTrack)));
        deltaRParticles = GetDeltaR(std::get<kTrackEta>(selectedParticleInformation.at(std::get<kMatchIndex>(selectedTrackInformation.at(iTrack)))), std::get<kTrackPhi>(selectedParticleInformation.at(std::get<kMatchIndex>(selectedTrackInformation.at(iTrack)))), std::get<kTrackEta>(selectedParticleInformation.at(std::get<kMatchIndex>(selectedTrackInformation.at(jTrack)))), std::get<kTrackPhi>(selectedParticleInformation.at(std::get<kMatchIndex>(selectedTrackInformation.at(jTrack)))));

        // Find the pair acceptance correction for the track pair regardless if it is disabled in the main code or not
        trackPairCorrectionFlag = fTrackPairEfficiencyCorrector->GetDisableCorrection();
        fTrackPairEfficiencyCorrector->SetDisableCorrection(false);
        std::tie(trackPairEfficiencyCorrection, trackPairEfficiencyError) = fTrackPairEfficiencyCorrector->GetTrackPairEfficiencyCorrection(deltaRTracks, centrality, std::get<kTrackPt>(selectedTrackInformation.at(iTrack)), std::get<kTrackPt>(selectedTrackInformation.at(jTrack)), jetPt);
        fTrackPairEfficiencyCorrector->SetDisableCorrection(trackPairCorrectionFlag);

        // Fill the matched deltaR values to the response matrix
        fillerParticleDeltaRResponseMatrix[0] = deltaRTracks;
        fillerParticleDeltaRResponseMatrix[1] = deltaRParticles;
        fillerParticleDeltaRResponseMatrix[2] = jetPt;
        fillerParticleDeltaRResponseMatrix[3] = std::get<kTrackPt>(selectedTrackInformation.at(jTrack));  // Since the vector is sorted, the smaller pT is at higher index
        fillerParticleDeltaRResponseMatrix[4] = centrality;
        fHistograms->fhParticleDeltaRResponse->Fill(fillerParticleDeltaRResponseMatrix, fTotalEventWeight * std::get<kTrackEfficiencyCorrection>(selectedTrackInformation.at(iTrack)) * std::get<kTrackEfficiencyCorrection>(selectedTrackInformation.at(jTrack)) * trackPairEfficiencyCorrection);

        // Fill the matched pT1*pT2 values to the response matrix
        trackMomentumProduct = std::get<kTrackPt>(selectedTrackInformation.at(iTrack)) * std::get<kTrackPt>(selectedTrackInformation.at(jTrack));
        particleMomentumProduct = std::get<kTrackPt>(selectedParticleInformation.at(std::get<kMatchIndex>(selectedTrackInformation.at(iTrack)))) * std::get<kTrackPt>(selectedParticleInformation.at(std::get<kMatchIndex>(selectedTrackInformation.at(jTrack))));

        // Use the defined weight exponents for the track momentum product
        weightIndex = 0;
        for(Double_t currentExponent : fWeightExponent){

          trackMomentumProduct = TMath::Power(trackMomentumProduct, currentExponent);
          particleMomentumProduct = TMath::Power(particleMomentumProduct, currentExponent);

          fillerParticlePtResponseMatrix[0] = trackMomentumProduct;
          fillerParticlePtResponseMatrix[1] = particleMomentumProduct;
          fillerParticlePtResponseMatrix[2] = jetPt;
          fillerParticlePtResponseMatrix[3] = std::get<kTrackPt>(selectedTrackInformation.at(jTrack));  // Since the vector is sorted, the smaller pT is at higher index
          fillerParticlePtResponseMatrix[4] = centrality;
          fillerParticlePtResponseMatrix[5] = trackMomentumProduct / particleMomentumProduct;
          fillerParticlePtResponseMatrix[6] = weightIndex;
          fHistograms->fhParticlePtResponse->Fill(fillerParticlePtResponseMatrix, fTotalEventWeight * std::get<kTrackEfficiencyCorrection>(selectedTrackInformation.at(iTrack)) * std::get<kTrackEfficiencyCorrection>(selectedTrackInformation.at(jTrack)) * trackPairEfficiencyCorrection);

          weightIndex++;
        } // Loop over defined weight exponents
      } // Inner track loop
    } // Outer track loop


  }  // Reconstructed jet loop

}

/*
 * Get the proper vz weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t vz = Vertex z position for the event
 *
 *   return: Multiplicative correction factor for vz
 */
Double_t EECAnalyzer::GetVzWeight(const Double_t vz) const{
  if(fIsRealData) return 1;  // No correction for real data
  return fVzWeightFunction->Eval(vz); // Weight for Monte Carlo
}

/*
 * Get the proper centrality weighting depending on analyzed system
 *
 *  Arguments:
 *   const Int_t hiBin = CMS hiBin
 *
 *   return: Multiplicative correction factor for the given CMS hiBin
 */
Double_t EECAnalyzer::GetCentralityWeight(const Int_t hiBin) const{
  if(fDataType != ForestReader::kPbPbMC) return 1;
  
  // No weighting for the most peripheral centrality bins. Different weight function for central and peripheral.
  if(hiBin < 60) return fCentralityWeightFunctionCentral->Eval(hiBin/2.0);
  return (hiBin < 194) ? fCentralityWeightFunctionPeripheral->Eval(hiBin/2.0) : 1;
}

/*
 *  Get the proper multiplicity weight for MC
 *
 *  Arguments:
 *   const Double_t multiplicity = Track multiplicity in the event
 *
 *   return: Multiplicative correction factor for the given multiplicity value
 */
Double_t EECAnalyzer::GetMultiplicityWeight(const Double_t multiplicity) const{
  if(fDataType != ForestReader::kPbPbMC) return 1;
  
  return fMultiplicityWeightFunction->Eval(multiplicity);
}

/*
 * Get the proper jet pT weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t jetPt = Jet pT for the weighted jet
 *
 *   return: Multiplicative correction factor for the jet pT
 */
Double_t EECAnalyzer::GetJetPtWeight(const Double_t jetPt) const{
  if(fJetPtWeightConfiguration == 0) return 1.0; // If jet weighting is disabled, apply no weight
  if(fIsRealData) return 1.0;  // No weight for data

  return fPtWeightFunction->Eval(jetPt);
}

/*
 * Get a smearing factor corresponding to worsening the smearing resolution in MC by 20 %
 * This is obtained by multiplying the MC smearing resolution by 0.666 and using this as additional
 * smearing for the data. Smearing factor depends on jet pT and centrality.
 *
 *  Arguments:
 *   Double_t jetPt = Jet pT
 *   const Double_t centrality = Centrality of the event
 *
 *  return: Additional smearing factor
 */
Double_t EECAnalyzer::GetSmearingFactor(Double_t jetPt, Double_t jetEta, const Double_t centrality) {
  
  // For all the jets above 500 GeV, use the resolution for 500 GeV jet
  if(jetPt > 500) jetPt = 500;
  
  // Find the correct centrality bin
  Int_t centralityBin = GetCentralityBin(centrality);
  
  // Set the parameters to the smearing function. pp and PbPb have different smearing function parameters
  if(fIsPpData){
    // Settings for pp
    // Determined using the macro constructJetPtClosures.C
    // Input file: ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noCorrelations_jetPtClosures_processed_2023-01-13.root
    fSmearingFunction->SetParameters(0.165659, -0.000775499, 2.70493e-06, -4.69846e-09, 3.15964e-12);

  } else {
    
    // Parameters for the smearing function
    // Determined using the macro constructJetPtClosures.C
    // Input file: PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtClosure_finalMcWeight_processed_2023-03-06.root
    Double_t resolutionFit[4][5] = {
      {0.424589, -0.00260826, 8.06713e-06, -1.17528e-08, 6.49404e-12},
      {0.415669, -0.00300066, 1.18443e-05, -2.30594e-08, 1.74213e-11},
      {0.358008, -0.00270528, 1.10689e-05, -2.15565e-08, 1.59491e-11},
      {0.237325, -0.00138461, 4.77259e-06, -7.80495e-09, 4.79538e-12}
    };
    
    for(int iParameter = 0; iParameter < 5; iParameter++){
      // Settings for PbPb
      
      fSmearingFunction->SetParameter(iParameter, resolutionFit[centralityBin][iParameter]);
    }
  }
  
  // Calculation for resolution worsening: we assume the jet energy resolution is a Gaussian distribution with some certain sigma, if you would like to add a Gaussian noise to make it worse, the sigma getting larger, then it obeys the random variable rule that X=Y+Z, where Y~N(y, sigmay) and Z~N(z,sigmaz), then X~N(y+z, sqrt(sigmay^2+sigmaz^2))). In this case, we assume that noise and the resolution are independent.
  // So let assume the sigmay is the jet energy resolution, then you want the sigmax = 1.2sigmay
  // which means that the sigmaz = sigmay * sqrt(1.2^2-1)
  
  // After the smearing function is set, read the value to return
  // Worsening resolution by 20%: 0.663
  // Worsening resolution by 10%: 0.458
  // Worsening resolution by 30%: 0.831

  // We want to worsen resolution in MC by the amount defined by JetMet group. The scaling factor is given by a JetMet manager
  return fSmearingFunction->Eval(jetPt)*fEnergyResolutionSmearingFinder->GetScalingFactor(jetEta);
  
}

/*
 * Check is a track passes the required subevent cut
 *
 *  Arguments:
 *   const Int_t subeventIndex = Subevent index for the track in consideration
 *
 *  return: true if subevent cut is passes, false if not
 */
Bool_t EECAnalyzer::PassSubeventCut(const Int_t subeventIndex) const{
  if(fSubeventCut == kSubeventAny) return true;
  if((fSubeventCut == kSubeventZero) && (subeventIndex == 0)) return true;
  if((fSubeventCut == kSubeventNonZero) && (subeventIndex > 0)) return true;
  return false;
}

/*
 * Check if the event passes all the track cuts
 *
 *  Arguments:
 *   ForestReader* eventReader = ForestReader containing the event information checked for event cuts
 *   const Bool_t fillHistograms = Flag for filling the event information histograms.
 *
 *   return = True if all event cuts are passed, false otherwise
 */
Bool_t EECAnalyzer::PassEventCuts(ForestReader* eventReader, const Bool_t fillHistograms){
  
  // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction. Only applied for data.
  if(eventReader->GetPrimaryVertexFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(EECHistograms::kPrimaryVertex);
  
  // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV. Only applied for PbPb data.
  if(eventReader->GetHfCoincidenceFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(EECHistograms::kHfCoincidence);
  
  // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible. Only applied for PbPb data.
  if(eventReader->GetClusterCompatibilityFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(EECHistograms::kClusterCompatibility);
  
  // Cut for HB/HE noise. Only applied for data.
  if(eventReader->GetHBHENoiseFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(EECHistograms::kHBHENoise);
  
  // Cut for beam scraping. Only applied for pp data.
  if(eventReader->GetBeamScrapingFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(EECHistograms::kBeamScraping);

  // Cut for pile-up.
  if(eventReader->GetPileupFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(EECHistograms::kPileup);
  
  // Cut for vertex z-position
  if(TMath::Abs(eventReader->GetVz()) > fVzCut) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(EECHistograms::kVzCut);
  
  return true;
  
}

/*
 * Check if a track passes all the track cuts
 *
 *  Arguments:
 *   ForestReader* trackReader = ForestReader from which the tracks are read
 *   const Int_t iTrack = Index of the checked track in reader
 *   TH1F* trackCutHistogram = Histogram to which the track cut performance is filled
 *   const Bool_t bypassFill = Pass filling the track cut histograms
 *
 *   return: True if all track cuts are passed, false otherwise
 */
Bool_t EECAnalyzer::PassTrackCuts(ForestReader* trackReader, const Int_t iTrack, TH1F* trackCutHistogram, const Bool_t bypassFill){
  
  // Only fill the track cut histograms for same event data
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kAllTracks);
  
  // Cuts specific to generator level MC tracks
  if(trackReader->GetTrackCharge(iTrack) == 0) return false;  // Require that the track is charged
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kMcCharge);
  
  if(!PassSubeventCut(trackReader->GetTrackSubevent(iTrack))) return false;  // Require desired subevent
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kMcSube);
  
  if(trackReader->GetTrackMCStatus(iTrack) != 1) return false;  // Require final state particles
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kMcStatus);
  
  Double_t trackPt = trackReader->GetTrackPt(iTrack);
  Double_t trackEta = trackReader->GetTrackEta(iTrack);
  Double_t trackEt = (trackReader->GetTrackEnergyEcal(iTrack)+trackReader->GetTrackEnergyHcal(iTrack))/TMath::CosH(trackEta);
  
  //  ==== Apply cuts for tracks and collect information on how much track are cut in each step ====
  
  // Cut for track pT
  if(trackPt <= fTrackMinPtCut) return false;                   // Minimum pT cut
  if(trackPt >= fTrackMaxPtCut) return false;                   // Maximum pT cut
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kPtCuts);
  
  // Cut for track eta
  if(TMath::Abs(trackEta) >= fTrackEtaCut) return false;          // Eta cut
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kEtaCut);
  
  // New cut for 2018 data based on track algorithm and MVA
  if(fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPbPbMC){
    if(trackReader->GetTrackAlgorithm(iTrack) == 6){
      if(trackReader->GetTrackMVA(iTrack) < 0.98){
        return false;
      }
    }
  }
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kTrackAlgorithm);
  
  // Cut for high purity
  if(!trackReader->GetTrackHighPurity(iTrack)) return false;     // High purity cut
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kHighPurity);
  
  // Cut for relative error for track pT
  if(trackReader->GetTrackPtError(iTrack)/trackPt >= fMaxTrackPtRelativeError) return false; // Cut for track pT relative error
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kPtError);
  
  // Cut for track distance from primary vertex
  if(TMath::Abs(trackReader->GetTrackVertexDistanceZ(iTrack)/trackReader->GetTrackVertexDistanceZError(iTrack)) >= fMaxTrackDistanceToVertex) return false; // Mysterious cut about track proximity to vertex in z-direction
  if(TMath::Abs(trackReader->GetTrackVertexDistanceXY(iTrack)/trackReader->GetTrackVertexDistanceXYError(iTrack)) >= fMaxTrackDistanceToVertex) return false; // Mysterious cut about track proximity to vertex in xy-direction
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kVertexDistance);

  // The final cuts are applied only for PbPb data
  if(fDataType != ForestReader::kPbPb && fDataType != ForestReader::kPbPbMC){
    if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kCaloSignal);
    if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kReconstructionQuality);
    return true;
  }
  
  // Cut for energy deposition in calorimeters for high pT tracks
  if(!(trackPt < fCalorimeterSignalLimitPt || (trackEt >= fHighPtEtFraction*trackPt))) return false;  // For high pT tracks, require signal also in calorimeters
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kCaloSignal);
  
  // Cuts for track reconstruction quality
  //if( trackReader->GetTrackChi2(iTrack) / (1.0*trackReader->GetNTrackDegreesOfFreedom(iTrack)) / (1.0*trackReader->GetNHitsTrackerLayer(iTrack)) >= fChi2QualityCut) return false; // Track reconstruction quality cut
  if(trackReader->GetTrackNormalizedChi2(iTrack) / (1.0*trackReader->GetNHitsTrackerLayer(iTrack)) >= fChi2QualityCut) return false; // Track reconstruction quality cut
  if(trackReader->GetNHitsTrack(iTrack) < fMinimumTrackHits) return false; // Cut for minimum number of hits per track
  if(fFillTrackHistograms && !bypassFill) trackCutHistogram->Fill(EECHistograms::kReconstructionQuality);
  
  // If passed all checks, return true
  return true;
}

/*
 * Check if a track passes all the track cuts
 *
 *  Arguments:
 *   UnfoldingForestReader* trackReader = UnfoldingForestReader from which the tracks are read
 *   const Int_t iTrack = Index of the checked track in reader
 *
 *   return: True if all track cuts are passed, false otherwise
 */
Bool_t EECAnalyzer::PassTrackCuts(UnfoldingForestReader* trackReader, const Int_t iTrack){
  
  Double_t trackPt = trackReader->GetTrackPt(iTrack);
  Double_t trackEta = trackReader->GetTrackEta(iTrack);
  Double_t trackEt = (trackReader->GetTrackEnergyEcal(iTrack)+trackReader->GetTrackEnergyHcal(iTrack))/TMath::CosH(trackEta);
  
  //  ==== Apply cuts for tracks and collect information on how much track are cut in each step ====
  
  // Cut for track pT
  if(trackPt <= fTrackMinPtCut) return false;                   // Minimum pT cut
  if(trackPt >= fTrackMaxPtCut) return false;                   // Maximum pT cut
  
  // Cut for track eta
  if(TMath::Abs(trackEta) >= fTrackEtaCut) return false;          // Eta cut
  
  // Cut for high purity
  if(!trackReader->GetTrackHighPurity(iTrack)) return false;     // High purity cut
  
  // Cut for relative error for track pT
  if(trackReader->GetTrackPtError(iTrack)/trackPt >= fMaxTrackPtRelativeError) return false; // Cut for track pT relative error
  
  // Cut for track distance from primary vertex
  if(TMath::Abs(trackReader->GetTrackVertexDistanceZ(iTrack)/trackReader->GetTrackVertexDistanceZError(iTrack)) >= fMaxTrackDistanceToVertex) return false; // Mysterious cut about track proximity to vertex in z-direction
  if(TMath::Abs(trackReader->GetTrackVertexDistanceXY(iTrack)/trackReader->GetTrackVertexDistanceXYError(iTrack)) >= fMaxTrackDistanceToVertex) return false; // Mysterious cut about track proximity to vertex in xy-direction
  
  // Cut for energy deposition in calorimeters for high pT tracks
  if(!(trackPt < fCalorimeterSignalLimitPt || (trackEt >= fHighPtEtFraction*trackPt))) return false;  // For high pT tracks, require signal also in calorimeters
  
  // Cuts for track reconstruction quality
  if(trackReader->GetTrackNormalizedChi2(iTrack) / (1.0*trackReader->GetNHitsTrackerLayer(iTrack)) >= fChi2QualityCut) return false; // Track reconstruction quality cut
  if(trackReader->GetNHitsTrack(iTrack) < fMinimumTrackHits) return false; // Cut for minimum number of hits per track
  
  // If passed all checks, return true
  return true;
}

/*
 * Check if a track passes the generator level particle selection criteria
 *
 *  Arguments:
 *   ForestReader *trackReader = ForestReader from which the generator level particles are read
 *   const Int_t iTrack = Index of the checked track in reader
 *
 *   return: True if all track cuts are passed, false otherwise
 */
Bool_t EECAnalyzer::PassGenParticleSelection(ForestReader* trackReader, const Int_t iTrack){

  // Cuts specific to generator level MC tracks
  if(trackReader->GetParticleCharge(iTrack) == 0) return false;  // Require that the track is charged

  if(!PassSubeventCut(trackReader->GetParticleSubevent(iTrack))) return false;  // Require desired subevent

  Double_t trackPt = trackReader->GetParticlePt(iTrack);

  // Cut for track pT
  if(trackPt <= fTrackMinPtCut) return false;       // Minimum track pT cut
  if(trackPt >= fTrackMaxPtCut) return false;       // Maximum track pT cut

  Double_t trackEta = trackReader->GetParticleEta(iTrack);

  // Cut for track eta
  if(TMath::Abs(trackEta) >= fTrackEtaCut) return false;          // Eta cut

  // If passed all checks, return true
  return true;

}

/*
 * Check if a track passes the generator level particle selection criteria
 *
 *  Arguments:
 *   UnfoldingForestReader *trackReader = UnfoldingForestReader from which the generator level particles are read
 *   const Int_t iTrack = Index of the checked track in reader
 *
 *   return: True if all track cuts are passed, false otherwise
 */
Bool_t EECAnalyzer::PassGenParticleSelection(UnfoldingForestReader* trackReader, const Int_t iTrack){
  
  // Cuts specific to generator level MC tracks
  if(trackReader->GetGenParticleCharge(iTrack) == 0) return false;  // Require that the track is charged
  
  if(!PassSubeventCut(trackReader->GetGenParticleSubevent(iTrack))) return false;  // Require desired subevent
  
  Double_t trackPt = trackReader->GetGenParticlePt(iTrack);
  
  // Cut for track pT
  if(trackPt <= fTrackMinPtCut) return false;       // Minimum track pT cut
  if(trackPt >= fTrackMaxPtCut) return false;       // Maximum track pT cut
  
  Double_t trackEta = trackReader->GetGenParticleEta(iTrack);
  
  // Cut for track eta
  if(TMath::Abs(trackEta) >= fTrackEtaCut) return false;          // Eta cut
  
  // If passed all checks, return true
  return true;
  
}


/*
 * Get the track efficiency correction for a given track
 *
 *  Arguments:
 *   const Int_t iTrack = Index of the track for which the efficiency correction is obtained
 *
 *   return: Multiplicative track efficiency correction
 */
Double_t EECAnalyzer::GetTrackEfficiencyCorrection(const Int_t iTrack){
  
  // No correction for generator level tracks
  if(fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenGen) return 1;
  if(fIsOOData) return 1.0; // Force no correction for OO data, as instructed by tracking group for now. This is because the tracking efficiency correction for OO data is not yet available.
  
  // Get track information
  Float_t trackPt = fTrackReader->GetTrackPt(iTrack);    // Track pT
  Float_t trackEta = fTrackReader->GetTrackEta(iTrack);  // Track eta
  Int_t hiBin = fTrackReader->GetHiBin();                // hiBin for 2018 track correction
  
  // Get the correction using the track and event information
  return GetTrackEfficiencyCorrection(trackPt, trackEta, hiBin);
  
}

/*
 * Get the track efficiency correction according to given information
 *
 *  Arguments:
 *   const Float_t trackPt = pT of the track
 *   const Float_t trackEta = Eta of the track
 *   const Int_t hiBin = CMS hiBin index. Basically centrality/2
 *
 *   return: Multiplicative track efficiency correction
 */
Double_t EECAnalyzer::GetTrackEfficiencyCorrection(const Float_t trackPt, const Float_t trackEta, const Int_t hiBin){
  
  // No correction for generator level tracks
  if(fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenGen) return 1;
  if(fDataType == ForestReader::kOO || fDataType == ForestReader::kOOMC) return 1;   // No correction for OO dataset, as instructed by tracking group for now. This is because the tracking efficiency correction for OO data is not yet available.


  
  // Weight factor only for 2017 pp MC as instructed be the tracking group
  double preWeight = 1.0;
  if(fDataType == ForestReader::kPpMC) preWeight = 0.979;
  
  // For PbPb2018 and pp2017, there is an efficiency table from which the correction comes
  return preWeight * fTrackEfficiencyCorrector->getCorrection(trackPt, trackEta, hiBin);
  
}

/*
 * Getter for EEC histograms
 */
EECHistograms* EECAnalyzer::GetHistograms() const{
  return fHistograms;
}

/*
 * Getter for centrality bin
 */
Int_t EECAnalyzer::GetCentralityBin(const Double_t centrality) const{
  
  // Find the correct centrality bin
  Int_t centralityBin = 0;
  for(int iCentrality = 1; iCentrality < fCard->GetNBin("CentralityBinEdges"); iCentrality++){
    if(centrality > fCard->Get("CentralityBinEdges",iCentrality)) centralityBin++;
  }

  return centralityBin;
}

/*
 * Get the track multiplicity in the current event
 */
Double_t EECAnalyzer::GetMultiplicity(){
  
  // Loop over all track in the event
  Int_t nTracks = fTrackReader->GetNTracks();
  Double_t trackMultiplicity = 0;
  Double_t trackEfficiencyCorrection = 0;
  
  // Disable subevent cut while determining the total multiplicity
  Int_t originalCut = fSubeventCut;
  fSubeventCut = kSubeventAny;
  
  for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
    
    // Check that all the track cuts are passed
    if(!PassTrackCuts(fTrackReader,iTrack,fHistograms->fhTrackCuts,true)) continue;
    
    // Get the efficiency correction
    trackEfficiencyCorrection = GetTrackEfficiencyCorrection(iTrack);
    
    //trackMultiplicity += 1;
    trackMultiplicity += trackEfficiencyCorrection;
    
  } // Track loop
  
  fSubeventCut = originalCut;
  
  return trackMultiplicity;
  
}

/*
 * Get the generator level multiplicity
 *
 * Argument:
 *  ForestReader* trackReader = Reader used to determine the multiplicity
 *  const Int_t subeventIndex = Subevents used to calculate the multiplicity
 *  const bool inMegaSkim = Multiplicity is asked for a mega skimmed file
 */
Int_t EECAnalyzer::GetGenMultiplicity(ForestReader* trackReader, const Int_t subeventIndex, const Bool_t isMegaSkim){

  // Find the number of particles
  Int_t nParticles = trackReader->GetNParticles();
  Int_t particleMultiplicity = 0;

  // Select the desired subevent cut, that might be different from general selection
  Int_t originalCut = fSubeventCut;
  if(subeventIndex >= 0){
    fSubeventCut = subeventIndex;
  }

  // Loop over all particles in the event
  for(Int_t iParticle = 0; iParticle < nParticles; iParticle++){

    if(isMegaSkim){

      // For mega skims, the kinematic selection is already done. We just need to pass the subevent cut
      if(!PassSubeventCut(trackReader->GetParticleSubevent(iParticle))) continue;

      // If we pass, increment the counter!
      particleMultiplicity++;

    } else {
  
      // Check that we are in the kinematic region of interest
      if(!PassGenParticleSelection(trackReader, iParticle)) continue;

      // If we are, increment the counter!
      particleMultiplicity++;

    }

  } // Particle loop

  // Set the subevent cut back to original value
  fSubeventCut = originalCut;

  // Return the determined multiplicity
  return particleMultiplicity;


}

/*
 * Get the analysis centrality bin corresponding to the given multiplicity value
 *
 *  Arguments:
 *   const Double_t multiplicity = Multiplicity in MC
 *
 *   return: Centrality corresponding to this value of multiplicity
 */
Double_t EECAnalyzer::GetCentralityFromMultiplicity(const Double_t multiplicity) const{
  
  // Centrality bin 0-10
  if(multiplicity > 2225) return 7;
  
  // Centrality bin 10-30
  if(multiplicity > 980) return 17;
  
  // Centrality bin 30-50
  if(multiplicity > 340) return 37;
  
  // Centrality bin 50-90
  return 57;
  
}

/*
 * Get deltaR between two objects
 *
 *  Arguments:
 *   const Double_t eta1 = Eta of the first object
 *   const Double_t phi1 = Phi of the first object
 *   const Double_t eta2 = Eta of the second object
 *   const Double_t phi2 = Phi of the second object
 *
 *  return: DeltaR between the two objects
 */
Double_t EECAnalyzer::GetDeltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2) const{
  
  Double_t deltaEta = eta1 - eta2;
  Double_t deltaPhi = phi1 - phi2;
  
  // Transform deltaPhi to interval [-pi,pi]
  while(deltaPhi > TMath::Pi()){deltaPhi += -2*TMath::Pi();}
  while(deltaPhi < -TMath::Pi()){deltaPhi += 2*TMath::Pi();}
  
  // Return the distance between the objects
  return TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
  
}

/*
 * Get the subevent combination type from two track subevents
 *
 *  Arguments:
 *   const Int_t subevent1 = Subevent index for the first track (0 = pythia, > 0 = hydjet)
 *   const Int_t subevent2 = Subevent index for the second track (0 = pythia, > 0 = hydjet)
 *
 * return: Subevent combination type: 0 = pythia-pythia, 1 = pythia-hydjet, 2 = hydjet-pythia, 3 = hydjet-hydjet
 */
Int_t EECAnalyzer::GetSubeventCombination(const Int_t subevent1, const Int_t subevent2) const{
  
  // For data or reconstructed tracks, just return -1. It will go to the underflow bin of the histogram
  if(subevent1 < 0 || subevent2 < 0) return -1;
  
  // If both tracks are from the pythia simulation, return flag the corresponding flag
  if(subevent1 == 0 && subevent2 == 0) return EECHistograms::kPythiaPythia;
  
  // If both tracks are from hydjet simulation, return the corresponding flag
  if(subevent1 > 0 && subevent2 > 0) return EECHistograms::kHydjetHydjet;
  
  // If the first track is from pythia and second from hydjet, return the corresponding flag
  if(subevent1 == 0 && subevent2 > 0) return EECHistograms::kPythiaHydjet;
  
  // The only option left is that the first track is from hydjet, and the second from pythia. Return the corresponding flag
  return EECHistograms::kHydjetPythia;
  
}

/*
 * Get the subevent index for a track
 *
 *  Arguments:
 *   const Int_t subevent = Subevent index for the track (0 = pythia, > 0 = hydjet)
 *
 * return: Subevent type: 0 = Pythia, 1 = Hydjet
 */
Int_t EECAnalyzer::GetSubeventIndex(const Int_t subevent) const{
  
  // For data or reconstructed tracks, just return -1. It will go to the underflow bin of the histogram
  if(subevent < 0) return -1;
  
  // If both tracks are from the pythia simulation, return 0
  if(subevent == 0) return 0;
  
  // If the particle is not from Pythia, it must be from Hydjet. Return 1.
  return 1;
  
}

/*
 * Get jet eta reflected around zero, avoiding overlapping jet cones
 *
 *  Arguments:
 *   const Double_t eta = Eta of the jet
 *
 * return: Eta of a reflected jet axis
 */
Double_t EECAnalyzer::GetReflectedEta(const Double_t eta) const{
  
  // If the jet eta is far enough from zero that the cones do not overlap, do a direct reflection
  if(TMath::Abs(eta) > fJetRadius) return -eta;
  
  // If eta is negative, add twice the jet radius to the eta to avoid overlapping cones
  if(eta < 0) return eta + 2*fJetRadius;
  
  // Now eta must be positive. Subtract twice the jet radius from the value to avoid overlapping cones
  return eta - 2*fJetRadius;
}

/*
 * Get perpendicular jet phi coordinate in a defined direction
 *
 *  Arguments:
 *   const Double_t phi = Eta of the jet
 *   const Int_t direction = 0: Positive phi direction, 1: Negative phi direction
 *
 * return: Perpendicular phi coordinate
 */
Double_t EECAnalyzer::GetPerpendicularPhi(const Double_t phi, const Int_t direction) const{

  // Make sure direction index is 0 or 1. Give a snarky comment if not.
  if(direction < 0 || direction > 1){
    std::cout << "EECAnalyzer::ERROR! Direction index in GetPerpendicularPhi out of bounds!" << std::endl;
    std::cout << "Fix the error in your code. Now it crashes because of your mistakes..." << std::endl;
    assert(0);
  }

  // Move the phi by 90 degrees to the defined direction
  Double_t adder[2] = {TMath::Pi()/2.0, -TMath::Pi()/2.0};
  Double_t perpendicularPhi = phi + adder[direction];

  // Make sure that the reflected phi index is in the range [-pi, pi]
  if(perpendicularPhi < TMath::Pi()) perpendicularPhi += 2*TMath::Pi();
  if(perpendicularPhi > TMath::Pi()) perpendicularPhi -= 2*TMath::Pi();

  // Return perpendicularized phi
  return perpendicularPhi;

}

/*
 * Transform the deltaR value to the unfolding axis
 *
 *  Arguments:
 *   const Double_t deltaR = DeltaR between the two particles
 *   const Double_t jetPt = pT of the jets the track pair is close to
 *
 * return: The new value transformed into the main unfolding axis
 */
Double_t EECAnalyzer::TransformToUnfoldingAxis(const Double_t deltaR, const Double_t jetPt) const{

  const Int_t nJetPtBinsEEC = fCard->GetNBin("JetPtBinEdgesEEC");
  const Double_t maxDeltaR = 0.8;
  Double_t transformedDeltaR = deltaR;
  for(Int_t iJetPt = 1; iJetPt < nJetPtBinsEEC+1; iJetPt++){
    if(jetPt >= fCard->Get("JetPtBinEdgesEEC",iJetPt)){
      transformedDeltaR += maxDeltaR;
    } else {
      return transformedDeltaR;
    }
  }

  // We should never reach this point. If we are here, just return error code -1
  return -1;
}

/*
 * Simple smearing for deltaR to see if it affects the final distributions
 *
 *  Arguments:
 *   const Double_t deltaR = Original DeltaR
 *
 *  return: Smeared deltaR value
 *
 */
Double_t EECAnalyzer::SimpleSmearDeltaR(const Double_t deltaR){

  double randomNumber = fRng->Rndm();

  // Go through the deltaR histogram bin by bin and perform smearing of the deltaR value
  if(deltaR < 0.00147203) return deltaR;

  // Low underflow bin
  if(deltaR < 0.00316074){
    if(randomNumber > 0.99) return 0.006;
    return deltaR;
  }

  // First underflow bin
  if(deltaR < 0.00509804){
    if(randomNumber > 0.91) return 0.006;
    return deltaR;
  }

  // Analysis bin 1
  if(deltaR < 0.00732051){
    if(randomNumber < 0.12) return 0.005;
    if(randomNumber > 0.88) return 0.008;
    return deltaR;
  }

  // Analysis bin 2
  if(deltaR < 0.00987013){
    if(randomNumber < 0.09) return 0.007;
    if(randomNumber > 0.89) return 0.01;
    return deltaR;
  }

  // Analysis bin 3
  if(deltaR < 0.0127951){
    if(randomNumber < 0.08) return 0.008;
    if(randomNumber > 0.91) return 0.013;
    return deltaR;
  }

  // Analysis bin 4
  if(deltaR < 0.0161506){
    if(randomNumber < 0.065) return 0.011;
    if(randomNumber > 0.925) return 0.017;
    return deltaR;
  }

  // Analysis bin 5
  if(deltaR < 0.02){
    if(randomNumber < 0.05) return 0.013;
    if(randomNumber > 0.93) return 0.021;
    return deltaR;
  }

  // Analysis bin 6
  if(deltaR < 0.0244161){
    if(randomNumber < 0.05) return 0.019;
    if(randomNumber > 0.94) return 0.025;
    return deltaR;
  }

  // Analysis bin 7
  if(deltaR < 0.0294822){
    if(randomNumber < 0.05) return 0.023;
    if(randomNumber > 0.95) return 0.03;
    return deltaR;
  }

  // Analysis bin 8
  if(deltaR < 0.0352941){
    if(randomNumber < 0.04) return 0.028;
    if(randomNumber > 0.95) return 0.036;
    return deltaR;
  }

  // Analysis bin 9
  if(deltaR < 0.0419615){
    if(randomNumber < 0.04) return 0.034;
    if(randomNumber > 0.96) return 0.042;
    return deltaR;
  }

  // Analysis bin 10
  if(deltaR < 0.0496104){
    if(randomNumber < 0.03) return 0.04;
    if(randomNumber > 0.96) return 0.05;
    return deltaR;
  }

  // Analysis bin 11
  if(deltaR < 0.0583852){
    if(randomNumber < 0.03) return 0.048;
    if(randomNumber > 0.96) return 0.06;
    return deltaR;
  }

  // Analysis bin 12
  if(deltaR < 0.0684517){
    if(randomNumber < 0.03) return 0.057;
    if(randomNumber > 0.97) return 0.07;
    return deltaR;
  }

  // Analysis bin 13
  if(deltaR < 0.08){
    if(randomNumber < 0.02) return 0.0666;
    if(randomNumber > 0.97) return 0.085;
    return deltaR;
  }

  // Analysis bin 14
  if(deltaR < 0.0932482){
    if(randomNumber < 0.02) return 0.07;
    if(randomNumber > 0.97) return 0.094;
    return deltaR;
  }

  // Analysis bin 15
  if(deltaR < 0.108447){
    if(randomNumber < 0.02) return 0.09;
    if(randomNumber > 0.98) return 0.11;
    return deltaR;
  }

  // Analysis bin 16
  if(deltaR < 0.125882){
    if(randomNumber < 0.02) return 0.1;
    if(randomNumber > 0.98) return 0.13;
    return deltaR;
  }

  // Analysis bin 17
  if(deltaR < 0.145885){
    if(randomNumber < 0.02) return 0.12;
    if(randomNumber > 0.98) return 0.15;
    return deltaR;
  }

  // Analysis bin 18
  if(deltaR < 0.168831){
    if(randomNumber < 0.01) return 0.14;
    if(randomNumber > 0.98) return 0.17;
    return deltaR;
  }

  // Analysis bin 19
  if(deltaR < 0.195156){
    if(randomNumber < 0.01) return 0.16;
    if(randomNumber > 0.98) return 0.2;
    return deltaR;
  }

  // Analysis bin 20
  if(deltaR < 0.225355){
    if(randomNumber < 0.01) return 0.19;
    if(randomNumber > 0.98) return 0.23;
    return deltaR;
  }

  // Analysis bin 21
  if(deltaR < 0.26){
    if(randomNumber < 0.01) return 0.22;
    if(randomNumber > 0.99) return 0.27;
    return deltaR;
  }

  // Analysis bin 22
  if(deltaR < 0.299745){
    if(randomNumber < 0.01) return 0.25;
    if(randomNumber > 0.99) return 0.3;
    return deltaR;
  }

  // Analysis bin 23
  if(deltaR < 0.34534){
    if(randomNumber < 0.01) return 0.29;
    if(randomNumber > 0.99) return 0.35;
    return deltaR;
  }

  // Analysis bin 24
  if(deltaR < 0.397647){
    if(randomNumber < 0.005) return 0.34;
    if(randomNumber > 0.995) return 0.4;
    return deltaR;
  }

  // Analysis bin 25
  if(deltaR < 0.457654){
    if(randomNumber < 0.005) return 0.39;
    if(randomNumber > 0.995) return 0.46;
    return deltaR;
  }

  return deltaR;

}

/*
 * Find matching events for the current event based on multiplicity
 *
 *  std::vector<int>& mixedEventIndices = Reference to a vector that will hold the indiced of matched events
 *  const Int_t nEventsToMatch = Number of matched events that are searched for
 *  const Double_t vz = Vertex z-position of the event
 *  const Int_t currentMultiplicity = Multiplicity of the event. HiBin in case of PbPb collisions
 *  const Int_t iEvent = Number of the current event
 *  const Bool_t binnedMixing = True: Use binned matching, False: Use unbinned matching
 */
void EECAnalyzer::FindMultiplicityMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, const Int_t currentMultiplicity, const Int_t eventNumber, const Int_t iEvent, const Bool_t binnedMixing){

  // If we are doing binned matching to find mixed events, use different method for that
  if(binnedMixing){
    FindBinnedMultiplicityMatchedEvents(mixedEventIndices, nEventsToMatch, vz, currentMultiplicity, eventNumber, iEvent);
    return;
  }

  // Event matching for Monte Carlo
  Int_t mixedEventIndex = fMixingStartIndex;     // Current index for mixed event file
  Bool_t allEventsWentThrough  = false;          // Flag telling if we have gone through all the events in mixing file without finding a matching event
  Double_t vzTolerance = 0.5;                    // Starting tolerance for vz position: 0.5 cm
  Int_t lowMultiplicityTolerance = 1;            // Multiplicity tolerance for events with low multiplicity start from 1
  Double_t multiplicityTolerance = 0.01;         // Percentage difference in multiplicity values that we allow between current and mixed events, starts from 1%.
  Int_t mixedEventFillIndex = 0;                 // Index telling how many mixed events have been found already
  Bool_t skipEvent = false;                      // Flag for skipping unwanted miked events
  Bool_t checkForDuplicates = false;             // Flag for checking duplicate events in the mixing list

  mixedEventIndices.clear();                     // Reset the previously found event index vector from previous events

  // Find events to mix until we have found defined number of matching events
  while(mixedEventFillIndex < nEventsToMatch){

    // By default, do not skip events
    skipEvent = false;

    // If we have checked all the events but not found an event to mix with, increase vz and hiBin tolerance
    if(allEventsWentThrough){
      if(fDebugLevel > 0){
        cout << "Could not find matching mixed events for event " << iEvent << endl;
        cout << "Increasing vz tolerance by 0.5 and multiplicity tolerance by 0.01" << endl;
      }
      
      multiplicityTolerance += 0.01;
      lowMultiplicityTolerance++;
      vzTolerance += 0.5;
      allEventsWentThrough = false;
      checkForDuplicates = true;
    }
    
    // Increment the counter for event index to be mixed with the current event
    mixedEventIndex++;
    
    // If we are out of bounds from the event in data file, reset the counter
    if(mixedEventIndex == fnEventsInMixingFile) {
      mixedEventIndex = -1;
      continue;
    }

    // If we come back to the first event index, we have gone through all the events without finding a matching mixed event
    if(mixedEventIndex == fMixingStartIndex) allEventsWentThrough = true;

    // Do not mix with events used in the previous iteration over the file
    if(checkForDuplicates){
      for(Int_t iMixedEvent : mixedEventIndices) {
        if(iMixedEvent == mixedEventIndex) skipEvent = true;
      }
    }

    // If we have already used the current mixed event, do not use the same event
    if(skipEvent) continue;

    // Match vz between the current event and mixed event
    if(TMath::Abs(fMixedEventVz.at(mixedEventIndex) - vz) > (vzTolerance + 1e-4)) continue;

    // Match multiplicity between the current event and mixed event
    if(currentMultiplicity > 100){
      if(TMath::Abs(fMixedEventMultiplicity.at(mixedEventIndex) - currentMultiplicity) > (currentMultiplicity * multiplicityTolerance)) continue;
    } else {
      if(TMath::Abs(fMixedEventMultiplicity.at(mixedEventIndex) - currentMultiplicity) > lowMultiplicityTolerance) continue;
    }

    // Avoid mixing an event with itself
    if(fMixedEventEventNumber.at(mixedEventIndex) == eventNumber) continue;

    // Mark that this event is good to use for mixing
    mixedEventIndices.push_back(mixedEventIndex);

    // Increment the index of event that we have found for mixing purposes
    mixedEventFillIndex++;
      
  } // While loop for searching for a good event to mix with

  // For the next event, randomize the starting event again
  fMixingStartIndex = fRng->Integer(fnEventsInMixingFile);  // Start mixing from random spot in file

}

/*
 * Find matching events for the current event based on hiBin
 *
 *  std::vector<int>& mixedEventIndices = Reference to a vector that will hold the indiced of matched events
 *  const Int_t nEventsToMatch = Number of matched events that are searched for
 *  const Double_t vz = Vertex z-position of the event
 *  const Int_t hiBin = HiBin index of the current event
 *  const Int_t iEvent = Number of the current event
 *  const Bool_t binnedMixing = True: Use binned matching, False: Use unbinned matching
 */ 
void EECAnalyzer::FindHiBinMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, const Int_t hiBin, const Int_t eventNumber, const Int_t iEvent, const Bool_t binnedMixing){

  // If we are doing binned matching to find mixed events, use different method for that
  if(binnedMixing){
    FindBinnedHiBinMatchedEvents(mixedEventIndices, nEventsToMatch, vz, hiBin, eventNumber, iEvent);
    return;
  }

  // Event matching for PbPb
  Int_t mixedEventIndex = fMixingStartIndex;     // Current index for mixed event file
  Bool_t allEventsWentThrough  = false;          // Flag telling if we have gone through all the events in mixing file without finding a matching event
  Double_t vzTolerance = 0.5;                    // Starting tolerance for vz position: 0.5 cm
  Int_t hiBinTolerance = 0;                      // Starting tolerance for hiBin: match the bin value
  Int_t mixedEventFillIndex = 0;                 // Index telling how many mixed events have been found already
  Bool_t skipEvent = false;                      // Flag for skipping unwanted miked events
  Bool_t checkForDuplicates = false;             // Flag for checking duplicate events in the mixing list

  mixedEventIndices.clear();  // Reset the previously found event index vector from previous events

  // Find events to mix until we have found defined number of matching events
  while(mixedEventFillIndex < nEventsToMatch){

    // By default, do not skip events
    skipEvent = false;

    // If we have checked all the events but not found an event to mix with, increase vz and hiBin tolerance
    if(allEventsWentThrough){
      if(fDebugLevel > 0){
        cout << "Could not find matching mixed events for event " << iEvent << endl;
        cout << "Increasing vz tolerance by 0.5 and hiBin tolerance by 1" << endl;
      }
      
      hiBinTolerance++;
      vzTolerance += 0.5;
      allEventsWentThrough = false;
      checkForDuplicates = true;
    }
    
    // Increment the counter for event index to be mixed with the current event
    mixedEventIndex++;
    
    // If we are out of bounds from the event in data file, reset the counter
    if(mixedEventIndex == fnEventsInMixingFile) {
      mixedEventIndex = -1;
      continue;
    }

    // If we come back to the first event index, we have gone through all the events without finding a matching mixed event
    if(mixedEventIndex == fMixingStartIndex) allEventsWentThrough = true;

    // Do not mix with events used in the previous iteration over the file
    if(checkForDuplicates){
      for(Int_t iMixedEvent : mixedEventIndices) {
        if(iMixedEvent == mixedEventIndex) skipEvent = true;
      }
    }

    // If we have already used the current mixed event, do not use the same event
    if(skipEvent) continue;

    // Match vz between the current event and mixed event
    if(TMath::Abs(fMixedEventVz.at(mixedEventIndex) - vz) > (vzTolerance + 1e-4)) continue;

    // Match hiBin between the current eent and mixed event
    if(TMath::Abs(fMixedEventHiBin.at(mixedEventIndex) - hiBin) > (hiBinTolerance + 1e-4)) continue;

    // Avoid mixing an event with itself
    if(fMixedEventEventNumber.at(mixedEventIndex) == eventNumber) continue;

    // Mark that this event is good to use for mixing
    mixedEventIndices.push_back(mixedEventIndex);

    // Increment the index of event that we have found for mixing purposes
    mixedEventFillIndex++;
      
  } // While loop for searching for a good event to mix with

  // For the next event, randomize the starting event again
  fMixingStartIndex = fRng->Integer(fnEventsInMixingFile);  // Start mixing from random spot in file

}

/*
 * Find matching events for the current event based on HF energy
 *
 *  std::vector<int>& mixedEventIndices = Reference to a vector that will hold the indiced of matched events
 *  const Int_t nEventsToMatch = Number of matched events that are searched for
 *  const Double_t vz = Vertex z-position of the event
 *  const Double_t hfEnergy = Energy in the HF calorimeters
 *  const Int_t eventNumber = Number of the current event
 *  const Int_t iEvent = Loop index of the current event
 *  const Bool_t binnedMixing = True: Use binned matching, False: Use unbinned matching
 */
void EECAnalyzer::FindHFEnergyMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, Double_t hfEnergy, const Int_t eventNumber, const Int_t iEvent, const Bool_t binnedMixing){

  // If we are doing binned matching to find mixed events, use different method for that
  if(binnedMixing){
    FindBinnedHFEnergyMatchedEvents(mixedEventIndices, nEventsToMatch, vz, hfEnergy, eventNumber, iEvent);
    return;
  }

  // Event matching for pPb
  Int_t mixedEventIndex = fMixingStartIndex;     // Current index for mixed event file
  Bool_t allEventsWentThrough  = false;          // Flag telling if we have gone through all the events in mixing file without finding a matching event
  Double_t vzTolerance = 0.5;                    // Starting tolerance for vz position: 0.5 cm
  Double_t lowHFenergyTolerance = 1;             // HF energy tolerance for events with low HF energy start from 1
  Double_t hfEnergyTolerance = 0.01;             // Percentage difference in HF energy values that we allow between current and mixed events, starts from 1%.
  Int_t mixedEventFillIndex = 0;                 // Index telling how many mixed events have been found already
  Bool_t skipEvent = false;                      // Flag for skipping unwanted miked events
  Bool_t checkForDuplicates = false;             // Flag for checking duplicate events in the mixing list

  mixedEventIndices.clear();  // Reset the previously found event index vector from previous events

  // Shift the hfEnergy in the data event to remove the contribution from the hard scattering
  if(hfEnergy < 2*fHFEnergyMatchingShift){
    hfEnergy /= 2;
  } else {
    hfEnergy -= fHFEnergyMatchingShift;
  }
  if(hfEnergy < 2) hfEnergy = 2; // TODO: There are very few events with energy below 5 units. This might need adjustment...

  // Find events to mix until we have found defined number of matching events
  while(mixedEventFillIndex < nEventsToMatch){

    // By default, do not skip events
    skipEvent = false;

    // If we have checked all the events but not found an event to mix with, increase vz and hiBin tolerance
    if(allEventsWentThrough){
      if(fDebugLevel > 0){
        cout << "Could not find matching mixed events for event " << iEvent << endl;
        cout << "Increasing vz tolerance by 0.5 and HF energy tolerance by 1%" << endl;
      }
      
      hfEnergyTolerance += 0.01;
      lowHFenergyTolerance++;
      vzTolerance += 0.5;
      allEventsWentThrough = false;
      checkForDuplicates = true;
    }
    
    // Increment the counter for event index to be mixed with the current event
    mixedEventIndex++;

    // If we are out of bounds from the event in data file, reset the counter
    if(mixedEventIndex == fnEventsInMixingFile) {
      mixedEventIndex = -1;
      continue;
    }

    // If we come back to the first event index, we have gone through all the events without finding a matching mixed event
    if(mixedEventIndex == fMixingStartIndex) allEventsWentThrough = true;

    // Do not mix with events used in the previous iteration over the file
    if(checkForDuplicates){
      for(Int_t iMixedEvent : mixedEventIndices) {
        if(iMixedEvent == mixedEventIndex) skipEvent = true;
      }
    }

    // If we have already used the current mixed event, do not use the same event
    if(skipEvent) continue;

    // Match vz between the current event and mixed event
    if(TMath::Abs(fMixedEventVz.at(mixedEventIndex) - vz) > (vzTolerance + 1e-4)) continue;

    // Match HF energy between the current event and mixed event
    if(hfEnergy > 100){
      if(TMath::Abs(fMixedEventHFEnergy.at(mixedEventIndex) - hfEnergy) > (hfEnergy * hfEnergyTolerance)) continue;
    } else {
      if(TMath::Abs(fMixedEventHFEnergy.at(mixedEventIndex) - hfEnergy) > lowHFenergyTolerance) continue;
    }

    // We are using minimum bias dataset for both mixing and analysis. Make sure not to match with the same event!
    if(fMixedEventEventNumber.at(mixedEventIndex) == eventNumber) continue;

    // Mark that this event is good to use for mixing
    mixedEventIndices.push_back(mixedEventIndex);

    // Increment the index of event that we have found for mixing purposes
    mixedEventFillIndex++;
      
  } // While loop for searching for a good event to mix with

  // For the next event, randomize the starting event again
  fMixingStartIndex = fRng->Integer(fnEventsInMixingFile);  // Start mixing from random spot in file

}

/*
 * Find matching events for the current event based on multiplicity bins
 *
 *  std::vector<int>& mixedEventIndices = Reference to a vector that will hold the indiced of matched events
 *  const Int_t nEventsToMatch = Number of matched events that are searched for
 *  const Double_t vz = Vertex z-position of the event
 *  const Int_t currentMultiplicity = Multiplicity of the event. HiBin in case of PbPb collisions
 *  const Int_t iEvent = Number of the current event
 *  
 */
void EECAnalyzer::FindBinnedMultiplicityMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, const Int_t currentMultiplicity, const Int_t eventNumber, const Int_t iEvent){


  // Event matching based on multiplicity bins
  Int_t mixedEventIndex = fMixingStartIndex;     // Current index for mixed event file
  Bool_t allEventsWentThrough  = false;          // Flag telling if we have gone through all the events in mixing file without finding a matching event
  Int_t vzTolerance = 0;                         // Starting tolerance for vz bin: exact match
  Int_t multiplicityTolerance = 0;               // Starting tolerance for multiplicity bin: exact match
  Int_t mixedEventFillIndex = 0;                 // Index telling how many mixed events have been found already
  Bool_t skipEvent = false;                      // Flag for skipping unwanted miked events
  Bool_t checkForDuplicates = false;             // Flag for checking duplicate events in the mixing list
  Int_t vzBin = FindMixingVzBin(vz);             // Mixing bin index for the vz in the current event
  Int_t multiplicityBin = FindMixingMultiplicityBin(currentMultiplicity); // Mixing bin index for the multiplicity in the current event

  mixedEventIndices.clear();                     // Reset the previously found event index vector from previous events

  // Find events to mix until we have found defined number of matching events
  while(mixedEventFillIndex < nEventsToMatch){

    // By default, do not skip events
    skipEvent = false;

    // If we have checked all the events but not found an event to mix with, increase vz and hiBin tolerance
    if(allEventsWentThrough){
      if(fDebugLevel > 0){
        cout << "Could not find matching mixed events for event " << iEvent << endl;
        cout << "Increasing vz tolerance by 0.5 and multiplicity tolerance by 0.01" << endl;
      }
      
      multiplicityTolerance++;
      vzTolerance++;
      allEventsWentThrough = false;
      checkForDuplicates = true;
    }
    
    // Increment the counter for event index to be mixed with the current event
    mixedEventIndex++;
    
    // If we are out of bounds from the event in data file, reset the counter
    if(mixedEventIndex == fnEventsInMixingFile) {
      mixedEventIndex = -1;
      continue;
    }

    // If we come back to the first event index, we have gone through all the events without finding a matching mixed event
    if(mixedEventIndex == fMixingStartIndex) allEventsWentThrough = true;

    // Do not mix with events used in the previous iteration over the file
    if(checkForDuplicates){
      for(Int_t iMixedEvent : mixedEventIndices) {
        if(iMixedEvent == mixedEventIndex) skipEvent = true;
      }
    }

    // If we have already used the current mixed event, do not use the same event
    if(skipEvent) continue;

    // Match vz bin between the current event and mixed event
    if(TMath::Abs(fMixedEventVzBin.at(mixedEventIndex) - vzBin) > (vzTolerance + 1e-4)) continue;

    // Match multiplicity bin between the current event and mixed event
    if(TMath::Abs(fMixedEventMultiplicityBin.at(mixedEventIndex) - multiplicityBin) > (multiplicityTolerance + 1e-4)) continue;

    // Avoid mixing an event with itself
    if(fMixedEventEventNumber.at(mixedEventIndex) == eventNumber) continue;

    // Mark that this event is good to use for mixing
    mixedEventIndices.push_back(mixedEventIndex);

    // Increment the index of event that we have found for mixing purposes
    mixedEventFillIndex++;
      
  } // While loop for searching for a good event to mix with

  // For the next event, randomize the starting event again
  fMixingStartIndex = fRng->Integer(fnEventsInMixingFile);  // Start mixing from random spot in file

}

/*
 * Find matching events for the current event based on hiBin bins
 *
 *  std::vector<int>& mixedEventIndices = Reference to a vector that will hold the indiced of matched events
 *  const Int_t nEventsToMatch = Number of matched events that are searched for
 *  const Double_t vz = Vertex z-position of the event
 *  const Int_t hiBin = HiBin index of the current event
 *  const Int_t iEvent = Number of the current event
 *  
 */ 
void EECAnalyzer::FindBinnedHiBinMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, const Int_t hiBin, const Int_t eventNumber, const Int_t iEvent){

  // Event matching for PbPb
  Int_t mixedEventIndex = fMixingStartIndex;     // Current index for mixed event file
  Bool_t allEventsWentThrough  = false;          // Flag telling if we have gone through all the events in mixing file without finding a matching event
  Double_t vzTolerance = 0;                      // Starting tolerance for vz bin: exact match
  Int_t hiBinTolerance = 0;                      // Starting tolerance for hiBin bin: exact match
  Int_t mixedEventFillIndex = 0;                 // Index telling how many mixed events have been found already
  Bool_t skipEvent = false;                      // Flag for skipping unwanted miked events
  Bool_t checkForDuplicates = false;             // Flag for checking duplicate events in the mixing list
  Int_t vzBin = FindMixingVzBin(vz);             // Mixing bin index for the vz in the current event
  Int_t hiBinBin = FindMixingHiBinBin(hiBin);    // Mixing bin index for the hiBin in the current event

  mixedEventIndices.clear();  // Reset the previously found event index vector from previous events

  // Find events to mix until we have found defined number of matching events
  while(mixedEventFillIndex < nEventsToMatch){

    // By default, do not skip events
    skipEvent = false;

    // If we have checked all the events but not found an event to mix with, increase vz and hiBin tolerance
    if(allEventsWentThrough){
      if(fDebugLevel > 0){
        cout << "Could not find matching mixed events for event " << iEvent << endl;
        cout << "Increasing vz tolerance by 0.5 and hiBin tolerance by 1" << endl;
      }
      
      hiBinTolerance++;
      vzTolerance++;
      allEventsWentThrough = false;
      checkForDuplicates = true;
    }
    
    // Increment the counter for event index to be mixed with the current event
    mixedEventIndex++;
    
    // If we are out of bounds from the event in data file, reset the counter
    if(mixedEventIndex == fnEventsInMixingFile) {
      mixedEventIndex = -1;
      continue;
    }

    // If we come back to the first event index, we have gone through all the events without finding a matching mixed event
    if(mixedEventIndex == fMixingStartIndex) allEventsWentThrough = true;

    // Do not mix with events used in the previous iteration over the file
    if(checkForDuplicates){
      for(Int_t iMixedEvent : mixedEventIndices) {
        if(iMixedEvent == mixedEventIndex) skipEvent = true;
      }
    }

    // If we have already used the current mixed event, do not use the same event
    if(skipEvent) continue;

    // Match vz bin between the current event and mixed event
    if(TMath::Abs(fMixedEventVzBin.at(mixedEventIndex) - vzBin) > (vzTolerance + 1e-4)) continue;

    // Match hiBin bin between the current eent and mixed event
    if(TMath::Abs(fMixedEventHiBinBin.at(mixedEventIndex) - hiBinBin) > (hiBinTolerance + 1e-4)) continue;

    // Avoid mixing an event with itself
    if(fMixedEventEventNumber.at(mixedEventIndex) == eventNumber) continue;

    // Mark that this event is good to use for mixing
    mixedEventIndices.push_back(mixedEventIndex);

    // Increment the index of event that we have found for mixing purposes
    mixedEventFillIndex++;
      
  } // While loop for searching for a good event to mix with

  // For the next event, randomize the starting event again
  fMixingStartIndex = fRng->Integer(fnEventsInMixingFile);  // Start mixing from random spot in file

}

/*
 * Find matching events for the current event based on HF energy bins
 *
 *  std::vector<int>& mixedEventIndices = Reference to a vector that will hold the indiced of matched events
 *  const Int_t nEventsToMatch = Number of matched events that are searched for
 *  const Double_t vz = Vertex z-position of the event
 *  const Double_t hfEnergy = Energy in the HF calorimeters
 *  const Int_t eventNumber = Number of the current event
 *  const Int_t iEvent = Loop index of the current event
 *
 */
void EECAnalyzer::FindBinnedHFEnergyMatchedEvents(std::vector<int>& mixedEventIndices, const Int_t nEventsToMatch, const Double_t vz, const Double_t hfEnergy, const Int_t eventNumber, const Int_t iEvent){

  // Event matching for pPb
  Int_t mixedEventIndex = fMixingStartIndex;     // Current index for mixed event file
  Bool_t allEventsWentThrough  = false;          // Flag telling if we have gone through all the events in mixing file without finding a matching event
  Int_t vzTolerance = 0;                         // Starting tolerance for vz bin: exact match
  Int_t hfEnergyTolerance = 0;                   // Starting tolerance for HF energy bin: exact match
  Int_t mixedEventFillIndex = 0;                 // Index telling how many mixed events have been found already
  Bool_t skipEvent = false;                      // Flag for skipping unwanted miked events
  Bool_t checkForDuplicates = false;             // Flag for checking duplicate events in the mixing list
  Int_t vzBin = FindMixingVzBin(vz);                   // Mixing bin index for the vz in the current event
  Int_t hfEnergyBin = FindMixingHFEnergyBin(hfEnergy); // Mixing bin index for the HF energy in the current event

  mixedEventIndices.clear();  // Reset the previously found event index vector from previous events

  // Find events to mix until we have found defined number of matching events
  while(mixedEventFillIndex < nEventsToMatch){

    // By default, do not skip events
    skipEvent = false;

    // If we have checked all the events but not found an event to mix with, increase vz and hiBin tolerance
    if(allEventsWentThrough){
      if(fDebugLevel > 0){
        cout << "Could not find matching mixed events for event " << iEvent << endl;
        cout << "Trying to match also with adjacent vz and HF energy bins" << endl;
      }
      
      hfEnergyTolerance ++;
      vzTolerance ++;
      allEventsWentThrough = false;
      checkForDuplicates = true;
    }
    
    // Increment the counter for event index to be mixed with the current event
    mixedEventIndex++;
    
    // If we are out of bounds from the event in data file, reset the counter
    if(mixedEventIndex == fnEventsInMixingFile) {
      mixedEventIndex = -1;
      continue;
    }

    // If we come back to the first event index, we have gone through all the events without finding a matching mixed event
    if(mixedEventIndex == fMixingStartIndex) allEventsWentThrough = true;

    // Do not mix with events used in the previous iteration over the file
    if(checkForDuplicates){
      for(Int_t iMixedEvent : mixedEventIndices) {
        if(iMixedEvent == mixedEventIndex) skipEvent = true;
      }
    }

    // If we have already used the current mixed event, do not use the same event
    if(skipEvent) continue;

    // Match vz bin between the current event and mixed event
    if(TMath::Abs(fMixedEventVzBin.at(mixedEventIndex) - vzBin) > (vzTolerance + 1e-4)) continue;

    // Match HF energy bin between the current event and mixed event
    if(TMath::Abs(fMixedEventHFEnergyBin.at(mixedEventIndex) - hfEnergyBin) > (hfEnergyTolerance + 1e-4)) continue;

    // We are using minimum bias dataset for both mixing and analysis. Make sure not to match with the same event!
    if(fMixedEventEventNumber.at(mixedEventIndex) == eventNumber) continue;

    // Mark that this event is good to use for mixing
    mixedEventIndices.push_back(mixedEventIndex);

    // Increment the index of event that we have found for mixing purposes
    mixedEventFillIndex++;
      
  } // While loop for searching for a good event to mix with

  // For the next event, randomize the starting event again
  fMixingStartIndex = fRng->Integer(fnEventsInMixingFile);  // Start mixing from random spot in file

}

/*
 * In case we are doing binless mixing, prepare vectors for vz and hiBin from the mixing file for faster event matching
 */
void EECAnalyzer::PrepareMixingVectors(){
  
  // Print out debug message
  if(fDebugLevel > 1) cout << "Preparing for binless mixing" << endl;
  
  // Start reading the file from a random point
  fMixingStartIndex = fRng->Integer(fnEventsInMixingFile);

  // Read vz and hiBin from each event in event mixing file to memory.
  // This way we avoid loading different mixed events in a loop several times
  fMixedEventVz.clear();                       // Clear the vectors  
  fMixedEventHiBin.clear();                    // for any
  fMixedEventMultiplicity.clear();             // possible
  fMixedEventHFEnergy.clear();                 // contents they
  fMixedEventEventNumber.clear();              // might have

  // Event matching based on HiBin
  if(fEventMatchForMixing == kMatchHiBin){
    for(Int_t iMixedEvent = 0; iMixedEvent < fnEventsInMixingFile; iMixedEvent++){
      fMixedEventReader->GetEvent(iMixedEvent);
      fMixedEventEventNumber.push_back(fMixedEventReader->GetEventNumber());
      if(fMegaSkimMode){
        // Event selection is already applied in mega skim mode. No need to check it again here.
        fMixedEventVz.push_back(fMixedEventReader->GetVz());
        fMixedEventHiBin.push_back(fMixedEventReader->GetHiBin());
      } else if(PassEventCuts(fMixedEventReader,false)){
        fMixedEventVz.push_back(fMixedEventReader->GetVz());
        fMixedEventHiBin.push_back(fMixedEventReader->GetHiBin());
      } else { // If event cuts not passed, input values such that events will never be mixed with these
        fMixedEventVz.push_back(kNonsenseValue);
        fMixedEventHiBin.push_back(kNonsenseValue);
      }
    } // Loop over all mixed events
  } // Event matching based on HiBin

  // Event matching based on background multiplicity
  else if(fEventMatchForMixing == kMatchMultiplicity){
    for(Int_t iMixedEvent = 0; iMixedEvent < fnEventsInMixingFile; iMixedEvent++){
      fMixedEventReader->GetEvent(iMixedEvent);
      fMixedEventEventNumber.push_back(fMixedEventReader->GetEventNumber());
      if(fMegaSkimMode){
        // Event selection is already applied in mega skim mode. No need to check it again here.
        fMixedEventVz.push_back(fMixedEventReader->GetVz());
        fMixedEventMultiplicity.push_back(GetGenMultiplicity(fMixedEventReader, kSubeventAny, true));
      } else if(PassEventCuts(fMixedEventReader,false)){
        fMixedEventVz.push_back(fMixedEventReader->GetVz());
        fMixedEventMultiplicity.push_back(GetGenMultiplicity(fMixedEventReader, kSubeventAny, false));
      } else { // If event cuts not passed, input values such that events will never be mixed with these
        fMixedEventVz.push_back(kNonsenseValue);
        fMixedEventMultiplicity.push_back(kNonsenseValue);
      }
    } // Loop over all mixed events
  } // Event matching based on background multiplicity

  // Event matching based on HF energy
  else if(fEventMatchForMixing == kMatchHFEnergy){

    for(Int_t iMixedEvent = 0; iMixedEvent < fnEventsInMixingFile; iMixedEvent++){
      fMixedEventReader->GetEvent(iMixedEvent);
      fMixedEventEventNumber.push_back(fMixedEventReader->GetEventNumber());
      if(fMegaSkimMode){
        // Event selection is already applied in mega skim mode. No need to check it again here.
        fMixedEventVz.push_back(fMixedEventReader->GetVz());
        fMixedEventHFEnergy.push_back(fMixedEventReader->GetHFSum());
      } else if(PassEventCuts(fMixedEventReader,false)){
        fMixedEventVz.push_back(fMixedEventReader->GetVz());
        fMixedEventHFEnergy.push_back(fMixedEventReader->GetHFSum());
      } else { // If event cuts not passed, input values such that events will never be mixed with these
        fMixedEventVz.push_back(kNonsenseValue);
        fMixedEventHFEnergy.push_back(kNonsenseValue);
      }
    } // Loop over all mixed events

  } // Event matching based on HF energy

  // If something else is selected, crash the code and complain to user
  else {
    throw std::invalid_argument(Form("EECAnalyzer::ERROR There are three different event matching methods implemented for unbinned mixing. %d = Match HiBin, %d = Match HF energy, %d = Match background multiplicity. In the card, you defined EventMatchForMixing = %d. I do not know what to do with that, so I am crashing your code instead. Please fix this in the configuration.", kMatchHiBin, kMatchHFEnergy, kMatchMultiplicity, fEventMatchForMixing));
  }
  
}

/*
 * In case we are doing binned mixing, prepare vectors for vz and hiBin from the mixing file for faster event matching
 */
void EECAnalyzer::PrepareBinnedMixingVectors(){
  
  // Print out debug message
  if(fDebugLevel > 1) cout << "Preparing for binned mixing" << endl;
  
  // Start reading the file from a random point
  fMixingStartIndex = fRng->Integer(fnEventsInMixingFile);

  // Read vz and hiBin from each event in event mixing file to memory.
  // This way we avoid loading different mixed events in a loop several times
  fMixedEventVzBin.clear();                       // Clear the vectors for 
  fMixedEventHFEnergyBin.clear();                 // any possible contents 
  fMixedEventEventNumber.clear();                 // they might have

  // Event matching based on HiBin
  if(fEventMatchForMixing == kMatchHiBin){

    for(Int_t iMixedEvent = 0; iMixedEvent < fnEventsInMixingFile; iMixedEvent++){
      fMixedEventReader->GetEvent(iMixedEvent);
      fMixedEventEventNumber.push_back(fMixedEventReader->GetEventNumber());

      if(fMegaSkimMode){
        // Event selection is already applied in mega skim mode. No need to check it again here.
        fMixedEventVzBin.push_back(FindMixingVzBin(fMixedEventReader->GetVz()));
        fMixedEventHiBin.push_back(FindMixingHiBinBin(fMixedEventReader->GetHiBin()));
      } else if(PassEventCuts(fMixedEventReader,false)){
        fMixedEventVzBin.push_back(FindMixingVzBin(fMixedEventReader->GetVz()));
        fMixedEventHiBin.push_back(FindMixingHiBinBin(fMixedEventReader->GetHiBin()));
      } else { // If event cuts not passed, input values such that events will never be mixed with these
        fMixedEventVzBin.push_back(kNonsenseBin);
        fMixedEventHiBin.push_back(kNonsenseBin);
      }
    } // Loop over all mixed events
  } // Event matching based on HiBin

  // For MC, fill the vz and generator level particle multiplicity arrays for more accurate event matching
  else if(fEventMatchForMixing == kMatchMultiplicity){
    for(Int_t iMixedEvent = 0; iMixedEvent < fnEventsInMixingFile; iMixedEvent++){
      fMixedEventReader->GetEvent(iMixedEvent);
      fMixedEventEventNumber.push_back(fMixedEventReader->GetEventNumber());

      if(fMegaSkimMode){
        // Event selection is already applied in mega skim mode. No need to check it again here.
        fMixedEventVzBin.push_back(FindMixingVzBin(fMixedEventReader->GetVz()));
        fMixedEventMultiplicityBin.push_back(FindMixingMultiplicityBin(GetGenMultiplicity(fMixedEventReader, kSubeventAny, true)));
      } else if(PassEventCuts(fMixedEventReader,false)){
        fMixedEventVzBin.push_back(FindMixingVzBin(fMixedEventReader->GetVz()));
        fMixedEventMultiplicityBin.push_back(FindMixingMultiplicityBin(GetGenMultiplicity(fMixedEventReader, kSubeventAny, false)));
      } else { // If event cuts not passed, input values such that events will never be mixed with these
        fMixedEventVzBin.push_back(kNonsenseBin);
        fMixedEventMultiplicityBin.push_back(kNonsenseBin);
      }
    } // Loop over all mixed events
  }

  // Event matching based on HF energy
  else if(fEventMatchForMixing == kMatchHFEnergy){

    for(Int_t iMixedEvent = 0; iMixedEvent < fnEventsInMixingFile; iMixedEvent++){
      fMixedEventReader->GetEvent(iMixedEvent);
      fMixedEventEventNumber.push_back(fMixedEventReader->GetEventNumber());

      if(fMegaSkimMode){
        // Event selection is already applied in mega skim mode. No need to check it again here.
        fMixedEventVzBin.push_back(FindMixingVzBin(fMixedEventReader->GetVz()));
        fMixedEventHFEnergyBin.push_back(FindMixingHFEnergyBin(fMixedEventReader->GetHFSum()));
      } else if(PassEventCuts(fMixedEventReader,false)){
        fMixedEventVzBin.push_back(FindMixingVzBin(fMixedEventReader->GetVz()));
        fMixedEventHFEnergyBin.push_back(FindMixingHFEnergyBin(fMixedEventReader->GetHFSum()));
      } else { // If event cuts not passed, input values such that events will never be mixed with these
        fMixedEventVzBin.push_back(kNonsenseBin);
        fMixedEventHFEnergyBin.push_back(kNonsenseBin);
      }
    } // Loop over all mixed events
  } // Event matching based on HF energy

  // If something else is selected, crash the code and complain to user
  else {
    throw std::invalid_argument(Form("EECAnalyzer::ERROR There are three different event matching methods implemented for binned mixing. %d = Match HiBin, %d = Match HF energy, %d = Match background multiplicity. In the card, you defined EventMatchForMixing = %d. I do not know what to do with that, so I am crashing your code instead. Please fix this in the configuration.", kMatchHiBin, kMatchHFEnergy, kMatchMultiplicity, fEventMatchForMixing));
  }
  
}

/*
 * Find a generic mixing bin based on example bin border and bin width
 * 
 * Arguments:
 *  Double_t value = Value for which a mixing bin is sought for
 *  Double_t anchor = Lower bin boundary of the 0 bin
 *  Double_t width = Bin width for mixing bins
 *
 * Return: Bin index determined from the input values 
 *
 */
Int_t EECAnalyzer::FindMixingBin(Double_t value, Double_t anchor, Double_t width){
  Double_t anchoredValue = value - anchor;
  if(anchoredValue < 0) anchoredValue -= width;
  int bin = anchoredValue / width;
  return bin;
}

/*
 * Find a generic mixing bin based from a vector that contains mixing bin borders
 * 
 * Arguments:
 *  Double_t value = Value for which a mixing bin is sought for
 *  std::vector<Double_t> borders = Lower bin boundary of the 0 bin
 *
 * Return: Bin index determined from the input values 
 *
 */
Int_t EECAnalyzer::FindMixingBin(Double_t value, std::vector<Double_t> borders){
  Int_t bin = -1;
  for(Double_t binBorder : borders){
    if(value < binBorder) return bin;
    bin++;
  }
  return bin;
}

// Find an example value from certain mixing bin
Double_t EECAnalyzer::FindExampleValue(Int_t iBin, std::vector<Double_t> borders){

  // If there is more than two values in the borders array, they are directly the bin borders. Find the center of i:th bin
  if(borders.size() > 2){
    return (borders.at(iBin) + borders.at(iBin+1)) / 2.0;
  }

  // If there are only two values, they are an example bin border (anchor) and bin width.
  // We will have to calculate the result from this information 
  double binAnchor = borders.at(0);
  double binWidth = borders.at(1);

  return binAnchor + iBin*binWidth + binWidth/2.0;

}

// Find mixing vz bin from a vz value
Int_t EECAnalyzer::FindMixingVzBin(Double_t vz){
  if(fVzMixingBinBorders.size() > 2) FindMixingBin(vz, fVzMixingBinBorders);
  return FindMixingBin(vz, fVzMixingBinBorders.at(0), fVzMixingBinBorders.at(1));
}

// Find mixing HF energy bin from hfEnergy value
Int_t EECAnalyzer::FindMixingHFEnergyBin(Double_t hfEnergy){
  if(fHFEnergyMixingBinBorders.size() > 2) FindMixingBin(hfEnergy, fHFEnergyMixingBinBorders);
  return FindMixingBin(hfEnergy, fHFEnergyMixingBinBorders.at(0), fHFEnergyMixingBinBorders.at(1));
}

// Find mixing HiBin bin from HiBin value
Int_t EECAnalyzer::FindMixingHiBinBin(Int_t hiBin){
  if(fHiBinMixingBinBorders.size() > 2) FindMixingBin(hiBin, fHiBinMixingBinBorders);
  return FindMixingBin(hiBin, fHiBinMixingBinBorders.at(0), fHiBinMixingBinBorders.at(1));
}

// Find mixing multiplicity bin from multiplicity value
Int_t EECAnalyzer::FindMixingMultiplicityBin(Int_t multiplicity){
  if(fMultiplicityMixingBinBorders.size() > 2) FindMixingBin(multiplicity, fMultiplicityMixingBinBorders);
  return FindMixingBin(multiplicity, fMultiplicityMixingBinBorders.at(0), fMultiplicityMixingBinBorders.at(1));
}

// Get vector with bin borders used in mixing
std::vector<Double_t> EECAnalyzer::GetMixingBinBorders(std::vector<Double_t> binBorderVector, std::vector<Double_t> valueVector){

  // If there are more than 2 elements in the bin border vector, these are directly the bin borders
  if(binBorderVector.size() > 2) return binBorderVector;

  // Otherwise, the two elements in the vector are interpreted as example bin border (anchor) and bin width
  double binAnchor = binBorderVector.at(0);
  double binWidth = binBorderVector.at(1);

  // Find the maximum and minimum values from the value vector that are not non-sensical
  Double_t minValue = 1e9;
  Double_t maxValue = -1e9;
  for(Double_t currentValue : valueVector){
    if(TMath::Abs(currentValue - kNonsenseValue) < 1e-6) continue;
    if(currentValue < minValue) minValue = currentValue;
    if(currentValue > maxValue) maxValue = currentValue;
  }

  // Find the maximum and minimum bin index for these values
  Int_t minBinIndex = FindMixingBin(minValue, binAnchor, binWidth);
  Int_t maxBinIndex = FindMixingBin(maxValue, binAnchor, binWidth);

  // Get the number of bins from the bin indices
  Int_t nMixingBins = maxBinIndex - minBinIndex + 1;

  // Find the first bin border
  Double_t firstBinBorder = binAnchor + (minBinIndex * binWidth);

  // Create the bin border vector and add all bin borders to it
  std::vector<Double_t> binBorders;

  for(int iBin = 0; iBin <= nMixingBins; iBin++){
    binBorders.push_back(firstBinBorder + iBin*binWidth);
  }

  return binBorders;
}

// Get vector with bin borders used in mixing
std::vector<Double_t> EECAnalyzer::GetBinnedMixingBinBorders(std::vector<Double_t> binBorderVector, std::vector<Int_t> binVector){

  // If there are more than 2 elements in the bin border vector, these are directly the bin borders
  if(binBorderVector.size() > 2) return binBorderVector;

  // Otherwise, the two elements in the vector are interpreted as example bin border (anchor) and bin width
  double binAnchor = binBorderVector.at(0);
  double binWidth = binBorderVector.at(1);

  // Find the maximum and minimum bin indices from the bin index vector that are not non-sensical
  Int_t minBinIndex = 1e9;
  Int_t maxBinIndex = -1e9;
  for(Int_t currentBin : binVector){
    if(currentBin == kNonsenseBin) continue;
    if(currentBin < minBinIndex) minBinIndex = currentBin;
    if(currentBin > maxBinIndex) maxBinIndex = currentBin;
  }

  // Get the number of bins from the bin indices
  Int_t nMixingBins = maxBinIndex - minBinIndex + 1;

  // Find the first bin border
  Double_t firstBinBorder = binAnchor + (minBinIndex * binWidth);

  // Create the bin border vector and add all bin borders to it
  std::vector<Double_t> binBorders;

  for(int iBin = 0; iBin <= nMixingBins; iBin++){
    binBorders.push_back(firstBinBorder + iBin*binWidth);
  }

  return binBorders;
}

#include "JDrawer.h"

/*
 * Macro for deriving weighting functions for MC to match vz and centrality in data
 */
void deriveMonteCarloWeightsRun2(){

  enum enumDataType {kData, kMC, knDataTypes};

  // Open files for MC and data
  TFile *dataFile = TFile::Open("data/eecAnalysis_akFlowJet_onlyJets_weightEventInfo_combinedTriggers_processed_2023-03-06.root");
  //TFile *dataFile = TFile::Open("data/ppData2017_highForest_pfJets_onlyJets_L2rel_wtaAxis_processed_2019-08-05.root"); // Note: old file from previous analysis
  TFile *mcFile = TFile::Open("data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_wtaAxis_onlyJets_noTrigger_ptHatWeight_processed_2023-03-06.root");
  //TFile *mcFile = TFile::Open("data/dijet_ppMC_GenGen_Pythia8_pfJets_wtaAxis_onlyJets_processed_2019-08-06.root"); // Note: old file from previous analysis

  // Read the histograms for vz and centrality
  TH1D *hVz[knDataTypes];  // 0 = Data, 1 = MC
  TH1D *hCentrality[knDataTypes]; // 0 = Data, 1 = MC

  hVz[kData] = (TH1D*) dataFile->Get("vertexZweighted");            // vz histogram for data (include event weight when combining triggers)
  hCentrality[kData] = (TH1D*) dataFile->Get("centralityWeighted"); // centrality histogram for data (include event weight when combining triggers)
  hVz[kMC] = (TH1D*) mcFile->Get("vertexZ");                        // vz histogram for MC
  hCentrality[kMC] = (TH1D*) mcFile->Get("centrality");             // centrality histogram for MC

  // Normalize all histograms to one
  for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
    hVz[iDataType]->Scale(1/hVz[iDataType]->Integral());
    hCentrality[iDataType]->Scale(1/hCentrality[iDataType]->Integral());
  }

  // Draw the histogram before weighting
  JDrawer *drawer = new JDrawer();

  drawer->DrawHistogram(hVz[kData],"v_{z} (cm)","counts", " ");
  hVz[kMC]->SetLineColor(kRed);
  hVz[kMC]->Draw("same");

  drawer->DrawHistogram(hCentrality[kData],"centrality","counts", " ");
  hCentrality[kMC]->SetLineColor(kRed);
  hCentrality[kMC]->Draw("same");

  // Do not draw over centrality
  drawer->CreateCanvas();

  // Calculate the ratio of the vz histograms and fit it with pol6
  TH1D *hVzRatio = (TH1D*) hVz[kData]->Clone();
  hVzRatio->Divide(hVz[kMC]);
  hVzRatio->Fit("pol6","","",-15,15);
  TF1* vzFit = hVzRatio->GetFunction("pol6");

  // Calculate the ratio of the centrality histograms and fit it with pol6
  TH1D *hCentralityRatioCentral = (TH1D*) hCentrality[kData]->Clone();
  TH1D *hCentralityRatioPeripheral = (TH1D*) hCentrality[kData]->Clone();
  hCentralityRatioCentral->Divide(hCentrality[kMC]);
  hCentralityRatioCentral->Fit("pol6","","",0,30);
  hCentralityRatioPeripheral->Divide(hCentrality[kMC]);
  hCentralityRatioPeripheral->Fit("pol6","","",30,90);
  TF1 *centralityFitCentral = hCentralityRatioCentral->GetFunction("pol6");
  TF1 *centralityFitPeripheral = hCentralityRatioPeripheral->GetFunction("pol6");

  cout << "Central at 30: " << centralityFitCentral->Eval(30) << endl;
  cout << "Peripheral at 30: " << centralityFitPeripheral->Eval(30) << endl;
  centralityFitPeripheral->SetParameter(0, centralityFitPeripheral->GetParameter(0) + (centralityFitCentral->Eval(30) - centralityFitPeripheral->Eval(30)));
  cout << "Adjusted peripheral at 30: " << centralityFitPeripheral->Eval(30) << endl;
  cout << "Adjusted parameter value: " << centralityFitPeripheral->GetParameter(0) << endl;

  // Draw the fitted ratios
  drawer->DrawHistogram(hVzRatio,"v_{z} (cm)","Data/MC ratio");
  drawer->DrawHistogram(hCentralityRatioCentral,"centrality","Data/MC ratio");
  drawer->DrawHistogram(hCentralityRatioPeripheral,"centrality","Data/MC ratio");

  // Check with a small simulation that we regain data v_z and centrality
  double randomNumberVz, weigthVz;
  double randomNumberCentrality, weigthCentrality;
  TH1D *vzReconstructed = (TH1D*) hVz[kData]->Clone("vzClone");
  vzReconstructed->Reset();
  TH1D *centralityReconstructed = (TH1D*) hCentrality[0]->Clone("centralityClone");
  centralityReconstructed->Reset();
  for(int i = 0; i < 1000000; i++){
    randomNumberVz = hVz[kMC]->GetRandom();
    weigthVz = vzFit->Eval(randomNumberVz);
    vzReconstructed->Fill(randomNumberVz,weigthVz);

    randomNumberCentrality = hCentrality[kMC]->GetRandom();
    weigthCentrality = randomNumberCentrality < 30 ? centralityFitCentral->Eval(randomNumberCentrality) : centralityFitPeripheral->Eval(randomNumberCentrality);
    centralityReconstructed->Fill(randomNumberCentrality,weigthCentrality);
  }

  vzReconstructed->Scale(1.0/vzReconstructed->Integral());
  drawer->DrawHistogram(hVz[kData],"v_{z} (cm)","counts", " ");
  vzReconstructed->SetLineColor(kMagenta);
  vzReconstructed->Draw("same");

  centralityReconstructed->Scale(1.0/centralityReconstructed->Integral());
  drawer->DrawHistogram(hCentrality[kData],"centrality","counts", " ");
  centralityReconstructed->SetLineColor(kMagenta);
  centralityReconstructed->Draw("same");
}

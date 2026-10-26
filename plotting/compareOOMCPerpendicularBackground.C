#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"
#include <TLine.h>
#include <algorithm>
#include <cmath>

/*
 * OO MC background closure plot, adapted from compareDifferentEECpairings.C.
 * True background: same-jet mixed-subevent + background-background pairs.
 * Perpendicular estimate: average(JP1,JP2) - P1P2 + average(P1P1,P2P2),
 * matching EECHistogramManager::SubtractBackground. Both curves retain the
 * same raw EEC scale; neither is independently normalized or fitted.
 */
void compareOOMCPerpendicularBackground(TString inputFileName = "projected/OOMC_GenGen_withSubevents.root", double minCentrality = 0., double maxCentrality = 30.){
  TFile* inputFile = TFile::Open(inputFileName, "READ");
  if(!inputFile || inputFile->IsZombie()){
    cout << "Cannot open " << inputFileName << endl;
    return;
  }
  EECCard* card = new EECCard(inputFile);
  if(!card->GetDataType().Contains("OO") || !card->GetDataType().Contains("MC")){
    cout << "This comparison requires OO MC with subevents." << endl;
    return;
  }

  // OO centrality is used directly, without the PbPb MC centrality shift.
  const std::pair<double,double> centrality = std::make_pair(minCentrality,maxCentrality);
  const std::pair<double,double> jetPt = std::make_pair(40.,80.);
  const double trackPt = 2.;
  const int iCentrality = card->FindBinIndexCentrality(centrality);
  const int iJetPt = card->FindBinIndexJetPtEEC(jetPt);
  const int iTrackPt = card->GetBinIndexTrackPtEEC(trackPt);
  if(iCentrality < 0 || iJetPt < 0 || iTrackPt < 0){
    cout << "The requested OO centrality, jet bin, or track threshold is unavailable." << endl;
    return;
  }
  EECHistogramManager* histograms = new EECHistogramManager(inputFile,card);
  auto get = [&](int pairing, int subevent, const char* name) -> TH1D* {
    TH1D* source = histograms->GetHistogramEnergyEnergyCorrelator(
      EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, pairing, subevent);
    if(!source){
      cout << "Missing histogram for pairing " << pairing << ", subevent " << subevent << endl;
      return nullptr;
    }
    TH1D* copy = (TH1D*)source->Clone(name);
    copy->SetDirectory(nullptr);
    return copy;
  };

  // Projected S1 already combines both mixed subevent orderings; use it once.
  TH1D* truth = get(EECHistograms::kSameJetPair, EECHistograms::kPythiaHydjet, "ooTrueBackground");
  TH1D* backgroundPairs = get(EECHistograms::kSameJetPair, EECHistograms::kHydjetHydjet, "ooBackgroundPairs");
  const int all = EECHistograms::knSubeventCombinations;
  TH1D* perpendicular = get(EECHistograms::kSignalPerpendicularConePair, all, "ooPerpendicularBackground");
  TH1D* jetSecond = get(EECHistograms::kSignalSecondPerpendicularConePair, all, "ooJetSecondPerpendicular");
  TH1D* betweenCones = get(EECHistograms::kPerpendicularPerpendicularConePair, all, "ooBetweenPerpendicularCones");
  TH1D* sameCone = get(EECHistograms::kPerpendicularConePair, all, "ooSamePerpendicularCone");
  TH1D* secondSameCone = get(EECHistograms::kSecondPerpendicularConePair, all, "ooSecondSamePerpendicularCone");
  if(!truth || !backgroundPairs || !perpendicular || !jetSecond || !betweenCones || !sameCone || !secondSameCone) return;

  truth->Add(backgroundPairs);
  perpendicular->Add(jetSecond);
  perpendicular->Scale(0.5);
  perpendicular->Add(betweenCones,-1.);
  perpendicular->Add(sameCone,0.5);
  perpendicular->Add(secondSameCone,0.5);
  TH1D* ratio = (TH1D*)perpendicular->Clone("ooPerpendicularToTrueBackground");
  ratio->SetDirectory(nullptr);
  ratio->Divide(truth);

  // Use the same Delta-r interval and split-canvas style as the reference macro.
  const double minDeltaR = 0.008;
  const double maxDeltaR = 0.39;
  double eecMin = 1.e100, eecMax = 0.;
  double ratioMin = 0.91, ratioMax = 1.09;
  const int firstBin = truth->FindBin(minDeltaR+1.e-6);
  const int lastBin = truth->FindBin(maxDeltaR-1.e-6);
  for(int bin = firstBin; bin <= lastBin; ++bin){
    for(TH1D* curve : {truth,perpendicular}){
      const double value = curve->GetBinContent(bin);
      if(value > 0.) eecMin = std::min(eecMin,value);
      eecMax = std::max(eecMax,value);
    }
    if(truth->GetBinContent(bin) > 0.){
      ratioMin = std::min(ratioMin,ratio->GetBinContent(bin)-ratio->GetBinError(bin));
      ratioMax = std::max(ratioMax,ratio->GetBinContent(bin)+ratio->GetBinError(bin));
    }
  }
  if(eecMax <= 0.){
    cout << "No background entries in the requested selection." << endl;
    return;
  }
  const double ratioPadding = 0.08*(ratioMax-ratioMin);
  ratioMin -= ratioPadding;
  ratioMax += ratioPadding;

  truth->SetMarkerStyle(kFullCircle);
  truth->SetMarkerColor(kBlack);
  truth->SetLineColor(kBlack);
  perpendicular->SetMarkerStyle(kOpenSquare);
  perpendicular->SetMarkerColor(kRed);
  perpendicular->SetLineColor(kRed);
  ratio->SetMarkerStyle(kOpenSquare);
  ratio->SetMarkerColor(kRed);
  ratio->SetLineColor(kRed);

  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.17);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.9);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogX(true);
  drawer->CreateSplitCanvas();
  drawer->SetLogY(true);
  truth->GetXaxis()->SetRangeUser(minDeltaR,maxDeltaR);
  truth->GetYaxis()->SetRangeUser(eecMin/2.,eecMax*2.);
  perpendicular->GetXaxis()->SetRangeUser(minDeltaR,maxDeltaR);
  drawer->DrawHistogramToUpperPad(truth,"#Deltar","EEC (A.U.)"," ","E1");
  perpendicular->Draw("E1 SAME");

  TLegend* legend = new TLegend(0.43,0.02,0.98,0.36);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->SetTextFont(62);
  legend->AddEntry((TObject*)nullptr,"Pythia+Hijing 5.36 TeV","");
  legend->AddEntry((TObject*)nullptr,Form("OO, Cent: %.0f-%.0f%%",minCentrality,maxCentrality),"");
  legend->AddEntry((TObject*)nullptr,"40 < jet p_{T} < 80 GeV","");
  legend->AddEntry((TObject*)nullptr,Form("p_{T}^{ch} > 2 GeV, n=%d",card->GetWeightExponent()),"");
  legend->AddEntry(truth,"True background","p");
  legend->AddEntry(perpendicular,"Perpendicular cone background","p");
  legend->Draw();

  drawer->SetLogY(false);
  drawer->SetGridY(true);
  ratio->GetXaxis()->SetRangeUser(minDeltaR,maxDeltaR);
  ratio->GetYaxis()->SetRangeUser(ratioMin,ratioMax);
  drawer->DrawHistogramToLowerPad(ratio,"#Deltar","#frac{Perp. cone bkg}{True background}"," ","E1");
  TLine* unity = new TLine(minDeltaR,1.,maxDeltaR,1.);
  unity->SetLineStyle(2);
  unity->Draw();
  drawer->SetGridY(false);

  TString output = Form("figures/OOMC_GenGen_perpendicularBackground_C%.0f-%.0f_J40-80_T2p0",minCentrality,maxCentrality);
  gPad->GetCanvas()->SaveAs(output+".pdf");
  cout << "Input: " << inputFileName << endl;
  cout << "True background integral = " << truth->Integral(firstBin,lastBin,"width")
       << ", perpendicular estimate = " << perpendicular->Integral(firstBin,lastBin,"width")
       << ", estimate/truth = " << perpendicular->Integral(firstBin,lastBin,"width")/truth->Integral(firstBin,lastBin,"width") << endl;
}

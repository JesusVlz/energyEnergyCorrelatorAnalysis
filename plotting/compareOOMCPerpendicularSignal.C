#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"
#include <TLine.h>
#include <algorithm>
#include <cmath>

/*
 * OO MC signal extraction plot, copied from compareOOMCPerpendicularBackground.C.
 * True signal: same-jet S0. Extracted signal: same-jet total minus B_perp.
 * Perpendicular estimate: average(JP1,JP2) - P1P2 + average(P1P1,P2P2),
 * matching EECHistogramManager::SubtractBackground. All three curves retain
 * the same raw EEC scale; none is independently normalized or fitted.
 */
void compareOOMCPerpendicularSignal(TString inputFileName = "projected/OOMC_GenGen_withSubevents.root", double minJetPt = 40., double maxJetPt = 80., double minCentrality = 0., double maxCentrality = 30.){
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
  const std::pair<double,double> jetPt = std::make_pair(minJetPt,maxJetPt);
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

  const int all = EECHistograms::knSubeventCombinations;
  TH1D* truth = get(EECHistograms::kSameJetPair, EECHistograms::kPythiaPythia, "ooTrueSignal");
  TH1D* total = get(EECHistograms::kSameJetPair, all, "ooSignalPlusBackground");
  TH1D* perpendicular = get(EECHistograms::kSignalPerpendicularConePair, all, "ooPerpendicularBackground");
  TH1D* jetSecond = get(EECHistograms::kSignalSecondPerpendicularConePair, all, "ooJetSecondPerpendicular");
  TH1D* betweenCones = get(EECHistograms::kPerpendicularPerpendicularConePair, all, "ooBetweenPerpendicularCones");
  TH1D* sameCone = get(EECHistograms::kPerpendicularConePair, all, "ooSamePerpendicularCone");
  TH1D* secondSameCone = get(EECHistograms::kSecondPerpendicularConePair, all, "ooSecondSamePerpendicularCone");
  if(!truth || !total || !perpendicular || !jetSecond || !betweenCones || !sameCone || !secondSameCone) return;

  perpendicular->Add(jetSecond);
  perpendicular->Scale(0.5);
  perpendicular->Add(betweenCones,-1.);
  perpendicular->Add(sameCone,0.5);
  perpendicular->Add(secondSameCone,0.5);
  TH1D* extracted = (TH1D*)total->Clone("ooPerpendicularExtractedSignal");
  extracted->SetDirectory(nullptr);
  extracted->Add(perpendicular,-1.);
  TH1D* totalRatio = (TH1D*)total->Clone("ooTotalToTrueSignal");
  TH1D* extractedRatio = (TH1D*)extracted->Clone("ooExtractedToTrueSignal");
  for(TH1D* ratio : {totalRatio,extractedRatio}){
    ratio->SetDirectory(nullptr);
    ratio->Divide(truth);
  }

  // Use the same Delta-r interval and split-canvas style as the reference macro.
  const double minDeltaR = 0.008;
  const double maxDeltaR = 0.39;
  double eecMin = 1.e100, eecMax = 0.;
  double ratioMin = 0.8, ratioMax = 1.2;
  int nonpositiveExtractedBins = 0;
  const int firstBin = truth->FindBin(minDeltaR+1.e-6);
  const int lastBin = truth->FindBin(maxDeltaR-1.e-6);
  for(int bin = firstBin; bin <= lastBin; ++bin){
    for(TH1D* curve : {truth,total,extracted}){
      const double value = curve->GetBinContent(bin);
      if(value > 0.) eecMin = std::min(eecMin,value);
      eecMax = std::max(eecMax,value);
    }
    if(extracted->GetBinContent(bin) <= 0.) ++nonpositiveExtractedBins;
    if(truth->GetBinContent(bin) > 0.){
      for(TH1D* ratio : {totalRatio,extractedRatio}){
        ratioMin = std::min(ratioMin,ratio->GetBinContent(bin)-ratio->GetBinError(bin));
        ratioMax = std::max(ratioMax,ratio->GetBinContent(bin)+ratio->GetBinError(bin));
      }
    }
  }
  if(eecMax <= 0.){
    cout << "No signal entries in the requested selection." << endl;
    return;
  }
  const double ratioPadding = 0.08*(ratioMax-ratioMin);
  ratioMin -= ratioPadding;
  ratioMax += ratioPadding;

  truth->SetMarkerStyle(kFullCircle);
  truth->SetMarkerColor(kBlack);
  truth->SetLineColor(kBlack);
  for(TH1D* curve : {total,totalRatio}){
    curve->SetMarkerStyle(kOpenSquare);
    curve->SetMarkerColor(kRed);
    curve->SetLineColor(kRed);
  }
  for(TH1D* curve : {extracted,extractedRatio}){
    curve->SetMarkerStyle(kOpenCross);
    curve->SetMarkerColor(kBlue);
    curve->SetLineColor(kBlue);
  }

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
  total->GetXaxis()->SetRangeUser(minDeltaR,maxDeltaR);
  extracted->GetXaxis()->SetRangeUser(minDeltaR,maxDeltaR);
  drawer->DrawHistogramToUpperPad(truth,"#Deltar","EEC (A.U.)"," ","E1");
  total->Draw("E1 SAME");
  extracted->Draw("E1 SAME");

  TLegend* legend = new TLegend(0.18,0.03,0.77,0.49);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->SetTextFont(62);
  legend->AddEntry((TObject*)nullptr,"Pythia+Hijing 5.36 TeV","");
  legend->AddEntry((TObject*)nullptr,Form("OO, Cent: %.0f-%.0f%%",minCentrality,maxCentrality),"");
  legend->AddEntry((TObject*)nullptr,Form("%.0f < jet p_{T} < %.0f GeV",minJetPt,maxJetPt),"");
  legend->AddEntry((TObject*)nullptr,Form("p_{T}^{ch} > 2 GeV, n=%d",card->GetWeightExponent()),"");
  legend->AddEntry(truth,"True signal","p");
  legend->AddEntry(total,"Signal + background","p");
  legend->AddEntry(extracted,"Extracted signal (perp. cone)","p");
  legend->Draw();

  drawer->SetLogY(false);
  drawer->SetGridY(true);
  for(TH1D* ratio : {totalRatio,extractedRatio}){
    ratio->GetXaxis()->SetRangeUser(minDeltaR,maxDeltaR);
    ratio->GetYaxis()->SetRangeUser(ratioMin,ratioMax);
  }
  drawer->DrawHistogramToLowerPad(totalRatio,"#Deltar","#frac{Color}{True signal}"," ","E1");
  extractedRatio->Draw("E1 SAME");
  TLine* unity = new TLine(minDeltaR,1.,maxDeltaR,1.);
  unity->SetLineStyle(2);
  unity->Draw();
  drawer->SetGridY(false);

  TString output = Form("figures/OOMC_GenGen_perpendicularSignal_C%.0f-%.0f_J%.0f-%.0f_T2p0",minCentrality,maxCentrality,minJetPt,maxJetPt);
  gPad->GetCanvas()->SaveAs(output+".pdf");
  cout << "Input: " << inputFileName << endl;
  cout << "True signal integral = " << truth->Integral(firstBin,lastBin,"width")
       << ", total = " << total->Integral(firstBin,lastBin,"width")
       << ", extracted signal = " << extracted->Integral(firstBin,lastBin,"width")
       << ", extracted/truth = " << extracted->Integral(firstBin,lastBin,"width")/truth->Integral(firstBin,lastBin,"width") << endl;
  cout << "Nonpositive extracted signal bins in displayed range: " << nonpositiveExtractedBins << endl;
}

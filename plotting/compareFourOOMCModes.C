#include "EECCard.h"
#include "EECHistogramManager.h"
#include "../src/EECHistograms.h"
R__LOAD_LIBRARY(plotting/DrawingClasses.so)

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorT.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace FourOOMCModes {
const std::array<const char*,4> kModes{{"RecoReco","RecoGen","GenReco","GenGen"}};
const std::array<const char*,4> kLabels{{
    "RecoReco: reco jets, reco tracks",
    "RecoGen: reco jets, gen particles*",
    "GenReco: gen jets, reco tracks",
    "GenGen: gen jets, gen particles*"}};
const std::array<double,4> kThresholds{{0.5,1.,2.,3.}};
const std::array<int,4> kColors{{kBlack,kBlue+1,kRed+1,kGreen+2}};
const std::array<int,4> kMarkers{{20,24,21,25}};
constexpr double kXLow=0.0008, kXHigh=0.8;
using Panel=std::array<std::unique_ptr<TH1D>,4>;

void Require(bool condition,const TString& message) {
  if(!condition) throw std::runtime_error(message.Data());
}
int Scalar(TFile& file,const char* name) {
  auto* value=dynamic_cast<TVectorT<float>*>(file.Get(Form("JCard/%s",name)));
  Require(value && value->GetNrows()==1,Form("Missing JCard/%s",name));
  return std::lround((*value)[value->GetLwb()]);
}
std::unique_ptr<TH1D> Shape(EECHistogramManager& manager,int mode,int selection) {
  auto* card=manager.GetCard();
  const int centrality=card->FindBinIndexCentrality(0.,30.);
  const int jet=card->FindBinIndexJetPtEEC(40.,80.);
  const double threshold=kThresholds[selection];
  const int track=card->GetBinIndexTrackPtEEC(threshold);
  Require(centrality>=0 && jet>=0 && track>=0 &&
      std::abs(card->GetLowBinBorderTrackPtEEC(track)-threshold)<1.e-6,
      "Requested centrality, jet bin, or track threshold is missing");
  // Select the same explicit bin indices for each mode. The default subevent
  // index selects the all-combinations total.
  auto* source=manager.GetHistogramEnergyEnergyCorrelator(
      EECHistogramManager::kEnergyEnergyCorrelator,centrality,jet,track,
      EECHistograms::kSameJetPair);
  Require(source,Form("Missing %s same-jet total C%d J%d T%d",kModes[mode],centrality,jet,track));
  std::unique_ptr<TH1D> shape(static_cast<TH1D*>(source->Clone(Form("fourOOMC_%s_T%d",kModes[mode],selection))));
  shape->SetDirectory(nullptr); shape->GetXaxis()->SetRange(0,0);
  const double area=shape->Integral("width");
  Require(std::isfinite(area) && area>0.,"Empty or nonfinite EEC area");
  shape->Scale(1./area);
  Require(std::abs(shape->Integral("width")-1.)<1.e-10,"EEC area normalization failed");
  for(int bin=0;bin<shape->GetNcells();++bin)
    Require(std::isfinite(shape->GetBinContent(bin)) && std::isfinite(shape->GetBinError(bin)),
        "Nonfinite EEC content or statistical uncertainty");
  std::cout<<"FOUR_OOMC_SHAPE mode="<<kModes[mode]<<" track>"<<threshold
      <<" source="<<source->GetName()<<" raw_area="<<area<<" normalized_area="<<shape->Integral("width")<<'\n';
  return shape;
}
void DrawPanel(TPad& pad,Panel& shapes,int selection,double yLow,double yHigh) {
  pad.cd();pad.SetLogx();pad.SetLogy();pad.SetTicks(1,1);
  pad.SetLeftMargin(0.14);pad.SetRightMargin(0.03);
  pad.SetBottomMargin(0.13);pad.SetTopMargin(0.025);
  auto* frame=pad.DrawFrame(kXLow,yLow,kXHigh,yHigh);
  frame->GetXaxis()->SetTitle("#Deltar");frame->GetYaxis()->SetTitle("Unit-area EEC");
  frame->GetXaxis()->SetTitleOffset(1.0);frame->GetYaxis()->SetTitleOffset(1.07);
  for(auto* axis:{frame->GetXaxis(),frame->GetYaxis()}) {
    axis->SetLabelFont(42);axis->SetTitleFont(42);
    axis->SetLabelSize(0.047);axis->SetTitleSize(0.054);axis->CenterTitle();
  }
  for(int mode=0;mode<4;++mode) {
    auto& shape=*shapes[mode];shape.SetLineColor(kColors[mode]);shape.SetMarkerColor(kColors[mode]);
    shape.SetMarkerStyle(kMarkers[mode]);shape.SetMarkerSize(0.65);shape.SetLineWidth(1);
    shape.Draw("E1 SAME");
  }
  TLatex text;text.SetNDC();text.SetTextFont(62);text.SetTextSize(0.048);
  text.DrawLatex(0.62,0.91,Form("p_{T}^{ch} > %.1f GeV",kThresholds[selection]));
  TLegend legend(0.18,0.16,0.83,0.42);
  legend.SetFillStyle(0);legend.SetBorderSize(0);legend.SetTextFont(42);legend.SetTextSize(0.038);
  legend.SetMargin(0.08);
  for(int mode=0;mode<4;++mode)legend.AddEntry(shapes[mode].get(),kLabels[mode],"lp");
  legend.DrawClone();pad.RedrawAxis();
}
} // namespace FourOOMCModes

// Run from repository root with the CMSSW ROOT runtime:
// root -l -b -q 'plotting/compareFourOOMCModes.C()'
// Diagnostic totals only: neither reconstructed-track mode has truth-origin
// labels. Generator-particle modes use a forest restricted to pT >= 5 GeV;
// therefore the four nominal thresholds select the same generator population.
void compareFourOOMCModes(const char* inputDirectory=".",const char* outputDirectory="figures") {
  using namespace FourOOMCModes;
  try {
    gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetEndErrorSize(2);
    std::array<Panel,4> panels;
    double smallest=std::numeric_limits<double>::max(),largest=0.;
    for(int mode=0;mode<4;++mode) {
      const TString name=Form("%s/merged_OOMC_%s_20261022_processed.root",inputDirectory,kModes[mode]);
      std::unique_ptr<TFile> input(TFile::Open(name,"READ"));
      Require(input && !input->IsZombie() && !input->TestBit(TFile::kRecovered),Form("Cannot open %s",name.Data()));
      Require(Scalar(*input,"DataType")==10 && Scalar(*input,"McCorrelationType")==mode,"Unexpected OOMC mode in input");
      // The large manager owns its card. All shapes are detached clones.
      auto manager=std::make_unique<EECHistogramManager>(input.get());
      for(int selection=0;selection<4;++selection) {
        panels[selection][mode]=FourOOMCModes::Shape(*manager,mode,selection);
        const auto& shape=*panels[selection][mode];
        for(int bin=1;bin<=shape.GetNbinsX();++bin) {
          const double x=shape.GetBinCenter(bin),value=shape.GetBinContent(bin);
          if(x<kXLow || x>kXHigh || value<=0.)continue;
          smallest=std::min(smallest,value);largest=std::max(largest,value+shape.GetBinError(bin));
        }
      }
    }
    Require(largest>0. && smallest<std::numeric_limits<double>::max(),"No positive EEC values");
    Require(gSystem->mkdir(outputDirectory,true)==0 || !gSystem->AccessPathName(outputDirectory),"Cannot create output directory");
    TCanvas canvas("OOMCFourModes","Four OOMC correlation modes",1500,1200);canvas.SetFillColor(kWhite);
    const double yLow=smallest*0.45,yHigh=largest*4.;
    std::array<std::unique_ptr<TPad>,4> pads;
    for(int selection=0;selection<4;++selection) {
      const int column=selection%2,row=selection/2;
      const double bottom=row==0?0.525:0.09,top=row==0?0.96:0.525;
      pads[selection]=std::make_unique<TPad>(Form("fourOOMCPanel%d",selection),"",0.5*column,bottom,0.5*(column+1),top);
      canvas.cd();pads[selection]->Draw();DrawPanel(*pads[selection],panels[selection],selection,yLow,yHigh);
    }
    canvas.cd();TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextAlign(22);
    text.SetTextSize(0.024);text.DrawLatex(0.5,0.985,"Pythia+Hijing O#minusO simulation: 0#minus30%, 40 < p_{T}^{jet} < 80 GeV");
    text.SetTextSize(0.0185);text.DrawLatex(0.5,0.067,"*Generator-particle input retains p_{T} #geq 5 GeV: nominal 0.5, 1, 2, 3 GeV cuts give identical gen selections.");
    text.SetTextSize(0.017);text.DrawLatex(0.5,0.043,"Same-jet total pairs; each curve normalized to unit area over the full histogram range; statistical uncertainties.");
    text.DrawLatex(0.5,0.020,"Jet p_{T} selection is reconstructed for RecoReco/RecoGen and generator level for GenReco/GenGen; populations differ.");
    canvas.Modified();canvas.Update();
    const TString stem=Form("%s/OOMC_fourModes_C0-30_J40-80",outputDirectory);
    canvas.Print(stem+".pdf");canvas.Print(stem+".png");
    std::cout<<"FOUR_OOMC_COMPARISON_PASS panels=4 modes=4 output="<<stem<<std::endl;
  } catch(const std::exception& error) {
    std::cerr<<"FOUR_OOMC_COMPARISON_FAILED "<<error.what()<<std::endl;gSystem->Exit(1);
  }
}

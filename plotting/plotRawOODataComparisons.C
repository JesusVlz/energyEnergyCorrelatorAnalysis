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
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace RawOODataComparisons {
constexpr double kDrawLow = 0.005, kDrawHigh = 0.8;
const std::array<std::pair<double,double>,4> kJetBins{{{40,80},{80,120},{120,160},{160,200}}};
const std::array<double,4> kTrackCuts{{0.5,1.,1.5,2.}};
const std::array<int,4> kColors{{kBlack,kBlue+1,kRed+1,kGreen+2}};
const std::array<int,4> kMarkers{{21,20,33,34}};

void Require(bool condition, const TString& message) {
  if (!condition) throw std::runtime_error(message.Data());
}
struct CentralityPlots {
  double low, high;
  std::array<std::unique_ptr<TH1D>,4> jets, tracks;
};
std::unique_ptr<TH1D> ReadShape(EECHistogramManager& manager, int centrality,
    std::pair<double,double> jets, double trackCut, const TString& name) {
  auto* card = manager.GetCard();
  const int jetBin = card->FindBinIndexJetPtEEC(jets);
  const int trackBin = card->GetBinIndexTrackPtEEC(trackCut);
  Require(jetBin >= 0 && jetBin < card->GetNJetPtBinsEEC(), "Requested jet bin is unavailable");
  Require(trackBin >= 0 && trackBin < card->GetNTrackPtBinsEEC() &&
      std::abs(card->GetLowBinBorderTrackPtEEC(trackBin)-trackCut)<1.e-6,
      "Requested track threshold is unavailable");
  // Resolve each stored centrality and jet selection explicitly.
  auto* source = manager.GetHistogramEnergyEnergyCorrelator(
      EECHistogramManager::kEnergyEnergyCorrelator, centrality, jetBin, trackBin,
      EECHistograms::kSameJetPair);
  Require(source, Form("Missing raw data histogram C%d J%d T%d",centrality,jetBin,trackBin));
  std::unique_ptr<TH1D> shape(static_cast<TH1D*>(source->Clone(name)));
  shape->SetDirectory(nullptr);
  for (int bin=0; bin<shape->GetNcells(); ++bin)
    Require(std::isfinite(shape->GetBinContent(bin)) && std::isfinite(shape->GetBinError(bin)),
        "Nonfinite histogram content or statistical uncertainty");
  // Match the raw comparison panels in EECDrawer.cxx. Projections already
  // contain bin-width scaling; normalize their full regular-bin width integral.
  shape->GetXaxis()->SetRange(0,0);
  const double integral = shape->Integral("width");
  Require(std::isfinite(integral) && integral>0, "Empty raw EEC normalization interval");
  shape->Scale(1./integral);
  Require(std::abs(shape->Integral("width")-1.)<1.e-10, "EEC area normalization failed");
  std::cout << "RAW_OO_DATA_SHAPE C=" << centrality << " jet=" << jets.first << "-" << jets.second
      << " track>" << trackCut << " source=" << source->GetName()
      << " raw_area=" << integral << " normalized_area=" << shape->Integral("width") << '\n';
  return shape;
}
void DrawPanel(TPad& pad, CentralityPlots& plots, bool jetScan, double yLow, double yHigh) {
  pad.cd(); pad.SetLogx(); pad.SetLogy(); pad.SetTicks(1,1);
  pad.SetLeftMargin(0.145); pad.SetRightMargin(0.035);
  pad.SetTopMargin(0.04); pad.SetBottomMargin(0.145);
  auto& shapes = jetScan ? plots.jets : plots.tracks;
  for (int i=0; i<4; ++i) {
    auto& h = *shapes[i];
    h.SetTitle(""); h.SetStats(false);
    h.SetMarkerColor(kColors[i]); h.SetLineColor(kColors[i]);
    h.SetMarkerStyle(kMarkers[i]); h.SetMarkerSize(i>=2 ? 1.15 : 0.9); h.SetLineWidth(1);
    h.GetXaxis()->SetRangeUser(kDrawLow,kDrawHigh); h.SetMinimum(yLow); h.SetMaximum(yHigh);
    h.GetXaxis()->SetTitle("#Deltar"); h.GetYaxis()->SetTitle("EEC");
    h.GetXaxis()->CenterTitle(); h.GetYaxis()->CenterTitle();
    for (auto* axis : {h.GetXaxis(),h.GetYaxis()}) {
      axis->SetLabelFont(42); axis->SetTitleFont(42);
      axis->SetLabelSize(0.043); axis->SetTitleSize(0.052);
    }
    h.GetXaxis()->SetTitleOffset(1.10); h.GetYaxis()->SetTitleOffset(1.22);
    h.Draw(i==0 ? "E1" : "E1 SAME");
  }
  TLatex label; label.SetNDC(); label.SetTextFont(62); label.SetTextSize(0.038);
  label.DrawLatex(0.255,0.60,"O#minusO  #sqrt{s_{NN}} = 5.36 TeV");
  label.DrawLatex(0.255,0.54,Form("Cent: %.0f#minus%.0f%%",plots.low,plots.high));
  label.DrawLatex(0.255,0.48,jetScan ? "p_{T}^{ch} > 2.0 GeV" : "40 < p_{T}^{jet} < 80 GeV");
  TLegend legend(0.215,0.155,0.83,0.44);
  legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42);
  legend.SetTextSize(0.038); legend.SetMargin(0.085);
  for (int i=0; i<4; ++i) {
    TString entry = jetScan ? Form("%.0f < p_{T}^{jet} < %.0f GeV",kJetBins[i].first,kJetBins[i].second)
        : Form("p_{T}^{ch} > %.1f GeV",kTrackCuts[i]);
    legend.AddEntry(shapes[i].get(),entry,"p");
  }
  legend.DrawClone(); pad.RedrawAxis();
}
} // namespace RawOODataComparisons

// Run from the repository root:
// root -l -b -q 'plotting/plotRawOODataComparisons.C()'
// Data only, raw same-jet pairs, full-range unit-area normalization as in
// EECDrawer::DrawEnergyEnergyCorrelators. No subtraction or unfolding.
void plotRawOODataComparisons(
    const char* inputName="projected/merged_eecOO2025_processed_121026.root",
    const char* outputDirectory="figures") {
  using namespace RawOODataComparisons;
  try {
    gROOT->SetBatch(true); gStyle->SetOptStat(0); gStyle->SetEndErrorSize(3);
    std::unique_ptr<TFile> input(TFile::Open(inputName,"READ"));
    Require(input && !input->IsZombie() && !input->TestBit(TFile::kRecovered),
        Form("Cannot open intact data file %s",inputName));
    // The manager owns its EECCard; allocate the large manager on the heap.
    auto manager = std::make_unique<EECHistogramManager>(input.get());
    auto* card = manager->GetCard();
    Require(card->GetDataType().Contains("OO") && !card->GetDataType().Contains("MC"),
        Form("Expected OO data, found %s",card->GetDataType().Data()));
    std::vector<CentralityPlots> all;
    double smallest=std::numeric_limits<double>::max(), largest=0.;
    for (int centrality=0; centrality<card->GetNCentralityBins(); ++centrality) {
      CentralityPlots plots;
      plots.low=card->GetLowBinBorderCentrality(centrality);
      plots.high=card->GetHighBinBorderCentrality(centrality);
      for (int i=0; i<4; ++i) {
        plots.jets[i]=ReadShape(*manager,centrality,kJetBins[i],2.,Form("rawOO_C%d_jetScan%d",centrality,i));
        plots.tracks[i]=ReadShape(*manager,centrality,kJetBins[0],kTrackCuts[i],Form("rawOO_C%d_trackScan%d",centrality,i));
        for (const auto* shape : {plots.jets[i].get(),plots.tracks[i].get()})
          for (int bin=1; bin<=shape->GetNbinsX(); ++bin) {
            const double x=shape->GetBinCenter(bin), value=shape->GetBinContent(bin);
            if (x<kDrawLow || x>kDrawHigh || value<=0.) continue;
            smallest=std::min(smallest,value); largest=std::max(largest,value+shape->GetBinError(bin));
          }
      }
      all.push_back(std::move(plots));
    }
    Require(largest>0. && smallest<std::numeric_limits<double>::max(),"No positive data to draw");
    const double yLow=smallest*0.35, yHigh=largest*3.; // Common axes for all centralities.
    Require(gSystem->mkdir(outputDirectory,true)==0 || !gSystem->AccessPathName(outputDirectory),
        "Cannot create figure directory");
    const TString book=Form("%s/rawEEC_OOData_allCentralities.pdf",outputDirectory);
    // Finish individual PDFs before opening ROOT's shared multipage PDF context.
    for (int pass=0; pass<2; ++pass) {
    for (size_t index=0; index<all.size(); ++index) {
      auto& plots=all[index];
      TCanvas canvas(Form("rawOOData_C%zu",index),"Raw EEC in OO data",1440,700);
      canvas.SetFillColor(kWhite);
      TPad left(Form("rawOOLeft_C%zu",index),"Jet momentum comparison",0.,0.055,0.5,0.925);
      TPad right(Form("rawOORight_C%zu",index),"Track threshold comparison",0.5,0.055,1.,0.925);
      left.Draw(); right.Draw();
      DrawPanel(left,plots,true,yLow,yHigh); DrawPanel(right,plots,false,yLow,yHigh);
      canvas.cd();
      TLatex heading; heading.SetNDC(); heading.SetTextFont(42); heading.SetTextAlign(22);
      heading.SetTextSize(0.034);
      heading.DrawLatex(0.5,0.965,"Raw energy-energy correlators in O#minusO data");
      heading.SetTextSize(0.023);
      heading.DrawLatex(0.5,0.025,"Each curve normalized to unit area over the full histogram range; statistical uncertainties");
      canvas.Modified(); canvas.Update();
      const TString stem=Form("%s/rawEEC_OOData_C%.0f-%.0f",outputDirectory,plots.low,plots.high);
      if (pass==0) {
        canvas.Print(stem+".pdf"); canvas.Print(stem+".png");
        std::cout << "RAW_OO_DATA_FIGURE_PASS centrality=" << plots.low << "-" << plots.high
            << " pdf=" << stem << ".pdf png=" << stem << ".png\n";
      } else {
        if (index==0) canvas.Print(book+"[");
        canvas.Print(book);
        if (index+1==all.size()) canvas.Print(book+"]");
      }
    }
    }
    std::cout << "RAW_OO_DATA_ALL_PASS centralities=" << all.size()
        << " input=" << inputName << " output=" << outputDirectory << std::endl;
  } catch (const std::exception& error) {
    std::cerr << "RAW_OO_DATA_PLOT_FAILED " << error.what() << std::endl;
    gSystem->Exit(1);
  }
}

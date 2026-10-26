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

namespace OOMCGenRecoTrackSelections {
const std::array<double,4> kThresholds{{0.5,1.,2.,3.}};
void Require(bool condition,const TString& message) {
  if(!condition) throw std::runtime_error(message.Data());
}
int Scalar(TFile& input,const char* name) {
  auto* value=dynamic_cast<TVectorT<float>*>(input.Get(Form("JCard/%s",name)));
  Require(value && value->GetNrows()==1,Form("Missing JCard/%s",name));
  return std::lround((*value)[value->GetLwb()]);
}
void Draw(TCanvas& canvas,TH1D& shape,double threshold,double low,double high) {
  auto* pad=new TPad(Form("pad_%s",canvas.GetName()),"",0.,0.14,1.,0.82);
  pad->Draw(); pad->cd(); pad->SetLogx(); pad->SetLogy(); pad->SetTicks(1,1);
  pad->SetLeftMargin(0.13); pad->SetRightMargin(0.04);
  pad->SetTopMargin(0.035); pad->SetBottomMargin(0.135);
  auto* frame=pad->DrawFrame(0.0008,low,0.8,high);
  frame->GetXaxis()->SetTitle("#Deltar"); frame->GetYaxis()->SetTitle("EEC");
  for(auto* axis:{frame->GetXaxis(),frame->GetYaxis()}) {
    axis->SetLabelFont(42); axis->SetTitleFont(42);
    axis->SetLabelSize(0.048); axis->SetTitleSize(0.058); axis->CenterTitle();
  }
  frame->GetXaxis()->SetTitleOffset(1.05); frame->GetYaxis()->SetTitleOffset(0.95);
  shape.SetLineColor(kBlack); shape.SetMarkerColor(kBlack);
  shape.SetMarkerStyle(20); shape.SetMarkerSize(0.65); shape.SetLineWidth(1);
  shape.Draw("E1 SAME"); pad->RedrawAxis(); canvas.cd();
  TLatex text; text.SetNDC(); text.SetTextAlign(22); text.SetTextFont(42);
  text.SetTextSize(0.031);
  text.DrawLatex(0.5,0.965,"Raw EEC in Pythia+Hijing O#minusO simulation");
  text.SetTextSize(0.028);
  text.DrawLatex(0.5,0.918,"GenReco: generator jets, reconstructed tracks; #sqrt{s_{NN}} = 5.36 TeV");
  text.SetTextSize(0.026);
  text.DrawLatex(0.5,0.868,Form("Cent: 0#minus30%%    40 < p_{T}^{jet} < 80 GeV    p_{T}^{ch} > %.1f GeV",threshold));
  TLegend legend(0.33,0.068,0.68,0.125);
  legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42); legend.SetTextSize(0.025);
  legend.AddEntry(&shape,"All combinations","lp"); legend.DrawClone();
  text.SetTextColor(kRed+1); text.SetTextSize(0.022);
  text.DrawLatex(0.5,0.044,"Pythia/Hijing decomposition unavailable: reconstructed-track subevent labels were not stored");
  text.SetTextColor(kBlack); text.SetTextSize(0.019);
  text.DrawLatex(0.5,0.017,"Raw same-jet pairs normalized to unit area over the full histogram range; statistical uncertainties");
  canvas.Modified(); canvas.Update();
}
} // namespace OOMCGenRecoTrackSelections

// Run from the repository root:
// root -l -b -q 'plotting/plotOOMCGenRecoTrackSelections.C()'
// This input resolves the low-pT selections, but has no reconstructed-track
// truth-origin labels. Draw only the available total; never invent components.
void plotOOMCGenRecoTrackSelections(
    const char* inputName="merged_OOMC_GenReco_20261022_processed.root",
    const char* outputDirectory="figures") {
  using namespace OOMCGenRecoTrackSelections;
  try {
    gROOT->SetBatch(true); gStyle->SetOptStat(0); gStyle->SetEndErrorSize(2);
    std::unique_ptr<TFile> input(TFile::Open(inputName,"READ"));
    Require(input && !input->IsZombie() && !input->TestBit(TFile::kRecovered),"Cannot read GenReco input");
    Require(Scalar(*input,"DataType")==10 && Scalar(*input,"McCorrelationType")==2,"Expected OOMC GenReco input");
    auto manager=std::make_unique<EECHistogramManager>(input.get());
    auto* card=manager->GetCard();
    const int centrality=card->FindBinIndexCentrality(0.,30.);
    const int jet=card->FindBinIndexJetPtEEC(40.,80.);
    Require(centrality>=0 && jet>=0,"Requested 0-30% / 40-80 GeV bin is missing");
    std::array<std::unique_ptr<TH1D>,4> shapes;
    double smallest=std::numeric_limits<double>::max(), largest=0., previous=std::numeric_limits<double>::max();
    for(size_t i=0;i<kThresholds.size();++i) {
      const int track=card->GetBinIndexTrackPtEEC(kThresholds[i]);
      Require(track>=0 && std::abs(card->GetLowBinBorderTrackPtEEC(track)-kThresholds[i])<1.e-6,"Missing track threshold");
      auto* source=manager->GetHistogramEnergyEnergyCorrelator(
          EECHistogramManager::kEnergyEnergyCorrelator,centrality,jet,track,EECHistograms::kSameJetPair);
      Require(source,"Missing raw GenReco EEC projection");
      shapes[i].reset(static_cast<TH1D*>(source->Clone(Form("OOMC_GenReco_T%zu",i))));
      shapes[i]->SetDirectory(nullptr); shapes[i]->GetXaxis()->SetRange(0,0);
      const double area=shapes[i]->Integral("width");
      Require(std::isfinite(area) && area>0. && area<previous,"Track selections do not have distinct decreasing raw yields");
      previous=area; shapes[i]->Scale(1./area);
      Require(std::abs(shapes[i]->Integral("width")-1.)<1.e-10,"Unit-area normalization failed");
      for(int bin=0;bin<shapes[i]->GetNcells();++bin) {
        const double value=shapes[i]->GetBinContent(bin), error=shapes[i]->GetBinError(bin);
        Require(std::isfinite(value) && std::isfinite(error),"Nonfinite histogram content/error");
        const double x=shapes[i]->GetBinCenter(bin);
        if(bin<1 || bin>shapes[i]->GetNbinsX() || x<0.0008 || x>0.8 || value<=0.) continue;
        smallest=std::min(smallest,value); largest=std::max(largest,value+error);
      }
      std::cout<<"OOMC_GENRECO_SELECTION_PASS track>"<<kThresholds[i]<<" source="<<source->GetName()
          <<" raw_area="<<area<<" normalized_area="<<shapes[i]->Integral("width")<<std::endl;
    }
    Require(largest>0. && smallest<std::numeric_limits<double>::max(),"No positive EEC values");
    Require(gSystem->mkdir(outputDirectory,true)==0 || !gSystem->AccessPathName(outputDirectory),"Cannot create figure directory");
    const TString prefix=Form("%s/OOMC_GenReco_totalEEC_C0-30_J40-80",outputDirectory);
    const TString book=prefix+"_allTracks.pdf";
    for(int pass=0;pass<2;++pass) {
      for(size_t i=0;i<kThresholds.size();++i) {
        TCanvas canvas(Form("OOMC_GenReco_T%zu_pass%d",i,pass),"OOMC GenReco track selection",1100,850);
        Draw(canvas,*shapes[i],kThresholds[i],smallest*0.35,largest*2.5);
        TString cut=Form("%.1f",kThresholds[i]); cut.ReplaceAll(".","p");
        const TString stem=prefix+"_T"+cut;
        if(pass==0) {canvas.Print(stem+".pdf"); canvas.Print(stem+".png");}
        else {if(i==0)canvas.Print(book+"[");canvas.Print(book);if(i+1==kThresholds.size())canvas.Print(book+"]");}
      }
    }
    std::cout<<"OOMC_GENRECO_PLOTS_PASS input="<<inputName<<" plots=4 centrality=0-30 jet=40-80 decomposition=unavailable"<<std::endl;
  } catch(const std::exception& error) {
    std::cerr<<"OOMC_GENRECO_PLOT_FAILED "<<error.what()<<std::endl;gSystem->Exit(1);
  }
}

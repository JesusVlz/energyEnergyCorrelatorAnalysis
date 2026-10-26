#include "EECCard.h"
R__LOAD_LIBRARY(plotting/DrawingClasses.so)

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <THnBase.h>
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
#include <vector>

namespace OOMCSubevents {
constexpr double kCentralityLow=0., kCentralityHigh=30.;
constexpr double kJetLow=40., kJetHigh=80.;
constexpr double kDrawLow=0.0008, kDrawHigh=0.8;
const std::array<double,4> kTrackCuts{{0.5,1.,2.,3.}};
const std::array<int,4> kColors{{kBlack,kBlue+1,kRed+1,kGreen+2}};
const std::array<const char*,4> kLabels{{
    "All combinations","Pythia#minusPythia","Pythia#minusHijing","Hijing#minusHijing"}};
const std::array<const char*,4> kNames{{"All","PythiaPythia","PythiaHijing","HijingHijing"}};
struct CutPlots {
  double trackCut;
  bool inputLimited=false;
  std::array<std::unique_ptr<TH1D>,4> histograms;
};
void Require(bool condition,const TString& message) {
  if(!condition) throw std::runtime_error(message.Data());
}
int Scalar(TFile& file,const char* name) {
  auto* value=dynamic_cast<TVectorT<float>*>(file.Get(Form("JCard/%s",name)));
  Require(value && value->GetNrows()==1,Form("Missing scalar JCard/%s",name));
  return std::lround((*value)[value->GetLwb()]);
}
void SelectExactRange(TAxis& axis,double low,double high,const char* label) {
  const int first=axis.FindFixBin(low+1.e-6), last=axis.FindFixBin(high-1.e-6);
  Require(first>=1 && last<=axis.GetNbins() && first<=last &&
      std::abs(axis.GetBinLowEdge(first)-low)<1.e-6 &&
      std::abs(axis.GetBinUpEdge(last)-high)<1.e-6,
      Form("Cannot select exact %s range %.1f-%.1f from the stored axis",label,low,high));
  axis.SetRange(first,last);
}
std::unique_ptr<TH1D> Project(THnBase& sparse,const TString& name,int firstSubevent,int lastSubevent) {
  sparse.GetAxis(5)->SetRange(firstSubevent,lastSubevent);
  std::unique_ptr<TH1D> result(sparse.Projection(0,"E"));
  Require(result!=nullptr,"Cannot project raw OOMC EEC");
  result->SetName(name); result->SetDirectory(nullptr);
  // Follow EECHistogramManager::FindHistogram's density convention.
  result->Scale(1.,"width");
  for(int bin=0;bin<result->GetNcells();++bin)
    Require(std::isfinite(result->GetBinContent(bin)) && std::isfinite(result->GetBinError(bin)),
        Form("Nonfinite projected content/error in %s",name.Data()));
  return result;
}
CutPlots MakeCut(THnBase& sparse,double threshold,int index) {
  auto* trackAxis=sparse.GetAxis(2);
  const int firstTrack=trackAxis->FindFixBin(threshold+1.e-6);
  Require(firstTrack>=1 && firstTrack<=trackAxis->GetNbins() &&
      std::abs(trackAxis->GetBinLowEdge(firstTrack)-threshold)<1.e-6,
      "Requested particle threshold is unavailable");
  // Axis 2 stores the lower-pT particle in each pair: integrate above threshold,
  // including its overflow, exactly as the standard EEC projection macro does.
  trackAxis->SetRange(firstTrack,trackAxis->GetNbins()+1);
  CutPlots plots; plots.trackCut=threshold;
  const std::array<int,4> first{{0,1,2,4}}, last{{0,1,3,4}};
  for(int component=0;component<4;++component)
    plots.histograms[component]=Project(sparse,Form("OOMC_C0to30_J40to80_T%d_%s",index,kNames[component]),
        first[component],last[component]);
  // Raw ordered mixed bins 1 and 2 are combined once. Processed S1/S2 must not
  // be summed: the existing manager already puts both orderings in either one.
  auto& total=*plots.histograms[0];
  for(int bin=0;bin<total.GetNcells();++bin) {
    double sum=0., sumError2=0.;
    for(int component=1;component<4;++component) {
      sum+=plots.histograms[component]->GetBinContent(bin);
      sumError2+=std::pow(plots.histograms[component]->GetBinError(bin),2);
    }
    const double value=total.GetBinContent(bin), error2=std::pow(total.GetBinError(bin),2);
    Require(std::abs(value-sum)<=1.e-9*std::max(1.,std::abs(value)),
        Form("Subevent content closure failed at track %.1f, bin %d",threshold,bin));
    Require(std::abs(error2-sumError2)<=1.e-9*std::max(1.,error2),
        Form("Subevent variance closure failed at track %.1f, bin %d",threshold,bin));
  }
  const double integral=total.Integral("width");
  Require(std::isfinite(integral) && integral>0.,"Empty selected OOMC distribution");
  for(auto& histogram:plots.histograms) histogram->Scale(1./integral);
  Require(std::abs(total.Integral("width")-1.)<1.e-10,"Total EEC area normalization failed");
  std::cout<<"OOMC_SUBEVENT_CLOSURE_PASS track>"<<threshold<<" raw_total="<<integral;
  double componentArea=0.;
  for(int component=1;component<4;++component) {
    const double fraction=plots.histograms[component]->Integral("width");
    Require(fraction>0.,Form("Empty %s component at track %.1f",kNames[component],threshold));
    componentArea+=fraction;
    std::cout<<" "<<kNames[component]<<"_fraction="<<fraction;
  }
  Require(std::abs(componentArea-1.)<1.e-10,"Component area fractions do not sum to one");
  std::cout<<std::endl;
  sparse.GetAxis(5)->SetRange(0,0);
  return plots;
}
void Draw(CutPlots& plots,TCanvas& canvas,double yLow,double yHigh,const TString& mode) {
  canvas.SetFillColor(kWhite);
  auto* pad=new TPad(Form("OOMC_plot_%s",canvas.GetName()),"",0.,0.14,1.,0.84);
  pad->Draw(); pad->cd();
  pad->SetLeftMargin(0.12); pad->SetRightMargin(0.035);
  pad->SetTopMargin(0.025); pad->SetBottomMargin(0.135);
  pad->SetTicks(1,1); pad->SetLogx(); pad->SetLogy();
  auto* frame=pad->DrawFrame(kDrawLow,yLow,kDrawHigh,yHigh);
  frame->SetTitle(""); frame->GetXaxis()->SetTitle("#Deltar"); frame->GetYaxis()->SetTitle("EEC");
  for(auto* axis:{frame->GetXaxis(),frame->GetYaxis()}) {
    axis->SetLabelFont(42); axis->SetTitleFont(42);
    axis->SetLabelSize(0.048); axis->SetTitleSize(0.058); axis->CenterTitle();
  }
  frame->GetXaxis()->SetTitleOffset(1.02); frame->GetYaxis()->SetTitleOffset(0.92);
  for(int component=0;component<4;++component) {
    auto& histogram=*plots.histograms[component];
    histogram.SetLineColor(kColors[component]); histogram.SetMarkerColor(kColors[component]);
    histogram.SetLineWidth(1); histogram.SetMarkerStyle(1); histogram.SetMarkerSize(0.5);
    histogram.Draw("E1 SAME");
  }
  pad->RedrawAxis(); canvas.cd();
  TLatex text; text.SetNDC(); text.SetTextAlign(22); text.SetTextFont(42);
  text.SetTextSize(0.031);
  text.DrawLatex(0.5,0.965,"Expected background in Pythia+Hijing O#minusO simulation");
  text.SetTextSize(0.027);
  text.DrawLatex(0.5,0.915,Form("Pythia+Hijing %s,  #sqrt{s_{NN}} = 5.36 TeV",mode.Data()));
  text.SetTextSize(0.026);
  if(plots.inputLimited)
    text.DrawLatex(0.5,0.867,"Cent: 0#minus30%    40 < p_{T}^{jet} < 80 GeV    stored gen. p_{T}^{ch} #geq 5 GeV");
  else
    text.DrawLatex(0.5,0.867,Form("Cent: 0#minus30%%     40 < p_{T}^{jet} < 80 GeV     p_{T}^{ch} > %.1f GeV",plots.trackCut));
  TLegend legend(0.105,0.046,0.985,0.108);
  legend.SetNColumns(4); legend.SetBorderSize(0); legend.SetFillStyle(0);
  legend.SetTextFont(42); legend.SetTextSize(0.022); legend.SetMargin(0.22);
  for(int component=0;component<4;++component)
    legend.AddEntry(plots.histograms[component].get(),kLabels[component],"l");
  legend.DrawClone();
  text.SetTextSize(0.018);
  text.DrawLatex(0.5,0.019,plots.inputLimited
      ? "Stored generator-particle range only; common total normalization; thresholds below 5 GeV cannot be resolved"
      : "Raw same-jet pairs; all components share the total EEC unit-area normalization; statistical uncertainties");
  canvas.Modified(); canvas.Update();
}
} // namespace OOMCSubevents

// Run from the repository root:
// root -l -b -q 'plotting/plotOOMCSubeventDecomposition.C()'
// Defaults to RecoGen (reconstructed jets, generator-level charged particles),
// matching the reference. Third argument true plots the verified stored-gen-pT>=5
// diagnostic only, without pretending to resolve the lower requested thresholds.
// Legacy Hydjet enum names encode the Hijing component
// for this OOMC input; use explicit Hijing labels in the resulting figures.
void plotOOMCSubeventDecomposition(
    const char* inputName="merged_OOMC_RecoGen.root",
    const char* outputDirectory="figures", bool useStoredGeneratorRange=false) {
  using namespace OOMCSubevents;
  try {
    gROOT->SetBatch(true); gStyle->SetOptStat(0); gStyle->SetEndErrorSize(2);
    std::unique_ptr<TFile> input(TFile::Open(inputName,"READ"));
    Require(input && !input->IsZombie() && !input->TestBit(TFile::kRecovered),"Cannot open intact OOMC file");
    Require(Scalar(*input,"DataType")==10,"Expected OOMC DataType 10");
    const int correlation=Scalar(*input,"McCorrelationType");
    Require(correlation==1 || correlation==3,"Truth subevents require RecoGen or GenGen MC");
    const TString mode=correlation==1 ? "RecoGen" : "GenGen";
    EECCard card(input.get());
    Require(card.FindBinIndexCentrality(0.,30.)>=0,"Input card has no exact 0-30% centrality bin");
    auto* sparse=dynamic_cast<THnBase*>(input->Get("energyEnergyCorrelator"));
    Require(sparse && sparse->GetNdimensions()==7,"Expected raw 7D OOMC EEC histogram");
    for(int axis=0;axis<7;++axis) sparse->GetAxis(axis)->SetRange(0,0);
    SelectExactRange(*sparse->GetAxis(1),kJetLow,kJetHigh,"jet pT");
    SelectExactRange(*sparse->GetAxis(3),kCentralityLow,kCentralityHigh,"centrality");
    Require(sparse->GetAxis(5)->GetNbins()==4,"Unexpected raw subevent axis");
    sparse->GetAxis(4)->SetRange(1,1); // Same-jet pairing type 0.
    const int weight=card.FindWeightExponentIndex(1.);
    Require(weight>=1 && weight<=sparse->GetAxis(6)->GetNbins(),"Nominal energy exponent 1 unavailable");
    sparse->GetAxis(6)->SetRange(weight,weight);
    if(!useStoredGeneratorRange) {
      // Reject a generator skim that cannot resolve the requested low-pT cuts.
      sparse->GetAxis(2)->SetRange(sparse->GetAxis(2)->FindFixBin(kTrackCuts.front()+1.e-6),
          sparse->GetAxis(2)->FindFixBin(kTrackCuts.back()-1.e-6));
      std::unique_ptr<TH1D> lowPtCheck(sparse->Projection(2,"E"));
      lowPtCheck->SetDirectory(nullptr);
      Require(lowPtCheck->Integral()>0.,
          "Input cannot resolve low-pT cuts: current OOMC generator skim starts at 5 GeV");
    } else {
      auto* track=dynamic_cast<THnBase*>(input->Get("track"));
      Require(track!=nullptr,"Cannot verify generator input particle range");
      for(int axis=0;axis<track->GetNdimensions();++axis)track->GetAxis(axis)->SetRange(0,0);
      std::unique_ptr<TH1D> spectrum(track->Projection(0,"E"));spectrum->SetDirectory(nullptr);
      double minimum=-1.;
      for(int bin=1;bin<=spectrum->GetNbinsX();++bin)
        if(spectrum->GetBinContent(bin)>0.) {minimum=spectrum->GetXaxis()->GetBinLowEdge(bin);break;}
      Require(std::abs(minimum-5.)<1.e-6,"This labeled diagnostic requires a verified 5 GeV generator input floor");
      std::cout<<"OOMC_STORED_GENERATOR_RANGE minimum_pt="<<minimum<<std::endl;
    }
    const std::vector<double> thresholds=useStoredGeneratorRange
        ? std::vector<double>{0.5} : std::vector<double>(kTrackCuts.begin(),kTrackCuts.end());
    std::vector<CutPlots> cuts;
    double smallest=std::numeric_limits<double>::max(), largest=0.;
    for(size_t index=0;index<thresholds.size();++index) {
      cuts.push_back(MakeCut(*sparse,thresholds[index],index));
      cuts.back().inputLimited=useStoredGeneratorRange;
      for(const auto& histogram:cuts.back().histograms)
        for(int bin=1;bin<=histogram->GetNbinsX();++bin) {
          const double x=histogram->GetBinCenter(bin), value=histogram->GetBinContent(bin);
          if(x<kDrawLow || x>kDrawHigh || value<=0.) continue;
          smallest=std::min(smallest,value);
          largest=std::max(largest,value+histogram->GetBinError(bin));
        }
    }
    Require(largest>0. && smallest<std::numeric_limits<double>::max(),"No positive OOMC data to draw");
    const double yLow=smallest*0.35, yHigh=largest*2.5;
    Require(gSystem->mkdir(outputDirectory,true)==0 || !gSystem->AccessPathName(outputDirectory),
        "Cannot create figure directory");
    TString prefix=Form("%s/OOMC_%s_subevents_C0-30_J40-80",outputDirectory,mode.Data());
    if(useStoredGeneratorRange)prefix+="_storedGenPt5";
    const TString book=prefix+"_allTracks.pdf";
    // Avoid interleaving individual PDFs with ROOT's shared multipage context.
    for(int pass=0;pass<(useStoredGeneratorRange ? 1 : 2);++pass) {
      for(size_t index=0;index<cuts.size();++index) {
        TCanvas canvas(Form("OOMC_%s_T%zu_pass%d",mode.Data(),index,pass),"OOMC subevent decomposition",1100,850);
        Draw(cuts[index],canvas,yLow,yHigh,mode);
        TString threshold=Form("%.1f",useStoredGeneratorRange ? 5. : cuts[index].trackCut); threshold.ReplaceAll(".","p");
        const TString stem=prefix+"_T"+threshold;
        if(pass==0) {
          canvas.Print(stem+".pdf"); canvas.Print(stem+".png");
          std::cout<<"OOMC_SUBEVENT_FIGURE_PASS particle_input_min="<<(useStoredGeneratorRange ? 5. : cuts[index].trackCut)<<" pdf="<<stem<<".pdf\n";
        } else {
          if(index==0) canvas.Print(book+"[");
          canvas.Print(book);
          if(index+1==cuts.size()) canvas.Print(book+"]");
        }
      }
    }
    std::cout<<"OOMC_SUBEVENT_ALL_PASS mode="<<mode<<" centrality=0-30 jet=40-80 stored_generator_range="<<useStoredGeneratorRange<<" input="
        <<inputName<<std::endl;
  } catch(const std::exception& error) {
    std::cerr<<"OOMC_SUBEVENT_PLOT_FAILED "<<error.what()<<std::endl;
    gSystem->Exit(1);
  }
}

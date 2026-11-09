#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TH3D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TPad.h>
#include <TStyle.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Run3Plot {
std::unique_ptr<TH1D> Project(TH3D* source,int stage,int low,int high,const std::string& name) {
  std::unique_ptr<TH1D> h(source->ProjectionX(name.c_str(),stage+1,stage+1,low,high,"e"));
  h->SetDirectory(nullptr); h->Rebin(5); h->Scale(1.,"width");
  return h;
}
void Draw(TH3D* source,const std::vector<int>& stages,int low,int high,const std::string& output,bool ratios) {
  TCanvas canvas("run3_qa_canvas","Run-3 jet QA",900,900);
  TPad upper("top","",0,ratios?.29:0,1,1),lower("bottom","",0,0,1,.29);
  upper.Draw(); upper.cd(); upper.SetLogy();
  TLegend legend(.34,.65,.88,.89); legend.SetBorderSize(0); legend.SetTextSize(.025);
  std::vector<std::unique_ptr<TH1D>> histograms;
  std::vector<std::unique_ptr<TGraph>> ratioGraphs;
  std::unique_ptr<TH1D> ratioFrame;
  int colors[]={1,2,4,8,6,7};
  for(size_t j=0;j<stages.size();++j) {
    auto h=Project(source,stages[j],low,high,"qa_stage_"+std::to_string(j));
    h->SetTitle(""); h->SetLineColor(colors[j]); h->SetLineWidth(2); h->SetMinimum(1e-3);
    h->GetXaxis()->SetRangeUser(0,600);
    h->GetYaxis()->SetTitle("Weighted jet count / GeV (common exposure)");
    legend.AddEntry(h.get(),source->GetYaxis()->GetBinLabel(stages[j]+1),"l");
    histograms.push_back(std::move(h));
  }
  double maximum=0;
  for(const auto& h:histograms) maximum=std::max(maximum,h->GetMaximum());
  if(maximum<=0) { std::cerr<<"No positive entries for "<<output<<"; no plot produced\n"; return; }
  histograms[0]->SetMaximum(maximum*10);
  for(size_t j=0;j<histograms.size();++j) histograms[j]->Draw(j?"HIST SAME":"HIST");
  legend.Draw();
  if(ratios) {
    canvas.cd(); lower.Draw(); lower.cd(); lower.SetBottomMargin(.28); lower.SetGridy();
    ratioFrame.reset(dynamic_cast<TH1D*>(histograms[0]->Clone("ratio_frame")));
    ratioFrame->SetDirectory(nullptr); ratioFrame->Reset();
    ratioFrame->SetMinimum(0); ratioFrame->SetMaximum(1.1);
    ratioFrame->GetYaxis()->SetTitle("Stage / no jet-quality cuts");
    ratioFrame->GetXaxis()->SetTitle("Final calibrated jet p_{T} (GeV)");
    ratioFrame->GetXaxis()->SetTitleSize(.10);ratioFrame->GetXaxis()->SetLabelSize(.08);
    ratioFrame->GetYaxis()->SetTitleSize(.075);ratioFrame->GetYaxis()->SetLabelSize(.075);ratioFrame->GetYaxis()->SetTitleOffset(.6);
    ratioFrame->Draw("AXIS");
    for(size_t j=1;j<histograms.size();++j) {
      auto ratio=std::make_unique<TGraph>();
      ratio->SetMarkerColor(colors[j]); ratio->SetMarkerStyle(20+j); ratio->SetMarkerSize(.7);
      for(int b=1;b<=histograms[0]->GetNbinsX();++b) {
        const double denominator=histograms[0]->GetBinContent(b);
        if(denominator<=0) continue; // Empty/unsupported bins are omitted, never drawn as zero efficiency.
        const double value=histograms[j]->GetBinContent(b)/denominator;
        if(std::isfinite(value)) ratio->SetPoint(ratio->GetN(),histograms[0]->GetBinCenter(b),value);
      }
      ratio->Draw("P SAME");ratioGraphs.push_back(std::move(ratio));
    }
  }
  canvas.SaveAs(output.c_str());
}
}

// centrality range is half-open; use -1,0 for pp and e.g. 50,80 for OO.
// Histograms and calibration/selection metadata must come from the updated analyzer.
void plotRun3JetQA(const char* input,const char* prefix="run3",double centralityLow=-1,double centralityHigh=100) {
  if(!(centralityLow<centralityHigh)) throw std::invalid_argument("Empty centrality range");
  std::unique_ptr<TFile> file(TFile::Open(input,"READ"));
  auto* selection=file && !file->IsZombie()?dynamic_cast<TH3D*>(file->Get("run3JetQA/selectionPtWeighted")):nullptr;
  auto* calibration=file?dynamic_cast<TH3D*>(file->Get("run3JetQA/calibrationPtWeighted")):nullptr;
  if(!selection || !calibration) throw std::runtime_error("Stage spectra absent: rerun eecAnalysis with Run3JetQA 1; final inclusiveJet alone is insufficient");
  const int low=selection->GetZaxis()->FindFixBin(centralityLow),high=selection->GetZaxis()->FindFixBin(std::nextafter(centralityHigh,centralityLow));
  if(low<1 || high>selection->GetNbinsZ()) throw std::out_of_range("Centrality range outside QA support");
  for(const char* key:{"JECState","JetIDTable","SelectionDefinition","CalibrationDefinition","VetoState"}) {
    auto* meta=dynamic_cast<TNamed*>(file->Get((std::string("run3JetQA/")+key).c_str()));
    if(meta) std::cout<<key<<": "<<meta->GetTitle()<<'\n';
  }
  gStyle->SetOptStat(0);
  Run3Plot::Draw(selection,{0,1,3,5},low,high,std::string(prefix)+"_selection_cumulative.png",true);
  Run3Plot::Draw(selection,{0,1,2,4},low,high,std::string(prefix)+"_selection_independent.png",true);
  std::vector<int> active;
  for(int level=0;level<4;++level) {
    auto p=Run3Plot::Project(calibration,level,low,high,"level_check_"+std::to_string(level));
    if(p->Integral(0,p->GetNbinsX()+1)!=0) active.push_back(level);
  }
  if(!active.empty()) Run3Plot::Draw(calibration,active,low,high,std::string(prefix)+"_calibration_fixed_jets.png",false);
  std::cout<<"Ratios are central values only; no independent/binomial error assumption. Do not normalize each selection separately.\n";
}

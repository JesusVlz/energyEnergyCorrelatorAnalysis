// Native ROOT view of a derived n=0 residual pair payload, never an EEC ratio.
#include <TH2.h>
#include "JDrawer.h"
#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>
#include <TNamed.h>
#include <TVectorT.h>
#include <TROOT.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>

void plotRun3PairEfficiency(const char* input,const char* prefix,const char* system,
                           const char* selection,int activity=0,int jet=-1,
                           const char* sampleLabel="Candidate; full training sample") {
  gROOT->SetBatch(true);TH1::AddDirectory(false);
  std::unique_ptr<TFile> file(TFile::Open(input,"READ"));
  if(!file || file->IsZombie() || !file->Get("Run3PairCorrectionVersion"))throw std::runtime_error("Missing Run-3 pair payload");
  auto* collision=dynamic_cast<TNamed*>(file->Get("Run3PairCollisionSystem"));
  auto* workingPoint=dynamic_cast<TNamed*>(file->Get("Run3TrackSelection"));
  if(!collision || !workingPoint || std::string(collision->GetTitle())!=system ||
     std::string(workingPoint->GetTitle())!=selection)throw std::runtime_error("Pair payload system/selection mismatch");
  auto* pt=dynamic_cast<TVectorT<float>*>(file->Get("JCard/TrackPairPtBinEdges"));
  auto* je=dynamic_cast<TVectorT<float>*>(file->Get("JCard/JetPtBinEdgesEEC"));
  auto* ce=dynamic_cast<TVectorT<float>*>(file->Get("JCard/CentralityBinEdges"));
  if(!pt || !je || !ce || activity<0 || activity>=ce->GetNrows()-1 || jet>=je->GetNrows()-1)throw std::runtime_error("Missing/incompatible pair binning");
  auto bin=[&](double lo,double hi){for(int i=0;i<pt->GetNrows()-1;++i)if(std::abs((*pt)[pt->GetLwb()+i]-lo)<1e-6 && std::abs((*pt)[pt->GetLwb()+i+1]-hi)<1e-6)return i;throw std::runtime_error("Exact 6-8 / 3-4 particle cells unavailable; recount exports");};
  const int high=bin(6,8),low=bin(3,4);
  TString key=Form("trackPairEfficiencyCorrection%s_C%dT%dA%d",jet<0?"":"CloseToJet",activity,high,low);
  if(jet>=0)key+=Form("J%d",jet);
  auto* h=dynamic_cast<TH1*>(file->Get(key));auto* mask=dynamic_cast<TH1*>(file->Get(key+"_valid"));
  if(!h || !mask)throw std::runtime_error("Missing pair efficiency/support mask");
  TGraphErrors graph;graph.SetName("normalized_pair_efficiency");
  double minimum=.5,maximum=1.5;
  std::ofstream csv(std::string(prefix)+".csv");csv<<std::setprecision(17)<<"deltaR_low,deltaR_high,epsilon,supported,training_uncertainty\n";
  for(int b=1;b<=h->GetNbinsX();++b){
    const double value=h->GetBinContent(b);const bool valid=mask->GetBinContent(b)>0 && std::isfinite(value) && value>0;
    csv<<h->GetXaxis()->GetBinLowEdge(b)<<','<<h->GetXaxis()->GetBinUpEdge(b)<<','<<value<<','<<valid<<",NOT_EVALUATED\n";
    if(!valid || h->GetXaxis()->GetBinLowEdge(b)<=0)continue;
    const int p=graph.GetN();graph.SetPoint(p,h->GetXaxis()->GetBinCenter(b),value);
    graph.SetPointError(p,h->GetXaxis()->GetBinWidth(b)/2.,0.);
    minimum=std::min(minimum,value*.9);maximum=std::max(maximum,value*1.1);
  }
  TFile numerical((std::string(prefix)+".root").c_str(),"RECREATE");
  TNamed("source_payload",input).Write();TNamed("sample_scope",sampleLabel).Write();
  TNamed("definition","epsilon = single-track-corrected n=0 reco pair yield / (alpha * gen pair yield); alpha from strict 0.1 < DeltaR < 0.4. No EEC momentum weights.").Write();
  TNamed("uncertainty","Finite-training uncertainty unavailable; horizontal bars show bin widths, no statistical vertical errors claimed.").Write();
  for(const char* suffix:{"","_valid","_genPairs","_recoPairs","_controlSingleOverGen"}){auto* object=file->Get(key+suffix);if(object)object->Write();}
  graph.Write();
  if(!graph.GetN()){TNamed("plot_status","NO_SUPPORTED_CELLS; no empty placeholder figure generated").Write();numerical.Close();std::cout<<"PAIR_CURVE_UNSUPPORTED "<<key<<std::endl;return;}
  std::unique_ptr<TH1> axis(static_cast<TH1*>(h->Clone("axis")));axis->Reset();axis->SetMinimum(minimum);axis->SetMaximum(maximum);
  axis->GetXaxis()->SetRangeUser(.004,.8);
  static auto* owner=new JDrawer();auto& drawer=*owner;drawer.Reset();drawer.SetCanvasSize(850,780);
  drawer.SetLeftMargin(.17);drawer.SetRightMargin(.06);drawer.SetTopMargin(.20);drawer.SetBottomMargin(.16);drawer.SetLogX(true);
  drawer.DrawHistogram(axis.get(),"#DeltaR","Track pair efficiency","","AXIS");
  graph.SetLineColor(kBlue);graph.SetMarkerColor(kBlue);graph.SetMarkerStyle(20);graph.SetMarkerSize(.55);graph.Draw("P SAME");
  TLine unity(.004,1,.8,1);unity.SetLineStyle(2);unity.Draw();
  TLatex label;label.SetNDC();label.SetTextFont(42);label.SetTextSize(.03);
  label.DrawLatex(.17,.96,"CMS Simulation Preliminary");label.DrawLatex(.17,.916,Form("%s MC, %s; %s",system,selection,sampleLabel));
  label.SetTextSize(.026);label.DrawLatex(.17,.875,"6 < p_{T,t} < 8 GeV; 3 < p_{T,a} < 4 GeV");
  label.SetTextSize(.022);label.DrawLatex(.17,.84,"Yield ratio normalized at 0.1 < #DeltaR < 0.4; map uncertainty pending");
  const TString activityLabel=std::string(system)=="pp"?"Inclusive pp":Form("OO centrality %g#minus%g%%",(*ce)[ce->GetLwb()+activity],(*ce)[ce->GetLwb()+activity+1]);
  label.DrawLatex(.22,.74,activityLabel);
  label.DrawLatex(.22,.70,jet<0?"Jet-independent selected-event pairs":Form("Common reco jets: %g#minus%g GeV",(*je)[je->GetLwb()+jet],(*je)[je->GetLwb()+jet+1]));
  label.DrawLatex(.22,.66,"Unsupported cells omitted; horizontal bars = bin widths");
  auto* canvas=gPad->GetCanvas();canvas->SaveAs((std::string(prefix)+".png").c_str());canvas->SaveAs((std::string(prefix)+".pdf").c_str());
  numerical.cd();canvas->Write("canvas");TNamed("plot_status","GENERATED_CANDIDATE_NOT_QUALIFIED").Write();numerical.Close();delete canvas;
  std::cout<<"PAIR_CURVE_SAVED "<<graph.GetN()<<" "<<key<<std::endl;
}

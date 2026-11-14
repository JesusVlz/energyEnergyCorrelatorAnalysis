// Run-3 adaptation of plotting/trackClosurePlotter.C (original is not modified).
// Preserves its corrected/uncorrected/generated comparison, per-event track
// normalization and JDrawer split-canvas convention. Adds a separate raw-scale
// EEC tracking diagnostic on reconstructed jets. No background subtraction,
// unfolding, pair-map derivation, or invented cross-covariance is performed.
// Reviewed interfaces: energyEnergyCorrelatorAnalysis 6f7464a8.
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"
R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorT.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Run3TrackingPlot {
using Hist = std::unique_ptr<TH1D>;
void Require(bool ok,const TString& message) {
  if(!ok) throw std::runtime_error(message.Data());
}
bool Near(double a,double b,double rel=1.e-5) {
  return std::isfinite(a) && std::isfinite(b) &&
    std::abs(a-b)<=rel*std::max({1.,std::abs(a),std::abs(b)});
}
Hist Copy(const TH1D* source,const TString& name) {
  Require(source,TString("Missing histogram: ")+name);
  Hist h(dynamic_cast<TH1D*>(source->Clone(name)));
  Require(bool(h),"Expected a TH1D");
  h->SetDirectory(nullptr); h->GetXaxis()->SetRange(0,0);
  for(int b=0;b<h->GetNcells();++b)
    Require(std::isfinite(h->GetBinContent(b)) && std::isfinite(h->GetBinError(b)),
            TString("Nonfinite histogram: ")+name);
  return h;
}
void Compatible(const TH1D* a,const TH1D* b) {
  Require(a && b && a->GetNbinsX()==b->GetNbinsX(),"Histogram bin counts differ");
  for(int i=1;i<=a->GetNbinsX()+1;++i)
    Require(Near(a->GetXaxis()->GetBinLowEdge(i),b->GetXaxis()->GetBinLowEdge(i),1.e-10),
            "Histogram bin edges differ");
}
double Scalar(TFile& f,const char* key) {
  auto* v=dynamic_cast<TVectorT<float>*>(f.Get(TString("JCard/")+key));
  Require(v && v->GetNrows()==1,TString("Missing scalar JCard/")+key);
  const double x=(*v)[v->GetLwb()];
  Require(std::isfinite(x),TString("Nonfinite JCard/")+key); return x;
}
void SameCard(TFile& a,TFile& b,const char* key) {
  auto* va=dynamic_cast<TVectorT<float>*>(a.Get(TString("JCard/")+key));
  auto* vb=dynamic_cast<TVectorT<float>*>(b.Get(TString("JCard/")+key));
  Require(va && vb && va->GetNrows()==vb->GetNrows(),TString("Missing/incompatible JCard/")+key);
  for(int i=0;i<va->GetNrows();++i)
    Require(Near((*va)[va->GetLwb()+i],(*vb)[vb->GetLwb()+i]),TString("Mismatched JCard/")+key);
}
void CheckInputs(TFile& reco,TFile& gen) {
  Require(std::lround(Scalar(reco,"DataType"))==10 && std::lround(Scalar(gen,"DataType"))==10,
          "This bounded macro requires OO MC (DataType=10)");
  Require(std::lround(Scalar(reco,"McCorrelationType"))==0 &&
          std::lround(Scalar(gen,"McCorrelationType"))==1,
          "Supply RecoReco first and RecoGen second (not GenGen)");
  for(const char* key : {"MatchJets","TriggerSelection","JetType","JetAxis","JetEtaCut",
      "MinJetPtCut","MaxJetPtCut","CutBadPhi","MinMaxTrackPtFraction","MaxMaxTrackPtFraction",
      "JetUncertainty","TrackEtaCut","MinTrackPtCut","MaxTrackPtCut","MaxTrackPtRelativeError",
      "VertexMaxDistance","SubeventCut","TrackEfficiencyVariation","JetPtWeight","ZVertexCut",
      "LowPtHatCut","HighPtHatCut","MultiplicityMode","JetRadius","WeightExponent",
      "CentralityBinEdges","JetPtBinEdgesEEC","TrackPtBinEdgesEEC","TrackPtBinEdges"})
    SameCard(reco,gen,key);
  Require(Near(Scalar(reco,"JetRadius"),0.4),"Configured for R=0.4 only");
  Require(std::lround(Scalar(reco,"SubeventCut"))==2,
          "Use all-particle totals (SubeventCut=2), not a Pythia-only denominator");
  Require(std::lround(Scalar(reco,"DisableTrackPairEfficiencyCorrection"))==1 &&
          std::lround(Scalar(gen,"DisableTrackPairEfficiencyCorrection"))==1,
          "This first-pass macro expects residual pair corrections OFF");
}
// Same per-event convention as the Run-2 tracking plotter. Fail if event flow
// would be silently excluded from that normalization.
double EventSum(EECHistogramManager& m,int c) {
  auto h=Copy(m.GetHistogramMultiplicity(c),"selected_event_multiplicity");
  Require(h->GetBinContent(0)==0. && h->GetBinContent(h->GetNbinsX()+1)==0. &&
          h->GetBinError(0)==0. && h->GetBinError(h->GetNbinsX()+1)==0.,
          "Nonempty multiplicity flow: establish event normalization before plotting");
  const double sum=h->Integral(1,h->GetNbinsX(),"width");
  Require(sum>0. && std::isfinite(sum),"No positive selected-event normalization");
  return sum;
}
void CheckJets(EECHistogramManager& a,EECHistogramManager& b,int ca,int cb) {
  auto ha=Copy(a.GetHistogramJetPt(ca),"reco_analysis_jets");
  auto hb=Copy(b.GetHistogramJetPt(cb),"gen_reference_analysis_jets");
  Compatible(ha.get(),hb.get());
  for(int k=0;k<ha->GetNcells();++k)
    Require(Near(ha->GetBinContent(k),hb->GetBinContent(k),1.e-4) &&
            Near(ha->GetBinError(k),hb->GetBinError(k),1.e-4),
            Form("RecoReco/RecoGen selected-jet histograms differ in bin %d; check inputs/cuts/weights",k));
}
// Rebin density histograms through bin-integrated sums. Never rebin EECs.
void RebinDensity(TH1D& h,int factor) {
  if(factor<=1) return;
  Require(h.GetNbinsX()%factor==0,"Requested rebin leaves a remainder");
  for(int b=1;b<=h.GetNbinsX();++b) {
    const double w=h.GetBinWidth(b);
    h.SetBinContent(b,h.GetBinContent(b)*w); h.SetBinError(b,h.GetBinError(b)*w);
  }
  h.Rebin(factor);
  for(int b=1;b<=h.GetNbinsX();++b) {
    const double w=h.GetBinWidth(b);
    h.SetBinContent(b,h.GetBinContent(b)/w); h.SetBinError(b,h.GetBinError(b)/w);
  }
}
struct RatioData {
  std::unique_ptr<TGraphErrors> graph;
  Hist mask;
  std::vector<double> value,error;
};
RatioData Ratio(const TH1D& num,const TH1D& den,const TString& name) {
  Compatible(&num,&den);
  RatioData r; r.graph=std::make_unique<TGraphErrors>(); r.graph->SetName(name);
  r.mask=Copy(&den,name+"_valid"); r.mask->Reset();
  r.value.assign(den.GetNcells(),std::numeric_limits<double>::quiet_NaN());
  r.error=r.value;
  for(int k=1;k<=den.GetNbinsX();++k) {
    const double d=den.GetBinContent(k),n=num.GetBinContent(k);
    if(!(d>0.)) continue;
    const double v=n/d,e=std::hypot(num.GetBinError(k)/d,v*den.GetBinError(k)/d);
    Require(std::isfinite(v) && std::isfinite(e),"Nonfinite ratio");
    const int p=r.graph->GetN(); r.graph->SetPoint(p,den.GetBinCenter(k),v);
    r.graph->SetPointError(p,0.5*den.GetBinWidth(k),e);
    r.mask->SetBinContent(k,1.); r.value[k]=v; r.error[k]=e;
  }
  Require(r.graph->GetN()>0,"No positive denominator bins"); return r;
}

// Drawing layout adapted from plotTrackClosure() in the original macro.
void Draw(TH1D& gen,TH1D& reco,TH1D* uncorrected,const TString& stem,
          const TString& xTitle,const TString& yTitle,bool logX,bool logY,
          double xMin,double xMax,const TString& selection,const TString& detail,
          const TString& provenance) {
  Compatible(&gen,&reco); if(uncorrected) Compatible(&gen,uncorrected);
  Require(gSystem->AccessPathName(stem+".root"),"Output exists; use a new output directory: "+stem);
  auto ratio=Ratio(reco,gen,"corrected_over_gen");
  std::unique_ptr<RatioData> rawRatio;
  if(uncorrected) rawRatio=std::make_unique<RatioData>(Ratio(*uncorrected,gen,"uncorrected_over_gen"));
  std::vector<TH1D*> hist{&gen,&reco}; if(uncorrected) hist.push_back(uncorrected);
  int colors[]={kBlack,kRed+1,kBlue+1},markers[]={21,20,24};
  double ymax=0.,ymin=std::numeric_limits<double>::max(); bool negative=false;
  for(size_t i=0;i<hist.size();++i) {
    auto& h=*hist[i]; h.SetTitle(""); h.SetMarkerColor(colors[i]);h.SetLineColor(colors[i]);
    h.SetMarkerStyle(markers[i]);h.SetMarkerSize(0.85);h.SetLineWidth(1);
    h.GetXaxis()->SetRangeUser(xMin,xMax);
    for(int k=1;k<=h.GetNbinsX();++k) {
      if(h.GetBinCenter(k)<xMin || h.GetBinCenter(k)>xMax)continue;
      const double y=h.GetBinContent(k),e=h.GetBinError(k);
      ymax=std::max(ymax,y+e); if(y>0.)ymin=std::min(ymin,y);
      if(y<0.)negative=true;
    }
  }
  Require(ymax>0. && std::isfinite(ymin),"No positive data in drawing range");
  if(negative)logY=false;
  // Keep one process-lifetime drawer, matching the legacy macro's ownership.
  // Its destructor deletes both canvas and owned pads; do not invoke it here.
  static JDrawer* legacyDrawer=new JDrawer();
  JDrawer& drawer=*legacyDrawer;
  drawer.SetDefaultAppearanceSplitCanvas();drawer.SetCanvasSize(800,800);
  drawer.SetLeftMargin(0.16);drawer.SetRightMargin(0.04);drawer.SetTopMargin(0.09);
  drawer.SetTitleOffsetY(1.55);drawer.SetLogX(logX);drawer.SetLogY(logY);
  drawer.CreateSplitCanvas();
  gen.SetMinimum(logY?ymin*0.25:0.); gen.SetMaximum(logY?ymax*25.:ymax*1.7);
  if(negative) { double min=0.;for(auto* h:hist)min=std::min(min,h->GetMinimum());gen.SetMinimum(min*1.2); }
  drawer.DrawHistogramToUpperPad(&gen,xTitle,yTitle,"","E1");
  reco.Draw("E1 SAME"); if(uncorrected)uncorrected->Draw("E1 SAME");
  TLatex text;text.SetNDC();text.SetTextFont(42);text.SetTextSize(0.032);
  text.DrawLatex(0.19,0.82,selection);text.DrawLatex(0.19,0.77,detail);
  const bool ptSpectrum=uncorrected && logX;
  TLegend legend(ptSpectrum?0.19:0.48,ptSpectrum?0.19:0.55,
                 ptSpectrum?0.70:0.92,ptSpectrum?0.38:0.72);
  legend.SetFillStyle(0);legend.SetBorderSize(0);
  legend.SetTextSize(0.034);legend.SetTextFont(42);
  legend.AddEntry(&gen,"Generated particles (RecoGen)","lp");
  legend.AddEntry(&reco,"Corrected tracks (RecoReco)","lp");
  if(uncorrected)legend.AddEntry(uncorrected,"Uncorrected tracks","lp");legend.Draw();
  auto frame=Copy(&gen,"ratio_frame");frame->Reset();frame->GetXaxis()->SetRangeUser(xMin,xMax);
  double rmin=0.8,rmax=1.2;
  for(auto* r : {&ratio,rawRatio.get()}) if(r)
    for(int k=1;k<=gen.GetNbinsX();++k) {
      if(gen.GetBinCenter(k)<xMin || gen.GetBinCenter(k)>xMax || !std::isfinite(r->value[k]))continue;
      rmin=std::min(rmin,r->value[k]-r->error[k]);rmax=std::max(rmax,r->value[k]+r->error[k]);
    }
  const double span=rmax-rmin;frame->SetMinimum(rmin-0.08*span);frame->SetMaximum(rmax+0.08*span);
  drawer.SetLogY(false);drawer.SetGridY(true);
  drawer.DrawHistogramToLowerPad(frame.get(),xTitle,"Reco / Gen","","AXIS");
  TLine unity(xMin,1.,xMax,1.);unity.SetLineStyle(2);unity.Draw();
  ratio.graph->SetMarkerColor(colors[1]);ratio.graph->SetLineColor(colors[1]);ratio.graph->SetMarkerStyle(20);
  ratio.graph->SetMarkerSize(0.7);ratio.graph->Draw("P SAME");
  if(rawRatio) {
    rawRatio->graph->SetMarkerColor(colors[2]);rawRatio->graph->SetLineColor(colors[2]);
    rawRatio->graph->SetMarkerStyle(24);rawRatio->graph->SetMarkerSize(0.7);rawRatio->graph->Draw("P SAME");
  }
  auto* canvas=gPad->GetCanvas();canvas->cd();
  text.SetTextFont(62);text.SetTextSize(0.031);text.DrawLatex(0.16,0.967,"CMS");
  text.SetTextFont(52);text.SetTextSize(0.025);text.DrawLatex(0.26,0.967,"Simulation Preliminary");
  text.SetTextFont(42);text.SetTextSize(0.023);text.SetTextAlign(31);
  text.DrawLatex(0.96,0.967,"Pythia+Hijing OO, 5.36 TeV");text.SetTextAlign(11);
  text.SetTextSize(0.019);text.DrawLatex(0.16,0.937,"Ratio errors: independent-histogram approximation; no event covariance");
  if(!uncorrected) {
    text.SetTextSize(0.017);text.DrawLatex(0.16,0.914,"Same-jet totals; pair correction OFF; no background subtraction or unfolding");
  }
  canvas->Modified();canvas->Update();canvas->Print(stem+".pdf");canvas->Print(stem+".png");
  TFile output(stem+".root","CREATE");Require(!output.IsZombie(),"Cannot create ROOT output");
  gen.Write("generated");reco.Write("reconstructed_corrected");
  if(uncorrected)uncorrected->Write("reconstructed_uncorrected");
  ratio.graph->Write();ratio.mask->Write();if(rawRatio){rawRatio->graph->Write();rawRatio->mask->Write();}
  TNamed("provenance",provenance.Data()).Write();
  TNamed("uncertainty_model","Independent histogram Sumw2 propagation; normalizers fixed; cross-file/event covariance not supplied").Write();
  canvas->Write("canvas");output.Close();
  std::ofstream csv((stem+".csv").Data());Require(bool(csv),"Cannot create CSV");csv<<std::setprecision(17);
  csv<<"low,high,gen,gen_error,reco,reco_error,ratio,ratio_error,valid";
  if(uncorrected)csv<<",uncorrected,uncorrected_error,uncorrected_ratio,uncorrected_ratio_error";csv<<'\n';
  for(int k=1;k<=gen.GetNbinsX();++k){
    csv<<gen.GetBinLowEdge(k)<<','<<gen.GetXaxis()->GetBinUpEdge(k)<<','<<gen.GetBinContent(k)<<','<<gen.GetBinError(k)
       <<','<<reco.GetBinContent(k)<<','<<reco.GetBinError(k)<<','<<ratio.value[k]<<','<<ratio.error[k]<<','<<ratio.mask->GetBinContent(k);
    if(uncorrected)csv<<','<<uncorrected->GetBinContent(k)<<','<<uncorrected->GetBinError(k)<<','<<rawRatio->value[k]<<','<<rawRatio->error[k];csv<<'\n';
  }
  std::cout<<"RUN3_TRACKING_PLOT_SAVED "<<stem<<'\n';
  // JDrawer owns ROOT pad objects through its canvas; delete only the canvas.
  delete canvas;
}
} // namespace Run3TrackingPlot

// Call from repository root. Files are projected, ONE energy power per file.
// Angular track plots use an exact 2-3 GeV slice by default, NOT the EEC >=2
// threshold. Single-track pT is displayed over 2-100 GeV. EEC reporting bins are
// exactly 40-80, 80-120, 120-160 GeV, displayed over 0.008 < Delta r < 0.39.
void trackClosurePlotterRun3(const char* recoFile,const char* genFile,
  const char* outputDirectory="figures/OO_MC_firstReview/tracking",double centralityLow=0.,
  double centralityHigh=30.,double trackThreshold=2.,const char* workingPoint="Nominal",
  bool drawTracks=true,bool drawEEC=true,double angularTrackHigh=3.) {
  using namespace Run3TrackingPlot;
  try {
    Require(drawTracks || drawEEC,"Select at least one plot family");
    Require(centralityLow<centralityHigh && trackThreshold>0.,"Invalid requested selection");
    const TString wp(workingPoint);
    Require(wp=="Nominal" || wp=="Loose" || wp=="Tight","Unknown working-point label");
    gROOT->SetBatch(true);gStyle->SetOptStat(0);gStyle->SetEndErrorSize(2);
    std::unique_ptr<TFile> rf(TFile::Open(recoFile,"READ")),gf(TFile::Open(genFile,"READ"));
    for(auto* f:{rf.get(),gf.get()})Require(f && !f->IsZombie() && !f->TestBit(TFile::kRecovered),"Input file absent, recovered or unreadable");
    CheckInputs(*rf,*gf);
    const double requiredVertex=wp=="Loose"?5.:(wp=="Tight"?2.:3.);
    Require(Near(Scalar(*rf,"VertexMaxDistance"),requiredVertex),"Working-point label disagrees with VertexMaxDistance");
    auto rm=std::make_unique<EECHistogramManager>(rf.get());
    auto gm=std::make_unique<EECHistogramManager>(gf.get());
    auto* rc=rm->GetCard();auto* gc=gm->GetCard();
    Require(rc->GetNWeightExponents()==1 && gc->GetNWeightExponents()==1,"Use separately projected n1/n2 files");
    const int power=rc->GetWeightExponent();Require(power==1 || power==2,"Expected n=1 or n=2");
    const int cr=rc->FindBinIndexCentrality(centralityLow,centralityHigh),cg=gc->FindBinIndexCentrality(centralityLow,centralityHigh);
    Require(cr>=0 && cg>=0,"Requested exact centrality interval is absent");
    const double nr=EventSum(*rm,cr),ng=EventSum(*gm,cg);
    Require(Near(nr,ng,1.e-4),"RecoReco/RecoGen event exposures differ; check campaign inputs before raw-scale closure");
    CheckJets(*rm,*gm,cr,cg);
    Require(gSystem->mkdir(outputDirectory,true)==0 || !gSystem->AccessPathName(outputDirectory),"Cannot create output directory");
    const TString tag=Form("%s_C%g-%g_n%d",workingPoint,centralityLow,centralityHigh,power);
    const TString selection=Form("%s; Cent: %g-%g%%; R = 0.4",workingPoint,centralityLow,centralityHigh);
    const TString provenance=Form("reco=%s; gen=%s; sumw_events_reco=%.17g; sumw_events_gen=%.17g; n=%d; common jet-spectrum check passed; event-ID equality not established by histograms",
                                  recoFile,genFile,nr,ng,power);
    std::ofstream log((TString(outputDirectory)+"/"+tag+"_input_checks.txt").Data());
    Require(bool(log),"Cannot write input check report");
    log<<provenance<<"\nTrack spectra: per-event; normalization uncertainty not propagated.\nEEC: unnormalized common exposure; no subtraction or unfolding; residual pair correction OFF.\n";
    std::cout<<"INPUT_CHECKS_PASSED events="<<nr<<" centrality="<<centralityLow<<'-'<<centralityHigh<<" n="<<power<<'\n';
    std::cout<<"NOTE: aggregate event/jet checks cannot prove identical physical event IDs; verify merge receipts.\n";
    int plots=0;
    if(drawTracks) {
      int tr=-1,tg=-1;
      for(int t=0;t<rc->GetNTrackPtBins();++t)
        if(Near(rc->GetLowBinBorderTrackPt(t),trackThreshold) && Near(rc->GetHighBinBorderTrackPt(t),angularTrackHigh))tr=t;
      for(int t=0;t<gc->GetNTrackPtBins();++t)
        if(Near(gc->GetLowBinBorderTrackPt(t),trackThreshold) && Near(gc->GetHighBinBorderTrackPt(t),angularTrackHigh))tg=t;
      Require(tr>=0 && tg>=0,"Exact angular track-pT slice unavailable");
      for(int v=0;v<3;++v) {
        auto get=[&](EECHistogramManager& m,int type,int c,int t)->TH1D*{
          if(v==0)return m.GetHistogramTrackPt(type,c);
          if(v==1)return m.GetHistogramTrackEta(type,c,t);
          return m.GetHistogramTrackPhi(type,c,t);
        };
        auto gen=Copy(get(*gm,0,cg,tg),"gen_tracks");auto reco=Copy(get(*rm,0,cr,tr),"corrected_tracks");
        auto raw=Copy(get(*rm,1,cr,tr),"uncorrected_tracks");
        for(auto* h:{gen.get(),reco.get(),raw.get()})RebinDensity(*h,v==0?10:2);
        gen->Scale(1./ng);reco->Scale(1./nr);raw->Scale(1./nr);
        const char* names[]={"trackPt","trackEta","trackPhi"};
        const char* x[]={"p_{T}^{ch} (GeV)","#eta","#varphi"};
        const char* y[]={"(1/#Sigmaw_{evt}) dN/dp_{T}","(1/#Sigmaw_{evt}) dN/d#eta","(1/#Sigmaw_{evt}) dN/d#varphi"};
        const double lo=v==0?trackThreshold:gen->GetXaxis()->GetXmin(),hi=v==0?100.:gen->GetXaxis()->GetXmax();
        TString detail=v==0?"Inclusive tracks; per-event normalization":Form("Inclusive tracks; %g < p_{T}^{ch} < %g GeV",trackThreshold,angularTrackHigh);
        Draw(*gen,*reco,raw.get(),TString(outputDirectory)+"/"+tag+"_"+names[v],x[v],y[v],v==0,v==0,lo,hi,selection,detail,provenance);
        ++plots;
      }
    }
    if(drawEEC) {
      const int tr=rc->GetBinIndexTrackPtEEC(trackThreshold),tg=gc->GetBinIndexTrackPtEEC(trackThreshold);
      Require(tr>=0 && tg>=0 && Near(rc->GetLowBinBorderTrackPtEEC(tr),trackThreshold) &&
              Near(gc->GetLowBinBorderTrackPtEEC(tg),trackThreshold),"Exact EEC track threshold unavailable");
      const double edges[]={40.,80.,120.,160.};
      for(int j=0;j<3;++j) {
        const int jr=rc->FindBinIndexJetPtEEC(edges[j],edges[j+1]),jg=gc->FindBinIndexJetPtEEC(edges[j],edges[j+1]);
        Require(jr>=0 && jg>=0,"An exact reporting jet-pT interval is absent");
        auto reco=Copy(rm->GetHistogramEnergyEnergyCorrelator(0,cr,jr,tr),"reco_eec_total");
        auto gen=Copy(gm->GetHistogramEnergyEnergyCorrelator(0,cg,jg,tg),"gen_eec_total");
        const TString detail=Form("%g < p_{T}^{jet,reco} < %g GeV; p_{T}^{ch} > %g GeV; n=%d",edges[j],edges[j+1],trackThreshold,power);
        Draw(*gen,*reco,nullptr,TString(outputDirectory)+"/"+tag+Form("_EEC_J%g-%g_T%g",edges[j],edges[j+1],trackThreshold),
             "#Deltar","EEC (common raw scale)",true,true,0.008,0.39,selection,detail,provenance+"; total same-jet pairs; residual pair correction OFF");
        ++plots;
      }
    }
    std::cout<<"RUN3_TRACKING_PLOTTING_PASS plots="<<plots<<" tag="<<tag<<" (software completion, not a physics closure verdict)\n";
  } catch(const std::exception& e) {
    std::cerr<<"RUN3_TRACKING_PLOTTING_FAIL "<<e.what()<<'\n';gSystem->Exit(1);
  }
}

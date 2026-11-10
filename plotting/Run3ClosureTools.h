#ifndef RUN3_CLOSURE_TOOLS_H
#define RUN3_CLOSURE_TOOLS_H
// Helpers for COPIES of the Run-2 closure macros; Run-2 files are not edited.
// Shared-signal errors are a pair-Sumw2 approximation, NOT full event covariance.
#include "EECCard.h"
#include "EECHistogramManager.h"
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Run3Closure {
inline void Require(bool ok, const TString& message) {
  if(!ok) throw std::runtime_error(message.Data());
}
inline bool Close(double a, double b, double tolerance=1.e-8) {
  return std::isfinite(a) && std::isfinite(b) &&
    std::abs(a-b)<=tolerance*std::max({1.,std::abs(a),std::abs(b)});
}
inline TFile* Open(const TString& name) {
  TFile* f=TFile::Open(name,"READ");
  Require(f && !f->IsZombie() && !f->TestBit(TFile::kRecovered),"Cannot open intact input: "+name);
  for(const char* key : {"DataType","McCorrelationType","JetRadius","WeightExponent","CentralityBinEdges",
                         "JetPtBinEdgesEEC","TrackPtBinEdgesEEC","BackgroundMethods"})
    Require(f->Get(TString("JCard/")+key),TString("Missing JCard/")+key);
  return f;
}
inline void CheckCard(EECCard& card) {
  Require(card.GetDataType().Contains("OO") && card.GetDataType().Contains("MC"),"Run3 closure copies require OO MC");
  const auto* v=card.GetCardEntry(EECCard::kMcCorrelationType);
  Require(v && v->GetNrows()==1,"Missing MC correlation mode");
  const int mode=std::lround((*v)[v->GetLwb()]);
  Require(mode==1 || mode==3,"Truth-origin closure requires RecoGen or GenGen with stored subevents");
  Require(card.GetNWeightExponents()==1,"Project each energy exponent into a separate input first");
  Require(card.GetWeightExponent()==1 || card.GetWeightExponent()==2,"Only n=1 or n=2 is supported");
  Require(std::abs(card.GetJetRadius()-0.4)<1.e-5,"This Run3 closure copy is configured for R=0.4");
}
inline TString Mode(EECCard& card) {
  const auto* v=card.GetCardEntry(EECCard::kMcCorrelationType);
  return std::lround((*v)[v->GetLwb()])==3 ? "GenGen" : "RecoGen";
}
inline int Track(EECCard& card,double threshold) {
  const int i=card.GetBinIndexTrackPtEEC(threshold);
  Require(i>=0 && i<card.GetNTrackPtBinsEEC() &&
          std::abs(card.GetLowBinBorderTrackPtEEC(i)-threshold)<1.e-5,
          Form("Missing exact track threshold %.6g GeV",threshold));
  return i;
}
inline void Compatible(const TH1D* a,const TH1D* b) {
  Require(a && b,"Required histogram is absent");
  Require(a->GetNbinsX()==b->GetNbinsX(),"Histogram bin counts differ");
  for(int k=1;k<=a->GetNbinsX()+1;++k)
    Require(Close(a->GetXaxis()->GetBinLowEdge(k),b->GetXaxis()->GetBinLowEdge(k),1.e-12),"Histogram bin edges differ");
}
inline TH1D* Clone(const TH1D* source,const TString& name) {
  Require(source,"Cannot clone a missing histogram");
  auto* h=dynamic_cast<TH1D*>(source->Clone(name));
  Require(h,"Expected TH1D"); h->SetDirectory(nullptr); h->GetXaxis()->SetRange(0,0);
  for(int k=0;k<h->GetNcells();++k)
    Require(std::isfinite(h->GetBinContent(k)) && std::isfinite(h->GetBinError(k)),"Nonfinite histogram value/error");
  return h;
}
inline TH1D* Read(EECHistogramManager& manager,int c,int j,int t,int pairing,int subevent,const TString& name) {
  auto* h=manager.GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator,c,j,t,pairing,subevent);
  Require(h,Form("Missing histogram C%d J%d T%d pairing%d subevent%d",c,j,t,pairing,subevent));
  return Clone(h,name);
}
inline TH1D* TrueBackground(EECHistogramManager& manager,int c,int j,int t) {
  auto* h=Read(manager,c,j,t,EECHistograms::kSameJetPair,EECHistograms::kPythiaHydjet,"run3_true_background");
  auto* hh=manager.GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator,c,j,t,
      EECHistograms::kSameJetPair,EECHistograms::kHydjetHydjet);
  Compatible(h,hh); h->Add(hh); return h; // Projected S1 already includes S2.
}
inline TH1D* Background(EECHistogramManager& manager,int c,int j,int t,int method,bool averageSameCone) {
  Require(method==0 || method==1,"method must be 0 (mixed) or 1 (perpendicular)");
  const int s=5*method, all=EECHistograms::knSubeventCombinations;
  auto get=[&](int p){return manager.GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator,c,j,t,p+s,all);};
  auto* b=Clone(get(EECHistograms::kSignalMixedConePair),"run3_background");
  auto add=[&](int p,double scale){auto* h=get(p);Compatible(b,h);b->Add(h,scale);};
  add(EECHistograms::kSignalSecondMixedConePair,1.); b->Scale(0.5);
  add(EECHistograms::kMixedMixedConePair,-1.);
  add(EECHistograms::kMixedConePair,averageSameCone?0.5:1.);
  if(averageSameCone) add(EECHistograms::kSecondMixedConePair,0.5);
  return b;
}
inline void Preflight(EECHistogramManager& manager,int c,int j,int t,int method=-1) {
  auto* card=manager.GetCard();
  Require(c>=0 && c<card->GetNCentralityBins() && j>=0 && j<card->GetNJetPtBinsEEC() && t>=0 && t<card->GetNTrackPtBinsEEC(),"Requested bin unavailable");
  std::unique_ptr<TH1D> total(Read(manager,c,j,t,EECHistograms::kSameJetPair,EECHistograms::knSubeventCombinations,"run3_check_total"));
  std::unique_ptr<TH1D> signal(Read(manager,c,j,t,EECHistograms::kSameJetPair,EECHistograms::kPythiaPythia,"run3_check_signal"));
  std::unique_ptr<TH1D> truth(TrueBackground(manager,c,j,t));
  Compatible(total.get(),signal.get()); Compatible(total.get(),truth.get());
  for(int k=0;k<total->GetNcells();++k)
    Require(Close(total->GetBinContent(k),signal->GetBinContent(k)+truth->GetBinContent(k)),"T=S0+S1+S3 identity fails");
  for(int m=0;m<2;++m) if(method<0 || method==m) {
    const int bit=m==0?4:2;
    Require((card->GetBackgroundMethods()&bit)!=0,"Requested background method not enabled in JCard");
    std::unique_ptr<TH1D> b(Background(manager,c,j,t,m,true)); Compatible(total.get(),b.get());
  }
}
inline void PositiveDenominator(const TH1D* h,double low=0.008,double high=0.39) {
  Require(h,"Missing ratio denominator");
  for(int k=1;k<=h->GetNbinsX();++k) {
    if(h->GetXaxis()->GetBinUpEdge(k)<=low || h->GetXaxis()->GetBinLowEdge(k)>=high) continue;
    Require(h->GetBinContent(k)>0. && std::isfinite(h->GetBinContent(k)),
      Form("No positive truth support in displayed bin %d; no closure plot is produced",k));
  }
}
// Under T=S+Btrue and independent S,Btrue,Best except for identical S in X:
// Var(R)=[Var(Btrue)+Var(Best)+(R-1)^2 Var(S)]/S^2.
// Uses the same stored marginal variances as the independent baseline.
inline void ApplySharedSignal(TH1D* ratio,const TH1D* extracted,const TH1D* signal,
                              const TH1D* total,const TH1D* trueBackground) {
  Compatible(ratio,extracted);Compatible(ratio,signal);Compatible(ratio,total);Compatible(ratio,trueBackground);
  for(int k=1;k<=ratio->GetNbinsX();++k) {
    const double den=signal->GetBinContent(k); if(!(den>0.)) continue;
    const double vs=std::pow(signal->GetBinError(k),2), vt=std::pow(total->GetBinError(k),2);
    const double vb=std::pow(trueBackground->GetBinError(k),2), vx=std::pow(extracted->GetBinError(k),2);
    Require(Close(vt,vs+vb),"Stored T variance is not Var(S)+Var(Btrue); approximation not applicable");
    double ve=vx-vt;
    Require(ve>=-1.e-10*std::max(1.,vx),"Cannot infer nonnegative estimated-background variance");
    ve=std::max(0.,ve);
    const double r=ratio->GetBinContent(k);
    Require(Close(r,extracted->GetBinContent(k)/den),"Shared-signal mode requires raw common-scale ratios");
    const double variance=(vb+ve+(r-1.)*(r-1.)*vs)/(den*den);
    Require(std::isfinite(variance) && variance>=0.,"Invalid shared-signal variance");
    ratio->SetBinError(k,std::sqrt(variance));
  }
}
inline void Header(TCanvas* canvas,bool correlated,bool averageSameCone,const TString& detail) {
  canvas->cd();
  TLatex title;title.SetNDC();title.SetTextFont(62);title.SetTextSize(0.032);
  title.DrawLatex(0.14,0.970,"CMS"); title.SetTextFont(52);title.SetTextSize(0.025);
  title.DrawLatex(0.24,0.970,"Simulation Preliminary");
  title.SetTextFont(42);title.SetTextSize(0.020);
  title.DrawLatex(0.14,0.938,correlated?"Ratio: shared-signal approximation (not full covariance)":"Ratio: legacy independent-histogram uncertainties");
  canvas->SetTitle(detail+(averageSameCone?"; averaged same-cone term":"; legacy single same-cone term"));
}
inline void Save(TCanvas* canvas,const TString& stem,const std::vector<TH1D*>& curves,
                 const std::vector<TH1D*>& ratios,const std::vector<TString>& labels,
                 const TString& input,bool correlated,bool averageSameCone,bool normalized) {
  Require(curves.size()==ratios.size() && curves.size()==labels.size(),"Export size mismatch");
  canvas->Modified();canvas->Update();canvas->SaveAs(stem+".pdf");canvas->SaveAs(stem+".png");
  TFile output(stem+".root","RECREATE");Require(!output.IsZombie(),"Cannot create numerical output");
  canvas->Write("canvas");
  TNamed("input_file",input.Data()).Write();
  TNamed("ratio_uncertainty",correlated?"SHARED_SIGNAL_APPROXIMATION_NOT_FULL_COVARIANCE":"LEGACY_UNCORRELATED_TH1_DIVIDE").Write();
  TNamed("background_definition",averageSameCone?"0.5*(JC1+JC2+C1C1+C2C2)-C1C2":"0.5*(JC1+JC2)-C1C2+C1C1 (legacy macro arithmetic)").Write();
  TNamed("normalization",normalized?"Independent unit-area shapes (legacy convention)":"Common raw input scale; no independent normalization").Write();
  TNamed("limitations","No measured event/shared-track/donor-reuse covariance or systematic uncertainties").Write();
  std::ofstream csv((stem+".csv").Data());Require(bool(csv),"Cannot create numerical CSV");
  csv<<"curve,bin,low,high,content,error,ratio,ratio_error\n"<<std::setprecision(17);
  for(size_t i=0;i<curves.size();++i) {
    curves[i]->Write(Form("curve%zu",i));ratios[i]->Write(Form("ratio%zu",i));TNamed(Form("label%zu",i),labels[i].Data()).Write();
    for(int b=1;b<=curves[i]->GetNbinsX();++b)
      csv<<i<<','<<b<<','<<curves[i]->GetXaxis()->GetBinLowEdge(b)<<','<<curves[i]->GetXaxis()->GetBinUpEdge(b)<<','
         <<curves[i]->GetBinContent(b)<<','<<curves[i]->GetBinError(b)<<','<<ratios[i]->GetBinContent(b)<<','<<ratios[i]->GetBinError(b)<<'\n';
  }
  output.Close();
}
} // namespace Run3Closure
#endif

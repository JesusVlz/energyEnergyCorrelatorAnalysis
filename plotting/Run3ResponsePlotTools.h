#ifndef RUN3_RESPONSE_PLOT_TOOLS_H
#define RUN3_RESPONSE_PLOT_TOOLS_H
// Additive Run-3 adaptations of plotJetPtResponseMatrix.C and
// trackParticleMatchingPlotter.C. Neither legacy file is changed.
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TNamed.h>
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
namespace Run3ResponsePlot {
inline void Require(bool ok,const TString& why){if(!ok)throw std::runtime_error(why.Data());}
inline bool Near(double a,double b){return std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<1.e-6*std::max({1.,std::abs(a),std::abs(b)});}
inline TString Unique(){static unsigned i=0;return Form("run3resp_%u",++i);}
template<class H> std::unique_ptr<H> Copy(const H* h){
  Require(h,"Required projected histogram is missing; inspect projection log and JCard");
  auto p=std::unique_ptr<H>(dynamic_cast<H*>(h->Clone(Unique())));
  Require(bool(p),"Unexpected histogram class");p->SetDirectory(nullptr);
  p->GetXaxis()->SetRange(0,0);if(p->GetDimension()>1)p->GetYaxis()->SetRange(0,0);
  for(int b=0;b<p->GetNcells();++b)Require(std::isfinite(p->GetBinContent(b))&&std::isfinite(p->GetBinError(b)),"Nonfinite input bin");
  return p;
}
inline double Scalar(TFile& f,const char* key){
  auto* v=dynamic_cast<TVectorT<float>*>(f.Get(TString("JCard/")+key));
  Require(v&&v->GetNrows()==1,TString("Missing scalar JCard/")+key);
  const double x=(*v)[v->GetLwb()];Require(std::isfinite(x),TString("Nonfinite JCard/")+key);return x;
}
struct Input {
  std::unique_ptr<TFile> file;std::unique_ptr<EECHistogramManager> manager;
  EECCard* card=nullptr;int centrality=0,track=0,power=1;TString selection,provenance,out;
  Input(const char* name,const char* output,double low,double high,double threshold,int mode):out(output){
    file.reset(TFile::Open(name,"READ"));Require(file&&!file->IsZombie()&&!file->TestBit(TFile::kRecovered),"Cannot open intact projected input");
    Require(Scalar(*file,"DataType")==10 && Scalar(*file,"McCorrelationType")==mode,"Wrong OO MC reconstruction mode for this plot family");
    for(const char* k:{"CentralityBinEdges","JetPtBinEdgesEEC","TrackPtBinEdgesEEC","WeightExponent","JetRadius"})
      Require(file->Get(TString("JCard/")+k),TString("Missing JCard/")+k);
    manager=std::make_unique<EECHistogramManager>(file.get());card=manager->GetCard();
    Require(Near(card->GetJetRadius(),.4),"This plotter requires R=0.4");
    Require(card->GetNWeightExponents()==1,"Use one separately projected energy power per file");
    power=card->GetWeightExponent();Require(power==1||power==2,"Expected n=1 or n=2");
    centrality=card->FindBinIndexCentrality(low,high);Require(low<high&&centrality>=0,"Exact activity interval missing (integrated is stored as -1,100)");
    track=card->GetBinIndexTrackPtEEC(threshold);
    Require(track>=0&&track<card->GetNTrackPtBinsEEC()&&Near(card->GetLowBinBorderTrackPtEEC(track),threshold),"Exact track threshold missing");
    Require(!gSystem->AccessPathName(output)||gSystem->mkdir(output,true)==0,"Cannot create output directory");
    selection=Near(low,-1)&&Near(high,100)?"OO integrated selected events":Form("OO centrality %g-%g%%",low,high);
    selection+=Form("; R=0.4; projected n=%d",power);
    provenance=Form("input=%s; centrality=[%g,%g); track=%g; n=%d; read-only; no unfolding, calibration fit or pair-map derivation",name,low,high,threshold,power);
  }
};
inline std::pair<int,int> Exact(const TAxis& a,double low,double high){
  const int l=a.FindFixBin(std::nextafter(low,high)),h=a.FindFixBin(std::nextafter(high,low));
  if(l<1||h>a.GetNbins()||l>h||!Near(a.GetBinLowEdge(l),low)||!Near(a.GetBinUpEdge(h),high))return {-1,-1};
  return {l,h};
}
inline void SameAxis(const TAxis& a,const TAxis& b){
  Require(a.GetNbins()==b.GetNbins(),"Axis dimensions differ");
  for(int i=1;i<=a.GetNbins()+1;++i)Require(Near(a.GetBinLowEdge(i),b.GetBinLowEdge(i)),"Axis boundaries differ");
}
// P(y-bin | x-bin, matched, both axes in regular range), as in the legacy
// particle-response plotter. Unsupported slices are blank and explicitly masked.
inline std::unique_ptr<TH2D> ConditionalX(const TH2D& raw,TH2D& valid){
  auto h=Copy(&raw);h->Reset();valid.Reset();
  for(int x=1;x<=raw.GetNbinsX();++x){double sum=0.;
    for(int y=1;y<=raw.GetNbinsY();++y){Require(raw.GetBinContent(x,y)>=0.,"Negative cell: conditional probability not defined");sum+=raw.GetBinContent(x,y);}
    if(!(sum>0.))continue;
    for(int y=1;y<=raw.GetNbinsY();++y){h->SetBinContent(x,y,raw.GetBinContent(x,y)/sum);h->SetBinError(x,y,raw.GetBinError(x,y)/sum);valid.SetBinContent(x,y,1.);}
  }return h;
}
inline std::unique_ptr<TH2D> IndexView(const TH2D& raw){
  auto h=std::make_unique<TH2D>(Unique(),"",raw.GetNbinsX(),.5,raw.GetNbinsX()+.5,raw.GetNbinsY(),.5,raw.GetNbinsY()+.5);h->SetDirectory(nullptr);h->Sumw2();
  for(int i=0;i<raw.GetNcells();++i){h->SetBinContent(i,raw.GetBinContent(i));h->SetBinError(i,raw.GetBinError(i));}return h;
}
inline void PositiveRange(TAxis& a){
  if(a.GetXmin()>0.)return;
  int first=1;while(first<=a.GetNbins()&&a.GetBinLowEdge(first)<=0.)++first;
  Require(first<=a.GetNbins(),"No positive axis support for logarithmic display");a.SetRange(first,a.GetNbins());
  std::cout<<"LOG_DISPLAY_FIRST_BIN "<<first<<"; excluded nonpositive-edge bins retained in ROOT/CSV\n";
}
inline JDrawer& Drawer(bool twoD,bool logX,bool logY,bool logZ){
  // Process-lifetime owner, as in the legacy macros. Delete canvas only.
  static auto* d=new JDrawer();d->Reset();d->SetCanvasSize(850,780);
  d->SetLeftMargin(.18);d->SetRightMargin(twoD?.19:.06);d->SetTopMargin(.20);d->SetBottomMargin(.17);
  d->SetTitleSizeX(.037);d->SetTitleSizeY(.037);d->SetTitleSizeZ(.031);
  d->SetLabelSizeX(.033);d->SetLabelSizeY(.033);d->SetLabelSizeZ(.030);d->SetTitleOffsetZ(1.35);
  d->SetTitleOffsetY(1.5);d->SetTitleOffsetX(1.2);d->SetLogX(logX);d->SetLogY(logY);d->SetLogZ(logZ);return *d;
}
inline void Header(const TString& selection,const TString& detail){
  TLatex t;t.SetNDC();t.SetTextFont(62);t.SetTextSize(.032);t.DrawLatex(.15,.963,"CMS");
  t.SetTextFont(52);t.SetTextSize(.025);t.DrawLatex(.25,.963,"Simulation Preliminary");
  t.SetTextFont(42);t.SetTextSize(.022);t.SetTextAlign(31);t.DrawLatex(.96,.963,"Pythia+Hijing, 5.36 TeV");t.SetTextAlign(11);
  t.SetTextSize(.024);t.DrawLatex(.15,.918,selection);t.SetTextSize(.020);t.DrawLatex(.15,.876,detail);
}
inline void Save(const TString& stem,TCanvas* canvas,const Input& in,const TString& description,
                 const std::vector<std::pair<TObject*,TString>>& objects){
  for(const char* ext:{".root",".png",".pdf",".csv"})Require(gSystem->AccessPathName(stem+ext),"Output exists; use a fresh directory: "+stem);
  TFile out(stem+".root","CREATE");Require(!out.IsZombie(),"Cannot create ROOT output");
  for(const auto& o:objects)o.first->Write(o.second);
  TNamed("provenance",in.provenance.Data()).Write();TNamed("definition",description.Data()).Write();
  in.card->Write(&out);out.cd();canvas->Write("canvas");out.Close();
  canvas->Modified();canvas->Update();canvas->Print(stem+".pdf");canvas->Print(stem+".png");
  std::cout<<"RUN3_RESPONSE_PLOT_SAVED "<<stem<<'\n';
}
inline void Matrix(const TH2D* source,Input& in,const TString& name,const TString& xt,const TString& yt,
                   const TString& detail,bool conditional=false,bool logAxes=false,bool indices=false){
  auto raw=Copy(source);Require(raw->Integral(1,raw->GetNbinsX(),1,raw->GetNbinsY())>0.,"Empty response matrix: "+name);
  auto mask=Copy(raw.get());mask->Reset();
  auto h=conditional?ConditionalX(*raw,*mask):Copy(raw.get());
  if(indices){h=IndexView(*h);mask=IndexView(*mask);}
  bool positive=true;for(int k=0;k<h->GetNcells();++k)if(h->GetBinContent(k)<0.)positive=false;
  h->SetTitle("");h->GetZaxis()->SetTitle(conditional?"Conditional fraction":"Stored weighted sum / bin");
  if(conditional){double zmin=1.;for(int k=0;k<h->GetNcells();++k)if(h->GetBinContent(k)>0.)zmin=std::min(zmin,h->GetBinContent(k));h->SetMinimum(zmin*.5);h->SetMaximum(1.);}
  if(logAxes){PositiveRange(*h->GetXaxis());PositiveRange(*h->GetYaxis());}
  auto& d=Drawer(true,logAxes,logAxes,positive);d.DrawHistogram(h.get(),xt,yt,"","COLZ");Header(in.selection,detail);
  if(conditional){TLatex t;t.SetNDC();t.SetTextSize(.018);t.DrawLatex(.15,.833,"Each reco slice normalized over gen bins; unsupported slices blank");}
  auto* canvas=gPad->GetCanvas();const TString stem=in.out+"/"+name;
  const TString definition=detail+(conditional?"; fixed-X regular-Y normalization; normalizer fixed in errors; not full covariance; unsupported slices blank":"; raw bin sums, no normalization");
  Save(stem,canvas,in,definition,{{raw.get(),"raw"},{h.get(),"display"},{mask.get(),"conditional_support"}});
  std::ofstream csv((stem+".csv").Data());Require(bool(csv),"Cannot write matrix CSV");csv<<std::setprecision(17)<<"ix,iy,xlow,xhigh,ylow,yhigh,raw,raw_error,display,display_error,supported,is_flow\n";
  const double nan=std::numeric_limits<double>::quiet_NaN();double flow=0.;
  for(int x=0;x<=raw->GetNbinsX()+1;++x)for(int y=0;y<=raw->GetNbinsY()+1;++y){
    const bool edge=x==0||y==0||x==raw->GetNbinsX()+1||y==raw->GetNbinsY()+1;
    const bool ok=!conditional||mask->GetBinContent(x,y)>0.;if(edge)flow+=raw->GetBinContent(x,y);
    csv<<x<<','<<y<<','<<raw->GetXaxis()->GetBinLowEdge(x)<<','<<raw->GetXaxis()->GetBinUpEdge(x)<<','<<raw->GetYaxis()->GetBinLowEdge(y)<<','<<raw->GetYaxis()->GetBinUpEdge(y)<<','<<raw->GetBinContent(x,y)<<','<<raw->GetBinError(x,y)<<','<<(ok?h->GetBinContent(x,y):nan)<<','<<(ok?h->GetBinError(x,y):nan)<<','<<ok<<','<<edge<<'\n';
  }
  std::cout<<"MATRIX_FLOW "<<name<<" sumw="<<flow<<" (flow retained, excluded from conditional normalization)\n";delete canvas;
}
inline void Spectrum(const TH1D* source,Input& in,const TString& name,const TString& xt,const TString& detail,
                     bool probability=false,bool logX=false){
  auto raw=Copy(source),h=Copy(source);const double sum=h->Integral(1,h->GetNbinsX());Require(sum>0.,"Empty spectrum: "+name);
  if(probability)for(int b=1;b<=h->GetNbinsX();++b)Require(h->GetBinContent(b)>=0.,"Signed spectrum cannot be labelled probability");
  if(probability)h->Scale(1./sum);h->SetTitle("");h->SetLineColor(kBlack);h->SetMarkerStyle(20);h->SetMarkerSize(.6);
  bool logY=true;double ymax=0.,ymin=0.;
  for(int b=1;b<=h->GetNbinsX();++b){if(h->GetBinContent(b)<0.)logY=false;ymax=std::max(ymax,h->GetBinContent(b)+h->GetBinError(b));ymin=std::min(ymin,h->GetBinContent(b)-h->GetBinError(b));}
  h->SetMaximum(ymax*2.);if(!logY)h->SetMinimum(ymin*1.2);if(logX)PositiveRange(*h->GetXaxis());auto& d=Drawer(false,logX,logY,false);
  d.DrawHistogram(h.get(),xt,probability?"Probability / stored bin":"Stored weighted sum / bin","","E1");Header(in.selection,detail);
  auto* canvas=gPad->GetCanvas();const TString stem=in.out+"/"+name;
  Save(stem,canvas,in,detail+(probability?"; regular-bin unit sum; normalizer fixed in errors":"; no area normalization"),{{raw.get(),"raw"},{h.get(),"display"}});
  std::ofstream csv((stem+".csv").Data());Require(bool(csv),"Cannot write spectrum CSV");csv<<std::setprecision(17)<<"bin,low,high,raw,raw_error,display,display_error,is_flow\n";
  for(int b=0;b<h->GetNcells();++b)csv<<b<<','<<h->GetBinLowEdge(b)<<','<<h->GetXaxis()->GetBinUpEdge(b)<<','<<raw->GetBinContent(b)<<','<<raw->GetBinError(b)<<','<<h->GetBinContent(b)<<','<<h->GetBinError(b)<<','<<(b==0||b==h->GetNbinsX()+1)<<'\n';
  std::cout<<"SPECTRUM_FLOW "<<name<<" sumw="<<raw->GetBinContent(0)+raw->GetBinContent(raw->GetNbinsX()+1)<<'\n';delete canvas;
}
inline int Slices(TH2D& response,Input& in,const TString& prefix){
  int n=0;const double edges[]={40,80,120,160};
  for(int j=0;j<3;++j)for(int fixGen=0;fixGen<2;++fixGen){
    auto bins=Exact(*(fixGen?response.GetYaxis():response.GetXaxis()),edges[j],edges[j+1]);
    if(bins.first<0){std::cout<<"UNSUPPORTED_NATIVE_INTERVAL "<<prefix<<' '<<edges[j]<<'-'<<edges[j+1]<<"; not relabelled or interpolated\n";continue;}
    std::unique_ptr<TH1D> h(fixGen?response.ProjectionX(Unique(),bins.first,bins.second,"e"):response.ProjectionY(Unique(),bins.first,bins.second,"e"));h->SetDirectory(nullptr);
    if(h->Integral(1,h->GetNbinsX())<=0.){std::cout<<"EMPTY_MATCHED_SLICE "<<edges[j]<<'-'<<edges[j+1]<<'\n';continue;}
    Spectrum(h.get(),in,prefix+Form("_%s%g-%g",fixGen?"fixedGen":"fixedReco",edges[j],edges[j+1]),fixGen?"Reconstructed jet p_{T} (GeV)":"Generated jet p_{T} (GeV)",Form("Matched jets; fixed %s p_{T}: %g-%g GeV; regular-bin unit sum",fixGen?"gen":"reco",edges[j],edges[j+1]),true,true);++n;
  }return n;
}
inline void PairPlots(const char* input,const char* out,double low,double high,double threshold){
  Input in(input,out,low,high,threshold,0);int n=0;const double edges[]={40,80,120,160};
  for(int j=0;j<3;++j){const int bin=in.card->FindBinIndexJetPtEEC(edges[j],edges[j+1]);Require(bin>=0,"Exact reporting jet bin absent");
    const TString tag=Form("J%g-%g_T%g",edges[j],edges[j+1],threshold);
    const TString detail=Form("Matched pairs; reco jets %g-%g GeV; min reco track p_{T}>%g GeV",edges[j],edges[j+1],threshold);
    Matrix(in.manager->GetHistogramTrackParticleDeltaRResponse(in.centrality,bin,in.track),in,"pairAngle_"+tag,"Reconstructed pair #Deltar","Generated pair #Deltar",detail+"; no n dependence",true,true);
    const TString x=Form("(p_{T,1}p_{T,2})^{%d}_{reco} (GeV^{%d})",in.power,2*in.power),y=Form("(p_{T,1}p_{T,2})^{%d}_{gen} (GeV^{%d})",in.power,2*in.power);
    Matrix(in.manager->GetHistogramTrackParticlePtResponse(in.centrality,bin,in.track),in,"pairMomentum_"+tag,x,y,detail,true,true);
    Spectrum(in.manager->GetHistogramTrackParticlePtClosure(in.centrality,bin,in.track),in,"pairMomentumRatio_"+tag,Form("[(p_{T,1}p_{T,2})_{reco}/(p_{T,1}p_{T,2})_{gen}]^{%d}",in.power),detail,true);n+=3;
  }std::cout<<"RUN3_RESPONSE_PLOTTING_PASS family=pairs plots="<<n<<" (not a pair-efficiency closure verdict)\n";
}
inline void JetPlots(const char* input,const char* out,double low,double high){
  Input in(input,out,low,high,2.,0);auto raw=Copy(in.manager->GetHistogramJetPtResponseMatrix(in.centrality));
  Matrix(raw.get(),in,"legacyJetResponse","Reconstructed jet p_{T} (GeV)","Generated jet p_{T} (GeV)","Diagnostic matched jets; legacy matrix bins retained",false,true);
  const int n=1+Slices(*raw,in,"legacyJetSlice");
  std::cout<<"RUN3_RESPONSE_PLOTTING_PASS family=jets plots="<<n<<"; native 40-80 support is checked, never assumed\n";
}
inline void AuditMarginal(const TH1D& total,const TH2D& response,bool reco,Input& in,const TString& name){
  SameAxis(*total.GetXaxis(),*(reco?response.GetXaxis():response.GetYaxis()));
  std::unique_ptr<TH1D> matched(reco?response.ProjectionX(Unique(),0,response.GetNbinsY()+1,"e"):response.ProjectionY(Unique(),0,response.GetNbinsX()+1,"e"));matched->SetDirectory(nullptr);
  std::ofstream csv((in.out+"/"+name+".csv").Data());Require(bool(csv),"Cannot write marginal audit");csv<<std::setprecision(17)<<"bin,total,matched_including_other_axis_flow,total_minus_matched,is_flow\n";
  for(int b=0;b<total.GetNcells();++b)csv<<b<<','<<total.GetBinContent(b)<<','<<matched->GetBinContent(b)<<','<<total.GetBinContent(b)-matched->GetBinContent(b)<<','<<(b==0||b==total.GetNbinsX()+1)<<'\n';
  // No errors assigned: totals and matched entries overlap.
}
inline void UnfoldingPlots(const char* input,const char* out,double low,double high,double threshold){
  Input in(input,out,low,high,threshold,1);int n=0;
  Require(Near(Scalar(*in.file,"JetPtWeight"),0.),"This marginal audit requires JetPtWeight=0");
  auto jet=Copy(in.manager->GetHistogramJetPtOneDimensionalUnfoldingResponse(in.centrality));
  Matrix(jet.get(),in,"nativeJetResponse","Reconstructed jet p_{T} (GeV)","Generated jet p_{T} (GeV)","Native jet response; full guard-bin support; matched entries",false,true);++n;
  n+=Slices(*jet,in,"nativeJetSlice");
  auto eec=Copy(in.manager->GetHistogramJetPtUnfoldingResponse(in.centrality,in.track));
  Matrix(eec.get(),in,"eecResponse","Reconstructed (jet p_{T}, #Deltar) bin index","Generated (jet p_{T}, #Deltar) bin index","EEC-weighted response; flattened indices, NOT GeV",false,false,true);++n;
  for(int reco=0;reco<2;++reco){
    auto h=Copy(reco?in.manager->GetHistogramJetPtUnfoldingMeasured(in.centrality,in.track):in.manager->GetHistogramJetPtUnfoldingTruth(in.centrality,in.track));
    Spectrum(h.get(),in,reco?"eecMeasured":"eecTruth","Stored flattened (jet p_{T}, #Deltar) coordinate","Response input marginal; includes unmatched contributions; not unfolded");++n;
    AuditMarginal(*h,*eec,reco,in,reco?"eecUnmatchedRecoAudit":"eecUnmatchedTruthAudit");
  }
  auto cov=Copy(in.manager->GetHistogramJetPtUnfoldingCovariance(0,in.centrality,in.track));
  Require(cov->GetNbinsX()==cov->GetNbinsY(),"Stored covariance is not square");
  Matrix(cov.get(),in,"storedCovariance","Reco (jet p_{T}, #Deltar) bin index","Reco (jet p_{T}, #Deltar) bin index","Stored within-jet covariance input; NOT full event covariance",false,false,true);++n;
  std::cout<<"RUN3_RESPONSE_PLOTTING_PASS family=unfoldingInputs plots="<<n<<"; no unfolding, iteration choice, or independent closure performed\n";
}
} // namespace Run3ResponsePlot
#endif

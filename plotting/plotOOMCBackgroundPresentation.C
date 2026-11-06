#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include <TH2.h>
#include "JDrawer.h"
#include "../src/EECHistograms.h"
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TNamed.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Raw EEC closure: no area matching, fitting, or independent normalization.
// S1 in the processed files already includes S2. Never add S2 again.
// B = (JC1 + JC2 + C1C1 + C2C2)/2 - C1C2, for either cone method.
namespace OOMCPresentation {
const double xmin = 0.008, xmax = 0.39;
const char* pairTokens[] = {"SignalConePair", "SignalReflectedConePair", "ReflectedConePair",
  "SignalMixedConePair", "MixedConePair", "SignalSecondMixedConePair", "MixedMixedConePair",
  "SecondMixedConePair", "SignalPerpendicularConePair", "PerpendicularConePair",
  "SignalSecondPerpendicularConePair", "PerpendicularPerpendicularConePair", "SecondPerpendicularConePair"};
struct Selection {
  TString flag, id, input, label;
  std::vector<int> cents;
  int jet = -1, track = -1;
  double jetLow = 0., jetHigh = 0., centLow = 0., centHigh = 0.;
};
struct Data {
  std::vector<std::unique_ptr<TH1D>> owned;
  std::map<std::string,TH1D*> h;
  TH1D* keep(TH1D* hnew, const std::string& key) {
    hnew->SetDirectory(nullptr); h[key] = hnew; owned.emplace_back(hnew); return hnew;
  }
  TH1D* clone(TH1D* source, const std::string& key, double scale=1.) {
    TH1D* result = keep(static_cast<TH1D*>(source->Clone(key.c_str())),key);
    result->Scale(scale); return result;
  }
};
struct Curve { TH1D* h; TString label; int color; int marker; };
struct Ratio { TH1D* numerator; TH1D* denominator; int color; int marker; };
struct PlotResult { bool noGrids=true, logY=false; int masked=0, positiveDenominators=0; };
std::string shellQuotePath(const TString& name);

std::string csv(const TString& value) {
  std::string out="\"";
  for(char c : std::string(value.Data())) { if(c=='\"') out+='\"'; out+=c; }
  return out+'\"';
}
TString centList(const Selection& s) {
  TString text; for(int i : s.cents) { if(!text.IsNull()) text+=';'; text+=Form("%d",i); } return text;
}
bool coreBin(const TH1D* h, int i) {
  return h->GetXaxis()->GetBinUpEdge(i)>xmin && h->GetXaxis()->GetBinLowEdge(i)<xmax;
}
bool populated(const TH1D* h) {
  for(int i=1;i<=h->GetNbinsX();++i)
    if(std::isfinite(h->GetBinContent(i)) && (h->GetBinContent(i)!=0. || h->GetBinError(i)!=0.)) return true;
  return false;
}
TH1D* read(TFile& file, const Selection& s, int pair, int subevent, Data& data, const std::string& key) {
  TH1D* sum=nullptr;
  for(int c : s.cents) {
    TString path=Form("energyEnergyCorrelator/%s/energyEnergyCorrelator%s_C%dT%dJ%d",
      pairTokens[pair],pairTokens[pair],c,s.track,s.jet);
    if(subevent != EECHistograms::knSubeventCombinations) path+=Form("S%d",subevent);
    TH1D* source=dynamic_cast<TH1D*>(file.Get(path));
    if(!source) throw std::runtime_error(std::string("Missing histogram: ")+file.GetName()+":"+path.Data());
    if(!sum) sum=data.clone(source,key);
    else {
      if(sum->GetNbinsX()!=source->GetNbinsX()) throw std::runtime_error("Mismatched input bin counts");
      for(int b=1;b<=sum->GetNbinsX()+1;++b)
        if(std::abs(sum->GetXaxis()->GetBinLowEdge(b)-source->GetXaxis()->GetBinLowEdge(b))>1.e-12)
          throw std::runtime_error("Mismatched input bin edges");
      sum->Add(source);
    }
  }
  return sum;
}
void load(TFile& file, const Selection& s, Data& d) {
  const int all=EECHistograms::knSubeventCombinations;
  read(file,s,EECHistograms::kSameJetPair,all,d,"T");
  read(file,s,EECHistograms::kSameJetPair,EECHistograms::kPythiaPythia,d,"PP");
  read(file,s,EECHistograms::kSameJetPair,EECHistograms::kPythiaHydjet,d,"PH");
  read(file,s,EECHistograms::kSameJetPair,EECHistograms::kHydjetHydjet,d,"HH");
  d.clone(d.h["PH"],"Btrue")->Add(d.h["HH"]);
  const int pairs[2][5]={{EECHistograms::kSignalMixedConePair,EECHistograms::kSignalSecondMixedConePair,EECHistograms::kMixedConePair,EECHistograms::kSecondMixedConePair,EECHistograms::kMixedMixedConePair},
    {EECHistograms::kSignalPerpendicularConePair,EECHistograms::kSignalSecondPerpendicularConePair,EECHistograms::kPerpendicularConePair,EECHistograms::kSecondPerpendicularConePair,EECHistograms::kPerpendicularPerpendicularConePair}};
  for(int m=0;m<2;++m) {
    const std::string method=m==0?"mix":"perp";
    for(int p=0;p<5;++p) read(file,s,pairs[m][p],all,d,method+std::to_string(p));
    TH1D* b=d.clone(d.h[method+"0"],"B"+method,0.5);
    for(int p=1;p<5;++p) b->Add(d.h[method+std::to_string(p)],p==4?-1.:0.5);
    d.clone(d.h["T"],"S"+method)->Add(b,-1.);
  }
}
void style(TH1D* h, int color, int marker) {
  h->SetLineColor(color); h->SetMarkerColor(color); h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.85); h->SetLineWidth(2); h->SetStats(false);
  h->GetXaxis()->SetRangeUser(xmin,xmax);
}
std::unique_ptr<TGraphErrors> graph(const Ratio& r, PlotResult& result) {
  std::unique_ptr<TGraphErrors> g(new TGraphErrors());
  g->SetLineColor(r.color); g->SetMarkerColor(r.color); g->SetMarkerStyle(r.marker); g->SetMarkerSize(0.8);
  for(int b=1;b<=r.denominator->GetNbinsX();++b) {
    if(!coreBin(r.denominator,b)) continue;
    const double den=r.denominator->GetBinContent(b), num=r.numerator->GetBinContent(b);
    if(!(den>0.) || !std::isfinite(den) || !std::isfinite(num)) { ++result.masked; continue; }
    const double en=r.numerator->GetBinError(b), ed=r.denominator->GetBinError(b);
    const double value=num/den, error=std::hypot(en/den,num*ed/(den*den));
    if(!std::isfinite(value) || !std::isfinite(error)) { ++result.masked; continue; }
    const int point=g->GetN(); g->SetPoint(point,r.denominator->GetBinCenter(b),value);
    g->SetPointError(point,0.,error); ++result.positiveDenominators;
  }
  return g;
}
// Every displayed ratio omits a nonpositive denominator.
// The conventional independent-error propagation is descriptive; input covariance
// between overlapping EEC components is not available in these projected files.
PlotResult draw(const Selection& s, const TString& family, const std::vector<Curve>& curves,
  const std::vector<Ratio>& ratios, const TString& ratioTitle, const TString& denominatorName,
  const TString& output) {
  PlotResult result;
  // Match existing macros: JDrawer owns pads as well as their parent canvas,
  // so its destructor cannot safely run after a split canvas is created.
  JDrawer& drawer=*new JDrawer();
  drawer.SetDefaultAppearanceSplitCanvas(); drawer.SetRelativeCanvasSize(1.1,1.1);
  drawer.SetLeftMargin(0.17); drawer.SetTopMargin(0.07);
  drawer.SetTitleOffsetY(1.9); drawer.SetTitleOffsetX(1.0);
  drawer.SetLogX(true); drawer.SetGridX(false); drawer.SetGridY(false);
  drawer.CreateSplitCanvas();
  TCanvas* canvas=drawer.GetUpperPad()->GetCanvas();
  canvas->SetName(s.id+"_"+family);
  double low=0., high=0., positive=1.e100;
  bool negative=false;
  for(const Curve& c : curves) {
    style(c.h,c.color,c.marker);
    for(int b=1;b<=c.h->GetNbinsX();++b) {
      if(!coreBin(c.h,b)) continue;
      const double y=c.h->GetBinContent(b), e=c.h->GetBinError(b);
      if(!std::isfinite(y) || !std::isfinite(e)) continue;
      if(y<0.) negative=true;
      if(y>0.) positive=std::min(positive,y);
      low=std::min(low,y-e); high=std::max(high,y+e);
    }
  }
  result.logY=!negative && positive<1.e99 && high>0.;
  std::unique_ptr<TH1D> frame(static_cast<TH1D*>(curves.front().h->Clone("oomc_upper_frame")));
  frame->SetDirectory(nullptr); frame->Reset("ICES"); frame->GetXaxis()->SetRangeUser(xmin,xmax);
  if(result.logY) frame->GetYaxis()->SetRangeUser(positive/2.,high*((family=="decomposition" || family=="cone_symmetry")?30.:10.));
  else { const double span=high>low?high-low:1.; frame->GetYaxis()->SetRangeUser(low-0.10*span,high+(family=="signal_extraction"?0.85:0.28)*span); }
  drawer.SetLogY(result.logY);
  drawer.SetTitleOffsetY(result.logY?1.9:1.5);
  drawer.DrawHistogramToUpperPad(frame.get(),"#Deltar","EEC (A.U.)"," ","AXIS");
  drawer.GetUpperPad()->GetListOfPrimitives()->Remove(frame.get());
  frame->DrawCopy("AXIS");
  for(const Curve& c : curves) c.h->DrawCopy("E1 SAME");
  if(negative) { TLine* zero=new TLine(xmin,0.,xmax,0.); zero->SetLineStyle(3); zero->SetLineColor(kGray+1); zero->Draw(); }
  TLatex header; header.SetNDC(); header.SetTextFont(62); header.SetTextSize(0.043);
  TString selection=Form("%.0f < jet p_{T} < %.0f GeV, p_{T}^{ch} > 2 GeV, n=1",s.jetLow,s.jetHigh);
  if(s.flag=="split") selection=Form("Cent: %.0f-%.0f%%, ",s.centLow,s.centHigh)+selection;
  header.DrawLatex(0.17,0.955,selection);
  header.SetTextFont(42); header.SetTextSize(0.039); header.DrawLatex(0.17,0.898,s.label);
  double legendLeft=0.48, legendBottom=0.51, legendRight=0.97;
  if(family=="background_closure" || family=="background_method_comparison") { legendLeft=0.19; legendBottom=0.58; legendRight=0.76; }
  if(family.Contains("signed_terms")) { legendLeft=0.19; legendRight=0.76; }
  if(family=="decomposition") legendBottom=0.62;
  if(family=="cone_symmetry") { legendLeft=0.19; legendBottom=0.66; legendRight=0.97; }
  TLegend legend(legendLeft,legendBottom,legendRight,0.87);
  legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextSize(curves.size()>6?0.036:0.041); legend.SetTextFont(62);
  if(curves.size()>6) { legend.SetNColumns(2); legend.SetHeader("Pythia+Hijing OO 5.36 TeV"); }
  else legend.AddEntry((TObject*)nullptr,"Pythia+Hijing OO 5.36 TeV","");
  for(const Curve& c : curves) legend.AddEntry(c.h,c.label,"lp");
  legend.DrawClone();
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  double ratioLow=0.8, ratioHigh=1.2;
  for(const Ratio& r : ratios) {
    graphs.emplace_back(graph(r,result));
    TGraphErrors* g=graphs.back().get();
    for(int i=0;i<g->GetN();++i) {
      double x,y; g->GetPoint(i,x,y);
      ratioLow=std::min(ratioLow,y-g->GetErrorY(i)); ratioHigh=std::max(ratioHigh,y+g->GetErrorY(i));
    }
  }
  const double padding=0.08*(ratioHigh-ratioLow);
  std::unique_ptr<TH1D> ratioFrame(static_cast<TH1D*>(frame->Clone("oomc_ratio_frame")));
  ratioFrame->SetDirectory(nullptr); ratioFrame->GetYaxis()->SetRangeUser(ratioLow-padding,ratioHigh+padding);
  drawer.SetLogY(false); drawer.SetTitleOffsetY(1.9); drawer.SetGridX(false); drawer.SetGridY(false);
  drawer.DrawHistogramToLowerPad(ratioFrame.get(),"#Deltar",ratioTitle," ","AXIS");
  drawer.GetLowerPad()->GetListOfPrimitives()->Remove(ratioFrame.get());
  ratioFrame->DrawCopy("AXIS");
  TLine unity(xmin,1.,xmax,1.); unity.SetLineStyle(2); unity.SetLineColor(kGray+2); unity.DrawClone();
  for(auto& g : graphs) g->DrawClone("P E1 SAME");
  TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.055);
  if(result.positiveDenominators==0) note.DrawLatex(0.20,0.84,"No positive "+denominatorName+" bins; ratios unavailable");
  else if(result.masked>0) note.DrawLatex(0.20,0.86,Form("%d ratio points masked (denominator #leq 0)",result.masked));
  for(TPad* pad : {drawer.GetUpperPad(),drawer.GetLowerPad()}) {
    pad->SetGridx(false); pad->SetGridy(false);
    result.noGrids=result.noGrids && !pad->GetGridx() && !pad->GetGridy(); pad->Modified();
  }
  canvas->SetGridx(false); canvas->SetGridy(false); canvas->Modified(); canvas->Update();
  canvas->SaveAs(output+".pdf");
  // Export the authoritative ROOT PDF to avoid ROOT's raster backend dropping
  // labels on sparse canvases. Poppler only rasterizes the C++ figure.
  const std::string raster="pdftoppm -cropbox -singlefile -r 144 -png "+shellQuotePath(output+".pdf")+" "+shellQuotePath(output);
  if(gSystem->Exec(raster.c_str())!=0 || gSystem->AccessPathName(output+".png",kReadPermission))
    throw std::runtime_error("Cannot export PNG preview: "+std::string(output.Data()));
  delete canvas;
  return result;
}
void row(std::ofstream& out, const Selection& s, const TString& status, const TString& family,
  const TString& pdf, const TString& png, const PlotResult& r=PlotResult()) {
  out<<csv(s.id)<<','<<csv(s.flag)<<','<<csv(s.input)<<','<<csv(centList(s))<<',';
  if(s.flag=="split") out<<s.centLow; out<<',';
  if(s.flag=="split") out<<s.centHigh; out<<',';
  out<<s.jet<<','<<s.jetLow<<','<<s.jetHigh<<','<<s.track<<",2,1,6,"<<csv(status)<<','<<csv(family)<<','
    <<(r.noGrids?"true":"false")<<",0,0,0,0,"<<(r.logY?"true":"false")<<','<<r.masked<<','<<r.positiveDenominators<<','
    <<csv(pdf)<<','<<csv(png)<<'\n';
}
// PDF assembly only: every source page has already been drawn by ROOT/JDrawer.
// Quoting protects spaces, apostrophes, and shell metacharacters in output paths.
std::string shellQuotePath(const TString& name) {
  const TString path=name.BeginsWith("/")?name:TString("./")+name;
  std::string quoted="'";
  for(char c : std::string(path.Data())) {
    if(c=='\'') quoted+="'\\''";
    else quoted+=c;
  }
  return quoted+"'";
}
void packAtlas(const std::vector<TString>& pdfs, const TString& output) {
  if(pdfs.empty()) throw std::runtime_error("No rendered PDF pages for atlas "+std::string(output.Data()));
  std::string command="pdfunite";
  for(const TString& pdf : pdfs) {
    if(gSystem->AccessPathName(pdf,kReadPermission)) throw std::runtime_error("Missing atlas source page: "+std::string(pdf.Data()));
    command+=" "+shellQuotePath(pdf);
  }
  command+=" "+shellQuotePath(output);
  const int status=gSystem->Exec(command.c_str());
  if(status!=0 || gSystem->AccessPathName(output,kReadPermission))
    throw std::runtime_error("pdfunite failed for "+std::string(output.Data())+" (status "+std::to_string(status)+")");
}
void saveHistograms(TFile& output, const Selection& s, Data& d) {
  output.cd(); TDirectory* dir=output.mkdir(s.id); dir->cd();
  for(const char* name : {"T","PP","PH","HH","Btrue","Bmix","Bperp","Smix","Sperp"}) d.h[name]->Write(name);
  TNamed("centrality_indices",centList(s).Data()).Write();
  TNamed("selection_flag",s.flag.Data()).Write(); TNamed("input_file",s.input.Data()).Write();
  output.cd();
}
void render(const Selection& s, Data& d, const TString& directory,
  std::ofstream& manifest, std::vector<std::pair<TString,TString>>& pages) {
  auto plot=[&](const TString& family,const std::vector<Curve>& c,const std::vector<Ratio>& r,
    const TString& title,const TString& denominator) {
    const TString base=directory+"/"+s.id+"_"+family;
    PlotResult status=draw(s,family,c,r,title,denominator,base);
    row(manifest,s,"plotted",family,base+".pdf",base+".png",status);
    pages.push_back({s.flag,base+".pdf"});
  };
  TH1D *t=d.h["T"], *pp=d.h["PP"], *ph=d.h["PH"], *hh=d.h["HH"], *bt=d.h["Btrue"];
  plot("decomposition",{{t,"Total",kBlack,20},{pp,"Pythia-Pythia (signal)",kBlue+1,24},
    {ph,"Pythia-Hijing",kRed+1,25},{hh,"Hijing-Hijing",kGreen+2,26},{bt,"True background",kMagenta+1,27}},
    {{pp,t,kBlue+1,24},{ph,t,kRed+1,25},{hh,t,kGreen+2,26},{bt,t,kMagenta+1,27}},"Component / total","total");
  plot("background_closure",{{bt,"True background",kBlack,20},{d.h["Bmix"],"Mixed cone background",kRed+1,25},
    {d.h["Bperp"],"Perpendicular cone background",kBlue+1,28}},
    {{d.h["Bmix"],bt,kRed+1,25},{d.h["Bperp"],bt,kBlue+1,28}},"Estimate / truth","true-background");
  plot("signal_extraction",{{pp,"Pythia-Pythia (signal)",kBlack,20},
    {d.h["Smix"],"Total - mixed background",kRed+1,25},{d.h["Sperp"],"Total - perp. background",kBlue+1,28}},
    {{d.h["Smix"],pp,kRed+1,25},{d.h["Sperp"],pp,kBlue+1,28}},"Signal / truth","true-signal");
  const int colors[]={kRed+1,kBlue+1,kGreen+2,kMagenta+1,kOrange+7};
  const int markers[]={24,25,26,27,28};
  for(int method=0;method<2;++method) {
    const std::string m=method==0?"mix":"perp";
    const char* labels[]={"+ 0.5 JC_{1}","+ 0.5 JC_{2}","+ 0.5 C_{1}C_{1}","+ 0.5 C_{2}C_{2}","- C_{1}C_{2}"};
    Data terms;
    std::vector<Curve> curves={{d.h["B"+m],method==0?"Mixed cone estimate":"Perp. cone estimate",kBlack,20}};
    std::vector<Ratio> ratios={{d.h["B"+m],bt,kBlack,20}};
    for(int p=0;p<5;++p) {
      TH1D* term=terms.clone(d.h[m+std::to_string(p)],"term"+std::to_string(p),p==4?-1.:0.5);
      curves.push_back({term,labels[p],colors[p],markers[p]}); ratios.push_back({term,bt,colors[p],markers[p]});
    }
    plot(method==0?"mixed_signed_terms":"perpendicular_signed_terms",curves,ratios,"Term / truth","true-background");
  }
  plot("cone_symmetry",{{d.h["mix0"],"Mixed: JC_{1}",kRed+1,24},{d.h["mix1"],"Mixed: JC_{2}",kRed+1,25},
    {d.h["mix2"],"Mixed: C_{1}C_{1}",kGreen+2,26},{d.h["mix3"],"Mixed: C_{2}C_{2}",kGreen+2,27},
    {d.h["perp0"],"Perp.: JC_{1}",kBlue+1,24},{d.h["perp1"],"Perp.: JC_{2}",kBlue+1,25},
    {d.h["perp2"],"Perp.: C_{1}C_{1}",kMagenta+1,26},{d.h["perp3"],"Perp.: C_{2}C_{2}",kMagenta+1,27}},
    {{d.h["mix1"],d.h["mix0"],kRed+1,25},{d.h["mix3"],d.h["mix2"],kGreen+2,27},
    {d.h["perp1"],d.h["perp0"],kBlue+1,25},{d.h["perp3"],d.h["perp2"],kMagenta+1,27}},"Cone 2 / cone 1","cone-1");
}
}

void plotOOMCBackgroundPresentation(
  TString inclusiveName="merged_eecAnalysis_OOMC_GenGen_processed.root",
  TString splitName="merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root",
  TString outputDirectory="figures/OOMC_GenGen_Run3_presentation_20261102") {
  using namespace OOMCPresentation;
  gROOT->SetBatch(true); gStyle->SetPadGridX(0); gStyle->SetPadGridY(0);
  gStyle->SetOptStat(0); TGaxis::SetExponentOffset(-0.10,0.,"y"); TH1::AddDirectory(false);
  if(gSystem->mkdir(outputDirectory,true)!=0 && gSystem->AccessPathName(outputDirectory))
    throw std::runtime_error("Cannot create output directory");
  std::unique_ptr<TFile> inclusive(TFile::Open(inclusiveName,"READ")), split(TFile::Open(splitName,"READ"));
  if(!inclusive || inclusive->IsZombie() || !split || split->IsZombie()) throw std::runtime_error("Cannot open OO MC input files");
  EECCard inclusiveCard(inclusive.get()), splitCard(split.get());
  for(EECCard* card : {&inclusiveCard,&splitCard}) {
    if(!card->GetDataType().Contains("OO") || !card->GetDataType().Contains("MC")) throw std::runtime_error("Expected OO MC inputs");
    if(card->GetNWeightExponents()!=1 || card->GetWeightExponent()!=1) throw std::runtime_error("Only n=1 inputs are supported");
    if(card->GetBackgroundMethods()!=6) throw std::runtime_error("Expected BackgroundMethods=6 (mixed and perpendicular)");
  }
  TFile numerical(outputDirectory+"/numerical_histograms.root","RECREATE");
  if(numerical.IsZombie()) throw std::runtime_error("Cannot create numerical ROOT output");
  TNamed("formula","B=0.5*(JC1+JC2+same1+same2)-cross; Btrue=S1+S3; PP=S0; S1 already includes S2").Write();
  TNamed("normalization","Common raw input scale; no independent normalization; inclusive sums all native centrality slots").Write();
  TNamed("ratio_errors","Independent numerator/denominator error propagation; covariance unavailable; denominator<=0 points omitted").Write();
  std::ofstream manifest((outputDirectory+"/run_manifest.csv").Data());
  if(!manifest) throw std::runtime_error("Cannot create run manifest");
  manifest<<"selection_id,selection_flag,input_file,centrality_indices,centrality_low,centrality_high,jet_index,jet_low,jet_high,track_index,track_threshold,exponent,background_methods,status,figure,no_grids,upper_grid_x,upper_grid_y,lower_grid_x,lower_grid_y,upper_log_y,masked_ratio_points,positive_ratio_denominators,pdf,png\n";
  std::vector<std::pair<TString,TString>> pages;
  int plotted=0, skipped=0;
  for(int kind=0;kind<2;++kind) {
    TFile& input=kind==0?*inclusive:*split;
    EECCard& card=kind==0?inclusiveCard:splitCard;
    const int track=card.GetBinIndexTrackPtEEC(2.);
    if(track<0 || std::abs(card.GetLowBinBorderTrackPtEEC(track)-2.)>1.e-6) throw std::runtime_error("Missing exact 2 GeV track threshold");
    std::vector<std::vector<int>> centGroups;
    if(kind==0) {
      std::vector<int> all;
      for(int c=0;c<card.GetNCentralityBins();++c) all.push_back(c);
      centGroups.push_back(all);
    } else for(int c=0;c<card.GetNCentralityBins();++c) centGroups.push_back({c});
    for(const auto& cents : centGroups) {
      for(const auto& jet : std::vector<std::pair<double,double>>{{40,80},{80,120},{120,160},{160,200}}) {
        Selection s; s.flag=kind==0?"inclusive":"split"; s.input=input.GetName(); s.cents=cents;
        s.jetLow=jet.first; s.jetHigh=jet.second; s.track=track; s.jet=card.FindBinIndexJetPtEEC(jet);
        if(s.jet<0 || s.jet>=card.GetNJetPtBinsEEC()) throw std::runtime_error("Missing requested finite jet bin");
        if(kind==0) { s.id="inclusive"; s.label="GenGen, R = 0.4, inclusive sample"; }
        else {
          s.centLow=card.GetLowBinBorderCentrality(cents.front()); s.centHigh=card.GetHighBinBorderCentrality(cents.front());
          s.id=Form("C%.0f-%.0f",s.centLow,s.centHigh); s.label="GenGen, R = 0.4";
        }
        s.id+=Form("_J%.0f-%.0f_T2_n1",jet.first,jet.second);
        Data data; load(input,s,data);
        if(!populated(data.h["T"])) { row(manifest,s,"skipped_empty","","",""); ++skipped; continue; }
        if(kind==1) {
          const bool native=(s.centLow==0. && s.centHigh==30.) || (s.centLow==30. && s.centHigh==50.) ||
            (s.centLow==50. && s.centHigh==90.) || (s.centLow==90. && s.centHigh==100.);
          if(!native) { row(manifest,s,"skipped_nonrequested_native_slot","","",""); ++skipped; continue; }
        }
        saveHistograms(numerical,s,data); render(s,data,outputDirectory,manifest,pages); ++plotted;
      }
    }
  }
  for(const TString& flag : {TString("inclusive"),TString("split"),TString("all")}) {
    std::vector<TString> pdfs;
    for(const auto& page : pages) if(flag=="all" || page.first==flag) pdfs.push_back(page.second);
    packAtlas(pdfs,outputDirectory+"/"+flag+"_atlas.pdf");
  }
  numerical.Write(); numerical.Close(); manifest.close();
  std::cout<<"OO presentation complete: "<<plotted<<" selections, "<<plotted*6<<" figures, "<<skipped<<" empty/nonrequested selections skipped. Output: "<<outputDirectory<<std::endl;
}

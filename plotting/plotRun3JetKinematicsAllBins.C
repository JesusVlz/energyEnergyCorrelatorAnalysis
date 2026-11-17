// Final-selected inclusive jets from the existing framework; no forest event loop.
// Reuse the reference Run-3 QA axes and visual conventions.
#include "plotRun3JetQA.C"
#include <THnSparse.h>
#include <TKey.h>
#include <array>
#include <iomanip>

namespace Run3AllBins {
const std::array<double,12> edges{{5,10,20,30,40,80,120,160,200,300,500,5020}};
double Card(TFile& file,const char* key) {
  auto* v=dynamic_cast<TVectorT<float>*>(file.Get((std::string("JCard/")+key).c_str()));
  if(!v || v->GetNrows()!=1)throw std::runtime_error(std::string("Missing scalar JCard/")+key);
  return (*v)[v->GetLwb()];
}
void Add(TH1D& h,int b,double w,double variance) {
  h.SetBinContent(b,h.GetBinContent(b)+w);
  h.SetBinError(b,std::sqrt(h.GetBinError(b)*h.GetBinError(b)+variance));
}
double Variance(const TH1D& h) {
  double v=0;for(int b=0;b<=h.GetNbinsX()+1;++b)v+=std::pow(h.GetBinError(b),2);return v;
}
void Check(double a,double b,const char* what) {
  if(std::abs(a-b)>1e-9*std::max({1.,std::abs(a),std::abs(b)}))throw std::runtime_error(std::string("Conservation failed: ")+what);
}
void Draw(TH1D& counts,TFile& out,const std::string& prefix,const std::string& heading,
          const std::string& scope,const std::string& level,const std::string& interval,
          const std::string& role,const std::string& status,const std::string& x,const std::string& y,bool pt,bool overview=false) {
  out.cd();counts.Write();
  std::unique_ptr<TH1D> h(static_cast<TH1D*>(counts.Clone((std::string(counts.GetName())+"_density").c_str())));
  h->Scale(1.,"width");Run3Plot::Axes(*h,false,false);
  h->GetXaxis()->SetTitle(x.c_str());h->GetYaxis()->SetTitle(y.c_str());
  const int color=level=="Reco"?kBlue:kRed;
  h->SetLineColor(color);h->SetMarkerColor(color);h->SetLineWidth(2);h->SetMarkerStyle(kOpenCircle);h->SetMarkerSize(.95);
  double bottom=0,top=0,positive=std::numeric_limits<double>::max();bool negative=false;
  for(int b=1;b<=h->GetNbinsX();++b){
    double w=h->GetBinContent(b),e=h->GetBinError(b);bottom=std::min(bottom,w-e);top=std::max(top,w+e);
    if(w<0)negative=true;if(w>0)positive=std::min(positive,w);
  }
  bool log=pt && !negative && positive<std::numeric_limits<double>::max();
  h->SetMinimum(log?positive/12.:(bottom<0?bottom*1.15:0));h->SetMaximum(top>0?top*(log?20.:3.2):1.);
  TCanvas canvas(counts.GetName(),"Run-3 selected jet QA",850,850);canvas.SetCanvasSize(850,850);
  canvas.SetLeftMargin(.17);canvas.SetRightMargin(.04);canvas.SetTopMargin(.055);canvas.SetBottomMargin(.13);canvas.SetTicks(1,1);canvas.SetLogy(log);canvas.SetLogx(overview);
  h->Draw("E1");
  TLegend legend(.21,.69,.95,.94);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextFont(42);legend.SetTextSize(.026);legend.SetMargin(.075);
  for(const auto& text:std::vector<std::string>{"CMS Simulation Preliminary",heading,scope,level+" jets, R = 0.4, WTA, |#eta| < 1.6",interval,role})legend.AddEntry((TObject*)nullptr,text.c_str(),"");
  legend.AddEntry(h.get(),level=="Gen"?"Final selection; no generator-jet JEC":"Final selection; active reco JEC + quality cuts","lp");
  if(status!="POPULATED")legend.AddEntry((TObject*)nullptr,status=="OUTSIDE_SELECTION"?"Outside QA selection: p_{T} #geq 10 GeV":(status=="EMPTY_BIN"?"No weighted entries in this sample":"Signed weights cancel; uncertainties retained"),"");
  legend.Draw();canvas.RedrawAxis();canvas.Modified();canvas.Update();
  std::string stem=prefix+"_"+counts.GetName();canvas.SaveAs((stem+".png").c_str());canvas.SaveAs((stem+".pdf").c_str());
  out.cd();h->Write();canvas.Write((std::string(counts.GetName())+"_canvas").c_str());
}
}

void plotRun3JetKinematicsAllBins(const char* input,const char* prefix,const char* jetLevel,
                                double activityLow=-1,double activityHigh=100) {
  using namespace Run3AllBins;
  gROOT->SetBatch(true);TH1::AddDirectory(false);gStyle->SetOptStat(0);gStyle->SetEndErrorSize(2);
  std::unique_ptr<TFile> file(TFile::Open(input,"READ"));if(!file || file->IsZombie())throw std::runtime_error("Unreadable input");
  const std::string level=jetLevel;if(level!="Reco" && level!="Gen")throw std::runtime_error("Specify Reco or Gen");
  if(Card(*file,"McCorrelationType")!=(level=="Gen"?3:0) || Card(*file,"MinJetPtCut")!=10 || Card(*file,"MaxJetPtCut")!=5020 || Card(*file,"JetAxis")!=1 || Card(*file,"MatchJets")!=0 || Card(*file,"JetUncertainty")!=0 || Card(*file,"JetPtWeight")!=0 || std::abs(Card(*file,"JetEtaCut")-1.6)>1e-6)
    throw std::runtime_error("Input is not the audited v6 ordinary-jet profile");
  auto* cardEdges=dynamic_cast<TVectorT<float>*>(file->Get("JCard/JetPtBinEdgesEEC"));
  if(!cardEdges || cardEdges->GetNrows()!=12)throw std::runtime_error("Missing 11-bin profile");
  for(int i=0;i<12;++i)Check((*cardEdges)[i+cardEdges->GetLwb()],edges[i],"card boundaries");
  auto* source=dynamic_cast<THnSparseF*>(file->Get("inclusiveJet"));if(!source || source->GetNdimensions()!=6 || !source->GetCalculateErrors())throw std::runtime_error("Missing inclusiveJet/Sumw2");
  auto* pt=source->GetAxis(0);auto* activity=source->GetAxis(4);
  // This entry point deliberately accepts only the audited 5 GeV source bins.
  // Its overflow retains eta/phi, but cannot resolve pT above 500 GeV.
  if(pt->GetNbins()!=100 || pt->GetXmin()!=0 || pt->GetXmax()!=500)throw std::runtime_error("Unexpected native pT binning");
  const int first=activity->FindFixBin(activityLow),last=activity->FindFixBin(std::nextafter(activityHigh,activityLow));
  if(first<1 || last>activity->GetNbins() || activity->GetBinLowEdge(first)!=activityLow || activity->GetBinUpEdge(last)!=activityHigh)throw std::runtime_error("Non-aligned activity range");
  std::array<std::array<std::unique_ptr<TH1D>,3>,11> h;
  const char* observables[]={"Pt","Eta","Phi"};const int axes[]={0,2,1};
  for(int j=0;j<11;++j)for(int o=0;o<3;++o){
    const std::string name=std::string("jet")+observables[o]+Form("_%.0f_%.0f",edges[j],edges[j+1]);
    auto* ax=source->GetAxis(axes[o]);
    h[j][o]=std::make_unique<TH1D>(name.c_str(),"",o==0?(j==10?1:int((edges[j+1]-edges[j])/5)):ax->GetNbins(),o==0?edges[j]:ax->GetXmin(),o==0?edges[j+1]:ax->GetXmax());h[j][o]->Sumw2();
  }
  std::array<double,11> sum{},var{};std::array<long long,11> occupied{};
  double total=0,totalVar=0;int coord[6];
  for(Long64_t b=0;b<source->GetNbins();++b){
    double w=source->GetBinContent(b,coord),v=source->GetBinError2(b);if(coord[4]<first || coord[4]>last)continue;
    if(w==0 && v==0)continue;
    if(!std::isfinite(w) || !std::isfinite(v) || v<0)throw std::runtime_error("Non-finite native bin");
    int j=-1;
    if(coord[0]==101)j=10;
    else for(int k=0;k<10;++k)if(pt->GetBinLowEdge(coord[0])>=edges[k] && pt->GetBinUpEdge(coord[0])<=edges[k+1]){j=k;break;}
    if(j<0 || j==0)throw std::runtime_error("Native content outside ordinary selected pT range");
    sum[j]+=w;var[j]+=v;++occupied[j];total+=w;totalVar+=v;
    for(int o=0;o<3;++o){int dest=o==0?(j==10?1:h[j][o]->FindFixBin(pt->GetBinCenter(coord[0]))):coord[axes[o]];Add(*h[j][o],dest,w,v);}
  }
  const std::string outPath=std::string(prefix)+"_numerical.root";
  TFile out(outPath.c_str(),"RECREATE");if(out.IsZombie())throw std::runtime_error("Cannot create numerical ROOT");
  TNamed("input_file",input).Write();TNamed("source_histogram","inclusiveJet").Write();TNamed("jet_level",jetLevel).Write();
  TNamed("selection","Native final-selected jets; MinJetPtCut=10, MaxJetPtCut=5020, MatchJets=0, JetAxis=1; detector event veto retained. Gen has no jet JEC/TrackMax/JetID; Reco has configured JEC/TrackMax/JetID.").Write();
  TNamed("upper_bin_definition","500 <= pT <= 5020: aggregate native pT overflow, bounded by inclusiveJet producer selection. Cannot resolve sub-bins or isolate exact 5020 endpoint. Not a half-open [500,5020) measurement.").Write();
  TNamed("normalization","Weighted jet yields divided by displayed bin width; no unit-area or event normalization. Counts and Sumw2 retained, including angular flow bins.").Write();
  TNamed("qualification",Card(*file,"MaxEventsPerFile")>0?"PILOT_DIAGNOSTIC_NOT_PHYSICS_QUALIFIED":"FULL_SAMPLE_JET_QA_DIAGNOSTIC").Write();
  if(auto* card=file->GetDirectory("JCard")){auto* dest=out.mkdir("JCard");TIter next(card->GetListOfKeys());while(auto* key=dynamic_cast<TKey*>(next())){dest->cd();std::unique_ptr<TObject> obj(key->ReadObj());obj->Write();}out.cd();}
  std::ofstream csv(std::string(prefix)+"_coverage.csv");csv<<std::setprecision(17)<<"level,jet_low,jet_high,role,observable,status,sumw,sumw2,occupied_native_bins,underflow,overflow,upper_endpoint\n";
  const bool pp=Card(*file,"DataType")==12;
  const std::string heading=std::string(pp?"pp MC 5.36 TeV":"OO MC 5.36 TeV")+(Card(*file,"MaxEventsPerFile")>0?" (pilot)":"");
  const std::string scope=pp?(Card(*file,"Run3JetID")==0?"Inclusive reference; Jet ID disabled":"Inclusive reference; Jet ID enabled"):(activityLow<0?"Integrated native selected population":Form("Centrality %.0f#minus%.0f%%",activityLow,activityHigh));
  TH1D all("jetPt_allStored","",11,edges.data());all.Sumw2();
  double check=0,checkVar=0;
  for(int j=0;j<11;++j){
    const std::string role=j>=1 && j<=6?"Reporting bin":"Additional stored interval";
    const std::string status=j==0?"OUTSIDE_SELECTION":(occupied[j]==0?"EMPTY_BIN":(sum[j]==0?"SIGNED_CANCELLATION":"POPULATED"));
    std::string interval=Form("%.0f #leq p_{T} %s %.0f GeV",edges[j],j==10?"#leq":"<",edges[j+1]);if(j==10)interval+=" (overflow aggregate)";
    for(int o=0;o<3;++o){
      Check(h[j][o]->Integral(0,h[j][o]->GetNbinsX()+1),sum[j],"per-bin sumw");Check(Variance(*h[j][o]),var[j],"per-bin Sumw2");
      csv<<level<<','<<edges[j]<<','<<edges[j+1]<<','<<role<<','<<observables[o]<<','<<status<<','<<sum[j]<<','<<var[j]<<','<<occupied[j]<<','<<h[j][o]->GetBinContent(0)<<','<<h[j][o]->GetBinContent(h[j][o]->GetNbinsX()+1)<<','<<(j==10?"inclusive_native_cut":"exclusive")<<'\n';
      Draw(*h[j][o],out,prefix,heading,scope,level,interval,role,status,o==0?"Jet p_{T} (GeV)":(o==1?"Jet #eta":"Jet #phi"),o==0?"dN^{jet}/dp_{T} (GeV^{-1})":(o==1?"dN^{jet}/d#eta":"dN^{jet}/d#phi"),o==0);
    }
    all.SetBinContent(j+1,sum[j]);all.SetBinError(j+1,std::sqrt(var[j]));check+=sum[j];checkVar+=var[j];
  }
  Check(check,total,"all interval sumw");Check(checkVar,totalVar,"all interval Sumw2");
  Draw(all,out,prefix,heading,scope,level,"10 #leq p_{T} #leq 5020 GeV","All 11 stored intervals; six reporting bins: 10#minus160 GeV","POPULATED","Jet p_{T} (GeV)","dN^{jet}/dp_{T} (GeV^{-1})",true,true);
  out.Close();std::cout<<"VALIDATED_ALL_BINS "<<level<<" sumw="<<std::setprecision(17)<<total<<" sumw2="<<totalVar<<"\n";
}

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TH2D.h>
#include <fstream>
#include <TH3D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorT.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
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
std::string SelectionLabel(int stage,bool jetID=true) {
  if(!jetID && stage==5)return "TrackMax+VetoMap (Jet ID off)";
  const char* labels[]={"No Selection","TrackMax","JetID","TrackMax+JetID",
                        "JetID+VetoMap","TrackMax+JetID+VetoMap","VetoMap"};
  return labels[stage];
}
std::string CalibrationLabel(int stage) {
  const char* labels[]={"Raw","L1FastJet","L2Relative","L2Residual"};
  return labels[stage];
}
void Axes(TH1D& h,bool lower,bool split) {
  h.SetStats(false); h.SetTitle("");
  for(auto* axis:{h.GetXaxis(),h.GetYaxis()}) {
    axis->SetTitleFont(42); axis->SetLabelFont(42); axis->CenterTitle(true);
    axis->SetLabelSize(lower?.095:(split?.044:.035));
    axis->SetTitleSize(lower?.115:(split?.055:.045));
  }
  h.GetXaxis()->SetTitleOffset(lower?1.05:1.1);
  h.GetYaxis()->SetTitleOffset(lower?.68:(split?1.35:1.6));
  h.GetYaxis()->SetNdivisions(lower?505:510);
  h.GetXaxis()->SetNdivisions(506);
  if(split && !lower) { h.GetXaxis()->SetLabelSize(0); h.GetXaxis()->SetTitleSize(0); }
}
void Draw(TH3D* source,const std::vector<int>& stages,int low,int high,const std::string& output,
          bool ratios,const std::string& heading,const std::string& selection,const std::string& kinematics,const std::string& description,
          TDirectory* numerical,bool jetID=true) {
  const std::string base=output.substr(0,output.size()-4);
  const std::string name=base.substr(base.find_last_of('/')+1);
  TCanvas canvas(name.c_str(),"Run-3 jet QA",850,850); canvas.SetCanvasSize(850,850);
  TPad upper("top","",0,ratios?.30:0,1,1),lower("bottom","",0,0,1,.30);
  upper.SetLeftMargin(.17); upper.SetRightMargin(.04); upper.SetTopMargin(.055);
  upper.SetBottomMargin(ratios?.014:.13); upper.SetTicks(1,1); upper.SetLogy();
  upper.Draw(); upper.cd();
  std::vector<std::unique_ptr<TH1D>> histograms;
  std::vector<std::unique_ptr<TGraph>> ratioGraphs;
  const int colors[]={kBlue,kRed,kBlack,kGreen+2};
  double maximum=0,smallest=std::numeric_limits<double>::max();
  for(size_t j=0;j<stages.size();++j) {
    auto h=Project(source,stages[j],low,high,"stage_"+std::to_string(stages[j]));
    h->SetLineColor(colors[j]); h->SetMarkerColor(colors[j]);
    h->SetLineWidth(2); h->SetLineStyle(1); h->SetMarkerStyle(kOpenCircle); h->SetMarkerSize(.95);
    for(int b=1;b<=h->GetNbinsX();++b) {
      if(h->GetBinCenter(b)>450) continue;
      const double y=h->GetBinContent(b);
      maximum=std::max(maximum,y+h->GetBinError(b));
      if(y>0) smallest=std::min(smallest,y);
    }
    histograms.push_back(std::move(h));
  }
  if(maximum<=0) throw std::runtime_error("No positive entries for "+output);
  auto frame=std::unique_ptr<TH1D>(static_cast<TH1D*>(histograms[0]->Clone("spectrum_frame")));
  frame->Reset("ICES"); frame->GetXaxis()->SetRangeUser(0,450);
  frame->SetMinimum(smallest/12.); frame->SetMaximum(maximum*3.);
  frame->GetXaxis()->SetTitle(ratios?"Jet p_{T} (GeV)":"Jet p_{T} at calibration level (GeV)");
  frame->GetYaxis()->SetTitle("dN^{jet}/dp_{T} (GeV^{-1})"); Axes(*frame,false,ratios);
  frame->DrawCopy("AXIS");
  for(const auto& h:histograms) h->DrawCopy("E1 SAME");

  // The upper-right region is clear of the steeply falling spectrum and its sparse tails.
  TLegend legend(.46,ratios?.53:.60,.95,ratios?.92:.93);
  legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(42);
  legend.SetTextSize(ratios?.037:.029); legend.SetMargin(.085);
  legend.AddEntry((TObject*)nullptr,"CMS Simulation Preliminary","");
  legend.AddEntry((TObject*)nullptr,heading.c_str(),"");
  legend.AddEntry((TObject*)nullptr,selection.c_str(),"");
  legend.AddEntry((TObject*)nullptr,kinematics.c_str(),"");
  legend.AddEntry((TObject*)nullptr,description.c_str(),"");
  for(size_t j=0;j<stages.size();++j) {
    const auto label=ratios?SelectionLabel(stages[j],jetID):CalibrationLabel(stages[j]);
    legend.AddEntry(histograms[j].get(),label.c_str(),"lp");
  }
  legend.DrawClone();

  if(ratios) {
    double ratioLow=.80,ratioHigh=1.12;
    for(size_t j=1;j<histograms.size();++j) {
      auto ratio=std::make_unique<TGraph>(); ratio->SetName(("ratio_stage_"+std::to_string(stages[j])).c_str());
      ratio->SetMarkerColor(colors[j]); ratio->SetLineColor(colors[j]);
      ratio->SetMarkerStyle(kOpenCircle); ratio->SetMarkerSize(.95); ratio->SetLineStyle(1); ratio->SetLineWidth(2);
      for(int b=1;b<=histograms[0]->GetNbinsX();++b) {
        const double denominator=histograms[0]->GetBinContent(b);
        if(denominator<=0) continue; // Unsupported bins are omitted, never drawn as zero efficiency.
        const double value=histograms[j]->GetBinContent(b)/denominator;
        const double x=histograms[0]->GetBinCenter(b);
        if(std::isfinite(value)) {
          ratio->SetPoint(ratio->GetN(),x,value);
          if(x<=450) { ratioLow=std::min(ratioLow,value-.025); ratioHigh=std::max(ratioHigh,value+.025); }
        }
      }
      ratioGraphs.push_back(std::move(ratio));
    }
    canvas.cd(); lower.SetLeftMargin(.17); lower.SetRightMargin(.04);
    lower.SetTopMargin(.015); lower.SetBottomMargin(.30); lower.SetTicks(1,1);
    lower.Draw(); lower.cd();
    auto ratioFrame=std::unique_ptr<TH1D>(static_cast<TH1D*>(frame->Clone("ratio_frame")));
    ratioFrame->Reset("ICES"); ratioFrame->SetMinimum(ratioLow); ratioFrame->SetMaximum(ratioHigh);
    ratioFrame->GetYaxis()->SetTitle("#frac{X}{No Selection}");
    ratioFrame->GetXaxis()->SetTitle("Jet p_{T} (GeV)"); Axes(*ratioFrame,true,true);
    ratioFrame->DrawCopy("AXIS");
    for(double y=std::ceil(ratioLow*10)/10.;y<ratioHigh;y+=.1) {
      TLine line(0,y,450,y); line.SetLineColor(kGray+1); line.SetLineStyle(3); line.DrawClone();
    }
    TLine unity(0,1,450,1); unity.SetLineColor(kGray+1); unity.DrawClone();
    for(const auto& ratio:ratioGraphs) ratio->DrawClone("P SAME");
    lower.RedrawAxis();
  }
  upper.RedrawAxis(); canvas.Modified(); canvas.Update();
  canvas.SaveAs(output.c_str()); canvas.SaveAs((base+".pdf").c_str());
  auto* saved=numerical->mkdir(name.c_str()); saved->cd();
  TNamed("source_histogram",source->GetName()).Write();
  TNamed("selection",selection.c_str()).Write(); TNamed("plot_definition",description.c_str()).Write();
  for(size_t j=0;j<histograms.size();++j) {
    histograms[j]->Write();
    const auto label=ratios?SelectionLabel(stages[j],jetID):CalibrationLabel(stages[j]);
    TNamed(("label_stage_"+std::to_string(stages[j])).c_str(),label.c_str()).Write();
  }
  for(const auto& ratio:ratioGraphs) ratio->Write();
  canvas.Write("canvas");
}
}

namespace Run3Plot {
void Kinematics(TFile* input,TDirectory* numerical,const std::string& prefix,int low,int high,
                const std::string& heading,const std::string& selection) {
  const double bounds[]={10,20,30,40,80,120,160};
  for(const char* coordinate:{"Eta","Phi"}) {
    auto* source=dynamic_cast<TH3D*>(input->Get((std::string("run3JetQA/selectedJet")+coordinate+"Weighted").c_str()));
    if(!source){std::cout<<"MISSING_INPUT selected jet "<<coordinate<<"; requires new native QA producer\n";continue;}
    for(int j=0;j<6;++j){
      const int first=source->GetYaxis()->FindFixBin(bounds[j]);
      const int last=source->GetYaxis()->FindFixBin(std::nextafter(bounds[j+1],bounds[j]));
      if(source->GetYaxis()->GetBinLowEdge(first)!=bounds[j] || source->GetYaxis()->GetBinUpEdge(last)!=bounds[j+1])
        throw std::runtime_error("Jet kinematics lacks exact requested pT boundaries");
      const std::string name=std::string("jet")+coordinate+Form("_%.0f_%.0f",bounds[j],bounds[j+1]);
      std::unique_ptr<TH1D> h(source->ProjectionX(name.c_str(),first,last,low,high,"e"));h->SetDirectory(nullptr);
      const double yield=h->Integral(0,h->GetNbinsX()+1);
      std::cout<<"JET_KINEMATICS_COUNT "<<selection<<' '<<coordinate<<' '<<bounds[j]<<' '<<bounds[j+1]<<' '<<yield<<'\n';
      numerical->cd();h->Write((name+"_counts").c_str());
      if(h->GetEntries()==0){std::cout<<"EMPTY no figure for "<<name<<'\n';continue;}
      h->Scale(1.,"width");Axes(*h,false,false);
      h->GetYaxis()->SetTitle(coordinate==std::string("Eta")?"dN^{jet}/d#eta":"dN^{jet}/d#phi");
      double envelope=0;for(int b=1;b<=h->GetNbinsX();++b)envelope=std::max(envelope,h->GetBinContent(b)+h->GetBinError(b));
      h->SetMarkerStyle(20);h->SetMarkerSize(.7);h->SetLineColor(kBlack);h->SetMinimum(0);h->SetMaximum(envelope*1.6);
      TCanvas canvas(name.c_str(),"Selected jet kinematics",850,850);canvas.SetLeftMargin(.17);canvas.SetRightMargin(.04);canvas.SetBottomMargin(.13);canvas.SetTopMargin(.06);canvas.SetTicks(1,1);
      h->Draw("E1");TLegend legend(.19,.73,.95,.93);legend.SetFillStyle(0);legend.SetBorderSize(0);legend.SetTextFont(42);legend.SetTextSize(.027);
      legend.AddEntry((TObject*)nullptr,"CMS Simulation Preliminary","");legend.AddEntry((TObject*)nullptr,heading.c_str(),"");
      legend.AddEntry((TObject*)nullptr,selection.c_str(),"");legend.AddEntry((TObject*)nullptr,Form("%.0f #leq p_{T}^{jet} < %.0f GeV, WTA axis",bounds[j],bounds[j+1]),"");legend.Draw();
      canvas.SaveAs((prefix+"_"+name+".png").c_str());canvas.SaveAs((prefix+"_"+name+".pdf").c_str());
      numerical->cd();h->Write();canvas.Write((name+"_canvas").c_str());
    }
  }
}
}

namespace Run3Plot {
void Acceptance(TFile* input,TDirectory* numerical,const std::string& prefix,double activityLow,double activityHigh,
                const std::string& heading,const std::string& selection) {
  auto* before=dynamic_cast<TH3D*>(input->Get("run3JetQA/signalVetoEtaPhiActivityBefore"));
  auto* after=dynamic_cast<TH3D*>(input->Get("run3JetQA/signalVetoEtaPhiActivityAfter"));
  auto* events=dynamic_cast<TH2D*>(input->Get("run3JetQA/signalVetoEventsActivity"));
  auto* veto=dynamic_cast<TH2*>(input->Get("run3JetQA/configuredJetVetoMap"));
  if(!before || !after || !events || !veto) {std::cout<<"MISSING_INPUT activity-resolved native veto QA\n";return;}
  const int first=before->GetZaxis()->FindFixBin(activityLow);
  const int last=before->GetZaxis()->FindFixBin(std::nextafter(activityHigh,activityLow));
  if(first<1 || last>before->GetNbinsZ() || before->GetZaxis()->GetBinLowEdge(first)!=activityLow || before->GetZaxis()->GetBinUpEdge(last)!=activityHigh)
    throw std::runtime_error("Requested activity does not align with veto occupancy axis");
  double denominator=0,kept=0,rejected=0;
  for(int b=first;b<=last;++b){denominator+=events->GetBinContent(1,b);kept+=events->GetBinContent(2,b);rejected+=events->GetBinContent(3,b);}
  if(denominator!=kept+rejected)throw std::runtime_error("Veto event denominator does not equal kept plus rejected");
  std::ofstream table(prefix+"_veto_events.csv");
  table<<"activity_low,activity_high,tested_events,kept_events,rejected_events,rejected_fraction,status\n";
  table<<activityLow<<','<<activityHigh<<','<<denominator<<','<<kept<<','<<rejected<<',';
  if(denominator>0)table<<rejected/denominator;table<<','<<(denominator>0?"PILOT_OR_PRODUCTION_COUNTS":"EMPTY_DENOMINATOR")<<'\n';
  numerical->cd();events->Write("veto_event_counts_by_physical_activity");
  TNamed("veto_denominator","Unweighted signal events tested after preceding event cuts, including events with zero veto-eligible jets; tested = kept + rejected. Detector coordinates, before whole-event veto.").Write();
  for(int stage=0;stage<3;++stage){
    std::unique_ptr<TH2> h;
    const std::string name=stage==0?"veto_occupancy_before":(stage==1?"veto_occupancy_after":"applied_veto_map");
    if(stage==2)h.reset(static_cast<TH2*>(veto->Clone(name.c_str())));
    else{
      auto* source=stage==0?before:after;source->GetZaxis()->SetRange(first,last);
      h.reset(dynamic_cast<TH2*>(source->Project3D("yxe")));source->GetZaxis()->SetRange(0,0);
      if(!h)throw std::runtime_error("Cannot project veto occupancy");h->SetName(name.c_str());
    }
    h->SetDirectory(nullptr);h->SetStats(false);h->SetTitle("");
    h->GetXaxis()->SetTitle("Detector jet #eta");h->GetYaxis()->SetTitle("Detector jet #phi");
    TCanvas canvas(name.c_str(),"Detector acceptance",900,850);canvas.SetLeftMargin(.13);canvas.SetRightMargin(.16);canvas.SetBottomMargin(.12);canvas.SetTopMargin(.20);
    h->Draw("COLZ");TLegend legend(.14,.81,.94,.99);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextFont(42);legend.SetTextSize(.026);
    legend.AddEntry((TObject*)nullptr,"CMS Simulation Preliminary","");legend.AddEntry((TObject*)nullptr,heading.c_str(),"");legend.AddEntry((TObject*)nullptr,selection.c_str(),"");
    legend.AddEntry((TObject*)nullptr,stage==0?"Veto-eligible jets before event veto":(stage==1?"Veto-eligible jets in kept events":"Actual configured veto map"),"");legend.Draw();
    canvas.SaveAs((prefix+"_"+name+".png").c_str());canvas.SaveAs((prefix+"_"+name+".pdf").c_str());
    numerical->cd();h->Write();canvas.Write((name+"_canvas").c_str());
  }
  std::cout<<"VETO_EVENTS "<<selection<<" tested="<<denominator<<" kept="<<kept<<" rejected="<<rejected<<'\n';
}
}

// Centrality range is half-open; use -1,0 for pp. A range of -1,100 is integrated.
// Read raw merged analyzer outputs; no separately projected ROOT input is needed.
void plotRun3JetQA(const char* input,const char* prefix="run3",double centralityLow=-1,double centralityHigh=100) {
  if(!(centralityLow<centralityHigh)) throw std::invalid_argument("Empty centrality range");
  gROOT->SetBatch(true); TH1::AddDirectory(false);
  gStyle->SetOptStat(0); gStyle->SetPaperSize(20,20); gStyle->SetEndErrorSize(2);
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
  auto* type=dynamic_cast<TVectorT<float>*>(file->Get("JCard/DataType"));
  const int dataType=type?std::lround((*type)[type->GetLwb()]):-1;
  const bool oo=dataType==9 || dataType==10;
  std::string heading=dataType==10?"Pythia+Hijing OO 5.36 TeV":
    (dataType==12?"Pythia pp 5.36 TeV":(oo?"OO 5.36 TeV":"pp 5.36 TeV"));
  auto* eventCap=dynamic_cast<TVectorT<float>*>(file->Get("JCard/MaxEventsPerFile"));
  if(eventCap && (*eventCap)[eventCap->GetLwb()]>0)heading+=" (pilot)";
  std::string centrality="Reco";
  if(oo) centrality+=centralityLow<0?", centrality integrated":Form(", Cent: %.0f#minus%.0f%%",centralityLow,centralityHigh);
  auto* minimumPt=dynamic_cast<TVectorT<float>*>(file->Get("JCard/MinJetPtCut"));
  if(!minimumPt) throw std::runtime_error("Missing jet-pT selection metadata");
  const std::string kinematics=Form("Jet p_{T} #geq %.0f GeV",(*minimumPt)[minimumPt->GetLwb()]);
  const std::string numericalPath=std::string(prefix)+"_numerical.root";
  std::unique_ptr<TFile> numerical(TFile::Open(numericalPath.c_str(),"RECREATE"));
  if(!numerical || numerical->IsZombie()) throw std::runtime_error("Cannot create "+numericalPath);
  TNamed("input_file",input).Write();
  TNamed("common_jet_pt_selection",kinematics.c_str()).Write();
  TNamed("normalization","Rebin by five and divide by bin width. Weighted yields at common exposure; no separate unit-area normalization.").Write();
  TNamed("legend_definitions","No Selection retains common event/kinematic cuts. TrackMax uses trackMax/rawPt. VetoMap is the configured whole-event veto. Cumulative labels list all applied quality cuts.").Write();
  TNamed("ratio_uncertainties","Central values only, unchanged from original jet QA plots. No independent/binomial or EEC rho=1 error assumption.").Write();
  auto* idCard=dynamic_cast<TVectorT<float>*>(file->Get("JCard/Run3JetID"));
  const bool jetID=!idCard || (*idCard)[idCard->GetLwb()]!=0;
  if(!jetID)centrality+=", Jet ID disabled";
  Run3Plot::Draw(selection,jetID?std::vector<int>{0,1,3,5}:std::vector<int>{0,1,5},low,high,std::string(prefix)+"_selection_cumulative.png",true,heading,centrality,kinematics,"Cumulative selections",numerical.get(),jetID);
  const bool independentVeto=selection->GetYaxis()->GetNbins()>=7;
  Run3Plot::Draw(selection,jetID?std::vector<int>{0,1,2,independentVeto?6:4}:std::vector<int>{0,1,6},low,high,std::string(prefix)+"_selection_independent.png",true,heading,centrality,kinematics,independentVeto?"Independent selections":"JetID+Veto is conditional; independent veto unavailable",numerical.get(),jetID);
  std::vector<int> active;
  for(int level=0;level<4;++level) {
    auto p=Run3Plot::Project(calibration,level,low,high,"level_check_"+std::to_string(level));
    if(p->Integral(0,p->GetNbinsX()+1)!=0) active.push_back(level);
  }
  if(!active.empty()) Run3Plot::Draw(calibration,active,low,high,std::string(prefix)+"_calibration_fixed_jets.png",false,heading,centrality,kinematics,"Fixed selected jets",numerical.get());
  Run3Plot::Kinematics(file.get(),numerical.get(),prefix,low,high,heading,centrality);
  Run3Plot::Acceptance(file.get(),numerical.get(),prefix,centralityLow,centralityHigh,heading,centrality);
  numerical->Write(); numerical->Close();
  std::cout<<"Ratios are central values only; no independent/binomial error assumption. Do not normalize each selection separately.\n";
}

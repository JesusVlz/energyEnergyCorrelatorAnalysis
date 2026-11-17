// Read-only views of the existing tracking payloads; no table derivation or averaging.
#include <TH2.h>
#include "JDrawer.h"
#include <TFile.h>
#include <TH2D.h>
#include <TH3.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <vector>

void plotRun3TrackingTables(const char* input,const char* output,const char* system,const char* workingPoint) {
  gROOT->SetBatch(true);gStyle->SetOptStat(0);TH1::AddDirectory(false);
  const bool oo=std::string(system)=="OO";
  if(!oo && std::string(system)!="pp")throw std::runtime_error("Expected OO or pp");
  std::unique_ptr<TFile> file(TFile::Open(input,"READ"));
  if(!file || file->IsZombie())throw std::runtime_error("Cannot open tracking payload");
  if(gSystem->mkdir(output,true)!=0 && gSystem->AccessPathName(output))throw std::runtime_error("Cannot create output");
  const TString prefix=TString(output)+"/";
  TFile numerical(prefix+"tracking_tables.root","RECREATE");
  TNamed("source_payload",input).Write();
  TNamed("interpretation","Existing payload bin contents/errors, without averaging, filling or refitting. Not a new MC efficiency derivation or closure.").Write();
  TNamed("lookup_policy",oo?"Runtime uses physical eta and event hiBin; existing 6.1 GeV pT cap retained.":"Runtime uses physical eta and pT; pp inclusive.").Write();
  std::ofstream csv((prefix+"tracking_table_bins.csv").Data());
  csv<<std::setprecision(17)<<"component,eta_low,eta_high,pt_low,pt_high,hiBin_low,hiBin_high,value,stored_error\n";
  const char* components[]={"Eff","Fak","Sec","Mul"};
  const char* titles[]={"Tracking efficiency","Fake fraction","Secondary fraction","Multiple-reconstruction rate"};
  int figures=0;
  for(int k=0;k<4;++k){
    const TString name=Form("h%s_%dD",components[k],oo?3:2);
    auto* source=dynamic_cast<TH1*>(file->Get(name));
    if(!source || source->GetDimension()!=(oo?3:2))throw std::runtime_error("Missing expected table component");
    auto* h3=dynamic_cast<TH3*>(source);auto* h2=dynamic_cast<TH2*>(source);
    if(oo && !h3)throw std::runtime_error("Expected TH3 payload");
    if(!oo && !h2)throw std::runtime_error("Expected TH2 payload");
    numerical.cd();source->Write(name);
    std::vector<double> x,y;
    for(int i=1;i<=source->GetNbinsX()+1;++i)x.push_back(source->GetXaxis()->GetBinLowEdge(i));
    for(int i=1;i<=source->GetNbinsY()+1;++i)y.push_back(source->GetYaxis()->GetBinLowEdge(i));
    for(int z=1;z<=(oo?source->GetNbinsZ():1);++z){
      const double lo=oo?source->GetZaxis()->GetBinLowEdge(z):-1;
      const double hi=oo?source->GetZaxis()->GetBinUpEdge(z):100;
      const TString stem=Form("%s_%s_%s_hiBin%g_%g",system,workingPoint,components[k],lo,hi);
      TH2D view(stem,"",x.size()-1,x.data(),y.size()-1,y.data());view.Sumw2();
      for(int ix=1;ix<=source->GetNbinsX();++ix)for(int iy=1;iy<=source->GetNbinsY();++iy){
        int b=oo?h3->GetBin(ix,iy,z):h2->GetBin(ix,iy);
        const double value=source->GetBinContent(b),error=source->GetBinError(b);
        view.SetBinContent(ix,iy,value);view.SetBinError(ix,iy,error);
        csv<<components[k]<<','<<x[ix-1]<<','<<x[ix]<<','<<y[iy-1]<<','<<y[iy]<<','<<lo<<','<<hi<<','<<value<<','<<error<<'\n';
      }
      numerical.cd();view.Write();
      // Restrict display only. Numerical output retains every original bin.
      view.GetXaxis()->SetRangeUser(-2.4,2.4);view.GetYaxis()->SetRangeUser(.5,500.);
      view.GetZaxis()->SetTitle(titles[k]);
      // Existing JDrawer owns overlapping canvas/pad pointers. Match the
      // process-lifetime drawer convention used by Run3ResponsePlotTools.
      static auto* owner=new JDrawer();auto& drawer=*owner;drawer.Reset();
      drawer.SetCanvasSize(850,780);drawer.SetLeftMargin(.15);drawer.SetRightMargin(.20);
      drawer.SetTopMargin(.18);drawer.SetBottomMargin(.14);drawer.SetLogY(true);
      drawer.DrawHistogram(&view,"Track #eta","Track p_{T} (GeV)","","COLZ");
      TLatex label;label.SetNDC();label.SetTextFont(42);label.SetTextSize(.031);
      label.DrawLatex(.15,.96,"CMS Simulation Preliminary");
      label.DrawLatex(.15,.915,Form("%s MC tracking payload; %s",system,workingPoint));
      label.SetTextSize(.025);
      label.DrawLatex(.15,.875,oo?Form("Physical hiBin [%g, %g); stored table values",lo,hi):"Inclusive pp; stored table values");
      if(oo)label.DrawLatex(.15,.84,"Runtime lookup retains the existing 6.1 GeV cap");
      auto* canvas=gPad->GetCanvas();canvas->SaveAs(prefix+stem+".png");canvas->SaveAs(prefix+stem+".pdf");
      numerical.cd();canvas->Write(stem+"_canvas");delete canvas;++figures;
    }
  }
  numerical.Close();std::cout<<"TRACKING_TABLE_VIEWS_SAVED "<<figures<<std::endl;
}

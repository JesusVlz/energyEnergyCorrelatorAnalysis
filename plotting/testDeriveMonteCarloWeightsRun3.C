// Synthetic ROOT smoke test. This is a SOFTWARE test, not physics closure.
// Run beside deriveMonteCarloWeightsRun3.C and the repository's JDrawer.h:
// root -l -b -q 'testDeriveMonteCarloWeightsRun3.C()'
#include "deriveMonteCarloWeightsRun3.C"

namespace Run3WeightSmokeTest {
void WriteSample(const TString& path, bool oo, bool data) {
  TFile file(path,"RECREATE");
  if (file.IsZombie()) throw std::runtime_error("Cannot create synthetic input.");
  TH1D vz("vertexZ","",80,-20,20), vzw("vertexZweighted","",80,-20,20);
  TH1D cent("centrality","",202,-1,100), centw("centralityWeighted","",202,-1,100);
  for (TH1D* h : {&vz,&vzw,&cent,&centw}) { h->SetDirectory(nullptr); h->Sumw2(); }
  const double exposure=data ? 5.0 : 1.0;
  for (int i=1; i<=vz.GetNbinsX(); ++i) {
    const double x=vz.GetBinCenter(i),t=x/15.0;
    if (std::abs(x)>=15) continue;
    const double base=1.e5*std::exp(-0.5*std::pow((x-0.2)/4.2,2));
    const double weight=0.9-0.05*t+0.32*t*t+0.07*t*t*t*t;
    const double count=base*exposure*(data ? weight : 1.0);
    vz.SetBinContent(i,count); vz.SetBinError(i,std::sqrt(count));
    vzw.SetBinContent(i,data ? count : base*weight);
    vzw.SetBinError(i,data ? std::sqrt(count) : std::sqrt(base)*weight);
  }
  for (int i=1; i<=cent.GetNbinsX(); ++i) {
    const double x=cent.GetXaxis()->GetBinLowEdge(i);
    if (x<0 || x>=90) continue;
    const double t=(x-30)/(x<30 ? 30.0 : 60.0);
    const double weight=x<30 ? 1+0.4*t+0.15*t*t : 1-0.2*t+0.1*t*t;
    const double base=2.e4*(0.8+std::exp(-x/35.0));
    const double count=base*exposure*(data ? weight : 1.0);
    cent.SetBinContent(i,count); cent.SetBinError(i,std::sqrt(count));
    centw.SetBinContent(i,data ? count : base*weight);
    centw.SetBinError(i,data ? std::sqrt(count) : std::sqrt(base)*weight);
  }
  vz.Write(); vzw.Write(); cent.Write(); centw.Write();
  TDirectory* card=file.mkdir("JCard"); card->cd();
  TVectorF v(1);
  v[0]=oo ? (data ? 9 : 10) : (data ? 11 : 12); v.Write("DataType");
  v[0]=15; v.Write("ZVertexCut");
  v[0]=0; v.Write("MultiplicityMode"); v.Write("JetPtWeight");
  v[0]=data ? 8 : 0; v.Write("TriggerSelection");
  TObjString("synthetic_software_test_not_physics").Write("GitHash");
  file.Close();
}

void CheckShape(TFile& result, const char* prefix) {
  auto* data=dynamic_cast<TH1*>(result.Get(TString(prefix)+"/data_normalized"));
  auto* preview=dynamic_cast<TH1*>(result.Get(TString(prefix)+"/mc_refit_preview"));
  Run3MCWeights::Require(data && preview,"Missing synthetic-test output histograms.");
  double maximum=0;
  for (int i=1; i<=data->GetNbinsX(); ++i)
    maximum=std::max(maximum,std::abs(preview->GetBinContent(i)/data->GetBinContent(i)-1));
  Run3MCWeights::Require(maximum<1e-7,std::string(prefix)+": synthetic polynomial shape was not recovered.");
}
}

void testDeriveMonteCarloWeightsRun3() {
  using namespace Run3WeightSmokeTest;
  try {
    gROOT->SetBatch(kTRUE);
    const TString base=TString::Format("%s/eec_weight_test_%d",gSystem->TempDirectory(),gSystem->GetPid());
    gSystem->mkdir(base,kTRUE);
    const TString od=base+"/oo_data.root",om=base+"/oo_mc.root";
    const TString pd=base+"/pp_data.root",pm=base+"/pp_mc.root";
    WriteSample(od,true,true); WriteSample(om,true,false);
    WriteSample(pd,false,true); WriteSample(pm,false,false);
    deriveMonteCarloWeightsRun3(od,om,pd,pm,base+"/output","both");
    TFile result(base+"/output/run3_mc_weights.root","READ");
    auto* complete=dynamic_cast<TParameter<int>*>(result.Get("complete"));
    Run3MCWeights::Require(complete && complete->GetVal()==1,"Macro did not finish.");
    CheckShape(result,"OO/vz"); CheckShape(result,"OO/centrality"); CheckShape(result,"pp/vz");
    auto* left=dynamic_cast<TF1*>(result.Get("OO/centrality/weight_centrality_0_30"));
    auto* right=dynamic_cast<TF1*>(result.Get("OO/centrality/weight_centrality_30_max"));
    Run3MCWeights::Require(left && right && std::abs(left->Eval(30)-right->Eval(30))<1e-8,
                          "Centrality continuity failed.");
    for (const char* figure : {"OO_eventWeights_stored.pdf","OO_eventWeights_refit.pdf",
                               "pp_vz_stored.png","pp_vz_refit.png"})
      Run3MCWeights::Require(!gSystem->AccessPathName(base+"/output/"+figure),"Missing figure output.");
    // A data-populated bin with no MC must fail derivation, not silently become zero weight.
    Run3MCWeights::Distribution bad;
    bad.data=new TH1D("bad_data","",10,-15,15);
    bad.before=new TH1D("bad_mc","",10,-15,15);
    bad.data->SetDirectory(nullptr); bad.before->SetDirectory(nullptr);
    bad.data->SetBinContent(1,1); bad.data->SetBinError(1,1);
    bool rejected=false; std::ostringstream sink;
    try { Run3MCWeights::Derive(bad,sink); } catch (const std::exception&) { rejected=true; }
    Run3MCWeights::Require(rejected,"Zero-MC support guard did not reject the input.");
    delete bad.data; delete bad.before;
    std::cout << "PASS: synthetic fits, export continuity, plots and zero-MC guard.\n"
              << "Test output: " << base << '\n';
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    gSystem->Exit(1);
  }
}

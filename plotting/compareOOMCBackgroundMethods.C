#include "plotOOMCBackgroundPresentation.C"

// Direct method comparison, on the same raw EEC scale as the closure figures.
// Run from the repository root with the CMSSW ROOT environment:
// root -l -b -q 'plotting/compareOOMCBackgroundMethods.C()'
void compareOOMCBackgroundMethods(
  TString inclusiveName="merged_eecAnalysis_OOMC_GenGen_processed.root",
  TString splitName="merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root",
  TString outputDirectory="figures/OOMC_GenGen_Run3_bkg_method_comparison_20261102") {
  using namespace OOMCPresentation;
  gROOT->SetBatch(true); gStyle->SetPadGridX(0); gStyle->SetPadGridY(0);
  gStyle->SetOptStat(0); TGaxis::SetExponentOffset(-0.10,0.,"y"); TH1::AddDirectory(false);
  if(gSystem->mkdir(outputDirectory,true)!=0 && gSystem->AccessPathName(outputDirectory))
    throw std::runtime_error("Cannot create method-comparison output directory");
  std::unique_ptr<TFile> inclusive(TFile::Open(inclusiveName,"READ")), split(TFile::Open(splitName,"READ"));
  if(!inclusive || inclusive->IsZombie() || !split || split->IsZombie())
    throw std::runtime_error("Cannot open OO MC inputs");
  EECCard inclusiveCard(inclusive.get()), splitCard(split.get());
  for(EECCard* card : {&inclusiveCard,&splitCard}) {
    if(!card->GetDataType().Contains("OO") || !card->GetDataType().Contains("MC"))
      throw std::runtime_error("Expected OO MC inputs");
    if(card->GetNWeightExponents()!=1 || card->GetWeightExponent()!=1)
      throw std::runtime_error("Expected n=1 inputs");
    if(card->GetBackgroundMethods()!=6)
      throw std::runtime_error("Expected mixed and perpendicular background methods");
  }
  TFile numerical(outputDirectory+"/numerical_histograms.root","RECREATE");
  if(numerical.IsZombie()) throw std::runtime_error("Cannot create numerical ROOT output");
  TNamed("formula","B=0.5*(JC1+JC2+same1+same2)-cross, independently for mixed and perpendicular cones").Write();
  TNamed("normalization","Common raw input scale; no independent normalization; inclusive sums all native centrality slots").Write();
  TNamed("ratio","Bmix/Bperp; nonpositive or nonfinite Bperp bins omitted").Write();
  TNamed("ratio_errors","Independent numerator/denominator error propagation; shared-event covariance unavailable").Write();
  std::ofstream manifest((outputDirectory+"/run_manifest.csv").Data());
  if(!manifest) throw std::runtime_error("Cannot create method-comparison manifest");
  manifest<<"selection_id,selection_flag,input_file,centrality_indices,centrality_low,centrality_high,jet_index,jet_low,jet_high,track_index,track_threshold,exponent,background_methods,status,figure,no_grids,upper_grid_x,upper_grid_y,lower_grid_x,lower_grid_y,upper_log_y,masked_ratio_points,positive_ratio_denominators,pdf,png\n";
  std::vector<std::pair<TString,TString>> pages;
  int plotted=0, skipped=0;
  for(int kind=0;kind<2;++kind) {
    TFile& input=kind==0?*inclusive:*split;
    EECCard& card=kind==0?inclusiveCard:splitCard;
    const int track=card.GetBinIndexTrackPtEEC(2.);
    if(track<0 || std::abs(card.GetLowBinBorderTrackPtEEC(track)-2.)>1.e-6)
      throw std::runtime_error("Missing exact 2 GeV track threshold");
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
          const bool requested=(s.centLow==0. && s.centHigh==30.) || (s.centLow==30. && s.centHigh==50.) ||
            (s.centLow==50. && s.centHigh==90.) || (s.centLow==90. && s.centHigh==100.);
          if(!requested) { row(manifest,s,"skipped_nonrequested_native_slot","","",""); ++skipped; continue; }
        }
        const TString family="background_method_comparison", base=outputDirectory+"/"+s.id+"_"+family;
        const Ratio ratio={data.h["Bmix"],data.h["Bperp"],kRed+1,25};
        // Reuse the established closure layout, including its upper-left legend.
        const PlotResult result=draw(s,family,
          {{data.h["Bmix"],"Mixed cone background",kRed+1,25},
           {data.h["Bperp"],"Perpendicular cone background",kBlue+1,28}},
          {ratio},"Mixed / perp.","perpendicular-background",base);
        if(!result.noGrids) throw std::runtime_error("Unexpected grid on comparison canvas");
        row(manifest,s,"plotted",family,base+".pdf",base+".png",result);
        saveHistograms(numerical,s,data);
        numerical.GetDirectory(s.id)->cd();
        PlotResult ratioResult;
        graph(ratio,ratioResult)->Write("Bmix_over_Bperp");
        numerical.cd();
        pages.push_back({s.flag,base+".pdf"}); ++plotted;
      }
    }
  }
  for(const TString& flag : {TString("inclusive"),TString("split"),TString("all")}) {
    std::vector<TString> pdfs;
    for(const auto& page : pages) if(flag=="all" || page.first==flag) pdfs.push_back(page.second);
    packAtlas(pdfs,outputDirectory+"/"+flag+"_atlas.pdf");
  }
  numerical.Write(); numerical.Close(); manifest.close();
  std::cout<<"OO background-method comparison complete: "<<plotted<<" figures, "<<skipped
    <<" empty/nonrequested selections skipped. Output: "<<outputDirectory<<std::endl;
}

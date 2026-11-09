#ifndef EEC_RUN3_JET_QA_H
#define EEC_RUN3_JET_QA_H
#include "ConfigurationCard.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3D.h>
#include <TNamed.h>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

// Diagnostics are filled before the event-veto early return. Selection spectra
// all use the same final calibrated pT and the same input-event exposure.
// Correction spectra use exactly the same final-selected jets at every level.
// Sumw2 is a jet-weight sum, not an event-cluster covariance estimator.
class Run3JetQA {
public:
  explicit Run3JetQA(ConfigurationCard* card)
    : fRequestedJEC(card->HasKey("Run3JEC") && card->Get("Run3JEC")==1),
      fPpL1(card->HasKey("Run3UsePpL1FastJet") && card->Get("Run3UsePpL1FastJet")==1),
      fOO(card->Get("DataType")==9 || card->Get("DataType")==10),
      fData(card->Get("DataType")==9 || card->Get("DataType")==11),
      fVeto(card->HasKey("ApplyJetVeto") && card->Get("ApplyJetVeto")==1),
      fTable(card->HasKey("Run3JetIDTable")?card->GetStr("Run3JetIDTable").Data():"MiniPOG20261108") {
    TDirectory::TContext context(nullptr);
    fSelection.reset(new TH3D("selectionPtWeighted",";Final calibrated jet p_{T} (GeV);Selection;Centrality (%)",500,0,1000,6,-.5,5.5,101,-1,100));
    fCounts.reset(new TH3D("selectionPtCounts",";Final calibrated jet p_{T} (GeV);Selection;Centrality (%)",500,0,1000,6,-.5,5.5,101,-1,100));
    fCalibration.reset(new TH3D("calibrationPtWeighted",";Jet p_{T} at this calibration level (GeV);Actual JEC level;Centrality (%)",500,0,1000,4,-.5,3.5,101,-1,100));
    fTrackFraction.reset(new TH3D("trackMaxOverRawPt",";trackMax/rawPt;Raw jet p_{T} (GeV);Centrality (%)",120,0,1.2,100,0,1000,101,-1,100));
    for(TH3D* h:{fSelection.get(),fCounts.get(),fCalibration.get(),fTrackFraction.get()}) { h->SetDirectory(nullptr); h->Sumw2(); }
    const char* selections[]={"No jet-quality selection","trackMax/rawPt only","Jet ID only","trackMax/rawPt AND Jet ID","Jet ID AND event veto","trackMax/rawPt AND Jet ID AND event veto"};
    for(int i=0;i<6;++i) {fSelection->GetYaxis()->SetBinLabel(i+1,selections[i]);fCounts->GetYaxis()->SetBinLabel(i+1,selections[i]);}
    const char* corrections[]={"rawpt","after L1FastJet","after L2Relative","after L2Residual (data only)"};
    for(int i=0;i<4;++i) fCalibration->GetYaxis()->SetBinLabel(i+1,corrections[i]);
    if(fVeto) {
      fMapPath=card->GetStr("JetVetoMapFile").Data(); fMapKey=card->GetStr("JetVetoMapHistogram").Data();
      std::unique_ptr<TFile> file(TFile::Open(fMapPath.c_str(),"READ"));
      auto* map=file && !file->IsZombie()?dynamic_cast<TH2*>(file->Get(fMapKey.c_str())):nullptr;
      if(!map) throw std::runtime_error("Run3JetQA cannot snapshot the configured veto histogram");
      fMap.reset(dynamic_cast<TH2*>(map->Clone("configuredJetVetoMap")));
      fMap->SetDirectory(nullptr);
    }
  }

  void Fill(double pt,double raw,double trackMax,double centrality,double weight,
            bool trackPass,bool idPass,bool vetoPass,const std::vector<double>& stages) {
    if(!std::isfinite(weight) || !std::isfinite(pt) || !std::isfinite(centrality))
      throw std::invalid_argument("Nonfinite Run-3 jet QA input");
    const bool selected[]={true,trackPass,idPass,trackPass&&idPass,
                           fVeto&&idPass&&vetoPass,fVeto&&trackPass&&idPass&&vetoPass};
    for(int i=0;i<6;++i) if(selected[i]) { fSelection->Fill(pt,i,centrality,weight); fCounts->Fill(pt,i,centrality); }
    if(raw>0 && std::isfinite(trackMax) && trackMax>=0) fTrackFraction->Fill(trackMax/raw,raw,centrality,weight);
    if(!(trackPass && idPass && vetoPass)) return;
    const size_t expected=fOO?2:(fRequestedJEC?(fData?(fPpL1?4:3):(fPpL1?3:2)):(fData?3:2));
    if(stages.size()!=expected) throw std::logic_error("Unexpected calibration stage count in Run3JetQA");
    fCalibration->Fill(stages[0],0.0,centrality,weight);
    size_t relative=1;
    if(!fOO && fRequestedJEC && fPpL1) {
      fCalibration->Fill(stages[1],1,centrality,weight);
      relative=2;
    }
    fCalibration->Fill(stages[relative],2,centrality,weight);
    if(!fOO && fData) fCalibration->Fill(stages.back(),3,centrality,weight);
  }

  void Write() const {
    TDirectory::TContext context;
    auto* directory=gDirectory->GetDirectory("run3JetQA");
    if(!directory) directory=gDirectory->mkdir("run3JetQA");
    directory->cd();
    fSelection->Write(); fCounts->Write(); fCalibration->Write(); fTrackFraction->Write();
    TNamed("JetIDTable",fTable.c_str()).Write();
    const std::string jecState = !fRequestedJEC ? "LEGACY calibration: Run-3 jet-pT JEC not active" :
      (fOO ? "Run-3 jet-pT JEC: OO L2Relative" :
       (fPpL1 ? "Run-3 jet-pT JEC: pp L1FastJet->L2Relative V2->[data residual]" :
                "Run-3 jet-pT JEC: pp L2Relative V2->[data residual]; L1FastJet disabled by card"));
    TNamed("JECState",jecState.c_str()).Write();
    TNamed("PpL1FastJetState",fPpL1?"ENABLED":"DISABLED").Write();
    TNamed("SelectionDefinition","Reco analysis-jet modes only; after trigger and standard event cuts; common analysis-axis eta and final-pT/matching cuts; before trackMax/ID/event-veto rejection. No selection means no jet-quality selection, not no event/kinematic selection. All curves share the same input exposure; do not unit-normalize separately.").Write();
    TNamed("CalibrationDefinition","Fixed jets passing trackMax/rawPt AND Jet ID and configured event selection, with final calibrated pT in the analysis range; identical jets at every active JEC level. No inactive or absent level is fabricated.").Write();
    TNamed("TrackMaxDefinition","MinMaxTrackPtFraction < trackMax/rawpt < MaxMaxTrackPtFraction; NOT trackMax/correctedPt. Raw and final-pT flow bins retained.").Write();
    TNamed("StatisticalDefinition","Sumw2 records individual jet weights. Jets share events; the plotting macro does not claim independent/binomial ratio errors. Use event-level resampling for physics uncertainties.").Write();
    TNamed("VetoState",fVeto?"ENABLED: whole-event mask; trackMax is NOT an eligibility condition for the veto scan":"DISABLED: veto stages are empty, not equivalent to zero loss").Write();
    if(fMap) {
      fMap->Write();
      TNamed("VetoMapFile",fMapPath.c_str()).Write(); TNamed("VetoMapKey",fMapKey.c_str()).Write();
      TNamed("VetoMapSnapshot","Payload bin contents, NOT jet occupancy; do not sum this object with hadd. Plot a single-job snapshot or the original pinned payload.").Write();
    }
  }
private:
  bool fRequestedJEC,fPpL1,fOO,fData,fVeto;
  std::string fTable,fMapPath,fMapKey;
  std::unique_ptr<TH3D> fSelection,fCounts,fCalibration,fTrackFraction;
  std::unique_ptr<TH2> fMap;
};
#endif

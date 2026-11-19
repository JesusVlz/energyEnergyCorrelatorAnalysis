#ifndef EEC_RUN3_JET_QA_H
#define EEC_RUN3_JET_QA_H
#include "ConfigurationCard.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3D.h>
#include <TNamed.h>
#include <TMath.h>
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
      fJetID(card->HasKey("Run3JetID") && card->Get("Run3JetID")==1),
      fTable(card->HasKey("Run3JetIDTable")?card->GetStr("Run3JetIDTable").Data():"MiniPOG20261108") {
    TDirectory::TContext context(nullptr);
    fSelection.reset(new TH3D("selectionPtWeighted",";Final calibrated jet p_{T} (GeV);Selection;Centrality (%)",500,0,1000,7,-.5,6.5,101,-1,100));
    fCounts.reset(new TH3D("selectionPtCounts",";Final calibrated jet p_{T} (GeV);Selection;Centrality (%)",500,0,1000,7,-.5,6.5,101,-1,100));
    fSelectionEta.reset(new TH3D("selectionEtaWeighted",";Analysis jet #eta;Selection;Centrality (%)",80,-2,2,7,-.5,6.5,101,-1,100));
    fSelectionEtaCounts.reset(new TH3D("selectionEtaCounts",";Analysis jet #eta;Selection;Centrality (%)",80,-2,2,7,-.5,6.5,101,-1,100));
    fSelectionPhi.reset(new TH3D("selectionPhiWeighted",";Analysis jet #phi;Selection;Centrality (%)",72,-TMath::Pi(),TMath::Pi(),7,-.5,6.5,101,-1,100));
    fSelectionPhiCounts.reset(new TH3D("selectionPhiCounts",";Analysis jet #phi;Selection;Centrality (%)",72,-TMath::Pi(),TMath::Pi(),7,-.5,6.5,101,-1,100));
    const std::vector<TH3D*> selectionHistograms={fSelection.get(),fCounts.get(),fSelectionEta.get(),fSelectionEtaCounts.get(),fSelectionPhi.get(),fSelectionPhiCounts.get()};
    fCalibration.reset(new TH3D("calibrationPtWeighted",";Jet p_{T} at this calibration level (GeV);Actual JEC level;Centrality (%)",500,0,1000,4,-.5,3.5,101,-1,100));
    fTrackFraction.reset(new TH3D("trackMaxOverRawPt",";trackMax/rawPt;Raw jet p_{T} (GeV);Centrality (%)",120,0,1.2,100,0,1000,101,-1,100));
    for(TH3D* h:selectionHistograms) { h->SetDirectory(nullptr); h->Sumw2(); }
    for(TH3D* h:{fCalibration.get(),fTrackFraction.get()}) { h->SetDirectory(nullptr); h->Sumw2(); }
    const char* selections[]={"No jet-quality selection","trackMax/rawPt only","Jet ID only","trackMax/rawPt AND Jet ID","Jet ID AND event veto","trackMax/rawPt AND Jet ID AND event veto","event veto only"};
    for(TH3D* h:selectionHistograms)
      for(int i=0;i<7;++i) h->GetYaxis()->SetBinLabel(i+1,selections[i]);
    if(!fJetID) {
      for(TH3D* h:selectionHistograms) {
        for(int stage:{2,3,4}) h->GetYaxis()->SetBinLabel(stage+1,"INACTIVE: Jet ID disabled");
        h->GetYaxis()->SetBinLabel(6,"trackMax/rawPt AND event veto (Jet ID disabled)");
      }
    }
    const char* corrections[]={"rawpt","after L1FastJet","after L2Relative","after L2Residual (data only)"};
    for(int i=0;i<4;++i) fCalibration->GetYaxis()->SetBinLabel(i+1,corrections[i]);
    std::vector<double> ptEdges,etaEdges,phiEdges,activityEdges;
    for(int i=0;i<card->GetNBin("JetPtBinEdgesEEC")+1;++i) ptEdges.push_back(card->Get("JetPtBinEdgesEEC",i));
    for(int i=0;i<=80;++i) etaEdges.push_back(-2.+.05*i);
    for(int i=0;i<=72;++i) phiEdges.push_back(-TMath::Pi()+2*TMath::Pi()*i/72.);
    for(int i=0;i<=101;++i) activityEdges.push_back(-1.+i);
    fEta.reset(new TH3D("selectedJetEtaWeighted",";WTA jet #eta;Final calibrated jet p_{T} (GeV);Centrality (%)",80,etaEdges.data(),ptEdges.size()-1,ptEdges.data(),101,activityEdges.data()));
    fPhi.reset(new TH3D("selectedJetPhiWeighted",";WTA jet #phi;Final calibrated jet p_{T} (GeV);Centrality (%)",72,phiEdges.data(),ptEdges.size()-1,ptEdges.data(),101,activityEdges.data()));
    const double activity[]={-1,0,30,50,90,100};
    std::vector<double> detectorEta;for(int i=0;i<=104;++i)detectorEta.push_back(-5.2+.1*i);
    fVetoBefore.reset(new TH3D("signalVetoEtaPhiActivityBefore",";Detector #eta;Detector #phi;Centrality (%)",104,detectorEta.data(),72,phiEdges.data(),5,activity));
    fVetoAfter.reset(new TH3D("signalVetoEtaPhiActivityAfter",";Detector #eta;Detector #phi;Centrality (%)",104,detectorEta.data(),72,phiEdges.data(),5,activity));
    const double states[]={-.5,.5,1.5,2.5};
    fVetoEvents.reset(new TH2D("signalVetoEventsActivity",";Event state;Centrality (%)",3,states,5,activity));
    for(int i=0;i<3;++i) fVetoEvents->GetXaxis()->SetBinLabel(i+1,i==0?"eligible":(i==1?"kept":"rejected"));
    fVetoEvents->SetDirectory(nullptr);fVetoEvents->Sumw2();
    for(TH3D* h:{fEta.get(),fPhi.get(),fVetoBefore.get(),fVetoAfter.get()}){h->SetDirectory(nullptr);h->Sumw2();}
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
            bool trackPass,bool idPass,bool vetoPass,const std::vector<double>& stages,double eta,double phi) {
    if(!std::isfinite(weight) || !std::isfinite(pt) || !std::isfinite(centrality))
      throw std::invalid_argument("Nonfinite Run-3 jet QA input");
    if(!std::isfinite(eta) || !std::isfinite(phi)) throw std::invalid_argument("Nonfinite jet axis in QA");
    double wrappedPhi=std::remainder(phi,2*TMath::Pi());if(wrappedPhi>=TMath::Pi())wrappedPhi-=2*TMath::Pi();
    const bool selected[]={true,trackPass,fJetID&&idPass,fJetID&&trackPass&&idPass,
                           fVeto&&fJetID&&idPass&&vetoPass,fVeto&&trackPass&&(!fJetID||idPass)&&vetoPass,fVeto&&vetoPass};
    for(int i=0;i<7;++i) if(selected[i]) {
      fSelection->Fill(pt,i,centrality,weight); fCounts->Fill(pt,i,centrality);
      fSelectionEta->Fill(eta,i,centrality,weight); fSelectionEtaCounts->Fill(eta,i,centrality);
      fSelectionPhi->Fill(wrappedPhi,i,centrality,weight); fSelectionPhiCounts->Fill(wrappedPhi,i,centrality);
    }
    if(raw>0 && std::isfinite(trackMax) && trackMax>=0) fTrackFraction->Fill(trackMax/raw,raw,centrality,weight);
    if(!(trackPass && (!fJetID || idPass) && vetoPass)) return;
    fEta->Fill(eta,pt,centrality,weight);fPhi->Fill(wrappedPhi,pt,centrality,weight);
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

  void FillVeto(double centrality,bool kept,const std::vector<std::pair<double,double>>& eligible) {
    if(!fVeto)throw std::logic_error("Veto QA called when disabled");
    fVetoEvents->Fill(0.,centrality);fVetoEvents->Fill(kept?1.:2.,centrality);
    for(const auto& jet:eligible){fVetoBefore->Fill(jet.first,jet.second,centrality);if(kept)fVetoAfter->Fill(jet.first,jet.second,centrality);}
  }

  void Write() const {
    TDirectory::TContext context;
    auto* directory=gDirectory->GetDirectory("run3JetQA");
    if(!directory) directory=gDirectory->mkdir("run3JetQA");
    directory->cd();
    fSelection->Write(); fCounts->Write(); fCalibration->Write(); fTrackFraction->Write();
    fSelectionEta->Write(); fSelectionEtaCounts->Write(); fSelectionPhi->Write(); fSelectionPhiCounts->Write();
    TNamed("SelectionKinematicsDefinition","selection{Pt,Eta,Phi}{Weighted,Counts}: identical seven selection stages, event exposure and card-defined kinematic/matching cuts. Eta/phi use the analysis axis (JetAxis=0: E-scheme, JetAxis=1: WTA); phi is wrapped to [-pi,pi). Eta/phi spectra integrate over the full selected jet-pT range; no additional pT-bin selection is stored. Weighted yields use the same event weight as selectionPtWeighted; Counts are unweighted. Coordinate underflow/overflow is retained.").Write();
    TNamed("VetoEventDenominator","All signal events reaching the veto after trigger/standard event cuts, including zero eligible veto jets; unweighted, exactly once per event. Tested = kept + rejected in each physical-activity bin.").Write();
    fEta->Write();fPhi->Write();fVetoBefore->Write();fVetoAfter->Write();fVetoEvents->Write();
    TNamed("KinematicsDefinition","Final-selected reconstructed jets; analysis WTA eta/phi, final calibrated pT; weighted yields at common event exposure. Veto occupancy uses detector coordinates and unweighted eligible jets/events before veto rejection.").Write();
    TNamed("JetIDTable",fTable.c_str()).Write();
    TNamed("JetIDState",fJetID?"ENABLED":"DISABLED by explicit pp MC card; also omitted from detector-veto eligibility; ID-only stages are inactive").Write();
    const std::string jecState = !fRequestedJEC ? "LEGACY calibration: Run-3 jet-pT JEC not active" :
      (fOO ? "Run-3 jet-pT JEC: OO L2Relative" :
       (fPpL1 ? "Run-3 jet-pT JEC: pp L1FastJet->L2Relative V2->[data residual]" :
                "Run-3 jet-pT JEC: pp L2Relative V2->[data residual]; L1FastJet disabled by card"));
    TNamed("JECState",jecState.c_str()).Write();
    TNamed("PpL1FastJetState",fPpL1?"ENABLED":"DISABLED").Write();
    TNamed("SelectionDefinition","Reco analysis-jet modes only; after trigger and standard event cuts; common analysis-axis eta and final-pT/matching cuts; before trackMax/ID/event-veto rejection. No selection means no jet-quality selection, not no event/kinematic selection. All curves share the same input exposure; do not unit-normalize separately.").Write();
    TNamed("CalibrationDefinition","Fixed jets passing trackMax/rawPt AND enabled Jet ID and configured event selection, with final calibrated pT in the analysis range; identical jets at every active JEC level. No inactive or absent level is fabricated.").Write();
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
  bool fRequestedJEC,fPpL1,fOO,fData,fVeto,fJetID;
  std::string fTable,fMapPath,fMapKey;
  std::unique_ptr<TH3D> fSelection,fCounts,fCalibration,fTrackFraction,fEta,fPhi,fVetoBefore,fVetoAfter;
  std::unique_ptr<TH3D> fSelectionEta,fSelectionEtaCounts,fSelectionPhi,fSelectionPhiCounts;
  std::unique_ptr<TH2D> fVetoEvents;
  std::unique_ptr<TH2> fMap;
};
#endif

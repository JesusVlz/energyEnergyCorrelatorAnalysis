// End-to-end C++/ROOT verification for MixingLookupTable, single-handle ForestReader binding,
// grouped single-open donor track prefetch, Run-3 jet veto companion alignment, and OOClosurePilot.
// Compiled with -fno-access-control in validation so deterministic TRandom3 seeds can be set.

#include "ConfigurationCard.h"
#include "EECAnalyzer.h"
#include "MixingLookupTable.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TROOT.h>
#include <TTree.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void Require(bool ok, const std::string& msg) {
  if (!ok) throw std::runtime_error(msg);
}

void SetNumeric(ConfigurationCard& card, const char* key, double value) {
  std::ostringstream line;
  line << key << " " << value;
  if (!card.HasKey(key)) {
    card.ReadInputLine(line.str().c_str());
  } else {
    auto* vec = card.GetVector(key);
    (*vec)[vec->GetLwb()] = value;
  }
}

void SetString(ConfigurationCard& card, const char* key, const std::string& value) {
  if (!card.HasKey(key)) {
    card.ReadInputLine((std::string(key) + " " + value).c_str());
  } else {
    card.fValueString[card.GetTVectorIndex(key)] = value.c_str();
  }
}

void WriteSyntheticVetoMap(const std::string& path) {
  TDirectory::TContext ctx;
  TFile file(path.c_str(), "RECREATE");
  TH2D map("veto_map", ";eta;phi", 104, -5.2, 5.2, 72, -TMath::Pi(), TMath::Pi());
  // Veto hot zone at bin center (eta=2.85, phi=0.25)
  map.SetBinContent(map.FindBin(2.85, 0.25), 100.0);
  map.Write();
  file.Close();
}

void WriteOODataForest(const std::string& path,
                       int nEvents,
                       ULong64_t firstEvt,
                       bool includeSignalJets,
                       int fileSalt) {
  TDirectory::TContext ctx;
  TFile file(path.c_str(), "RECREATE");

  file.mkdir("hiEvtAnalyzer")->cd();
  TTree hiTree("HiTree", "OO events");
  UInt_t run = 390000 + fileSalt;
  UInt_t lumi = 10 + fileSalt;
  ULong64_t evt = firstEvt;
  Float_t vz = 0.0f;
  Int_t hiBin = 20;
  Float_t hf = 45.0f, hfPlus = 22.5f, hfMinus = 22.5f;
  hiTree.Branch("run", &run, "run/i");
  hiTree.Branch("lumi", &lumi, "lumi/i");
  hiTree.Branch("evt", &evt, "evt/l");
  hiTree.Branch("vz", &vz, "vz/F");
  hiTree.Branch("hiBin", &hiBin, "hiBin/I");
  hiTree.Branch("hiHF_pf", &hf, "hiHF_pf/F");
  hiTree.Branch("hiHFPlus_pf", &hfPlus, "hiHFPlus_pf/F");
  hiTree.Branch("hiHFMinus_pf", &hfMinus, "hiHFMinus_pf/F");

  file.cd();
  file.mkdir("hltanalysis")->cd();
  TTree hltTree("HltTree", "HLT");
  Int_t mbTrig = 1;
  hltTree.Branch("HLT_MinimumBiasHF_OR_BptxAND_v1", &mbTrig, "HLT_MinimumBiasHF_OR_BptxAND_v1/I");

  file.cd();
  file.mkdir("skimanalysis")->cd();
  TTree skimTree("HltTree", "Skim");
  Int_t pv = 1, hfCoinc = 1, pileup = 1;
  skimTree.Branch("pprimaryVertexFilter", &pv, "pprimaryVertexFilter/I");
  skimTree.Branch("pOOprimaryVertexFilter", &pv, "pOOprimaryVertexFilter/I");
  skimTree.Branch("OOpfCoincFilterPF2Th4", &hfCoinc, "OOpfCoincFilterPF2Th4/I");
  skimTree.Branch("OOphfCoincFilterPF2Th4", &hfCoinc, "OOphfCoincFilterPF2Th4/I");
  skimTree.Branch("pileupVertexFilter", &pileup, "pileupVertexFilter/I");

  file.cd();
  file.mkdir("akCs4PFJetAnalyzer")->cd();
  TTree jetTree("t", "Jets");
  Int_t nref = 0;
  Float_t jtpt[8] = {}, rawpt[8] = {}, jteta[8] = {}, jtphi[8] = {};
  Float_t wtaEta[8] = {}, wtaPhi[8] = {}, trackMax[8] = {};
  Float_t nhf[8] = {}, nef[8] = {}, chf[8] = {}, muf[8] = {}, cef[8] = {};
  Int_t chm[8] = {}, cem[8] = {}, nhm[8] = {}, nem[8] = {}, mum[8] = {};
  jetTree.Branch("nref", &nref, "nref/I");
  jetTree.Branch("jtpt", jtpt, "jtpt[nref]/F");
  jetTree.Branch("rawpt", rawpt, "rawpt[nref]/F");
  jetTree.Branch("jteta", jteta, "jteta[nref]/F");
  jetTree.Branch("jtphi", jtphi, "jtphi[nref]/F");
  jetTree.Branch("WTAeta", wtaEta, "WTAeta[nref]/F");
  jetTree.Branch("WTAphi", wtaPhi, "WTAphi[nref]/F");
  jetTree.Branch("trackMax", trackMax, "trackMax[nref]/F");
  jetTree.Branch("jtPfNHF", nhf, "jtPfNHF[nref]/F");
  jetTree.Branch("jtPfNEF", nef, "jtPfNEF[nref]/F");
  jetTree.Branch("jtPfCHF", chf, "jtPfCHF[nref]/F");
  jetTree.Branch("jtPfMUF", muf, "jtPfMUF[nref]/F");
  jetTree.Branch("jtPfCEF", cef, "jtPfCEF[nref]/F");
  jetTree.Branch("jtPfCHM", chm, "jtPfCHM[nref]/I");
  jetTree.Branch("jtPfCEM", cem, "jtPfCEM[nref]/I");
  jetTree.Branch("jtPfNHM", nhm, "jtPfNHM[nref]/I");
  jetTree.Branch("jtPfNEM", nem, "jtPfNEM[nref]/I");
  jetTree.Branch("jtPfMUM", mum, "jtPfMUM[nref]/I");

  file.cd();
  file.mkdir("PbPbTracks")->cd();
  TTree trkTree("trackTree", "Tracks");
  Int_t nTrk = 0;
  std::vector<float> trkPt, trkPtErr, trkEta, trkPhi, dz, dzErr, dxy, dxyErr;
  std::vector<bool> highPurity;
  std::vector<char> trkCharge;
  trkTree.Branch("nTrk", &nTrk, "nTrk/I");
  trkTree.Branch("trkPt", &trkPt);
  trkTree.Branch("trkPtError", &trkPtErr);
  trkTree.Branch("trkEta", &trkEta);
  trkTree.Branch("trkPhi", &trkPhi);
  trkTree.Branch("highPurity", &highPurity);
  trkTree.Branch("trkDzFirstVtx", &dz);
  trkTree.Branch("trkDzErrFirstVtx", &dzErr);
  trkTree.Branch("trkDxyFirstVtx", &dxy);
  trkTree.Branch("trkDxyErrFirstVtx", &dxyErr);
  trkTree.Branch("trkCharge", &trkCharge);

  const int hiBins[4] = {20, 20, 40, 40};
  const float vzVals[4] = {-1.2f, -1.0f, 1.1f, 1.3f};

  for (int i = 0; i < nEvents; ++i) {
    evt = firstEvt + static_cast<ULong64_t>(i);
    vz = vzVals[i % 4] + 0.02f * static_cast<float>(fileSalt % 5);
    hiBin = hiBins[i % 4];
    hf = 40.0f + 5.0f * static_cast<float>(i % 4);
    hfPlus = 0.5f * hf;
    hfMinus = 0.5f * hf;
    // Include one rejected event per file (if nEvents >= 6) to test event-cut & MB filtering
    mbTrig = (i == 5) ? 0 : 1;
    pv = (i == 6) ? 0 : 1;
    hfCoinc = 1;
    pileup = 1;

    if (includeSignalJets) {
      nref = 1;
      const float basePhi = -1.2f + 0.15f * static_cast<float>(i % 4);
      const float baseEta = 0.25f * static_cast<float>((i % 3) - 1);
      jtpt[0] = 115.0f + 5.0f * static_cast<float>(i);
      rawpt[0] = 110.0f + 5.0f * static_cast<float>(i);
      jteta[0] = baseEta;
      jtphi[0] = basePhi;
      wtaEta[0] = baseEta;
      wtaPhi[0] = basePhi;
      trackMax[0] = 18.0f;
      nhf[0] = 0.1f; nef[0] = 0.1f; chf[0] = 0.6f; muf[0] = 0.0f; cef[0] = 0.1f;
      chm[0] = 6; cem[0] = 0; nhm[0] = 2; nem[0] = 2; mum[0] = 0;
    } else {
      nref = 0;
    }

    trkPt.clear(); trkPtErr.clear(); trkEta.clear(); trkPhi.clear();
    dz.clear(); dzErr.clear(); dxy.clear(); dxyErr.clear();
    highPurity.clear(); trkCharge.clear();

    // Populate tracks near every jet candidate axis so both signal and mixed cones get pairs
    for (int t = 0; t < 4; ++t) {
      const float basePhi = -1.2f + 0.15f * static_cast<float>(t);
      const float baseEta = 0.25f * static_cast<float>((t % 3) - 1);
      for (int sub = 0; sub < 2; ++sub) {
        const float pt = 2.5f + 0.4f * static_cast<float>(sub) + 0.05f * static_cast<float>((fileSalt + i) % 5);
        trkPt.push_back(pt);
        trkPtErr.push_back(0.02f * pt);
        trkEta.push_back(baseEta + (sub == 0 ? 0.04f : -0.05f));
        trkPhi.push_back(basePhi + (sub == 0 ? 0.03f : -0.06f));
        dz.push_back(0.01f);
        dzErr.push_back(0.05f);
        dxy.push_back(0.01f);
        dxyErr.push_back(0.05f);
        highPurity.push_back(true);
        trkCharge.push_back(sub == 0 ? 1 : -1);
      }
    }
    nTrk = static_cast<Int_t>(trkPt.size());

    hiTree.Fill();
    hltTree.Fill();
    skimTree.Fill();
    jetTree.Fill();
    trkTree.Fill();
  }

  file.cd("hiEvtAnalyzer"); hiTree.Write();
  file.cd("hltanalysis"); hltTree.Write();
  file.cd("skimanalysis"); skimTree.Write();
  file.cd("akCs4PFJetAnalyzer"); jetTree.Write();
  file.cd("PbPbTracks"); trkTree.Write();
  file.Close();
}

void WriteOOMCForest(const std::string& path,
                     int nEvents,
                     ULong64_t firstEvt,
                     bool pureHijingDonor,
                     int fileSalt) {
  TDirectory::TContext ctx;
  TFile file(path.c_str(), "RECREATE");

  file.mkdir("hiEvtAnalyzer")->cd();
  TTree hiTree("HiTree", "OOMC events");
  UInt_t run = 1;
  UInt_t lumi = 100 + fileSalt;
  ULong64_t evt = firstEvt;
  Float_t vz = 0.0f, hf = 20.0f, hfPlus = 10.0f, hfMinus = 10.0f, pthat = 80.0f, weight = 1.0f;
  Int_t hiBin = 20;
  hiTree.Branch("run", &run, "run/i");
  hiTree.Branch("lumi", &lumi, "lumi/i");
  hiTree.Branch("evt", &evt, "evt/l");
  hiTree.Branch("vz", &vz, "vz/F");
  hiTree.Branch("hiBin", &hiBin, "hiBin/I");
  hiTree.Branch("hiHF_pf", &hf, "hiHF_pf/F");
  hiTree.Branch("hiHFPlus_pf", &hfPlus, "hiHFPlus_pf/F");
  hiTree.Branch("hiHFMinus_pf", &hfMinus, "hiHFMinus_pf/F");
  hiTree.Branch("pthat", &pthat, "pthat/F");
  hiTree.Branch("weight", &weight, "weight/F");

  file.cd();
  file.mkdir("hltanalysis")->cd();
  TTree hltTree("HltTree", "HLT");
  Int_t mbTrig = 1;
  hltTree.Branch("HLT_MinimumBiasHF_OR_BptxAND_v1", &mbTrig, "HLT_MinimumBiasHF_OR_BptxAND_v1/I");

  file.cd();
  file.mkdir("skimanalysis")->cd();
  TTree skimTree("HltTree", "Skim");
  Int_t pv = 1, hfCoinc = 1, pileup = 1;
  skimTree.Branch("pprimaryVertexFilter", &pv, "pprimaryVertexFilter/I");
  skimTree.Branch("pOOprimaryVertexFilter", &pv, "pOOprimaryVertexFilter/I");
  skimTree.Branch("OOpfCoincFilterPF2Th4", &hfCoinc, "OOpfCoincFilterPF2Th4/I");
  skimTree.Branch("OOphfCoincFilterPF2Th4", &hfCoinc, "OOphfCoincFilterPF2Th4/I");
  skimTree.Branch("pileupVertexFilter", &pileup, "pileupVertexFilter/I");

  file.cd();
  file.mkdir("akCs4PFJetAnalyzer")->cd();
  TTree jetTree("t", "MC Jets");
  Int_t nref = 1, ngen = 1;
  Float_t jtpt[4] = {}, rawpt[4] = {}, jteta[4] = {}, jtphi[4] = {};
  Float_t wtaEta[4] = {}, wtaPhi[4] = {}, trackMax[4] = {};
  Float_t refpt[4] = {}, refeta[4] = {}, refphi[4] = {};
  Int_t refparton[4] = {1};
  Float_t nhf[4] = {}, nef[4] = {}, chf[4] = {}, muf[4] = {}, cef[4] = {};
  Int_t chm[4] = {}, cem[4] = {}, nhm[4] = {}, nem[4] = {}, mum[4] = {};
  Float_t genpt[4] = {}, geneta[4] = {}, genphi[4] = {}, wtaGenEta[4] = {}, wtaGenPhi[4] = {};
  Int_t gensubid[4] = {0};

  jetTree.Branch("nref", &nref, "nref/I");
  jetTree.Branch("jtpt", jtpt, "jtpt[nref]/F");
  jetTree.Branch("rawpt", rawpt, "rawpt[nref]/F");
  jetTree.Branch("jteta", jteta, "jteta[nref]/F");
  jetTree.Branch("jtphi", jtphi, "jtphi[nref]/F");
  jetTree.Branch("WTAeta", wtaEta, "WTAeta[nref]/F");
  jetTree.Branch("WTAphi", wtaPhi, "WTAphi[nref]/F");
  jetTree.Branch("trackMax", trackMax, "trackMax[nref]/F");
  jetTree.Branch("refpt", refpt, "refpt[nref]/F");
  jetTree.Branch("refeta", refeta, "refeta[nref]/F");
  jetTree.Branch("refphi", refphi, "refphi[nref]/F");
  jetTree.Branch("refparton_flavorForB", refparton, "refparton_flavorForB[nref]/I");
  jetTree.Branch("jtPfNHF", nhf, "jtPfNHF[nref]/F");
  jetTree.Branch("jtPfNEF", nef, "jtPfNEF[nref]/F");
  jetTree.Branch("jtPfCHF", chf, "jtPfCHF[nref]/F");
  jetTree.Branch("jtPfMUF", muf, "jtPfMUF[nref]/F");
  jetTree.Branch("jtPfCEF", cef, "jtPfCEF[nref]/F");
  jetTree.Branch("jtPfCHM", chm, "jtPfCHM[nref]/I");
  jetTree.Branch("jtPfCEM", cem, "jtPfCEM[nref]/I");
  jetTree.Branch("jtPfNHM", nhm, "jtPfNHM[nref]/I");
  jetTree.Branch("jtPfNEM", nem, "jtPfNEM[nref]/I");
  jetTree.Branch("jtPfMUM", mum, "jtPfMUM[nref]/I");
  jetTree.Branch("ngen", &ngen, "ngen/I");
  jetTree.Branch("genpt", genpt, "genpt[ngen]/F");
  jetTree.Branch("geneta", geneta, "geneta[ngen]/F");
  jetTree.Branch("genphi", genphi, "genphi[ngen]/F");
  jetTree.Branch("WTAgeneta", wtaGenEta, "WTAgeneta[ngen]/F");
  jetTree.Branch("WTAgenphi", wtaGenPhi, "WTAgenphi[ngen]/F");
  jetTree.Branch("gensubid", gensubid, "gensubid[ngen]/I");

  file.cd();
  file.mkdir("PbPbTracks")->cd();
  TTree trkTree("trackTree", "MC Reco Tracks");
  Int_t nTrk = 0;
  std::vector<float> trkPt, trkPtErr, trkEta, trkPhi, dz, dzErr, dxy, dxyErr;
  std::vector<bool> highPurity;
  std::vector<char> trkCharge;
  trkTree.Branch("nTrk", &nTrk, "nTrk/I");
  trkTree.Branch("trkPt", &trkPt);
  trkTree.Branch("trkPtError", &trkPtErr);
  trkTree.Branch("trkEta", &trkEta);
  trkTree.Branch("trkPhi", &trkPhi);
  trkTree.Branch("highPurity", &highPurity);
  trkTree.Branch("trkDzFirstVtx", &dz);
  trkTree.Branch("trkDzErrFirstVtx", &dzErr);
  trkTree.Branch("trkDxyFirstVtx", &dxy);
  trkTree.Branch("trkDxyErrFirstVtx", &dxyErr);
  trkTree.Branch("trkCharge", &trkCharge);

  file.cd();
  file.mkdir("HiGenParticleAna")->cd();
  TTree genTree("hi", "Gen particles");
  std::vector<float> pt, eta, phi;
  std::vector<int> chg, sube, sta;
  genTree.Branch("pt", &pt);
  genTree.Branch("eta", &eta);
  genTree.Branch("phi", &phi);
  genTree.Branch("chg", &chg);
  genTree.Branch("sube", &sube);
  genTree.Branch("sta", &sta);

  for (int i = 0; i < nEvents; ++i) {
    evt = firstEvt + static_cast<ULong64_t>(i);
    vz = (i % 2 == 0) ? 0.5f : -0.5f;
    hiBin = 20;
    hf = 25.0f + static_cast<float>(i);
    hfPlus = 0.5f * hf;
    hfMinus = 0.5f * hf;

    // Put a vetoed reco jet at (eta=2.85, phi=0.25) in one donor event per file to verify Pass 1b veto companion
    const bool vetoThisEvent = pureHijingDonor && (i == 2);
    nref = 1;
    jtpt[0] = vetoThisEvent ? 30.0f : 120.0f;
    rawpt[0] = vetoThisEvent ? 28.0f : 115.0f;
    jteta[0] = vetoThisEvent ? 2.85f : 0.2f;
    jtphi[0] = vetoThisEvent ? 0.25f : -1.0f;
    wtaEta[0] = jteta[0];
    wtaPhi[0] = jtphi[0];
    trackMax[0] = 15.0f;
    refpt[0] = 120.0f; refeta[0] = 0.2f; refphi[0] = -1.0f; refparton[0] = 1;
    nhf[0] = 0.1f; nef[0] = 0.1f; chf[0] = 0.6f; muf[0] = 0.0f; cef[0] = 0.1f;
    chm[0] = 5; cem[0] = 0; nhm[0] = 2; nem[0] = 2; mum[0] = 0;

    ngen = pureHijingDonor ? 0 : 1;
    genpt[0] = 125.0f;
    geneta[0] = 0.2f;
    genphi[0] = -1.0f;
    wtaGenEta[0] = 0.2f;
    wtaGenPhi[0] = -1.0f;
    gensubid[0] = 0;

    trkPt.clear(); trkPtErr.clear(); trkEta.clear(); trkPhi.clear();
    dz.clear(); dzErr.clear(); dxy.clear(); dxyErr.clear();
    highPurity.clear(); trkCharge.clear();
    for (int sub = 0; sub < 2; ++sub) {
      const float rpt = 4.0f + static_cast<float>(sub);
      trkPt.push_back(rpt);
      trkPtErr.push_back(0.02f * rpt);
      trkEta.push_back(0.2f + (sub == 0 ? 0.03f : -0.04f));
      trkPhi.push_back(-1.0f + (sub == 0 ? 0.03f : -0.04f));
      dz.push_back(0.01f);
      dzErr.push_back(0.05f);
      dxy.push_back(0.01f);
      dxyErr.push_back(0.05f);
      highPurity.push_back(true);
      trkCharge.push_back(sub == 0 ? 1 : -1);
    }
    nTrk = static_cast<Int_t>(trkPt.size());

    pt.clear(); eta.clear(); phi.clear(); chg.clear(); sube.clear(); sta.clear();
    // 6 background particles (sube=1) with unique phi per (fileSalt, i) so background fingerprints are distinct
    for (int k = 0; k < 6; ++k) {
      pt.push_back(2.5f + 0.2f * static_cast<float>(k));
      eta.push_back(0.2f + (k < 2 ? 0.03f * static_cast<float>(k + 1) : 0.6f));
      phi.push_back(-1.0f + (k < 2 ? 0.04f * static_cast<float>(k + 1) : (0.3f * k + 0.017f * (10 * fileSalt + i))));
      chg.push_back((k % 2 == 0) ? 1 : -1);
      sube.push_back(1);
      sta.push_back(1);
    }
    if (!pureHijingDonor) {
      // 2 signal PYTHIA particles (sube=0) inside the jet cone
      for (int s = 0; s < 2; ++s) {
        pt.push_back(8.0f + 2.0f * static_cast<float>(s));
        eta.push_back(0.2f + (s == 0 ? -0.02f : 0.05f));
        phi.push_back(-1.0f + (s == 0 ? -0.03f : 0.04f));
        chg.push_back(s == 0 ? 1 : -1);
        sube.push_back(0);
        sta.push_back(1);
      }
    }

    hiTree.Fill();
    hltTree.Fill();
    skimTree.Fill();
    jetTree.Fill();
    trkTree.Fill();
    genTree.Fill();
  }

  file.cd("hiEvtAnalyzer"); hiTree.Write();
  file.cd("hltanalysis"); hltTree.Write();
  file.cd("skimanalysis"); skimTree.Write();
  file.cd("akCs4PFJetAnalyzer"); jetTree.Write();
  file.cd("PbPbTracks"); trkTree.Write();
  file.cd("HiGenParticleAna"); genTree.Write();
  file.Close();
}

void CompareTH1Exact(TH1* a, TH1* b, const std::string& label) {
  Require(a != nullptr && b != nullptr, "Null TH1 for " + label);
  Require(a->GetNbinsX() == b->GetNbinsX(), "Bin count X mismatch for " + label);
  for (int i = 0; i <= a->GetNbinsX() + 1; ++i) {
    if (std::abs(a->GetBinContent(i) - b->GetBinContent(i)) > 1e-12 ||
        std::abs(a->GetBinError(i) - b->GetBinError(i)) > 1e-12) {
      std::ostringstream err;
      err << label << " bin " << i << " differs: " << a->GetBinContent(i) << " vs " << b->GetBinContent(i);
      throw std::runtime_error(err.str());
    }
  }
}

void CompareTH2Exact(TH2* a, TH2* b, const std::string& label) {
  Require(a != nullptr && b != nullptr, "Null TH2 for " + label);
  Require(a->GetNbinsX() == b->GetNbinsX() && a->GetNbinsY() == b->GetNbinsY(), "Bin count 2D mismatch for " + label);
  for (int ix = 0; ix <= a->GetNbinsX() + 1; ++ix) {
    for (int iy = 0; iy <= a->GetNbinsY() + 1; ++iy) {
      if (std::abs(a->GetBinContent(ix, iy) - b->GetBinContent(ix, iy)) > 1e-12) {
        std::ostringstream err;
        err << label << " bin (" << ix << "," << iy << ") differs: "
            << a->GetBinContent(ix, iy) << " vs " << b->GetBinContent(ix, iy);
        throw std::runtime_error(err.str());
      }
    }
  }
}

double CompareTHnSparseExact(THnSparse* a, THnSparse* b, const std::string& label) {
  Require(a != nullptr && b != nullptr, "Null THnSparse for " + label);
  Require(a->GetNbins() == b->GetNbins(),
          label + " filled sparse bin count differs: " + std::to_string(a->GetNbins()) +
          " vs " + std::to_string(b->GetNbins()));
  const int nDim = a->GetNdimensions();
  std::vector<Int_t> coord(nDim, 0);
  double totalAbsWeight = 0.0;
  for (Long64_t idx = 0; idx < a->GetNbins(); ++idx) {
    const double va = a->GetBinContent(idx, coord.data());
    const double ea = a->GetBinError(idx);
    const Long64_t bIdx = b->GetBin(coord.data(), false);
    Require(bIdx >= 0, label + " missing sparse bin in second histogram");
    const double vb = b->GetBinContent(bIdx);
    const double eb = b->GetBinError(bIdx);
    if (std::abs(va - vb) > 1e-12 || std::abs(ea - eb) > 1e-12) {
      std::ostringstream err;
      err << label << " sparse bin content/error differs: (" << va << " +- " << ea
          << ") vs (" << vb << " +- " << eb << ")";
      throw std::runtime_error(err.str());
    }
    totalAbsWeight += std::abs(va);
  }
  return totalAbsWeight;
}

std::string ReadWholeFile(const std::string& path) {
  std::ifstream in(path);
  Require(static_cast<bool>(in), "Cannot read file " + path);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

void MakeFixtures(const std::filesystem::path& dir) {
  std::filesystem::create_directories(dir);
  const ULong64_t base64 = (ULong64_t(1) << 43);

  std::ofstream dataList(dir / "oo_donors.txt");
  for (int f = 0; f < 10; ++f) {
    std::ostringstream name;
    name << "donor_0" << f << ".root";
    const auto p = (dir / name.str()).string();
    WriteOODataForest(p, 8, base64 + static_cast<ULong64_t>(1000 * (f + 1)), false, f);
    dataList << p << "\n";
  }
  dataList.close();

  WriteOODataForest((dir / "signal_zero_jets.root").string(), 4, base64 + 500000ULL, false, 90);
  WriteOODataForest((dir / "signal_jets.root").string(), 5, base64 + 600000ULL, true, 91);

  WriteSyntheticVetoMap((dir / "veto_map.root").string());
  std::ofstream mcList(dir / "oomc_donors.txt");
  for (int f = 0; f < 4; ++f) {
    std::ostringstream name;
    name << "mc_donor_0" << f << ".root";
    const auto p = (dir / name.str()).string();
    WriteOOMCForest(p, 4, base64 + static_cast<ULong64_t>(20000 * (f + 1)), true, f);
    mcList << p << "\n";
  }
  mcList.close();
  WriteOOMCForest((dir / "mc_signal.root").string(), 4, base64 + 900000ULL, false, 50);
}

void VerifyEquivalence(const std::filesystem::path& dir) {
  // 1. CanonicalizeMixingPath checks (local vs XRootD PFN vs /store/)
  Require(MixingLookupTable::CanonicalizeMixingPath("root://eosuser.cern.ch//eos/user/j/jvelazqu/pool_0.root") ==
          MixingLookupTable::CanonicalizeMixingPath("/eos/user/j/jvelazqu/pool_0.root"),
          "CanonicalizeMixingPath failed on /eos/user XRootD rewrite");
  Require(MixingLookupTable::CanonicalizeMixingPath("root://cmsxrootd.fnal.gov//store/hidata/OIRun2025/f1.root") ==
          "/store/hidata/OIRun2025/f1.root",
          "CanonicalizeMixingPath failed on /store/ LFN extraction");

  const std::string dataDonorList = (dir / "oo_donors.txt").string();
  const std::string dataSignal = (dir / "signal_jets.root").string();
  const std::string dataLookup = dataDonorList + ".lookup.tsv";

  // Build OO Data lookup table once
  {
    ConfigurationCard buildCard("cardEECOOInclusiveR04.input");
    SetNumeric(buildCard, "BackgroundMethods", 6);
    SetString(buildCard, "MixingFileList", dataDonorList);
    SetNumeric(buildCard, "ApplyJetVeto", 0);
    SetNumeric(buildCard, "Run3JetID", 0);
    SetNumeric(buildCard, "Run3JetQA", 0);
    EECAnalyzer builder({}, &buildCard, 1, 1, 0);
    builder.BuildAndWriteMixingLookupTable(dataDonorList, dataLookup, 0);
  }
  Require(gROOT->GetListOfFiles()->GetSize() == 0, "Leaked TFile handles after BuildAndWriteMixingLookupTable");

  // 2. Unbinned OO Data (EventMatchForMixing 0) with covariance enabled: Lookup Table vs On-the-Fly Fallback
  for (int binned : {0, 1}) {
    ConfigurationCard cardLookup("cardEECOOInclusiveR04.input");
    ConfigurationCard cardFallback("cardEECOOInclusiveR04.input");
    for (ConfigurationCard* c : {&cardLookup, &cardFallback}) {
      SetNumeric(*c, "BackgroundMethods", 6);
      SetNumeric(*c, "EventMatchForMixing", 0);
      SetNumeric(*c, "BinnedMixing", binned);
      SetNumeric(*c, "MixingMaxOpenFiles", 0);
      SetNumeric(*c, "SkipCovarianceMatrix", 0);
      SetNumeric(*c, "ApplyJetVeto", 0);
      SetNumeric(*c, "Run3JetID", 0);
      SetNumeric(*c, "Run3JetQA", 0);
      SetString(*c, "MixingFileList", dataDonorList);
    }
    SetString(cardLookup, "MixingLookupTable", dataLookup);
    SetString(cardFallback, "MixingLookupTable", "none");

    auto anaLookup = std::make_unique<EECAnalyzer>(std::vector<TString>{dataSignal.c_str()}, &cardLookup, 1, 1, 0);
    auto anaFallback = std::make_unique<EECAnalyzer>(std::vector<TString>{dataSignal.c_str()}, &cardFallback, 1, 1, 0);
    anaLookup->fRng->SetSeed(20261120 + binned);
    anaFallback->fRng->SetSeed(20261120 + binned);

    anaLookup->RunAnalysis();
    anaFallback->RunAnalysis();

    const std::string tag = binned ? "OO_Data_Binned" : "OO_Data_Unbinned";
    const double eecSum = CompareTHnSparseExact(
      anaLookup->GetHistograms()->fhEnergyEnergyCorrelator,
      anaFallback->GetHistograms()->fhEnergyEnergyCorrelator,
      tag + "_EEC");
    Require(eecSum > 0.0, tag + " produced zero EEC weight");
    const double covSum = CompareTHnSparseExact(
      anaLookup->GetHistograms()->fhJetPtUnfoldingCovariance,
      anaFallback->GetHistograms()->fhJetPtUnfoldingCovariance,
      tag + "_Covariance");
    Require(covSum > 0.0, tag + " produced zero covariance weight");
    CompareTH2Exact(
      anaLookup->GetHistograms()->fhMixedEventCounts,
      anaFallback->GetHistograms()->fhMixedEventCounts,
      tag + "_MixedEventCounts");

    anaLookup.reset();
    anaFallback.reset();
    Require(gROOT->GetListOfFiles()->GetSize() == 0, "Leaked TFile handles after " + tag);
  }

  // 3. OO MC GenGen (EventMatchForMixing 2) + ApplyJetVeto 1 (companion HighForestReader Pass 1b)
  const std::string mcDonorList = (dir / "oomc_donors.txt").string();
  const std::string mcSignal = (dir / "mc_signal.root").string();
  const std::string mcLookup = (dir / "oomc_veto.lookup.tsv").string();
  const std::string vetoMap = (dir / "veto_map.root").string();
  {
    auto configureMcVeto = [&](ConfigurationCard& c) {
      SetNumeric(c, "BackgroundMethods", 6);
      SetNumeric(c, "McCorrelationType", 3); // kGenGen
      SetNumeric(c, "EventMatchForMixing", 2); // kMatchMultiplicity
      SetNumeric(c, "BinnedMixing", 0);
      SetNumeric(c, "MixingMaxOpenFiles", 0);
      SetNumeric(c, "Run3JetID", 1);
      SetNumeric(c, "ApplyJetVeto", 1);
      SetString(c, "JetVetoMapFile", vetoMap);
      SetString(c, "JetVetoMapHistogram", "veto_map");
      SetString(c, "JetVetoPtConvention", "stored");
      SetNumeric(c, "JetVetoInputMinPt", 15);
      SetNumeric(c, "Run3JetQA", 0);
      SetString(c, "MixingFileList", mcDonorList);
    };

    ConfigurationCard buildCard("cardEECOOMCInclusiveR04.input");
    configureMcVeto(buildCard);
    {
      EECAnalyzer builder({}, &buildCard, 1, 1, 0);
      builder.BuildAndWriteMixingLookupTable(mcDonorList, mcLookup, 0);
    }

    ConfigurationCard cardLookup("cardEECOOMCInclusiveR04.input");
    ConfigurationCard cardFallback("cardEECOOMCInclusiveR04.input");
    configureMcVeto(cardLookup);
    configureMcVeto(cardFallback);
    SetString(cardLookup, "MixingLookupTable", mcLookup);
    SetString(cardFallback, "MixingLookupTable", "none");

    auto anaLookup = std::make_unique<EECAnalyzer>(std::vector<TString>{mcSignal.c_str()}, &cardLookup, 1, 1, 0);
    auto anaFallback = std::make_unique<EECAnalyzer>(std::vector<TString>{mcSignal.c_str()}, &cardFallback, 1, 1, 0);
    anaLookup->fRng->SetSeed(99991);
    anaFallback->fRng->SetSeed(99991);

    anaLookup->RunAnalysis();
    anaFallback->RunAnalysis();

    const double mcEecSum = CompareTHnSparseExact(
      anaLookup->GetHistograms()->fhEnergyEnergyCorrelator,
      anaFallback->GetHistograms()->fhEnergyEnergyCorrelator,
      "OOMC_GenGen_Veto_EEC");
    Require(mcEecSum > 0.0, "OOMC_GenGen_Veto produced zero EEC weight");
    CompareTH1Exact(
      anaLookup->GetHistograms()->fhJetVetoEventCounts,
      anaFallback->GetHistograms()->fhJetVetoEventCounts,
      "OOMC_GenGen_Veto_Counts");
    // 4 donor files * 4 events = 16 standard-cut donor events; 1 vetoed per file = 4 vetoed, 12 kept
    auto* vetoCounts = anaLookup->GetHistograms()->fhJetVetoEventCounts;
    Require(vetoCounts->GetBinContent(4) == 16 && vetoCounts->GetBinContent(5) == 12 && vetoCounts->GetBinContent(6) == 4,
            "Unexpected mixing veto counts in OOMC GenGen companion check");
  }
  Require(gROOT->GetListOfFiles()->GetSize() == 0, "Leaked TFile handles after OOMC GenGen veto check");

  // 4. OOClosurePilot byte-for-byte ledger equivalence (crossing donor files with global pool cap 10)
  const std::string pilotLookup = (dir / "oomc_pilot.lookup.tsv").string();
  {
    auto configurePilot = [&](ConfigurationCard& c, const std::string& prefix) {
      SetNumeric(c, "TriggerSelection", 0);
      SetNumeric(c, "BackgroundMethods", 6);
      SetNumeric(c, "McCorrelationType", 3);
      SetNumeric(c, "JetAxis", 1);
      SetNumeric(c, "SubeventCut", 2);
      SetNumeric(c, "EventMatchForMixing", 2);
      SetNumeric(c, "BinnedMixing", 0);
      SetNumeric(c, "MegaSkimMixing", 0);
      SetNumeric(c, "HFEnergyShift", 0);
      SetNumeric(c, "ApplyJetVeto", 0);
      SetNumeric(c, "Run3JetID", 0);
      SetNumeric(c, "Run3JetQA", 0);
      SetNumeric(c, "DisableTrackPairEfficiencyCorrection", 1);
      SetNumeric(c, "OOClosurePilot", 1);
      SetNumeric(c, "OOClosureMaxSignalEvents", 4);
      SetNumeric(c, "OOClosureMaxPoolEvents", 10);
      SetNumeric(c, "OOClosureSeed", 424242);
      SetNumeric(c, "OOClosurePartitionPolicy", 0);
      SetString(c, "OOClosureOutputPrefix", prefix);
      SetString(c, "MixingFileList", mcDonorList);
    };

    ConfigurationCard buildCard("cardEECOOMCInclusiveR04.input");
    configurePilot(buildCard, (dir / "pilot_build_unused").string());
    {
      EECAnalyzer builder({}, &buildCard, 1, 1, 0);
      builder.BuildAndWriteMixingLookupTable(mcDonorList, pilotLookup, 0);
    }

    const std::string prefLookup = (dir / "pilot_lookup").string();
    const std::string prefFallback = (dir / "pilot_fallback").string();
    ConfigurationCard cardLookup("cardEECOOMCInclusiveR04.input");
    ConfigurationCard cardFallback("cardEECOOMCInclusiveR04.input");
    configurePilot(cardLookup, prefLookup);
    configurePilot(cardFallback, prefFallback);
    SetString(cardLookup, "MixingLookupTable", pilotLookup);
    SetString(cardFallback, "MixingLookupTable", "none");

    auto anaLookup = std::make_unique<EECAnalyzer>(std::vector<TString>{mcSignal.c_str()}, &cardLookup, 1, 1, 0);
    auto anaFallback = std::make_unique<EECAnalyzer>(std::vector<TString>{mcSignal.c_str()}, &cardFallback, 1, 1, 0);
    anaLookup->RunAnalysis();
    anaFallback->RunAnalysis();

    CompareTHnSparseExact(
      anaLookup->GetHistograms()->fhEnergyEnergyCorrelator,
      anaFallback->GetHistograms()->fhEnergyEnergyCorrelator,
      "OOClosurePilot_EEC");

    for (const char* suffix : {".partners.csv", ".pool.csv", ".jets.csv", ".contributions.csv"}) {
      const std::string lhs = ReadWholeFile(prefLookup + suffix);
      const std::string rhs = ReadWholeFile(prefFallback + suffix);
      Require(lhs == rhs, std::string("OOClosurePilot ledger mismatch for ") + suffix);
    }
  }
  Require(gROOT->GetListOfFiles()->GetSize() == 0, "Leaked TFile handles after OOClosurePilot check");

  // 5. Cross-lumi reused MC event numbers (real OO2025 MC pattern: evt=1..3 in lumi=1 and lumi=2)
  //    plus RecoGen lookup-table build destructor safety when readTrackTree == false.
  {
    const std::string lumi1File = (dir / "oomc_lumi1_same_evt.root").string();
    const std::string lumi2File = (dir / "oomc_lumi2_same_evt.root").string();
    // fileSalt=1 -> lumi=101, evt=1..3; fileSalt=2 -> lumi=102, evt=1..3 (identical evt numbers!)
    WriteOOMCForest(lumi1File, 3, 1ULL, false, 1);
    WriteOOMCForest(lumi2File, 3, 1ULL, true, 2);

    const std::string sameEvtPoolList = (dir / "oomc_same_evt_pool.txt").string();
    {
      std::ofstream out(sameEvtPoolList);
      out << lumi1File << "\n" << lumi2File << "\n";
    }
    const std::string recogenLookup = (dir / "oomc_recogen_same_evt.lookup.tsv").string();

    auto configureRecoGenSameEvt = [&](ConfigurationCard& c) {
      SetNumeric(c, "TriggerSelection", 0);
      SetNumeric(c, "BackgroundMethods", 6);
      SetNumeric(c, "McCorrelationType", 1); // kRecoGen: exercises fTrackReader == fJetReader when readTrackTree == false
      SetNumeric(c, "JetAxis", 1);
      SetNumeric(c, "SubeventCut", 2);
      SetNumeric(c, "EventMatchForMixing", 0); // HiBin matching
      SetNumeric(c, "BinnedMixing", 1);
      SetNumeric(c, "MegaSkimMixing", 0);
      SetNumeric(c, "ApplyJetVeto", 0);
      SetNumeric(c, "Run3JetID", 1);
      SetNumeric(c, "Run3JetQA", 0);
      SetNumeric(c, "DisableTrackPairEfficiencyCorrection", 1);
      SetNumeric(c, "MixingMaxOpenFiles", 2);
      SetString(c, "MixingFileList", sameEvtPoolList);
    };

    ConfigurationCard buildCard("cardEECOOMCInclusiveR04.input");
    configureRecoGenSameEvt(buildCard);
    {
      EECAnalyzer builder({}, &buildCard, 1, 1, 0);
      builder.BuildAndWriteMixingLookupTable(sameEvtPoolList, recogenLookup, 0);
    }

    ConfigurationCard cardLookup("cardEECOOMCInclusiveR04.input");
    ConfigurationCard cardFallback("cardEECOOMCInclusiveR04.input");
    configureRecoGenSameEvt(cardLookup);
    configureRecoGenSameEvt(cardFallback);
    SetString(cardLookup, "MixingLookupTable", recogenLookup);
    SetString(cardFallback, "MixingLookupTable", "none");

    // Note: lumi1File has evt=1 in vz=+0.5, evt=2 in vz=-0.5, evt=3 in vz=+0.5.
    // In vz=-0.5, lumi1File ONLY has evt=2 (1 event!), so its ONLY valid mixing partner in vz=-0.5
    // is lumi2File's evt=2 (run=1, lumi=102, evt=2). Comparing only evt would falsely reject it and fail!
    auto anaLookup = std::make_unique<EECAnalyzer>(std::vector<TString>{lumi1File.c_str()}, &cardLookup, 1, 1, 0);
    auto anaFallback = std::make_unique<EECAnalyzer>(std::vector<TString>{lumi1File.c_str()}, &cardFallback, 1, 1, 0);
    anaLookup->fRng->SetSeed(314159);
    anaFallback->fRng->SetSeed(314159);
    anaLookup->RunAnalysis();
    anaFallback->RunAnalysis();

    const double recogenSum = CompareTHnSparseExact(
      anaLookup->GetHistograms()->fhEnergyEnergyCorrelator,
      anaFallback->GetHistograms()->fhEnergyEnergyCorrelator,
      "OOMC_RecoGen_CrossLumiSameEvt_EEC");
    Require(recogenSum > 0.0, "OOMC_RecoGen_CrossLumiSameEvt produced zero EEC weight");
  }
  Require(gROOT->GetListOfFiles()->GetSize() == 0, "Leaked TFile handles after CrossLumiSameEvt check");

  std::cout << "PASS: CanonicalizeMixingPath, OO unbinned/binned lookup==fallback, OOMC GenGen+veto companion, OOClosurePilot ledgers, CrossLumiSameEvt+RecoGen, and zero leaked TFiles verified.\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    TH1::AddDirectory(kFALSE);
    if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " [--make-fixtures | --verify-equivalence | --all] <dir>\n";
      return 2;
    }
    const std::string mode = argv[1];
    const std::filesystem::path dir = argv[2];
    if (mode == "--make-fixtures") {
      MakeFixtures(dir);
    } else if (mode == "--verify-equivalence") {
      VerifyEquivalence(dir);
    } else if (mode == "--all") {
      MakeFixtures(dir);
      VerifyEquivalence(dir);
    } else {
      throw std::invalid_argument("Unknown mode: " + mode);
    }
    return 0;
  } catch (const std::exception& err) {
    std::cerr << "FAIL: " << err.what() << "\n";
    return 1;
  }
}

#include "trackingEfficiency2024pp.h"

// C++ includes
#include <iostream>
#include <stdexcept>

// ROOT includes
#include <TAxis.h>
#include <TFile.h>
#include <TH2.h>

TrkEff2024ppRef::TrkEff2024ppRef(
  bool isQuiet,
  const std::string& filePath
) :
  fIsQuiet(isQuiet),
  fInputFile(NULL),
  fEfficiency(NULL),
  fFakeRate(NULL),
  fSecondaryRate(NULL),
  fMultipleRecoRate(NULL)
{
  if(!fIsQuiet){
    std::cout
      << "Opening Run-3 ppRef tracking correction file: "
      << filePath
      << std::endl;
  }

  fInputFile = TFile::Open(
    filePath.c_str(),
    "READ"
  );

  if(fInputFile == NULL || fInputFile->IsZombie()){
    if(fInputFile != NULL){
      fInputFile->Close();
      delete fInputFile;
      fInputFile = NULL;
    }

    throw std::runtime_error(
      "Could not open Run-3 ppRef tracking correction file: " +
      filePath
    );
  }

  fEfficiency =
    dynamic_cast<TH2*>(
      fInputFile->Get("hEff_2D")
    );

  fFakeRate =
    dynamic_cast<TH2*>(
      fInputFile->Get("hFak_2D")
    );

  fSecondaryRate =
    dynamic_cast<TH2*>(
      fInputFile->Get("hSec_2D")
    );

  fMultipleRecoRate =
    dynamic_cast<TH2*>(
      fInputFile->Get("hMul_2D")
    );

  if(fEfficiency == NULL ||
     fFakeRate == NULL ||
     fSecondaryRate == NULL ||
     fMultipleRecoRate == NULL){

    fInputFile->Close();
    delete fInputFile;
    fInputFile = NULL;

    throw std::runtime_error(
      "The Run-3 ppRef tracking file does not contain all required "
      "histograms: hEff_2D, hFak_2D, hSec_2D, hMul_2D."
    );
  }
}

TrkEff2024ppRef::~TrkEff2024ppRef(){

  if(fInputFile != NULL){
    fInputFile->Close();
    delete fInputFile;
    fInputFile = NULL;
  }
}

bool TrkEff2024ppRef::CheckBounds(
  float pt,
  float eta
) const {

  if(fEfficiency == NULL){
    return false;
  }

  const TAxis* etaAxis =
    fEfficiency->GetXaxis();

  const TAxis* ptAxis =
    fEfficiency->GetYaxis();

  const double minimumEta =
    etaAxis->GetBinLowEdge(1);

  const double maximumEta =
    etaAxis->GetBinUpEdge(
      etaAxis->GetNbins()
    );

  const double minimumPt =
    ptAxis->GetBinLowEdge(1);

  const double maximumPt =
    ptAxis->GetBinUpEdge(
      ptAxis->GetNbins()
    );

  // Avoid underflow and overflow bins.
  if(eta < minimumEta || eta >= maximumEta){
    if(!fIsQuiet){
      std::cout
        << "TrkEff2024ppRef: eta = "
        << eta
        << " is outside ["
        << minimumEta
        << ", "
        << maximumEta
        << "). Returning zero."
        << std::endl;
    }

    return false;
  }

  if(pt < minimumPt || pt >= maximumPt){
    if(!fIsQuiet){
      std::cout
        << "TrkEff2024ppRef: pT = "
        << pt
        << " is outside ["
        << minimumPt
        << ", "
        << maximumPt
        << "). Returning zero."
        << std::endl;
    }

    return false;
  }

  return true;
}

float TrkEff2024ppRef::getCorrection(
  float pt,
  float eta,
  int hiBin
){

  // No centrality coordinate in this pp correction.
  (void)hiBin;

  if(!CheckBounds(pt, eta)){
    return 0;
  }

  const int efficiencyBin =
    fEfficiency->FindBin(
      eta,
      pt
    );

  const int fakeBin =
    fFakeRate->FindBin(
      eta,
      pt
    );

  const int secondaryBin =
    fSecondaryRate->FindBin(
      eta,
      pt
    );

  const int multipleRecoBin =
    fMultipleRecoRate->FindBin(
      eta,
      pt
    );

  const float efficiency =
    fEfficiency->GetBinContent(
      efficiencyBin
    );

  const float fakeRate =
    fFakeRate->GetBinContent(
      fakeBin
    );

  const float secondaryRate =
    fSecondaryRate->GetBinContent(
      secondaryBin
    );

  const float multipleRecoRate =
    fMultipleRecoRate->GetBinContent(
      multipleRecoBin
    );

  const float denominator =
    efficiency * (1.0f + multipleRecoRate);

  if(denominator <= 0.001f){
    if(!fIsQuiet){
      std::cout
        << "TrkEff2024ppRef: correction denominator is too small "
        << "for pT = "
        << pt
        << ", eta = "
        << eta
        << ". Returning zero."
        << std::endl;
    }

    return 0;
  }

  return
    (1.0f - fakeRate) *
    (1.0f - secondaryRate) /
    denominator;
}
#include "trackingEfficiency2025OO.h"

// ROOT includes
#include "TFile.h"
#include "TH3.h"
#include "TMath.h"

// C++ includes
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

  /**
   * Current nominal OO prescription:
   *
   * For tracks above 6.1 GeV, evaluate the correction using the table content
   * at pT = 6.1 GeV.
   */
  const float kMaximumLookupPt = 6.1;

}


/**
 * Check whether the track and event values are inside the supported range.
 */
bool TrkEff2025OO::checkBounds(float pt, float eta, int hiBin) {

  const bool inBounds =
    (
      TMath::Abs(eta) <= 2.4 &&
      hiBin >= 0 &&
      hiBin <= 199 &&
      pt > 0 &&
      pt <= 500
    );

  if (inBounds) {
    return true;
  }

  if (!isQuiet && !warnedAboutBounds) {

    std::cout
      << "TrkEff2025OO: correction requested outside the supported range "
      << "(|eta| <= 2.4, 0 <= hiBin <= 199, 0 < pT <= 500 GeV). "
      << "Returning correction 0. Further warnings are suppressed. "
      << "Received pT=" << pt
      << ", eta=" << eta
      << ", hiBin=" << hiBin
      << "."
      << std::endl;

    warnedAboutBounds = true;
  }

  return false;
}


/**
 * Calculate the multiplicative single-track correction.
 *
 * The nominal OO correction is:
 *
 *                  (1 - fake rate) * (1 - secondary rate)
 *   correction = ------------------------------------------------
 *                    efficiency * (1 + multiple rate)
 */
float TrkEff2025OO::getCorrection(float pt, float eta, int hiBin) {

  if (!checkBounds(pt, eta, hiBin)) {
    return 0;
  }

  /**
   * Follow the nominal OO table prescription:
   *
   * For tracks above 6.1 GeV, use the correction stored at pT = 6.1 GeV.
   */
  const float lookupPt =
    (pt > kMaximumLookupPt) ? kMaximumLookupPt : pt;

  /**
   * TH3::FindBin returns the global ROOT bin corresponding to:
   *
   *   X = eta
   *   Y = pT
   *   Z = hiBin
   */
  const double efficiencyValue =
    efficiency->GetBinContent(
      efficiency->FindBin(eta, lookupPt, hiBin)
    );

  const double fakeRateValue =
    fakeRate->GetBinContent(
      fakeRate->FindBin(eta, lookupPt, hiBin)
    );

  const double multipleRateValue =
    multipleRate->GetBinContent(
      multipleRate->FindBin(eta, lookupPt, hiBin)
    );

  const double secondaryRateValue =
    secondaryRate->GetBinContent(
      secondaryRate->FindBin(eta, lookupPt, hiBin)
    );

  /**
   * Protect against corrupted or undefined table values.
   */
  if (
    !std::isfinite(efficiencyValue) ||
    !std::isfinite(fakeRateValue) ||
    !std::isfinite(multipleRateValue) ||
    !std::isfinite(secondaryRateValue)
  ) {

    if (!isQuiet && !warnedAboutInvalidCorrection) {

      std::cout
        << "TrkEff2025OO: non-finite table content encountered. "
        << "Returning correction 0. Further warnings are suppressed."
        << std::endl;

      warnedAboutInvalidCorrection = true;
    }

    return 0;
  }

  /**
   * Protect against division by an efficiency that is zero or very close
   * to zero.
   */
  if (efficiencyValue <= 0.001) {

    if (!isQuiet && !warnedAboutLowEfficiency) {

      std::cout
        << "TrkEff2025OO: efficiency <= 0.001 encountered. "
        << "Returning correction 0. Further warnings are suppressed."
        << std::endl;

      warnedAboutLowEfficiency = true;
    }

    return 0;
  }

  const double denominator =
    efficiencyValue * (1.0 + multipleRateValue);

  if (denominator <= 0) {

    if (!isQuiet && !warnedAboutInvalidCorrection) {

      std::cout
        << "TrkEff2025OO: non-positive correction denominator encountered. "
        << "Returning correction 0. Further warnings are suppressed."
        << std::endl;

      warnedAboutInvalidCorrection = true;
    }

    return 0;
  }

  const double correction =
    (1.0 - fakeRateValue) *
    (1.0 - secondaryRateValue) /
    denominator;

  if (!std::isfinite(correction) || correction <= 0) {

    if (!isQuiet && !warnedAboutInvalidCorrection) {

      std::cout
        << "TrkEff2025OO: invalid correction factor encountered. "
        << "Returning correction 0. Further warnings are suppressed."
        << std::endl;

      warnedAboutInvalidCorrection = true;
    }

    return 0;
  }

  return static_cast<float>(correction);
}


/**
 * Constructor.
 *
 * Opens the nominal OO ROOT file and retrieves the four three-dimensional
 * correction histograms.
 */
TrkEff2025OO::TrkEff2025OO(
  bool isQuiet_,
  std::string fileName
) :
  isQuiet(isQuiet_),
  warnedAboutBounds(false),
  warnedAboutLowEfficiency(false),
  warnedAboutInvalidCorrection(false),
  correctionFile(NULL),
  efficiency(NULL),
  fakeRate(NULL),
  multipleRate(NULL),
  secondaryRate(NULL)
{
  if (fileName.empty()) {

    throw std::invalid_argument(
      "TrkEff2025OO: tracking correction file name is empty."
    );
  }

  correctionFile = TFile::Open(fileName.c_str(), "READ");

  if (
    correctionFile == NULL ||
    correctionFile->IsZombie() ||
    !correctionFile->IsOpen()
  ) {

    if (correctionFile != NULL) {
      correctionFile->Close();
      delete correctionFile;
      correctionFile = NULL;
    }

    throw std::runtime_error(
      std::string(
        "TrkEff2025OO: could not open tracking correction file: "
      ) + fileName
    );
  }

  /**
   * Histograms expected in:
   *
   * Eff_OO_2025_Hijing_MB_Centrality_fromHihfpf_NoPU_3D_
   * Nominal_Official_18Nov2025.root
   */
  efficiency =
    dynamic_cast<TH3*>(
      correctionFile->Get("hEff_3D")
    );

  fakeRate =
    dynamic_cast<TH3*>(
      correctionFile->Get("hFak_3D")
    );

  multipleRate =
    dynamic_cast<TH3*>(
      correctionFile->Get("hMul_3D")
    );

  secondaryRate =
    dynamic_cast<TH3*>(
      correctionFile->Get("hSec_3D")
    );

  /**
   * Crash immediately if the table does not have the expected structure.
   *
   * This is preferable to obtaining a null-pointer crash later in the
   * event loop.
   */
  if (
    efficiency == NULL ||
    fakeRate == NULL ||
    multipleRate == NULL ||
    secondaryRate == NULL
  ) {

    correctionFile->Close();
    delete correctionFile;
    correctionFile = NULL;

    throw std::runtime_error(
      std::string(
        "TrkEff2025OO: expected hEff_3D, hFak_3D, "
        "hMul_3D and hSec_3D in file: "
      ) + fileName
    );
  }

  if (!isQuiet) {

    std::cout
      << "TrkEff2025OO: using nominal OO tracking correction table: "
      << fileName
      << std::endl;

    std::cout
      << "TrkEff2025OO: pT lookups above 6.1 GeV "
      << "are capped at 6.1 GeV."
      << std::endl;
  }
}


/**
 * Destructor.
 */
TrkEff2025OO::~TrkEff2025OO() {

  if (correctionFile != NULL) {

    correctionFile->Close();
    delete correctionFile;
    correctionFile = NULL;
  }
}
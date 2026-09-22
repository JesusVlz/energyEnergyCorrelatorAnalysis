#ifndef TRKEFF2025OO
#define TRKEFF2025OO

#include <string>

#include "TrackingEfficiencyInterface.h"

class TFile;
class TH3;

/**
 * Tracking-efficiency correction provider for 2025 OO collisions.
 *
 * The correction is read from the nominal OO ROOT table as a function of:
 *
 *   - track eta
 *   - track pT
 *   - event hiBin
 *
 * The class implements TrackingEfficiencyInterface so that EECAnalyzer can use
 * it through the existing fTrackEfficiencyCorrector pointer.
 */
class TrkEff2025OO : public TrackingEfficiencyInterface {
public:

  /**
   * Constructor.
   *
   * @param isQuiet_ Suppress informational and warning messages when true.
   * @param fileName Complete path to the OO tracking-correction ROOT file.
   */
  TrkEff2025OO(bool isQuiet_ = false, std::string fileName = "");

  /**
   * Destructor.
   *
   * Closes and deletes the ROOT correction file.
   */
  virtual ~TrkEff2025OO();

  /**
   * Return the multiplicative single-track correction.
   *
   * @param pt Track transverse momentum in GeV.
   * @param eta Track pseudorapidity.
   * @param hiBin Event centrality index.
   *
   * @return Multiplicative tracking correction.
   */
  float getCorrection(float pt, float eta, int hiBin);

private:

  /**
   * Verify that the requested track and event values are inside the supported
   * analysis range.
   */
  bool checkBounds(float pt, float eta, int hiBin);

  // Control informational output.
  bool isQuiet;

  // Suppress repeated warnings during large event loops.
  bool warnedAboutBounds;
  bool warnedAboutLowEfficiency;
  bool warnedAboutInvalidCorrection;

  // ROOT file containing the correction histograms.
  TFile* correctionFile;

  // Three-dimensional correction components.
  //
  // Expected axis order:
  //   X = track eta
  //   Y = track pT
  //   Z = event hiBin
  TH3* efficiency;
  TH3* fakeRate;
  TH3* multipleRate;
  TH3* secondaryRate;
};

#endif
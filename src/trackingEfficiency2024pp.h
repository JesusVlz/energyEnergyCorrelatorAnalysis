#ifndef TRKEFF2024PP
#define TRKEFF2024PP

// C++ includes
#include <string>

// Own includes
#include "TrackingEfficiencyInterface.h"

// Forward declarations
class TFile;
class TH2;

class TrkEff2024ppRef : public TrackingEfficiencyInterface {

public:

  TrkEff2024ppRef(
    bool isQuiet,
    const std::string& filePath
  );

  ~TrkEff2024ppRef() override;

  // hiBin is unused because the ppRef correction is two-dimensional:
  // eta x pT.
  float getCorrection(
    float pt,
    float eta,
    int hiBin
  ) override;

private:

  bool CheckBounds(
    float pt,
    float eta
  ) const;

  bool fIsQuiet;

  TFile* fInputFile;

  TH2* fEfficiency;
  TH2* fFakeRate;
  TH2* fSecondaryRate;
  TH2* fMultipleRecoRate;
};

#endif
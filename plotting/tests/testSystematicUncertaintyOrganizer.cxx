/* Build from the repository root with:
 * g++ $(root-config --cflags) -Iplotting plotting/tests/testSystematicUncertaintyOrganizer.cxx
 *   plotting/SystematicUncertaintyOrganizer.cxx plotting/EECCard.cxx
 *   $(root-config --libs) -o /tmp/testSystematicUncertaintyOrganizer
 * (Join the compiler command onto one line.)
 */

#include "SystematicUncertaintyOrganizer.h"

#include <TMemFile.h>
#include <TVectorT.h>

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void Require(bool condition, const char* message){
  if(!condition) throw std::runtime_error(message);
}

void WriteSingleBinCard(TMemFile& file){
  file.mkdir("JCard")->cd();
  TVectorT<float> value(1, 1);
  value[1] = 0;
  const char* entries[] = {
    "DataType", "McCorrelationType", "FirstUnfoldedCentralityBin",
    "LastUnfoldedCentralityBin", "FirstUnfoldedJetPtBin", "LastUnfoldedJetPtBin",
    "FirstUnfoldedTrackPtBin", "LastUnfoldedTrackPtBin"
  };
  for(const char* entry : entries) value.Write(entry);
  file.cd();
}

void CheckMissingInput(){
  TMemFile file("missing_uncertainties.root", "RECREATE");
  WriteSingleBinCard(file);

  bool rejected = false;
  try {
    SystematicUncertaintyOrganizer organizer(&file);
  } catch(const std::runtime_error& error){
    const std::string message = error.what();
    rejected = message.find("no systematic uncertainty histograms") != std::string::npos &&
               message.find(file.GetName()) != std::string::npos &&
               message.find("C[0,0], J[0,0], T[0,0]") != std::string::npos;
  }
  Require(rejected, "Missing uncertainty histograms must report an explicit input error");
}

void CheckPartialInput(){
  TMemFile file("partial_uncertainties.root", "RECREATE");
  WriteSingleBinCard(file);

  const double edges[] = {0.01, 0.05, 0.2, 0.4};
  TH1D reference("systematicUncertainty_jetEnergyResolution_C0J0T0", "", 3, edges);
  reference.SetDirectory(nullptr);
  for(int bin = 1; bin <= reference.GetNbinsX(); ++bin){
    reference.SetBinContent(bin, 10. * bin);
    reference.SetBinError(bin, 0.1 * bin);
  }
  reference.Write();

  SystematicUncertaintyOrganizer organizer(&file);
  for(int source = 0; source < organizer.GetNUncertaintySources(); ++source){
    TH1D* histogram = organizer.GetSystematicUncertainty(0, 0, 0, source);
    Require(histogram != nullptr, "A missing source must receive a fallback histogram");
    Require(histogram->GetNbinsX() == reference.GetNbinsX(), "Fallback must retain the reference bin count");
    for(int bin = 1; bin <= reference.GetNbinsX(); ++bin){
      Require(histogram->GetBinLowEdge(bin) == reference.GetBinLowEdge(bin), "Fallback must retain bin edges");
      Require(histogram->GetBinWidth(bin) == reference.GetBinWidth(bin), "Fallback must retain bin widths");
      Require(histogram->GetBinContent(bin) == reference.GetBinContent(bin), "Fallback must retain central values");
      const double expectedError = source == SystematicUncertaintyOrganizer::kJetEnergyResolution ? reference.GetBinError(bin) : 0.;
      Require(histogram->GetBinError(bin) == expectedError, "Existing errors must remain and missing-source errors must be zero");
    }
  }
}

} // namespace

int main(){
  try {
    CheckMissingInput();
    CheckPartialInput();
  } catch(const std::exception& error){
    std::cerr << error.what() << std::endl;
    return 1;
  }
  std::cout << "SystematicUncertaintyOrganizer regression checks passed" << std::endl;
  return 0;
}

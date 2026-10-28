// Raw selected-event projections only. No unit-area normalization, subtraction,
// unfolding, radius ratio or closure inference is performed by this macro.
#include <TFile.h>
#include <THnBase.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TVectorT.h>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

void projectInclusiveRaw(const char* inputName, const char* outputName){
  TFile input(inputName, "READ");
  if(input.IsZombie()) throw std::runtime_error("Cannot read inclusive input");
  auto raw = dynamic_cast<THnBase*>(input.Get("inclusive/inclusiveRawEEC"));
  auto exponents = dynamic_cast<TVectorT<float>*>(input.Get("JCard/WeightExponent"));
  auto radiusVector = dynamic_cast<TVectorT<float>*>(input.Get("JCard/JetRadius"));
  auto mode = dynamic_cast<TVectorT<float>*>(input.Get("JCard/InclusiveMode"));
  auto trackSelection = dynamic_cast<TNamed*>(input.Get("inclusive/TrackSelection"));
  if(!raw || raw->GetNdimensions() != 6 || !exponents || !radiusVector || !mode ||
     (*mode)[mode->GetLwb()] != 1)
    throw std::runtime_error("Requires explicit InclusiveMode 1 raw diagnostic input");
  if(!trackSelection || std::string(trackSelection->GetTitle()) != "strict pT > 2 GeV; input threshold remains configured separately")
    throw std::runtime_error("Requires strict pT > 2 inclusive fill; earlier raw files must be rerun");
  const double radius = (*radiusVector)[radiusVector->GetLwb()];
  const int radiusCode = std::lround(radius*10);
  if((radiusCode != 4 && radiusCode != 8) || std::abs(radius-radiusCode/10.) > 1e-5) throw std::runtime_error("Unsupported radius");
  if(std::abs(raw->GetAxis(0)->GetXmax()-2*radius) > 1e-5)
    throw std::runtime_error("Raw angular support does not reach 2R");
  const double jetEdges[] = {40,80,120,160};
  auto jetAxis = raw->GetAxis(1);
  for(double edge : jetEdges){
    bool present = false;
    for(int bin = 1; bin <= jetAxis->GetNbins()+1; ++bin)
      if(std::abs(jetAxis->GetBinLowEdge(bin)-edge) < 1e-5) present = true;
    if(!present) throw std::runtime_error("Requested primary jet boundary is unavailable");
  }
  auto trackAxis = raw->GetAxis(2);
  const int firstTrack = trackAxis->FindFixBin(2.000001);
  if(std::abs(trackAxis->GetBinLowEdge(firstTrack)-2) > 1e-5)
    throw std::runtime_error("Track threshold 2 GeV cannot be resolved");
  // NEW prevents overwriting a previous result.
  TFile output(outputName, "NEW");
  if(output.IsZombie()) throw std::runtime_error("Output exists or cannot be created");
  raw->GetAxis(2)->SetRange(firstTrack,trackAxis->GetNbins()+1);
  raw->GetAxis(3)->SetRange(1,1); // same-signal-cone pair type = 0
  for(int n : {1,2}){
    int exponentIndex = -1;
    for(int i=0; i<exponents->GetNrows(); ++i)
      if((*exponents)[exponents->GetLwb()+i] == n) exponentIndex = i;
    if(exponentIndex < 0) throw std::runtime_error("Missing nominal energy power");
    raw->GetAxis(5)->SetRange(exponentIndex+1,exponentIndex+1);
    for(int j=0; j<3; ++j){
      jetAxis->SetRange(jetAxis->FindFixBin(jetEdges[j]+1e-5),
                        jetAxis->FindFixBin(jetEdges[j+1]-1e-5));
      std::unique_ptr<TH1D> spectrum(raw->Projection(0,"E"));
      spectrum->SetName(Form("INC_R%02d_J%.0f_%.0f_N%d_rawSelected",radiusCode,jetEdges[j],jetEdges[j+1],n));
      spectrum->SetTitle("Raw selected-event diagnostic; no single-track/pair/centrality correction; not final");
      spectrum->Write();
    }
  }
  if(auto exposure = input.Get("inclusive/inclusiveJetExposure")) exposure->Write();
  if(auto coverage = input.Get("inclusive/inclusiveEventCoverage")) coverage->Write();
  TObjString source(inputName); source.Write("source_input");
  TNamed state("correction_state", "RAW_DIAGNOSTIC; final corrections, normalization range, closure and final covariance BLOCKED/NOT_RUN");
  state.Write();
  trackSelection->Write("TrackSelection");
  if(auto card = input.GetDirectory("JCard")){
    auto destination=output.mkdir("JCard");
    TIter next(card->GetListOfKeys());
    while(auto key=next()){
      destination->cd();
      if(auto object=card->Get(key->GetName())) object->Write(key->GetName());
    }
  }
  output.Write();
}

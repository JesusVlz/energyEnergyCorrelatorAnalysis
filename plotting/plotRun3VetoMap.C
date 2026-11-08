#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <cmath>

// Plot the original payload (or one single-job run3JetQA/configuredJetVetoMap).
// An occupancy map is NOT a substitute. Do not use a hadd-summed map snapshot.
void plotRun3VetoMap(const char* fileName, const char* output="jetVetoMap.png", const char* key="jetvetomap_all") {
  std::unique_ptr<TFile> file(TFile::Open(fileName,"READ"));
  auto* source=file && !file->IsZombie()?dynamic_cast<TH2*>(file->Get(key)):nullptr;
  if(!source) throw std::runtime_error("Missing requested two-dimensional veto payload");
  std::unique_ptr<TH2> map(dynamic_cast<TH2*>(source->Clone("map_to_plot")));
  map->SetDirectory(nullptr);
  double maximum=1;
  for(int x=1;x<=map->GetNbinsX();++x) for(int y=1;y<=map->GetNbinsY();++y) {
    const double value=map->GetBinContent(x,y);
    if(!std::isfinite(value)) throw std::runtime_error("Nonfinite veto-map content");
    maximum=std::max(maximum,std::abs(value));
  }
  gStyle->SetOptStat(0);
  TCanvas canvas("veto_map_canvas","Veto payload",1000,760);
  canvas.SetLeftMargin(.14); canvas.SetRightMargin(.24);
  canvas.SetBottomMargin(.14); canvas.SetTopMargin(.06);
  for(auto* axis:{map->GetXaxis(),map->GetYaxis(),map->GetZaxis()}) {
    axis->SetTitleFont(42); axis->SetLabelFont(42);
    axis->SetTitleSize(.035); axis->SetLabelSize(.032);
    axis->SetLabelOffset(.008); axis->CenterTitle(false);
  }
  map->GetXaxis()->SetTitleOffset(1.15);
  map->GetYaxis()->SetTitleOffset(1.35);
  map->GetZaxis()->SetTitleOffset(1.45);
  map->SetTitle(""); map->GetXaxis()->SetTitle("Detector jet #eta"); map->GetYaxis()->SetTitle("Detector jet #phi (rad)");
  map->GetZaxis()->SetTitle("Veto-map value (nonzero rejects)");
  map->SetMinimum(-maximum); map->SetMaximum(maximum); map->Draw("COLZ");
  canvas.SaveAs(output);
}

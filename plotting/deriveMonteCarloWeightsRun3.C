// Run-3 event-weight derivation and production-weight QA.
// Put this file beside plotting/JDrawer.h. No libEECDrawing.so is required.
//
// root -l -b -q 'deriveMonteCarloWeightsRun3.C("OO_DATA.root","OO_MC.root","PP_DATA.root","PP_MC.root")'
// Use empty strings for BOTH members of a system to skip that system.
// Modes: "stored" (actual production histograms), "derive" (1D refit checks),
//        "both" (default). A refit check is NOT independent event-level closure.
// Source audit: JesusVlz/energyEnergyCorrelatorAnalysis @ 6d4800b8c57747e1db761fc4fecf72a38a1a00c3.
// This macro never edits the analyzer, event weights, or input ROOT files.

#ifndef EEC_DERIVE_MONTE_CARLO_WEIGHTS_RUN3_C
#define EEC_DERIVE_MONTE_CARLO_WEIGHTS_RUN3_C

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h> // Protect local input aliases/symlinks from RECREATE.

#include <TCanvas.h>
#include <TDecompSVD.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2.h>   // JDrawer.h also uses TH2F.
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TPad.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorD.h>
#include <TVectorF.h>

#include "JDrawer.h"

namespace Run3MCWeights {

// Plot-only choices. The fits always use the ORIGINAL fine bins.
constexpr double kVzMax = 15.0;
constexpr double kCentralitySplit = 30.0;
constexpr double kRatioMin = 0.50;
constexpr double kRatioMax = 1.50;
const char* const kOOLabel = "P+H"; // Keep the generator shorthand in the supplied example.
const char* const kPPLabel = "PYTHIA";

// BEGIN PORTABLE MATH
// Scaled monomials avoid fitting tiny coefficients in powers of c ~ 90.
// The two centrality pieces have a SHARED value at 30%, fitted jointly.
std::vector<double> Basis(double x, bool centrality, double cmax) {
  if (!centrality) {
    std::vector<double> b(7, 1.0);
    for (int j = 1; j < 7; ++j) b[j] = b[j-1] * (x/kVzMax);
    return b;
  }
  if (!(cmax > kCentralitySplit)) throw std::runtime_error("Centrality range must extend above 30%.");
  std::vector<double> b(13, 0.0);
  b[0] = 1.0;
  const bool left = x < kCentralitySplit;
  const double t = (x-kCentralitySplit)/(left ? kCentralitySplit : cmax-kCentralitySplit);
  double power = 1.0;
  for (int j = 1; j <= 6; ++j) {
    power *= t;
    b[left ? j : j+6] = power;
  }
  return b;
}

double Evaluate(const std::vector<double>& p, double x, bool centrality, double cmax) {
  const auto b = Basis(x, centrality, cmax);
  if (b.size() != p.size()) throw std::runtime_error("Wrong number of fit coefficients.");
  double value = 0.0;
  for (std::size_t j = 0; j < p.size(); ++j) value += p[j]*b[j];
  return value;
}

double Choose(int n, int k) {
  double answer = 1.0;
  for (int j = 1; j <= k; ++j) answer *= double(n-j+1)/j;
  return answer;
}

// Convert sum a_j ((x-origin)/scale)^j to ROOT pol6 parameters in x.
std::array<double,7> ToPhysicalPolynomial(const std::vector<double>& p,
                                         bool centrality, bool peripheral,
                                         double cmax) {
  if (p.size() != (centrality ? 13u : 7u)) throw std::runtime_error("Invalid fit parameter count.");
  const double origin = centrality ? kCentralitySplit : 0.0;
  const double scale = centrality ? (peripheral ? cmax-kCentralitySplit : kCentralitySplit) : kVzMax;
  std::array<double,7> result{};
  for (int j = 0; j <= 6; ++j) {
    const double a = p[j == 0 ? 0 : (centrality && peripheral ? j+6 : j)];
    for (int k = 0; k <= j; ++k)
      result[k] += a * Choose(j,k) * std::pow(-origin,j-k) / std::pow(scale,j);
  }
  return result;
}
// END PORTABLE MATH

TString Unique(const TString& stem) {
  static unsigned long counter = 0;
  return TString::Format("%s_%lu", stem.Data(), ++counter);
}

void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

std::unique_ptr<TFile> OpenInput(const char* filename) {
  std::unique_ptr<TFile> f(TFile::Open(filename, "READ"));
  Require(f && !f->IsZombie(), std::string("Cannot open input: ")+filename);
  return f;
}

TH1* InputHistogram(TFile& file, const char* key, bool required=true) {
  TH1* h = dynamic_cast<TH1*>(file.Get(key));
  if (!h && !required) return nullptr;
  Require(h && h->GetDimension()==1,
          std::string(file.GetName())+": missing one-dimensional histogram '"+key+"'. "
          "Use an analyzer/processed ROOT file with event-information histograms, not a HiForest.");
  return h;
}

bool CardNumber(TFile& file, const char* key, double& value) {
  const TString path = TString("JCard/")+key;
  TObject* object = file.Get(path);
  if (auto* v = dynamic_cast<TVectorF*>(object)) {
    if (v->GetNrows()>0) { value=(*v)[v->GetLwb()]; return true; }
  }
  if (auto* v = dynamic_cast<TVectorD*>(object)) {
    if (v->GetNrows()>0) { value=(*v)[v->GetLwb()]; return true; }
  }
  return false;
}

void AuditCard(TFile& file, int expectedType, std::ostream& log) {
  log << "Input: " << file.GetName() << '\n';
  double value = 0;
  if (CardNumber(file, "DataType", value)) {
    log << "  DataType = " << value << '\n';
    Require(std::abs(value-expectedType)<1e-5,
            std::string(file.GetName())+": wrong JCard/DataType; expected "+std::to_string(expectedType));
  } else {
    log << "  WARNING: JCard/DataType absent; verify the supplied system manually.\n";
  }
  for (const char* key : {"TriggerSelection", "ZVertexCut", "MultiplicityMode", "JetPtWeight"}) {
    if (CardNumber(file,key,value)) log << "  " << key << " = " << value << '\n';
  }
  if (CardNumber(file,"MultiplicityMode",value))
    Require(value==0.0, "MultiplicityMode is not supported: derive centrality weights using hiBin, not inferred multiplicity.");
  if (CardNumber(file,"ZVertexCut",value))
    Require(std::abs(value-kVzMax)<1e-5, "ZVertexCut differs from 15 cm; change kVzMax explicitly for a different selection.");
  if (auto* s=dynamic_cast<TObjString*>(file.Get("JCard/GitHash")))
    log << "  producer GitHash = " << s->GetString() << '\n';
  if (auto* s=dynamic_cast<TObjString*>(file.Get("JCard/ProjectionGitHash")))
    log << "  projection GitHash = " << s->GetString() << '\n';
}

// Copy entire native bins within [low,high); exclude flow bins explicitly.
// Do not mutate file-owned histograms or rebin before determining fit coordinates.
TH1D* Crop(const TH1& input, const TString& name, double low, double high) {
  const TAxis* axis = input.GetXaxis();
  const int first = axis->FindFixBin(low+1e-8);
  const int last = axis->FindFixBin(high-1e-8);
  Require(first>=1 && last<=input.GetNbinsX() && first<=last,
          std::string(input.GetName())+": requested range lies outside histogram support.");
  Require(std::abs(axis->GetBinLowEdge(first)-low)<1e-7 &&
          std::abs(axis->GetBinUpEdge(last)-high)<1e-7,
          "Range boundaries must coincide with native histogram bin edges.");
  std::vector<double> edges;
  for (int i=first; i<=last+1; ++i) edges.push_back(axis->GetBinLowEdge(i));
  auto* out = new TH1D(name, "", last-first+1, edges.data());
  out->SetDirectory(nullptr);
  out->Sumw2();
  for (int i=first; i<=last; ++i) {
    const double y=input.GetBinContent(i), e=input.GetBinError(i);
    Require(std::isfinite(y) && y>=0 && std::isfinite(e) && e>=0,
            std::string(input.GetName())+": invalid/negative histogram bin; signed-weight samples are not supported here.");
    out->SetBinContent(i-first+1,y);
    out->SetBinError(i-first+1,e);
  }
  out->SetEntries(input.GetEntries());
  return out;
}

void SameBinning(const TH1& a, const TH1& b) {
  Require(a.GetNbinsX()==b.GetNbinsX(), "Data/MC bin counts differ; provide matching native binning.");
  for (int i=1; i<=a.GetNbinsX()+1; ++i)
    Require(std::abs(a.GetXaxis()->GetBinLowEdge(i)-b.GetXaxis()->GetBinLowEdge(i))<1e-8,
            "Data/MC bin edges differ; do not interpolate event-count histograms.");
}

void Normalize(TH1& h) {
  const double area=h.Integral(1,h.GetNbinsX());
  Require(std::isfinite(area) && area>0, std::string(h.GetName())+": empty or invalid histogram.");
  h.Scale(1.0/area); // Probability per bin, NOT density per unit x.
}

TH1D* Copy(const TH1D& h, const TString& name) {
  auto* result=static_cast<TH1D*>(h.Clone(name));
  result->SetDirectory(nullptr);
  return result;
}

struct Distribution {
  TString system, observable;
  bool centrality=false;
  double cmax=90;
  TH1D *data=nullptr, *before=nullptr, *stored=nullptr;
};

Distribution LoadDistribution(TFile& dataFile, TFile& mcFile, const TString& system,
                              bool centrality, double cmax, bool requireStored,
                              std::ostream& log) {
  Distribution d;
  d.system=system; d.centrality=centrality; d.cmax=cmax;
  d.observable=centrality ? "centrality" : "vz";
  const char* rawKey=centrality ? "centrality" : "vertexZ";
  const char* weightedKey=centrality ? "centralityWeighted" : "vertexZweighted";
  const double lo=centrality ? 0 : -kVzMax, hi=centrality ? cmax : kVzMax;
  // Data weighted QA retains any trigger-combination event weights, as in Jussi's original macro.
  d.data=Crop(*InputHistogram(dataFile,weightedKey),Unique("data"),lo,hi);
  d.before=Crop(*InputHistogram(mcFile,rawKey),Unique("mc_before"),lo,hi);
  TH1* storedSource=InputHistogram(mcFile,weightedKey,requireStored);
  if (storedSource) d.stored=Crop(*storedSource,Unique("mc_stored"),lo,hi);
  SameBinning(*d.data,*d.before);
  if (d.stored) SameBinning(*d.before,*d.stored);
  if (centrality) {
    for (int i=1; i<=d.before->GetNbinsX(); ++i) {
      const double edge=d.before->GetXaxis()->GetBinLowEdge(i);
      Require(std::abs(d.before->GetBinWidth(i)-0.5)<1e-8 &&
              std::abs(2*edge-std::round(2*edge))<1e-8,
              "Centrality needs native 0.5%-wide bins aligned with hiBin/2. "
              "A rebinned histogram cannot recover the discrete event coordinate.");
    }
  }
  log << "\n[" << system << '/' << d.observable << "] range [" << lo << ',' << hi << ")\n"
      << "  data key: " << weightedKey << ", MC baseline key: " << rawKey << '\n'
      << "  data integral = " << d.data->Integral() << ", MC baseline integral = " << d.before->Integral() << '\n';
  if (d.stored) {
    double difference=0;
    for (int i=1; i<=d.before->GetNbinsX(); ++i)
      difference=std::max(difference,std::abs(d.before->GetBinContent(i)-d.stored->GetBinContent(i)));
    log << "  MC stored-weighted integral = " << d.stored->Integral()
        << ", max absolute raw/stored bin difference = " << difference << '\n';
    if (difference==0)
      log << "  WARNING: raw and stored-weighted CONTENTS ARE IDENTICAL. This file may predate the active weights; "
             "do not assume it includes current source-code weights.\n";
  }
  Normalize(*d.data); Normalize(*d.before);
  if (d.stored) Normalize(*d.stored);
  return d;
}

// Ratios are represented by graphs so bins with a zero denominator are omitted,
// rather than plotted as a misleading zero. Errors are diagonal display errors.
TGraphErrors* RatioGraph(const TH1& numerator, const TH1& denominator,
                        bool centralityCoordinate=false) {
  SameBinning(numerator,denominator);
  auto* graph=new TGraphErrors();
  graph->SetName(Unique("ratio"));
  int point=0;
  for (int i=1; i<=numerator.GetNbinsX(); ++i) {
    const double a=numerator.GetBinContent(i), b=denominator.GetBinContent(i);
    if (b<=0) continue;
    const double ea=numerator.GetBinError(i), eb=denominator.GetBinError(i);
    const double x=centralityCoordinate ? numerator.GetXaxis()->GetBinLowEdge(i) : numerator.GetBinCenter(i);
    graph->SetPoint(point,x,a/b);
    graph->SetPointError(point,0,std::hypot(ea/b,a*eb/(b*b)));
    ++point;
  }
  return graph;
}

struct Fit {
  std::vector<double> p;
  TMatrixD covariance;
  double chi2=0;
  int ndof=0;
  TF1 *vertex=nullptr, *central=nullptr, *peripheral=nullptr;
  TH1D* after=nullptr;
};

TF1* MakePolynomial(const TString& name, const std::array<double,7>& p,
                    double lo, double hi) {
  auto* f=new TF1(Unique(name), "pol6",lo,hi);
  f->SetParameters(p.data());
  f->SetNpx(1000);
  return f;
}

Fit Derive(const Distribution& d, std::ostream& log) {
  std::vector<double> x,y,error;
  int skipped=0;
  for (int i=1; i<=d.data->GetNbinsX(); ++i) {
    const double a=d.data->GetBinContent(i), b=d.before->GetBinContent(i);
    Require(!(a>0 && b<=0), "Data populate a bin with zero baseline MC: a multiplicative weight cannot supply missing MC support.");
    if (a<=0 || b<=0) { ++skipped; continue; }
    const double e=std::hypot(d.data->GetBinError(i)/b,a*d.before->GetBinError(i)/(b*b));
    if (!(e>0) || !std::isfinite(e)) { ++skipped; continue; }
    // The filled c is hiBin/2.0, at the LOWER EDGE of each native 0.5% bin.
    x.push_back(d.centrality ? d.data->GetXaxis()->GetBinLowEdge(i) : d.data->GetBinCenter(i));
    y.push_back(a/b); error.push_back(e);
  }
  const int npar=d.centrality ? 13 : 7;
  Require(int(x.size())>npar, "Too few populated bins to fit pol6; increase MC statistics.");
  TMatrixD design(x.size(),npar);
  TVectorD rhs(x.size());
  for (int i=0; i<int(x.size()); ++i) {
    const auto b=Basis(x[i],d.centrality,d.cmax);
    rhs[i]=y[i]/error[i];
    for (int j=0; j<npar; ++j) design(i,j)=b[j]/error[i];
  }
  TDecompSVD svd(design);
  Require(svd.Decompose(), "SVD decomposition of the polynomial fit failed.");
  const TVectorD singular=svd.GetSig();
  Require(singular.GetNrows()==npar && singular[npar-1]>1e-12*singular[0],
          "Rank-deficient or ill-conditioned fit; inspect the available MC support.");
  Bool_t solved=false;
  TVectorD solution=svd.Solve(rhs,solved);
  Require(solved && solution.GetNrows()>=npar, "Weighted least-squares fit failed.");
  Fit result;
  result.p.resize(npar);
  for (int j=0; j<npar; ++j) {
    result.p[j]=solution[j];
    Require(std::isfinite(solution[j]), "Non-finite fit coefficient.");
  }
  result.ndof=int(x.size())-npar;
  for (int i=0; i<int(x.size()); ++i) {
    const double pull=(y[i]-Evaluate(result.p,x[i],d.centrality,d.cmax))/error[i];
    result.chi2+=pull*pull;
  }
  result.covariance.ResizeTo(npar,npar);
  const TMatrixD& v=svd.GetV();
  for (int i=0; i<npar; ++i)
    for (int j=0; j<npar; ++j) {
      result.covariance(i,j)=0;
      for (int k=0; k<npar; ++k)
        result.covariance(i,j)+=v(i,k)*v(j,k)/(singular[k]*singular[k]);
    }
  const double low=d.centrality ? 0 : -kVzMax, high=d.centrality ? d.cmax : kVzMax;
  double minWeight=std::numeric_limits<double>::max(),maxWeight=0;
  // Check between bins too: negative polynomial excursions are never clipped away silently.
  for (int i=0; i<=10000; ++i) {
    const double at=low+(high-low)*i/10000.;
    const double weight=Evaluate(result.p,at,d.centrality,d.cmax);
    Require(std::isfinite(weight) && weight>0, "Candidate polynomial is non-positive within its domain; no deployable fit is written.");
    minWeight=std::min(minWeight,weight); maxWeight=std::max(maxWeight,weight);
  }
  if (d.centrality) {
    result.central=MakePolynomial("weight_central",ToPhysicalPolynomial(result.p,true,false,d.cmax),0,kCentralitySplit);
    result.peripheral=MakePolynomial("weight_peripheral",ToPhysicalPolynomial(result.p,true,true,d.cmax),kCentralitySplit,d.cmax);
    Require(std::abs(result.central->Eval(kCentralitySplit)-result.peripheral->Eval(kCentralitySplit))<1e-6,
            "Exported centrality polynomials lost continuity numerically.");
  } else {
    result.vertex=MakePolynomial("weight_vz",ToPhysicalPolynomial(result.p,false,false,d.cmax),-kVzMax,kVzMax);
  }
  // Verify the exported physical-x pol6 reproduces the well-conditioned fit.
  for (int i=0; i<=1000; ++i) {
    const double at=low+(high-low)*i/1000.;
    const double expected=Evaluate(result.p,at,d.centrality,d.cmax);
    TF1* exported=d.centrality ? (at<kCentralitySplit ? result.central : result.peripheral) : result.vertex;
    const double value=exported->Eval(at);
    Require(std::isfinite(value) && value>0 && std::abs(value-expected)<1e-7*std::max(1.,std::abs(expected)),
            "Exported physical-x polynomial does not reproduce the scaled fit accurately.");
  }
  result.after=Copy(*d.before,Unique("mc_refit_preview"));
  for (int i=1; i<=result.after->GetNbinsX(); ++i) {
    const double at=d.centrality ? result.after->GetXaxis()->GetBinLowEdge(i) : result.after->GetBinCenter(i);
    const double w=Evaluate(result.p,at,d.centrality,d.cmax);
    result.after->SetBinContent(i,result.after->GetBinContent(i)*w);
    result.after->SetBinError(i,result.after->GetBinError(i)*w);
  }
  Normalize(*result.after);
  log << "  candidate fit: diagonal-error chi2/ndof = " << result.chi2 << '/' << result.ndof
      << ", skipped zero-content/error bins = " << skipped << '\n'
      << "  positive weight range on dense grid = [" << minWeight << ',' << maxWeight << "]\n"
      << "  CAUTION: same-sample 1D check; does not test the product wVz*wCentrality.\n";
  return result;
}

TH1D* ForDisplay(const TH1D& input, bool centrality, int rebin) {
  Require(rebin>0 && input.GetNbinsX()%rebin==0,
          "Plot rebin factor must divide the number of selected bins exactly.");
  TH1D* rebinned=Copy(input,Unique("display_native"));
  if (rebin>1) rebinned->Rebin(rebin);
  if (!centrality) return rebinned;
  // Display c as hiBin without changing the fit coordinate or probabilities.
  std::vector<double> edges;
  for (int i=1; i<=rebinned->GetNbinsX()+1; ++i)
    edges.push_back(2*rebinned->GetXaxis()->GetBinLowEdge(i));
  auto* result=new TH1D(Unique("display_hiBin"),"",rebinned->GetNbinsX(),edges.data());
  result->SetDirectory(nullptr); result->Sumw2();
  for (int i=1; i<=rebinned->GetNbinsX(); ++i) {
    result->SetBinContent(i,rebinned->GetBinContent(i));
    result->SetBinError(i,rebinned->GetBinError(i));
  }
  delete rebinned;
  return result;
}

void StyleHistogram(TH1& h, int color, int marker) {
  h.SetLineColor(color); h.SetMarkerColor(color); h.SetMarkerStyle(marker);
  h.SetMarkerSize(0.85); h.SetLineWidth(1); h.SetStats(false);
}

void StyleGraph(TGraphErrors& g, int color, int marker) {
  g.SetLineColor(color); g.SetMarkerColor(color); g.SetMarkerStyle(marker);
  g.SetMarkerSize(0.85); g.SetLineWidth(1);
}

void SaveCanvas(TCanvas& c, const TString& directory, const TString& stem) {
  c.Modified(); c.Update();
  c.SaveAs(directory+"/"+stem+".pdf");
  c.SaveAs(directory+"/"+stem+".png");
  c.Write(stem);
}

TCanvas* DrawComparison(const Distribution& d, const TH1D& after, bool refit,
                       int rebin, const TString& outputDir, std::ostream& log) {
  TH1D* hd=ForDisplay(*d.data,d.centrality,rebin);
  TH1D* hb=ForDisplay(*d.before,d.centrality,rebin);
  TH1D* ha=ForDisplay(after,d.centrality,rebin);
  StyleHistogram(*hd,kRed+1,21); StyleHistogram(*hb,kBlue,20); StyleHistogram(*ha,kGreen+2,22);
  double ymin=std::numeric_limits<double>::max(),ymax=0;
  for (TH1D* h : {hd,hb,ha})
    for (int i=1; i<=h->GetNbinsX(); ++i) {
      const double value=h->GetBinContent(i);
      if (value>0) ymin=std::min(ymin,value);
      ymax=std::max(ymax,value+h->GetBinError(i));
    }
  hd->SetMinimum(0.45*ymin); hd->SetMaximum(ymax*(d.centrality ? 12.0 : 3.5));

  // Use the SAME header-only drawer as the existing EEC macros.
  // Keep the drawer and plot objects alive for interactive ROOT sessions.
  auto* drawer=new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetCanvasSize(700,650);
  drawer->SetSplitRatio(0.32);
  drawer->SetMargins(0.17,0.045,0.09,0.12);
  drawer->CreateSplitCanvas();
  drawer->SetLogY(true);
  const char* xtitle=d.centrality ? "hiBin" : "v_{z} (cm)";
  drawer->DrawHistogramToUpperPad(hd,xtitle,"Normalized distribution","","E1");
  hd->GetXaxis()->SetLabelSize(0); hd->GetXaxis()->SetTitleSize(0);
  hb->Draw("E1 SAME"); ha->Draw("E1 SAME");
  auto* legend=new TLegend(d.centrality ? 0.46 : 0.22,
                          d.centrality ? 0.63 : 0.14,
                          d.centrality ? 0.93 : 0.79,
                          d.centrality ? 0.83 : 0.34);
  legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetTextFont(42); legend->SetTextSize(0.041);
  const TString mcLabel=d.system=="OO" ? kOOLabel : kPPLabel;
  legend->AddEntry(hd,d.system=="OO" ? "OO Data" : "pp Data","lp");
  legend->AddEntry(hb,mcLabel+" (before reweight)","lp");
  legend->AddEntry(ha,mcLabel+(refit ? " (after refit)" : " (stored weights)"),"lp");
  legend->Draw();
  auto* label=new TLatex(); label->SetNDC(); label->SetTextFont(42); label->SetTextSize(0.040);
  label->DrawLatex(0.17,0.925,"#bf{CMS} #it{Preliminary}");
  label->SetTextAlign(31); label->DrawLatex(0.955,0.925,"5.36 TeV");
  if (refit) {
    label->SetTextAlign(13); label->SetTextSize(0.031);
    label->DrawLatex(0.20,0.83,"1D refit check");
  }
  drawer->GetUpperPad()->RedrawAxis();

  auto* ratioBefore=RatioGraph(*hb,*hd);
  auto* ratioAfter=RatioGraph(*ha,*hd);
  StyleGraph(*ratioBefore,kBlue,20); StyleGraph(*ratioAfter,kGreen+2,22);
  TH1D* frame=Copy(*hd,Unique("ratio_frame")); frame->Reset();
  frame->SetMinimum(kRatioMin); frame->SetMaximum(kRatioMax);
  drawer->SetLogY(false);
  drawer->DrawHistogramToLowerPad(frame,xtitle,"#frac{MC}{Data}","","AXIS");
  frame->GetYaxis()->SetNdivisions(505);
  const double xlow=frame->GetXaxis()->GetXmin(),xhigh=frame->GetXaxis()->GetXmax();
  auto* unity=new TLine(xlow,1,xhigh,1); unity->SetLineStyle(2); unity->SetLineColor(kGray+2); unity->Draw();
  ratioBefore->Draw("P SAME"); ratioAfter->Draw("P SAME");
  drawer->GetLowerPad()->RedrawAxis();
  int clipped=0;
  for (TGraphErrors* graph : {ratioBefore,ratioAfter})
    for (int i=0; i<graph->GetN(); ++i) {
      double x,y; graph->GetPoint(i,x,y);
      if (y<kRatioMin || y>kRatioMax) ++clipped;
    }
  if (clipped) log << "  display note: " << clipped << " ratio points outside [" << kRatioMin << ',' << kRatioMax
                   << "]; full ratios are saved in ROOT.\n";
  const TString stem=d.system+"_"+(d.centrality ? "hiBin" : "vz")+(refit ? "_refit" : "_stored");
  TCanvas* canvas=drawer->GetUpperPad()->GetCanvas();
  canvas->SetName(Unique(stem));
  SaveCanvas(*canvas,outputDir,stem);
  ratioBefore->Write(stem+"_ratio_before"); ratioAfter->Write(stem+"_ratio_after");
  return canvas;
}

void DrawOOCombined(TCanvas* centrality, TCanvas* vz, const TString& outDir, bool refit) {
  if (!centrality || !vz) return;
  auto* canvas=new TCanvas(Unique("OO_eventWeights"),"OO event weights",1400,650);
  canvas->Divide(2,1,0,0);
  canvas->cd(1); centrality->DrawClonePad();
  canvas->cd(2); vz->DrawClonePad();
  canvas->cd();
  SaveCanvas(*canvas,outDir,refit ? "OO_eventWeights_refit" : "OO_eventWeights_stored");
}

void ReportAgreement(const Distribution& d, const TH1D& mc, const char* stage, std::ostream& log) {
  double rms=0,totalVariation=0,covered=0,largest=0;
  for (int i=1; i<=d.data->GetNbinsX(); ++i) {
    const double a=d.data->GetBinContent(i),b=mc.GetBinContent(i);
    totalVariation+=std::abs(a-b)/2;
    if (a<=0) continue;
    const double residual=b/a-1;
    rms+=a*residual*residual; covered+=a; largest=std::max(largest,std::abs(residual));
  }
  log << "  " << stage << ": data-weighted relative RMS [%] = " << 100*std::sqrt(rms/covered)
      << ", max relative deviation [%] = " << 100*largest
      << ", total variation [%] = " << 100*totalVariation << '\n';
}

void WriteFit(const Distribution& d, const Fit& fit, std::ostream& parameters) {
  fit.after->Write("mc_refit_preview");
  TVectorD p(fit.p.size());
  for (int j=0; j<int(fit.p.size()); ++j) p[j]=fit.p[j];
  p.Write("scaled_basis_parameters"); fit.covariance.Write("scaled_basis_covariance_diagonal_errors");
  TParameter<double>("fit_chi2_diagonal_errors",fit.chi2).Write();
  TParameter<int>("fit_ndof",fit.ndof).Write();
  parameters << "\n// " << d.system << " / " << d.observable << ": candidate REPLACEMENT, not a residual multiplier.\n";
  const auto writeOne=[&](TF1* f,const char* key,const char* member) {
    if (!f) return;
    f->Write(key);
    parameters << "// pol6, physical x; domain [" << f->GetXmin() << ',' << f->GetXmax() << "]\n";
    parameters << member << "->SetParameters(";
    for (int j=0; j<7; ++j) parameters << (j ? ", " : "") << f->GetParameter(j);
    parameters << ");\n";
  };
  writeOne(fit.vertex,"weight_vz","fVzWeightFunction");
  writeOne(fit.central,"weight_centrality_0_30","fCentralityWeightFunctionCentral");
  writeOne(fit.peripheral,"weight_centrality_30_max","fCentralityWeightFunctionPeripheral");
  if (d.centrality)
    parameters << "// Evaluate with c = hiBin/2.0. Require 0 <= hiBin < " << 2*d.cmax
               << "; outside the fitted range use unity, not polynomial extrapolation.\n";
}

} // namespace Run3MCWeights

void deriveMonteCarloWeightsRun3(const char* ooData="", const char* ooMC="",
                                const char* ppData="", const char* ppMC="",
                                const char* outputDirectory="figures/mcWeightsRun3",
                                const char* mode="both", double maxCentrality=90.0,
                                int vertexPlotRebin=1, int centralityPlotRebin=1) {
  using namespace Run3MCWeights;
  Require(ooData && ooMC && ppData && ppMC && outputDirectory && mode,
          "Use empty strings, not null pointers, for unused input arguments.");
  if (!std::strlen(ooData) && !std::strlen(ooMC) && !std::strlen(ppData) && !std::strlen(ppMC)) {
    std::cout << "Run from plotting/:\n"
      << "  root -l -b -q 'deriveMonteCarloWeightsRun3.C(\"OO_DATA.root\",\"OO_MC.root\",\"PP_DATA.root\",\"PP_MC.root\")'\n"
      << "Modes in argument 6: stored, derive, both. Use \"\",\"\" to skip a system.\n"
      << "Inputs must be analyzer/processed event-information histograms, not HiForest files.\n";
    return;
  }
  const TString selectedMode=mode;
  Require(selectedMode=="stored" || selectedMode=="derive" || selectedMode=="both", "Mode must be stored, derive, or both.");
  const bool stored=selectedMode!="derive", derive=selectedMode!="stored";
  Require(bool(std::strlen(ooData))==bool(std::strlen(ooMC)), "Supply both OO data and OO MC, or neither.");
  Require(bool(std::strlen(ppData))==bool(std::strlen(ppMC)), "Supply both pp data and pp MC, or neither.");
  Require(maxCentrality>kCentralitySplit && maxCentrality<=90 &&
          std::abs(maxCentrality*2-std::round(maxCentrality*2))<1e-8,
          "maxCentrality must be on the 0.5% grid, above 30% and at most 90%.");
  TString outDir=outputDirectory;
  gSystem->ExpandPathName(outDir);
  Require(!outDir.IsNull(), "Output directory must not be empty.");
  gSystem->mkdir(outDir,kTRUE);
  Require(!gSystem->AccessPathName(outDir), "Cannot create output directory.");
  const TString rootName=outDir+"/run3_mc_weights.root";
  // Never RECREATE an input file by an accidental identical output path.
  struct stat outputInfo{};
  const bool outputExists=(::stat(rootName.Data(),&outputInfo)==0);
  for (const char* input : {ooData,ooMC,ppData,ppMC}) {
    if (!std::strlen(input)) continue;
    TString path=input; gSystem->ExpandPathName(path);
    Require(path!=rootName,"Output ROOT file coincides with an input file.");
    struct stat inputInfo{};
    if (outputExists && ::stat(path.Data(),&inputInfo)==0)
      Require(inputInfo.st_dev!=outputInfo.st_dev || inputInfo.st_ino!=outputInfo.st_ino,
              "Output ROOT file aliases an input file; refusing to overwrite it.");
  }
  std::ofstream log((outDir+"/weight_audit.txt").Data());
  std::ofstream parameters((outDir+"/candidate_parameters.txt").Data());
  Require(bool(log) && bool(parameters), "Cannot open output text files.");
  log << std::setprecision(17); parameters << std::setprecision(17);
  log << "Run-3 event-weight QA and derivation; mode=" << mode << '\n'
      << "Source audited at 6d4800b8c57747e1db761fc4fecf72a38a1a00c3.\n"
      << "MC baseline = vertexZ/centrality. Stored after = vertexZweighted/centralityWeighted.\n"
      << "At that source revision Run-3 generator/pThat normalization is unity.\n"
      << "The stored weighted histograms contain the FULL event-weight product.\n"
      << "1D refits do not validate the combined product, generator normalization, or independent closure.\n"
      << "Fit/ratio errors neglect unit-area normalization correlations; refit preview errors hold fitted weights fixed.\n"
      << "Vertex refits use bin-center evaluation. Centrality refits use discrete hiBin/2.0.\n"
      << "Centrality plotting/fitting range is [0," << maxCentrality << ")%; vertex plots retain the full selected-event centrality acceptance.\n"
      << "WARNING: do not combine disjoint pThat samples using these shape weights as cross-section weights.\n\n";
  parameters << "// Candidate fits only. These lines do not modify the analyzer automatically.\n"
             << "// Apply to the matching Run-3 constructor branch only AFTER validation.\n"
             << "// Keep the generator/pThat normalization policy separate.\n"
             << "// Outside |vz| <= 15 cm or the stated centrality support, do not extrapolate.\n";
  std::unique_ptr<TFile> output(TFile::Open(rootName,"RECREATE"));
  Require(output && !output->IsZombie(), "Cannot create output ROOT file.");
  TNamed("mode",mode).Write();
  TNamed("scope","stored=production product; refit=same-sample single-marginal preview, not event-level closure").Write();
  TParameter<int>("complete",0).Write();

  for (int system=0; system<2; ++system) {
    const char* dataName=system==0 ? ooData : ppData;
    const char* mcName=system==0 ? ooMC : ppMC;
    if (!std::strlen(dataName)) continue;
    const TString label=system==0 ? "OO" : "pp";
    auto dataFile=OpenInput(dataName),mcFile=OpenInput(mcName);
    AuditCard(*dataFile,system==0 ? 9 : 11,log);
    AuditCard(*mcFile,system==0 ? 10 : 12,log);
    output->cd();
    TDirectory* sysDir=output->mkdir(label);
    Require(sysDir,"Cannot create system output directory."); sysDir->cd();
    TNamed("data_file",dataName).Write(); TNamed("mc_file",mcName).Write();
    TCanvas *storedVz=nullptr,*storedCent=nullptr,*refitVz=nullptr,*refitCent=nullptr;
    for (int observable=0; observable<(system==0 ? 2 : 1); ++observable) {
      const bool centrality=observable==1;
      Distribution d=LoadDistribution(*dataFile,*mcFile,label,centrality,maxCentrality,stored,log);
      sysDir->cd();
      TDirectory* obsDir=sysDir->mkdir(d.observable); obsDir->cd();
      d.data->Write("data_normalized"); d.before->Write("mc_before_normalized");
      if (d.stored) d.stored->Write("mc_stored_normalized");
      RatioGraph(*d.data,*d.before,centrality)->Write("data_over_mc_for_fit");
      ReportAgreement(d,*d.before,"before",log);
      const int rebin=centrality ? centralityPlotRebin : vertexPlotRebin;
      if (stored) {
        ReportAgreement(d,*d.stored,"stored production weights",log);
        TCanvas* c=DrawComparison(d,*d.stored,false,rebin,outDir,log);
        if (centrality) storedCent=c; else storedVz=c;
      }
      if (derive) {
        Fit fit=Derive(d,log);
        WriteFit(d,fit,parameters);
        ReportAgreement(d,*fit.after,"single-marginal refit preview",log);
        TCanvas* c=DrawComparison(d,*fit.after,true,rebin,outDir,log);
        if (centrality) refitCent=c; else refitVz=c;
      }
    }
    sysDir->cd();
    if (system==0) {
      if (stored) DrawOOCombined(storedCent,storedVz,outDir,false);
      if (derive) DrawOOCombined(refitCent,refitVz,outDir,true);
    }
  }
  output->cd();
  TParameter<int>("complete",1).Write("complete",TObject::kOverwrite);
  output->Write(); output->Close();
  log << "\nCOMPLETE. No analyzer configuration or input files were changed.\n";
  std::cout << "Run-3 weight outputs: " << outDir << '\n'
            << "  run3_mc_weights.root, weight_audit.txt, candidate_parameters.txt, PDF/PNG figures\n"
            << "  Stored plots show the input file's weights; refit plots are one-dimensional previews.\n";
}

#endif // EEC_DERIVE_MONTE_CARLO_WEIGHTS_RUN3_C

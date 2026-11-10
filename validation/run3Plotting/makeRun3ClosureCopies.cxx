// Reproducible, bounded copy-and-patch of the two legacy macros. No Python.
// Inputs are unchanged. Each replacement must match the expected source.
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
namespace fs=std::filesystem;
std::string read(const fs::path& p){std::ifstream f(p);if(!f)throw std::runtime_error("Cannot read "+p.string());return {std::istreambuf_iterator<char>(f),{}};}
void write(const fs::path&p,const std::string&s){std::ofstream f(p);if(!f)throw std::runtime_error("Cannot write "+p.string());f<<s;}
void replace(std::string&s,const std::string&a,const std::string&b,int expected=1){size_t pos=0;int n=0;while((pos=s.find(a,pos))!=s.npos){s.replace(pos,a.size(),b);pos+=b.size();++n;}if(n!=expected)throw std::runtime_error("Unexpected matches "+std::to_string(n)+" for "+a.substr(0,70));}
void block(std::string&s,const std::string&a,const std::string&z,const std::string&b){auto i=s.find(a),j=s.find(z,i);if(i==s.npos||j==s.npos)throw std::runtime_error("Missing block "+a);s.replace(i,j-i,b);}
int main(int argc,char**argv){try{
 const fs::path base=argc>1?argv[1]:".";const fs::path p=base/"plotting";
 std::string s=read(p/"compareDifferentEECpairings.C");
 replace(s,"#include \"../src/EECHistograms.h\"","#include \"../src/EECHistograms.h\"\n#include \"Run3ClosureTools.h\"");
 replace(s,"void compareDifferentEECpairings(){",R"PATCH(// Copy of compareDifferentEECpairings.C at 792d76c0; original remains unchanged.
// plotKind: signal, background, signalWithTotal, or methods. The default is
// raw-scale closure (the legacy closure presets), NOT normalized pPb shapes.
// averageSameCone=false preserves the actual legacy C1C1 arithmetic. Set true
// explicitly for the two-cone-averaged definition; never mix definitions.
void compareDifferentEECpairingsRun3(
  TString fileName="merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root",
  TString plotKind="signal", bool correlatedRatio=false,
  double minJetPt=80., double maxJetPt=120.,
  double minCentrality=0., double maxCentrality=30., double trackThreshold=2.,
  TString outputDirectory="figures/run3LegacyClosure", bool averageSameCone=false,
  bool normalizeShapes=false, double ratioMinimum=0.8, double ratioMaximum=1.2,
  bool extendRatioRange=false){
  using Run3Closure::Require;
  Require(plotKind=="signal" || plotKind=="background" || plotKind=="signalWithTotal" || plotKind=="methods","Unsupported plotKind");
  Require(!correlatedRatio || ((plotKind=="signal" || plotKind=="signalWithTotal") && !normalizeShapes),
    "Shared-signal approximation is only defined for raw-scale SIGNAL closure, not background/normalized ratios");
  Require(minJetPt<maxJetPt && minCentrality<maxCentrality && ratioMinimum<ratioMaximum,"Invalid selection/range");
  Require(gSystem->mkdir(outputDirectory,true)==0 || !gSystem->AccessPathName(outputDirectory),"Cannot create output directory");
)PATCH");
 block(s,"  // Studied file","  EECCard* card",R"PATCH(  // Explicit Run-3 file; public-figure presets must not override the input.
  TFile* inputFile=Run3Closure::Open(fileName);

)PATCH");
 replace(s,"  bool isPbPb = card->GetDataType().Contains(\"PbPb\");",R"PATCH(  Run3Closure::CheckCard(*card);
  const bool hasCentrality=true; // OO bins, no PbPb MC centrality shift.)PATCH");
 replace(s,"isPbPb","hasCentrality",13);
 block(s,"  // Define the bins that are compared","  // Enumeration for different pairing types",R"PATCH(  // Native OO selections are resolved by physical boundaries, never by PbPb indices.
  std::vector<std::pair<double,double>> comparedCentralityBin{{minCentrality,maxCentrality}};
  std::vector<std::pair<double,double>> comparedJetPtBin{{minJetPt,maxJetPt}};
  std::vector<double> comparedTrackPtBin{trackThreshold};
  int adjustCentralityPercentage=0;
  TString systemString="Pythia+Hijing OO 5.36 TeV";

)PATCH");
 replace(s,"manualLegend.push_back(\"Mixed cone, plain HFSum\");","manualLegend.push_back(\"Mixed cone signal\");");
 replace(s,"bool normalizeDistributions = true;","bool normalizeDistributions = normalizeShapes;");
 block(s,"  // If we are dealing with MC, shift the centrality","  // Override the selected configuration", "  // OO uses the stored centrality directly; never apply the PbPb +4% shift.\n\n");
 block(s,"  // Override the selected configuration","  // ====================================================\n  //                Drawing configuration",R"PATCH(  // Retain legacy local pairing enums and rendering; select the Run-3 closure.
  if(plotKind=="background") {
    comparedEnergyEnergyCorrelatorPairings={{kSameJetPair,kAllBackground},
      {kCompiledBackgroundMixed,knSubeventCombinations},{kCompiledBackgroundPerpendicularCone,knSubeventCombinations}};
    manualLegend={"True background","Mixed cone estimate","Perpendicular cone estimate"};
    legendPosition=kPaperClosure;
  } else if(plotKind=="methods") {
    comparedEnergyEnergyCorrelatorPairings={{kCompiledBackgroundPerpendicularCone,knSubeventCombinations},
      {kCompiledBackgroundMixed,knSubeventCombinations}};
    manualLegend={"Perpendicular cone background","Mixed cone background"};
    legendPosition=kPaperClosure;
  } else {
    legendPosition=kBottomLeft;
    if(plotKind=="signalWithTotal") {
      comparedEnergyEnergyCorrelatorPairings={{kSameJetPair,kPythiaPythia},
        {kSameJetPair,knSubeventCombinations},{kMixedBackgroundSubtracted,knSubeventCombinations},
        {kPerpendicularBackgroundSubtracted,knSubeventCombinations}};
      manualLegend={"True signal","Signal + background","Mixed cone signal","Perpendicular cone signal"};
    }
  }

)PATCH");
 replace(s,"TString saveComment = \"_plainHFMatch\";","TString saveComment = \"_Run3_\"+plotKind+(correlatedRatio?\"_sharedSignalApprox\":\"_uncorrelated\");");
 replace(s,"const char* figureFormat = \"png\";","const char* figureFormat = \"pdf\";");
 replace(s,"std::pair<double, double> ratioZoom = std::make_pair(0.85, 1.15);","std::pair<double, double> ratioZoom = std::make_pair(ratioMinimum, ratioMaximum);");
 replace(s,"if(card->GetBinIndexTrackPtEEC(trackPtBin) < 0){","if(Run3Closure::Track(*card,trackPtBin) < 0){");
 const std::string manager="  EECHistogramManager* histograms = new EECHistogramManager(inputFile, card);";
 replace(s,manager,manager+R"PATCH(
  for(const auto& c:comparedCentralityBin) for(const auto& j:comparedJetPtBin) for(double t:comparedTrackPtBin)
    Run3Closure::Preflight(*histograms,card->FindBinIndexCentrality(c),card->FindBinIndexJetPtEEC(j),Run3Closure::Track(*card,t));
)PATCH");
 const std::string load="hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = histograms->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType, iSubevent);";
 replace(s,load,load+R"PATCH(
            if(auto* original=hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent])
              hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent]=
                Run3Closure::Clone(original,Form("run3_input_C%dJ%dT%dP%dS%d",iCentrality,iJetPt,iTrackPt,iPairingType,iSubevent));
)PATCH");
 const std::string compiled="hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]);";
 replace(s,compiled,"// Default reproduces legacy arithmetic; optional averaged C1C1/C2C2 is explicit.\n            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(averageSameCone ? helperHistogram : hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]);");
 replace(s,"backgroundHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]);","backgroundHistogram->Add(averageSameCone ? helperHistogram : hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]);");
 const std::string div="hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Divide(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][comparedEnergyEnergyCorrelatorPairings.at(0).first][comparedEnergyEnergyCorrelatorPairings.at(0).second]);";
 replace(s,div,R"PATCH(auto* denominator=hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][comparedEnergyEnergyCorrelatorPairings.at(0).first][comparedEnergyEnergyCorrelatorPairings.at(0).second];
          Run3Closure::PositiveDenominator(denominator,drawingRange.first,drawingRange.second);
          )PATCH"+div+R"PATCH(
          if(correlatedRatio && pairingType!=comparedEnergyEnergyCorrelatorPairings.front()) {
            std::unique_ptr<TH1D> trueBackground(Run3Closure::TrueBackground(*histograms,iCentrality,iJetPt,iTrackPt));
            Run3Closure::ApplySharedSignal(
              hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second],
              hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second],denominator,
              hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSameJetPair][knSubeventCombinations],trueBackground.get());
          }
)PATCH");
 replace(s,"// Create a new canvas for the plot\n        drawer->CreateSplitCanvas();",R"PATCH(// Optional range extension is explicit: the unmodified legacy zoom remains default.
        if(extendRatioRange) {
          double low=ratioMinimum,high=ratioMaximum;
          for(size_t curve=1;curve<comparedEnergyEnergyCorrelatorPairings.size();++curve) {
            auto pairing=comparedEnergyEnergyCorrelatorPairings[curve];
            auto* h=hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairing.first][pairing.second];
            for(int b=1;b<=h->GetNbinsX();++b) {
              if(h->GetXaxis()->GetBinUpEdge(b)<=drawingRange.first || h->GetXaxis()->GetBinLowEdge(b)>=drawingRange.second) continue;
              low=std::min(low,h->GetBinContent(b)-h->GetBinError(b));
              high=std::max(high,h->GetBinContent(b)+h->GetBinError(b));
            }
          }
          const double pad=.08*(high-low);ratioZoom={low-pad,high+pad};
        }
        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();)PATCH");
 const std::string legend="legend->AddEntry((TObject*) 0, Form(\"%s %.0f-%.0f%%\", systemString.Data(), centralityBin.first-4*adjustCentralityPercentage, centralityBin.second-4*adjustCentralityPercentage), \"\");";
 replace(s,legend,R"PATCH(legend->AddEntry((TObject*)0,systemString,"");
          legend->AddEntry((TObject*)0,Form("Cent: %.0f-%.0f%%, %s",centralityBin.first,centralityBin.second,Run3Closure::Mode(*card).Data()),"");)PATCH");
 replace(s,"legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.06);legend->SetTextFont(62);","legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.048);legend->SetTextFont(62);");
 replace(s,"Form(\"p_{T}^{ch} > %.0f GeV, %s\", trackPtBin, energyWeightString.Data())","Form(\"p_{T}^{ch} > %.1f GeV, %s\", trackPtBin, energyWeightString.Data())");
 const std::string save="gPad->GetCanvas()->SaveAs(Form(\"figures/eecPairingTypeComparison%s%s%s%s.%s\", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));";
 replace(s,save,R"PATCH(TCanvas* canvas=gPad->GetCanvas();
          Run3Closure::Header(canvas,correlatedRatio,averageSameCone,plotKind);
          TString stem=outputDirectory+"/"+Form("C%.0f-%.0f_J%.0f-%.0f_T%.1f_n%d_",centralityBin.first,centralityBin.second,jetPtBin.first,jetPtBin.second,trackPtBin,weightExponent)+plotKind;
          stem+=(averageSameCone?"_averagedCones":"_legacyC1");
          stem+=(normalizeDistributions?"_unitArea":"_raw");
          stem+=(correlatedRatio?"_sharedSignalApprox":"_uncorrelated");
          std::vector<TH1D*> curves,ratios;
          for(auto pairing:comparedEnergyEnergyCorrelatorPairings) {
            curves.push_back(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairing.first][pairing.second]);
            ratios.push_back(hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairing.first][pairing.second]);
          }
          Run3Closure::Save(canvas,stem,curves,ratios,manualLegend,fileName,correlatedRatio,averageSameCone,normalizeDistributions);
          std::cout<<"RUN3_CLOSURE_SAVED "<<stem<<std::endl;)PATCH");
 write(p/"compareDifferentEECpairingsRun3.C",s);

 // The second copy retains its independent shape normalization by default.
 s=read(p/"validateBackgroundSubtraction.C");
 replace(s,"#include \"JDrawer.h\"","#include \"JDrawer.h\"\n#include \"Run3ClosureTools.h\"");
 replace(s,"void validateBackgroundSubtraction(){",R"PATCH(// Copy of validateBackgroundSubtraction.C at 792d76c0; original untouched.
// normalizeShapes=true retains the original normalized-shape comparison.
// For raw closure set normalizeShapes=false; only then may shared-signal
// approximation be requested. backgroundMethod: 0 mixed, 1 perpendicular,
// -1 read an already-processed signal (no shared-signal mode without provenance).
void validateBackgroundSubtractionRun3(
  TString inputFileName="merged_eecAnalysis_OOMC_GenGen_Csplit_processed.root",
  int backgroundMethod=0, bool correlatedRatio=false, bool normalizeShapes=true,
  double minJetPt=80.,double maxJetPt=120.,double minCentrality=0.,double maxCentrality=30.,
  double trackThreshold=2.,TString outputDirectory="figures/run3LegacyValidation",bool averageSameCone=false){
  using Run3Closure::Require;
  Require(backgroundMethod>=-1 && backgroundMethod<=1,"Invalid background method");
  Require(!correlatedRatio || (!normalizeShapes && backgroundMethod>=0),"Shared-signal mode needs raw scale and explicitly reconstructed background");
  Require(gSystem->mkdir(outputDirectory,true)==0 || !gSystem->AccessPathName(outputDirectory),"Cannot create output directory");
)PATCH");
 block(s,"  // Open the input file","  EECCard* card", "  TFile* inputFile=Run3Closure::Open(inputFileName);\n\n");
 replace(s,"  const int weightExponent = card->GetWeightExponent();","  Run3Closure::CheckCard(*card);\n  const int weightExponent = card->GetWeightExponent();");
 block(s,"  // Select which histograms are fitted","  // Select the types",R"PATCH(  int firstDrawnCentralityBin=card->FindBinIndexCentrality(minCentrality,maxCentrality);
  int lastDrawnCentralityBin=firstDrawnCentralityBin;
  int firstDrawnJetPtBinEEC=card->FindBinIndexJetPtEEC(minJetPt,maxJetPt);
  int lastDrawnJetPtBinEEC=firstDrawnJetPtBinEEC;
  int firstDrawnTrackPtBinEEC=Run3Closure::Track(*card,trackThreshold);
  int lastDrawnTrackPtBinEEC=firstDrawnTrackPtBinEEC;
  Require(firstDrawnCentralityBin>=0 && firstDrawnJetPtBinEEC>=0,"Requested OO bin is missing");

)PATCH");
 replace(s,"const bool saveFigures = false;","const bool saveFigures = true;");
 replace(s,"histograms = new EECHistogramManager(inputFile,card);",R"PATCH(histograms = new EECHistogramManager(inputFile,card);
  if(backgroundMethod>=0) Run3Closure::Preflight(*histograms,firstDrawnCentralityBin,firstDrawnJetPtBinEEC,firstDrawnTrackPtBinEEC,backgroundMethod);)PATCH");
 std::string lhs="hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]";
 std::string truth=lhs+"[0] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaPythia);";
 replace(s,truth,truth+"\n          "+lhs+"[0]=Run3Closure::Clone("+lhs+"[0],\"run3_validation_truth\");");
 std::string extract=lhs+"[1] = histograms->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorSignal);";
 replace(s,extract,R"PATCH(if(backgroundMethod<0) {
            )PATCH"+extract+"\n            "+lhs+"[1]=Run3Closure::Clone("+lhs+"[1],\"run3_validation_extracted\");\n          } else {\n            "+lhs+R"PATCH([1]=Run3Closure::Read(*histograms,iCentrality,iJetPt,iTrackPt,EECHistograms::kSameJetPair,EECHistograms::knSubeventCombinations,"run3_validation_extracted");
            std::unique_ptr<TH1D> background(Run3Closure::Background(*histograms,iCentrality,iJetPt,iTrackPt,backgroundMethod,averageSameCone));
            )PATCH"+lhs+"[1]->Add(background.get(),-1.);\n          }");
 replace(s,"          for(int iSignalType = 0; iSignalType < 2; iSignalType++){\n            "+lhs+"[iSignalType]->Scale", "          if(normalizeShapes) for(int iSignalType = 0; iSignalType < 2; iSignalType++){\n            Require("+lhs+"[iSignalType]->Integral(lowSignalBin,highSignalBin,\"width\")>0.,\"Nonpositive shape normalization\");\n            "+lhs+"[iSignalType]->Scale");
 std::string rlhs="hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]";
 std::string divide=rlhs+"->Divide("+lhs+"[0]);";
 replace(s,divide,"Run3Closure::PositiveDenominator("+lhs+"[0],drawingRange.first,drawingRange.second);\n          "+divide+R"PATCH(
          if(correlatedRatio) {
            std::unique_ptr<TH1D> total(Run3Closure::Read(*histograms,iCentrality,iJetPt,iTrackPt,EECHistograms::kSameJetPair,EECHistograms::knSubeventCombinations,"run3_validation_total"));
            std::unique_ptr<TH1D> truthBackground(Run3Closure::TrueBackground(*histograms,iCentrality,iJetPt,iTrackPt));
            Run3Closure::ApplySharedSignal()PATCH"+rlhs+","+lhs+"[1],"+lhs+"[0],total.get(),truthBackground.get());\n          }");
 replace(s,"Form(\"%s 5.02 TeV\",histograms->GetCard()->GetAlternativeDataType(false).Data())","Form(\"Pythia+Hijing OO 5.36 TeV, %s\",Run3Closure::Mode(*card).Data())");
 std::string save2="gPad->GetCanvas()->SaveAs(Form(\"figures/%sSignalValidityCheck%s%s%s%s.%s\", histograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));";
 replace(s,save2,R"PATCH(TCanvas* canvas=gPad->GetCanvas();
            Run3Closure::Header(canvas,correlatedRatio,averageSameCone,normalizeShapes?"Unit-area shape closure":"Raw-scale signal closure");
            TString stem=outputDirectory+"/"+Form("C%.0f-%.0f_J%.0f-%.0f_T%.1f_n%d_method%d",minCentrality,maxCentrality,minJetPt,maxJetPt,trackThreshold,weightExponent,backgroundMethod);
            stem+=(normalizeShapes?"_unitArea":"_raw");stem+=(averageSameCone?"_averagedCones":"_legacyC1");
            stem+=(correlatedRatio?"_sharedSignalApprox":"_uncorrelated");
            auto* truth=)PATCH"+lhs+"[0];\n            auto* extracted="+lhs+"[1];\n            auto* ratio="+rlhs+R"PATCH(;
            std::unique_ptr<TH1D> unity(Run3Closure::Clone(truth,"truth_reference_ratio")); unity->Divide(truth);
            Run3Closure::Save(canvas,stem,{truth,extracted},{unity.get(),ratio},{"True signal","Extracted signal"},inputFileName,correlatedRatio,averageSameCone,normalizeShapes);)PATCH");
 write(p/"validateBackgroundSubtractionRun3.C",s);
 std::cout<<"Created Run-3 copies; legacy sources unchanged.\n";
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}

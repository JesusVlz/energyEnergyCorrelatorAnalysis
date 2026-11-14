// Native response-input diagnostics only: does not execute unfolding.
#include "Run3ResponsePlotTools.h"
R__LOAD_LIBRARY(plotting/DrawingClasses.so)
void plotUnfoldingInputsRun3(const char* input,const char* output="figures/run3UnfoldingInputs",double centralityLow=0.,double centralityHigh=30.,double trackThreshold=2.){
  try {gROOT->SetBatch(true);gStyle->SetOptStat(0);Run3ResponsePlot::UnfoldingPlots(input,output,centralityLow,centralityHigh,trackThreshold);}
  catch(const std::exception& e){std::cerr<<"RUN3_RESPONSE_PLOTTING_FAIL "<<e.what()<<'\n';gSystem->Exit(1);}
}

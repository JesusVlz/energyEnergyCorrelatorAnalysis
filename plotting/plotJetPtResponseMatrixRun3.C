// Run-3 adaptation of plotJetPtResponseMatrix.C; legacy file unchanged.
#include "Run3ResponsePlotTools.h"
R__LOAD_LIBRARY(plotting/DrawingClasses.so)
void plotJetPtResponseMatrixRun3(const char* input,const char* output="figures/run3JetResponse",double centralityLow=0.,double centralityHigh=30.){
  try {gROOT->SetBatch(true);gStyle->SetOptStat(0);Run3ResponsePlot::JetPlots(input,output,centralityLow,centralityHigh);}
  catch(const std::exception& e){std::cerr<<"RUN3_RESPONSE_PLOTTING_FAIL "<<e.what()<<'\n';gSystem->Exit(1);}
}

// Run-3 adaptation of trackParticleMatchingPlotter.C; legacy file unchanged.
#include "Run3ResponsePlotTools.h"
R__LOAD_LIBRARY(plotting/DrawingClasses.so)
void trackParticleMatchingPlotterRun3(const char* input,const char* output="figures/run3PairResponse",double centralityLow=0.,double centralityHigh=30.,double trackThreshold=2.){
  try {gROOT->SetBatch(true);gStyle->SetOptStat(0);Run3ResponsePlot::PairPlots(input,output,centralityLow,centralityHigh,trackThreshold);}
  catch(const std::exception& e){std::cerr<<"RUN3_RESPONSE_PLOTTING_FAIL "<<e.what()<<'\n';gSystem->Exit(1);}
}

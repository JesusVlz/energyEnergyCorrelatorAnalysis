// Reuse the plotEEChistograms copy so both jet-comparison PDFs have the same style.
#include "plotEEChistogramsJetPtComparison.C"

void plotOOMCJetPtComparison(
    const char* inputName = "projected/OOMC_GenGen_withSubevents.root",
    const char* outputDirectory = "figures", double trackThreshold = 2., double minCentrality = 0., double maxCentrality = 30.){
  plotEEChistogramsJetPtComparison(inputName,trackThreshold,outputDirectory,minCentrality,maxCentrality);
}

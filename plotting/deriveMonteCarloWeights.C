// Run-3 entry point. Original Run-2 implementation: deriveMonteCarloWeightsRun2.C.
#include "deriveMonteCarloWeightsRun3.C"

void deriveMonteCarloWeights(const char* ooData="", const char* ooMC="",
                            const char* ppData="", const char* ppMC="",
                            const char* outputDirectory="figures/mcWeightsRun3",
                            const char* mode="both", double maxCentrality=90.0,
                            int vertexPlotRebin=1, int centralityPlotRebin=1) {
  deriveMonteCarloWeightsRun3(ooData,ooMC,ppData,ppMC,outputDirectory,mode,
                             maxCentrality,vertexPlotRebin,centralityPlotRebin);
}

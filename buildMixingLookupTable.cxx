// Standalone utility to prebuild a lightweight donor metadata lookup table (.lookup.tsv)
// for event-mixing pools so analysis jobs open zero donor ROOT files at startup
// and open only the sharded/matched donor files once for track prefetch.

#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TStopwatch.h>
#include <TString.h>

#include "src/ConfigurationCard.h"
#include "src/EECAnalyzer.h"

int main(int argc, char** argv) {
  if (argc < 3 || argc > 6) {
    std::cout << "Usage:\n"
              << "  " << argv[0]
              << " <mixingFileList> <configurationCard> [outputLookupFile] [localMixing=0] [maxPoolEvents=0]\n\n"
              << "Defaults:\n"
              << "  outputLookupFile : <mixingFileList>.lookup.tsv\n"
              << "  localMixing      : 0 (preserve paths/URLs; 2 strips root://<redirector>/ for local mounts)\n"
              << "  maxPoolEvents    : 0 (all events in mixingFileList; or OOClosureMaxPoolEvents if enabled)\n";
    return 1;
  }

  const std::string mixingListPath = argv[1];
  const std::string cardPath = argv[2];
  const std::string outputLookupPath =
    (argc >= 4 && std::string(argv[3]) != "") ? std::string(argv[3]) : (mixingListPath + ".lookup.tsv");
  const int localMixingMode = (argc >= 5) ? std::atoi(argv[4]) : 0;
  const int maxPoolEvents = (argc >= 6) ? std::atoi(argv[5]) : 0;

  try {
    TStopwatch timer;
    timer.Start();

    ConfigurationCard card(cardPath.c_str());
    const int runLocalArg = (localMixingMode >= 2) ? 2 : 1;
    EECAnalyzer analyzer(std::vector<TString>{}, &card, runLocalArg, 1, 0);
    analyzer.BuildAndWriteMixingLookupTable(mixingListPath, outputLookupPath, maxPoolEvents);

    timer.Stop();
    std::cout << "Lookup table build completed in:\n";
    timer.Print();
  } catch (const std::exception& err) {
    std::cerr << "Error building mixing lookup table: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

// C++ includes
#include <iostream>   // Input/output stream. Needed for cout.
#include <fstream>    // File stream for input/output to/from files
#include <stdlib.h>   // Standard utility libraries
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <vector>     // C++ vector class
#include <sstream>    // Libraries for checking boolean input
#include <string>     // Libraries for checking boolean input
#include <iomanip>    // Libraries for checking boolean input
#include <algorithm>  // Libraries for checking boolean input
#include <cctype>     // Libraries for checking boolean input

// Includes from Root
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TStopwatch.h>

// Own includes
#include "src/EECAnalyzer.h"
#include "src/ConfigurationCard.h"
#include "src/EECHistograms.h"

using namespace std;

/*
 * File list reader
 *
 *  Arguments:
 *    std::vector<TString> &fileNameVector = Vector filled with filenames found in the file
 *    TString fileNameFile = Text file containing one analysis file name in each line
 *    int debug = Level of debug messages shown
 *    int locationIndex = Where to find analysis files: 0 = Purdue EOS, 1 = CERN EOS, 2 = Vanderbilt T2, 3 = FNAL XRootD redirector, 4 = Nebraska T2
 *    bool runLocal = True: Local run mode. False: Crab run mode
 */
void ReadFileList(std::vector<TString> &fileNameVector, TString fileNameFile, int debug, int locationIndex, bool runLocal)
{
  
  // Possible location for the input files
  const char* fileLocation[] = {
    "root://eos.cms.rcac.purdue.edu/",
    "root://eoscms.cern.ch/",
    "root://xrootd-vanderbilt.sites.opensciencegrid.org/",
    "root://cmsxrootd.fnal.gov/",
    "root://xrootd-local.unl.edu:1094/"
  };
  const int numberOfFileLocations = sizeof(fileLocation) / sizeof(fileLocation[0]);
  // Local lists do not use this index; retain the negative ACCRE convention.
  if(!runLocal && (locationIndex < 0 || locationIndex >= numberOfFileLocations)){
    std::cerr << "Error: invalid file location index " << locationIndex
              << "; expected 0.." << numberOfFileLocations - 1 << " for CRAB running." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  // Set up the file names file for reading
  ifstream file_stream(fileNameFile);
  std::string line;
  fileNameVector.clear();
  if( debug > 0 ) std::cout << "Open file " << fileNameFile.Data() << " to extract files to run over" << std::endl;
  
  // Open the file names file for reading
  if( file_stream.is_open() ) {
    if( debug > 0) std::cout << "Opened " << fileNameFile.Data() << " for reading" << std::endl;
    int lineNumber = 0;
    
    // Loop over the lines in the file
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      if( debug > 0) std::cout << lineNumber << ": " << line << std::endl;
      TString lineString(line);
      
      // Put all non-empty lines to file names vector
      if( lineString.CompareTo("", TString::kExact) != 0 ) {
        
        if(runLocal){
          // For local running, it is assumed that the file name is directly the centents of the line
          fileNameVector.push_back(lineString);
          
        } else {
          // For crab running, the line will have format ["file1", "file2", ... , "fileN"]
          TObjArray* fileNameArray = lineString.Tokenize(" ");  // Tokenize the string from every ' ' character
          int numberOfFiles = fileNameArray->GetEntries();
          TObjString* currentFileNameObject;
          TString currentFileName;
          for(int i = 0; i < numberOfFiles; i++){   // Loop over all the files in the array
            currentFileNameObject = (TObjString *)fileNameArray->At(i);
            currentFileName = currentFileNameObject->String();
            
            // Strip unwanted characters
            currentFileName.Remove(TString::kBoth,'['); // Remove possible parantheses
            currentFileName.Remove(TString::kBoth,']'); // Remove possible parantheses
            currentFileName.Remove(TString::kBoth,','); // Remove commas
            currentFileName.Remove(TString::kBoth,'"'); // Remove quotation marks
            
            // Keep explicit URLs unchanged; ROOT selects the appropriate protocol.
            // Bare file names still use the selected XRootD endpoint.
            if(!currentFileName.Contains("://")){
              currentFileName.Prepend(fileLocation[locationIndex]);
            }
            fileNameVector.push_back(currentFileName);
          }
        }
        
      } // Empty line if
      
      
      lineNumber++;
    } // Loop over lines in the file
    
  // If cannot read the file, give error and end program
  } else {
    std::cout << "Error, could not open " << fileNameFile.Data() << " for reading" << std::endl;
    assert(0);
  }
}

/*
 *  Convert string to boolean value
 */
bool checkBool(string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b;
  is >> std::boolalpha >> b;
  return b;
}

/*
 *  Main program
 *
 *  Command line arguments:
 *  argv[1] = List of files to be analyzed, given in text file. For crab analysis a job ID instead.
 *  argv[2] = Card file with binning and cut information for the analysis
 *  argv[3] = .root file to which the histograms are written
 *  argv[4] = Index for the EOS location from where the input files are searched
 *  argv[5] = True: Search input files from local machine. False (default): Search input files from grid with xrootd
 *  argc[6] = Index for the used mixing list for CRAB running
 */
int main(int argc, char **argv) {
  
  //==== Read arguments =====
  if ( argc<5 ) {
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    cout<<"+ Usage of the macro: " << endl;
    cout<<"+  "<<argv[0]<<" [fileNameFile] [configurationCard] [outputFileName] [fileLocation] <runLocal> <mixingListIndex>"<<endl;
    cout<<"+  fileNameFile: Text file containing the list of files used in the analysis. Also a single .root file can be given. For crab analysis a job id should be given here." <<endl;
    cout<<"+  configurationCard: Card file with binning and cut information for the analysis." <<endl;
    cout<<"+  outputFileName: .root file to which the histograms are written." <<endl;
    cout<<"+  fileLocation: Where to find analysis files: 0 = Purdue EOS, 1 = CERN EOS, 2 = Vanderbilt T2, 3 = FNAL XRootD redirector, 4 = Nebraska T2. -1 with runLocal=true indicates running slurm jobs on ACCRE." << endl;
    cout<<"+  Complete URLs in CRAB input lists are preserved; fileLocation prefixes only bare file names." << endl;
    cout<<"+  runLocal: True: Search input files from local machine. False (default): Search input files from grid with xrootd." << endl;
    cout<<"+  mixingListIndex: Mixing file list index for CRAB running. If -1, used list is determined from CRAB job ID (default = -1)" << endl;
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    cout << endl << endl;
    exit(1);
  }
  
  // First, check if we are supposed to run locally or on crab
  bool runLocal = false;
  if(argc >= 6) runLocal = checkBool(argv[5]);
  
  // Find the file list name depending on if we run locally or on crab
  TString fileNameFile;
  int crabJobId;
  if(runLocal){
    fileNameFile = argv[1];
    crabJobId = 1;
  } else {
    crabJobId = atoi(argv[1]);
    fileNameFile = Form("job_input_file_list_%d.txt",crabJobId);
  }
  
  // Read the other command line arguments
  const char* cardName = argv[2];
  TString outputFileName = argv[3];
  const int fileSearchIndex = atoi(argv[4]);

  int mixingListIndex = -1;
  if(argc >= 7) mixingListIndex = atoi(argv[6]);

  // Set the mixing list ID to negative crab job ID value to signal that we are using CRAB job ID to determine the mixing list.
  if(mixingListIndex == -1) mixingListIndex = crabJobId * -1;
  
  // The git hash here will be replaced by the latest commit hash by makeEECAnalysisTar.sh script
  const char* gitHash = "GITHASHHERE";
  
  // Read the card
  ConfigurationCard *configurationCard = new ConfigurationCard(cardName);
  configurationCard->SetGitHash(gitHash);
  int debugLevel = configurationCard->Get("DebugLevel");
  if(debugLevel > 0){
    configurationCard->PrintOut();
    cout << endl;
  }
  
  // Read the file names used for the analysis to a vector
  std::vector<TString> fileNameVector;
  fileNameVector.clear();

  // Check if a root file is given, add it directly to the vector. Otherwise assume we are reading a file list
  if(fileNameFile.EndsWith(".root")){
    fileNameVector.push_back(fileNameFile);
  } else {
    ReadFileList(fileNameVector,fileNameFile,debugLevel,fileSearchIndex,runLocal);
  }

  // Negative file search index combined with true flag for local run marks slurm running
  int localRunIndex = runLocal;
  if(fileSearchIndex < 0 && runLocal){
    localRunIndex = 2;
  }
  
  // Variable for histograms in the analysis
  EECHistograms* histograms;

  // Time the analysis
  TStopwatch* analysisTimer = new TStopwatch();
  analysisTimer->Start();
  
  // Run the analysis over the list of files
  EECAnalyzer* eecAnalysis = new EECAnalyzer(fileNameVector, configurationCard, localRunIndex, mixingListIndex);
  eecAnalysis->RunAnalysis();
  histograms = eecAnalysis->GetHistograms();

  // Print the time it took to run the analysis
  analysisTimer->Stop();
  std::cout << "Analysis was compeleted in: " << std::endl;
  analysisTimer->Print();
  
  // Write the histograms and card to file
  TFile* outputFile = new TFile(outputFileName, "RECREATE");
  histograms->Write();
  configurationCard->WriteCard(outputFile);
  outputFile->Close();
  
  // After writing to the file, delete all created objects
  delete analysisTimer;
  delete configurationCard;
  delete eecAnalysis;
  delete outputFile;
  
}


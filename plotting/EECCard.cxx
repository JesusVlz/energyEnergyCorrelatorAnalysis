/*
 * Implementation of the EECCard class
 */

// Own includes
#include "EECCard.h"

/*
 * Contructor with input file
 *
 *  TFile *inFile = Input file
 */
EECCard::EECCard(TFile *inFile):
  fInputFile(inFile),
  fCardDirectory("JCard"),
  fDataType(-1),
  fDataTypeString(""),
  fAlternativeDataTypeString(""),
  fBackgroundLegacyMode(false),
  fGitHash(0),
  fProjectionGitHash(0),
  fProcessGitHash(0),
  fUnfoldingGitHash(0)
{

  // Read the vectors from the input file
  fInputFile->cd(fCardDirectory.Data());
  ReadVectors();
  fDataType = (*fCardEntries[kDataType])[1];
  
  // Read the Monte Carlo type from the vector: 0 = RecoReco, 1 = RecoGen, 2 = GenReco, 3 = GenGen
  if(fCardEntries[kMcCorrelationType]) {
    fMonteCarloType = (*fCardEntries[kMcCorrelationType])[1];
  } else {
    fMonteCarloType = 0;
  }
  
  FindDataTypeString();
  
}

/*
 * Destructor
 */
EECCard::~EECCard(){

}

/*
 * Reader for all the vectors from the input file
 */
void EECCard::ReadVectors(){
  
  // Read the git hash
  fGitHash = (TObjString*) gDirectory->Get("GitHash");
  fProjectionGitHash = (TObjString*) gDirectory->Get("ProjectionGitHash");
  fProcessGitHash = (TObjString*) gDirectory->Get("ProcessGitHash");
  fUnfoldingGitHash = (TObjString*) gDirectory->Get("UnfoldingGitHash");

  // Read the TVectorT<float>:s
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    fCardEntries[iEntry] = (TVectorT<float>*) gDirectory->Get(fCardEntryNames[iEntry]);
  }

  // Naming for variable defining which background methods are used has changed over time
  // Define a legacy mode in order to be able to function with the old naming convention
  if((fCardEntries[kBackgroundMethods]) == NULL){
    fCardEntryNames[kBackgroundMethods] = "DoReflectedCone";
    fCardEntries[kBackgroundMethods] = (TVectorT<float>*) gDirectory->Get(fCardEntryNames[kBackgroundMethods]);
    fBackgroundLegacyMode = true;
  }
  
  // Read the file names
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    fFileNames[iFileName] = (TObjString*) gDirectory->Get(fFileNameSaveName[iFileName]);
  }

}

/*
 * Construct a data type string based on information on the card
 */
void EECCard::FindDataTypeString(){
  
  // Define the different data types corresponding to certain indices
  TString dataTypes[11] = {"pp","PbPb","pp MC","PbPb MC","pPb p #rightarrow -#eta","pPb p #rightarrow +#eta", "pPb p #rightarrow -#eta 5TeV", "pPb MC p #rightarrow -#eta", "pPb MC p #rightarrow +#eta", "OO","OO MC"}; 
  TString alternativeDataTypes[11] = {"pp","PbPb","Pythia8","Pythia+Hydjet","pPb p #rightarrow -#eta","pPb p #rightarrow +#eta", "pPb", "Pythia+EPOS p #rightarrow -#eta","Pythia+EPOS p #rightarrow +#eta", "OO","OO MC"};
  if(fDataType < 0 || fDataType > 10){
    fDataTypeString = "Unknown";
    fAlternativeDataTypeString = "Unknown";
    fDataTypeStringWithoutMCType = "Unknown";
    return;
  }
  
  // Remember the data type without the MC type
  fDataTypeStringWithoutMCType = alternativeDataTypes[fDataType];
  
  // Define Monte Carlo types and add them to MC productions, which are data types 2 and 3
  TString monteCarloString[4] = {" RecoReco"," RecoGen"," GenReco"," GenGen"};
  if(fMonteCarloType < 0 || fMonteCarloType > 3){
    fDataTypeString = "Unknown";
    fAlternativeDataTypeString = "Unknown";
    return;
  }
  
  for(int iDataType = 2; iDataType <=3; iDataType++){
    dataTypes[iDataType].Append(monteCarloString[fMonteCarloType]);
    alternativeDataTypes[iDataType].Append(monteCarloString[fMonteCarloType]);
  }
  
  // Remember the constructed data type string
  fDataTypeString = dataTypes[fDataType];
  fAlternativeDataTypeString = alternativeDataTypes[fDataType];
}

/*
 *  Getter for data type string
 */
TString EECCard::GetDataType() const{
  return fDataTypeString;
}

/*
 *  Getter for alternative data type string
 */
TString EECCard::GetAlternativeDataType(const bool includeMCtype) const{
  if(includeMCtype) return fAlternativeDataTypeString;
  return fDataTypeStringWithoutMCType;
}

// Getter for the first unfolded centrality bin index
int EECCard::GetFirstUnfoldedCentralityBin() const{
  return (*fCardEntries[kFirstUnfoldedCentralityBin])[1];
}

// Getter for the last unfolded centrality bin index
int EECCard::GetLastUnfoldedCentralityBin() const{
  return (*fCardEntries[kLastUnfoldedCentralityBin])[1];
}  

// Getter for the first unfolded track pT bin index
int EECCard::GetFirstUnfoldedTrackPtBin() const{
  return (*fCardEntries[kFirstUnfoldedTrackPtBin])[1];
}

// Getter for the last unfolded track pT bin index
int EECCard::GetLastUnfoldedTrackPtBin() const{
  return (*fCardEntries[kLastUnfoldedTrackPtBin])[1];
} 

// Getter for the first unfolded jet pT bin index
int EECCard::GetFirstUnfoldedJetPtBin() const{
  return (*fCardEntries[kFirstUnfoldedJetPtBin])[1];
} 

// Getter for the last unfolded jet pT bin index
int EECCard::GetLastUnfoldedJetPtBin() const{
  return (*fCardEntries[kLastUnfoldedJetPtBin])[1];
} 

/*
 *  Getter for subevent cut
 */
int EECCard::GetSubeventCut() const{
  return (*fCardEntries[kSubeventCut])[1];
}

/*
 *  Getter for jet type
 */
int EECCard::GetJetType() const{
  return (*fCardEntries[kJetType])[1];
}

// Getter for the minimum jet pT cut
double EECCard::GetJetPtCut() const{
  return (*fCardEntries[kMinPtCut])[1];
}

// Getter for the jet radios used in energy-energy correlator analysis
double EECCard::GetJetRadius() const{
  return (*fCardEntries[kJetRadius])[1];   // valor escalar guardado en la posición 0
}

/*
 *  Getter for information which background methods are included in the file
 */
int EECCard::GetBackgroundMethods() const{

  // In legacy mode, we need to decode the information the old DoReflectedCone variable included
  if(fBackgroundLegacyMode){
    int doReflectedCone = (*fCardEntries[kBackgroundMethods])[1];
    if(doReflectedCone <= 1) return doReflectedCone;
    if(doReflectedCone > 1) return 5;
  }

  // For newer files, read the actual value
  return (*fCardEntries[kBackgroundMethods])[1];
}

/*
 * Get the number of bins for internal index
 * If no vector is found in the index, return 0.
 */
int EECCard::GetNBins(const int index) const{
  if(fCardEntries[index]) return fCardEntries[index]->GetNoElements()-1;
  return 0;
}

// Get the number of centrality bins
int EECCard::GetNCentralityBins() const{
  return GetNBins(kCentralityBinEdges);
}

// Get the number of track pT bins
int EECCard::GetNTrackPtBins() const{
  return GetNBins(kTrackPtBinEdges);
}

// Get the number of jet pT bins in energy-energy correlator analysis
int EECCard::GetNJetPtBinsEEC() const{
  return GetNBins(kJetPtBinEdgesEEC);
}

// Get the number of track pT bins in energy-energy correlator analysis
int EECCard::GetNTrackPtBinsEEC() const{
  return GetNBins(kTrackPtBinEdgesEEC);
}

/*
 * Get a bin index based on a given value.
 * If value is out of bounds, return -1
 */
int EECCard::GetBinIndex(const int index, const double value) const{
  
  // Find the number of bins in the array
  int nBins = GetNBins(index);
  
  // If the number given is smaller than the first value, it is out of bounds. Return -1 in this case.
  if(value < (*fCardEntries[index])[1]) return -1;
  
  // Find the bin in which the given value is and return it
  for(int iBin = 2; iBin <= nBins+1; iBin++){
    if(value < (*fCardEntries[index])[iBin]) return iBin-2;
  }
  
  // If a value was not found, it must be larger than anything in the array. Return -1 in this case.
  return -1;
}

// Get the bin index for a given centrality value
int EECCard::GetBinIndexCentrality(const double value) const{
  return GetBinIndex(kCentralityBinEdges,value);
}

// Get the bin index for a given track pT value
int EECCard::GetBinIndexTrackPt(const double value) const{
  return GetBinIndex(kTrackPtBinEdges,value);
}

// Get the bin index for a given jet pT value in energy-energy correlator analysis
int EECCard::GetBinIndexJetPtEEC(const double value) const{
  return GetBinIndex(kJetPtBinEdgesEEC,value);
}

// Get the bin index for a given track pT value in energy-energy correlator analysis
int EECCard::GetBinIndexTrackPtEEC(const double value) const{
  return GetBinIndex(kTrackPtBinEdgesEEC,value);
}


/*
 * Find if a bin with given borders exist in an array defined by internal index
 *
 *  const int index = Internal index defining the array to search the bin from
 *  const double lowBorder = Low bin border of the searched bin
 *  const double highBorder = High bin border of the searched bin
 *
 *  return: Return the index of the corresponding bin. If the bin is not found, return -1
 */
int EECCard::FindBinIndex(const int index, const double lowBorder, const double highBorder) const{
  
  // Find the number of bins in the array
  int nBins = GetNBins(index);

  // Define a small number
  double epsilon = 0.001;
 
  // Go through all the bins and check if there exist a bin that matches with in given bins
  for(int iBin = 0; iBin < nBins; iBin++){
    if(TMath::Abs(GetLowBinBorder(index, iBin) - lowBorder) < epsilon){
      if(TMath::Abs(GetHighBinBorder(index, iBin) - highBorder) < epsilon){
        return iBin;
      }
    }
  }
  
  // If a match for bin borders was not found, the given bin is not in the array. Return -1 to show that.
  return -1;
}

// Find if a centrality bin with given borders exists and return its index
int EECCard::FindBinIndexCentrality(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kCentralityBinEdges,lowBorder,highBorder);
}

// Find if a centrality bin with given borders exists and return its index
int EECCard::FindBinIndexCentrality(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kCentralityBinEdges,binBorders.first,binBorders.second);
}

// Find if a track pT bin with given borders exists and return its index
int EECCard::FindBinIndexTrackPt(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kTrackPtBinEdges,lowBorder,highBorder);
}

// Find if a track pT bin with given borders exists and return its index
int EECCard::FindBinIndexTrackPt(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kTrackPtBinEdges,binBorders.first,binBorders.second);
}

// Find if a jet pT bin in energy-energy correlator analysis with given borders exists and return its index
int EECCard::FindBinIndexJetPtEEC(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kJetPtBinEdgesEEC,lowBorder,highBorder);
}

// Find if a jet pT bin in energy-energy correlator analysis with given borders exists and return its index
int EECCard::FindBinIndexJetPtEEC(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kJetPtBinEdgesEEC,binBorders.first,binBorders.second);
}

// Find if a track pT bin in energy-energy correlator analysis with given borders exists and return its index
int EECCard::FindBinIndexTrackPtEEC(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kTrackPtBinEdgesEEC,lowBorder,highBorder);
}

// Find if a track pT bin in energy-energy correlator analysis with given borders exists and return its index
int EECCard::FindBinIndexTrackPtEEC(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kTrackPtBinEdgesEEC,binBorders.first,binBorders.second);
}

// Get the low border of i:th bin from internal index
double EECCard::GetLowBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin > GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+1];
  return -1;
}

// Get the low border of i:th centrality bin
double EECCard::GetLowBinBorderCentrality(const int iBin) const{
  return GetLowBinBorder(kCentralityBinEdges,iBin);
}

// Get the low border of i:th track pT bin
double EECCard::GetLowBinBorderTrackPt(const int iBin) const{
  return GetLowBinBorder(kTrackPtBinEdges,iBin);
}

// Get the low border of i:th jet pT bin in energy-energy correlator analysis
double EECCard::GetLowBinBorderJetPtEEC(const int iBin) const{
  return GetLowBinBorder(kJetPtBinEdgesEEC,iBin);
}

// Get the low border of i:th track pT bin in energy-energy correlator analysis
double EECCard::GetLowBinBorderTrackPtEEC(const int iBin) const{
  return GetLowBinBorder(kTrackPtBinEdgesEEC,iBin);
}

// Get the high border of i:th bin from internal index
double EECCard::GetHighBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin >= GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+2];
  return -1;
}

// Get the high border of i:th centrality bin
double EECCard::GetHighBinBorderCentrality(const int iBin) const{
  return GetHighBinBorder(kCentralityBinEdges,iBin);
}

// Get the high border of i:th track pT bin
double EECCard::GetHighBinBorderTrackPt(const int iBin) const{
  return GetHighBinBorder(kTrackPtBinEdges,iBin);
}

// Get the high border of i:th jet pT bin in energy-energy correlator analysis
double EECCard::GetHighBinBorderJetPtEEC(const int iBin) const{
  return GetHighBinBorder(kJetPtBinEdgesEEC,iBin);
}

// Get the high border of i:th track pT bin in energy-energy correlator analysis
double EECCard::GetHighBinBorderTrackPtEEC(const int iBin) const{
  return GetHighBinBorder(kTrackPtBinEdgesEEC,iBin);
}

// Get the bin borders of the i:th centrality bin
std::pair<double,double> EECCard::GetBinBordersCentrality(const int iBin) const{
  return std::make_pair(GetLowBinBorderCentrality(iBin), GetHighBinBorderCentrality(iBin)); 
}

// Get the bin borders of the i:th track pT bin
std::pair<double,double> EECCard::GetBinBordersTrackPt(const int iBin) const{
  return std::make_pair(GetLowBinBorderTrackPt(iBin), GetHighBinBorderTrackPt(iBin)); 
}

// Get the bin borders of the i:th jet pT bin in energy-energy correlator analysis
std::pair<double,double> EECCard::GetBinBordersJetPtEEC(const int iBin) const{
  return std::make_pair(GetLowBinBorderJetPtEEC(iBin), GetHighBinBorderJetPtEEC(iBin)); 
}

// Get the bin borders of the i:th track pT bin in energy-energy correlator analysis
std::pair<double,double> EECCard::GetBinBordersTrackPtEEC(const int iBin) const{
  return std::make_pair(GetLowBinBorderTrackPtEEC(iBin), GetHighBinBorderTrackPtEEC(iBin)); 
}

/*
 *  Getter for the weight exponent used in energy-energy correlators
 */
int EECCard::GetWeightExponent(int index) const{

  // If the entry does not exist, the exponent must be 1 since this entry was added to the card to make a study for higher exponents
  if(fCardEntries[kWeightExponent]){

    // Do a sanity check for the input index. Return -1 for nonsensical indices
    if(index < 1) return -1;
    if(index > fCardEntries[kWeightExponent]->GetNoElements()) return -1;
    
    // If the index is sane, return the weight exponent from the said index
    return (*fCardEntries[kWeightExponent])[index];
  } else {
    return 1;
  }
}

/*
 * Find the index in card for the input weight exponent
 */
int EECCard::FindWeightExponentIndex(double weightExponent) const{

  // In very old files the weight exponent does not exist. Return 0 in this case.
  if(!fCardEntries[kWeightExponent]) return 0;

  // If an entry exists, loop over all the entries in the card and search for the input value
  double epsilon = 0.001;
  for(int iElement = 1; iElement <= fCardEntries[kWeightExponent]->GetNoElements(); iElement++){
    if(TMath::Abs(weightExponent - (*fCardEntries[kWeightExponent])[iElement]) < epsilon) return iElement;
  }

  // If we did not find an entry from the vector, return -1 to tell the entry does not exist
  return -1;

} 

/*
 * Get the number of weight exponents that are defined in the file
 */
int EECCard::GetNWeightExponents() const{

  // In very old files the weight exponent does not exist. Return 1 in this case.
  if(!fCardEntries[kWeightExponent]) return 1;

  // The number of the defined exponents is in the number of elements of the vector
  return fCardEntries[kWeightExponent]->GetNoElements();

}

/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  float entryContent = Content to be put into the vector
 */
void EECCard::AddOneDimensionalVector(int entryIndex, float entryContent){
  
  // Only allow addition to postprocessing vectors
  if(entryIndex <= kPtHatBinEdges && entryIndex != kWeightExponent) return;
  
  // Make a new one dimensional vector to the desired index with given content
  float contents[1] = {entryContent};
  fCardEntries[entryIndex] = new TVectorT<float>(1,1,contents);

}

/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  int dimension = Number of entries in the given array
 *  float *contents = Content to be put into the vector
 */
void EECCard::AddVector(int entryIndex, int dimension, double *contents){
  
  // Convert double pointer to float pointer
  float* convertedContents = new float[dimension];
  for(int i = 0; i < dimension; i++){
    convertedContents[i] = contents[i];
  }
  
  // Make a new vector to the desired index with given content
  fCardEntries[entryIndex] = new TVectorT<float>(1,dimension,convertedContents);
  
  // Delete the converted contents array
  delete[] convertedContents;
  
}

/*
 * Add file name to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added file name in file name array
 *  TString fileName = Added file name
 */
void EECCard::AddFileName(int entryIndex, TString fileName){
  
  // Make convert the string to TObjString and add it to the file name array
  fFileNames[entryIndex] = new TObjString(fileName.Data());
  
}

/*
 * Add git hash for the projecting code used to project the histograms to the card
 *
 * Arguments:
 *  const char* gitHash = Git hash to be added for projection
 */
void EECCard::AddProjectionGitHash(const char* gitHash){
  
  // Convert the const char* to TObjString and assign it to projection git hash
  fProjectionGitHash = new TObjString(gitHash);
  
}

/*
 * Add a git hash used to process the histograms in the file to the card
 *
 * Arguments:
 *  const char* gitHash = Git hash to be added for processing
 */
void EECCard::AddProcessGitHash(const char* gitHash){
  
  // Convert the const char* to TObjString and assign it to processing git hash
  fProcessGitHash = new TObjString(gitHash);
  
}

/*
 * Add a git hash used to unfold the energy-energy correlators to the card
 *
 * Arguments:
 *  const char* gitHash = Git hash to be added for unfolding
 */
void EECCard::AddUnfoldingGitHash(const char* gitHash){
  
  // Convert the const char* to TObjString and assign it to unfolding git hash
  fUnfoldingGitHash = new TObjString(gitHash);
  
}

/*
 * Print the contents of the card to the console
 */
void EECCard::Print() const{

  const char* gitHash = "Unknown";
  if(fGitHash != NULL) gitHash = fGitHash->String().Data();
  
  std::cout<<std::endl<<"========================= EECCard =========================="<<std::endl;
  std::cout << Form("%25s","GitHash"); //print keyword
  std::cout << ": " << gitHash << std::endl;
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry]){
      std::cout << Form("%25s",fCardEntryNames[iEntry]); //print keyword
      std::cout << " (dim = "<<fCardEntries[iEntry]->GetNoElements() << ") ";//print size of TVector
      for(int iElement = 1; iElement <= fCardEntries[iEntry]->GetNoElements(); iElement++){
        std::cout << (*fCardEntries[iEntry])[iElement] << " ";//TVector components
      }
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;
  
  if(fProjectionGitHash != NULL) std::cout << "Git hash for projections: " << fProjectionGitHash->String().Data() << std::endl;
  if(fUnfoldingGitHash != NULL) std::cout << "Git hash for unfolding: " << fUnfoldingGitHash->String().Data() << std::endl;
  
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    if(fFileNames[iFileName]){
      std::cout << "Used " << fFileNameType[iFileName] << " file: " << fFileNames[iFileName]->String().Data() << std::endl;
    }
  }
  
  if(fProcessGitHash != NULL) std::cout << "Git hash for histogram processing: " << fProcessGitHash->String().Data() << std::endl;
  
}

/*
 * Write the contents of the EECCard to a file
 */
void EECCard::Write(TDirectory* file){
  
  // Create a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write the git hash to the file
  if(fGitHash) fGitHash->Write("GitHash");
  
  // Write all the vectors to the file. Not all of these exist in older versions of cards, thus check if exists before writing.
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry])  fCardEntries[iEntry]->Write(fCardEntryNames[iEntry]);
  }
   
  // Write the git hash used for projections to the file
  if(fProjectionGitHash) fProjectionGitHash->Write("ProjectionGitHash");
  
  // Write all the data names to the file.
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    if(fFileNames[iFileName]) fFileNames[iFileName]->Write(fFileNameSaveName[iFileName]);
  }
  
  // Write the git hash used to unfold the energy-energy correlators
  if(fUnfoldingGitHash) fUnfoldingGitHash->Write("UnfoldingGitHash");

  // Write the git hash used to process the histograms
  if(fProcessGitHash) fProcessGitHash->Write("ProcessGitHash");
  
  // Return back to the main directory
  file->cd("../");
}

/*
 * Write the git hash used for unfolding energy-energy correlators to the file
 */
void EECCard::WriteUnfoldInfo(TDirectory* file){
  
  // Go to a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write the git hash used to unfold the histograms
  if(fUnfoldingGitHash) fUnfoldingGitHash->Write("UnfoldingGitHash");

  // Write the file name for the response matrix
  if(fFileNames[kResponseMatrixFile]) fFileNames[kResponseMatrixFile]->Write(fFileNameSaveName[kResponseMatrixFile]);

  // Write the unfold entries
  for(int iEntry = kFirstUnfoldedCentralityBin; iEntry <= kLastUnfoldedJetPtBin; iEntry++){
    if(fCardEntries[iEntry])  fCardEntries[iEntry]->Write(fCardEntryNames[iEntry]);
  }

  // Return back to the main directory
  file->cd("../");
  
}

/*
 * Write the git hash used for processing histograms to the file
 */
void EECCard::WriteProcessHash(TDirectory* file){
  
  // Go to a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write the git hash used to process the histograms
  if(fProcessGitHash) fProcessGitHash->Write("ProcessGitHash");
  
  // Return back to the main directory
  file->cd("../");
  
}

/*
 * Write the information related to background subtraction
 */
void EECCard::WriteBackground(TDirectory* file){
  
  // Go to a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");

  // Write the entries related to background subtraction
  for(int iEntry = kBackgroundMethod; iEntry <= kBackgroundSystematic; iEntry++){
    if(fCardEntries[iEntry])  fCardEntries[iEntry]->Write(fCardEntryNames[iEntry]);
  }

  // Return back to the main directory
  file->cd("../");
  
}

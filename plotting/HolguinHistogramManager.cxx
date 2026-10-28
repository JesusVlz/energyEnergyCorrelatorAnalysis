/*
 * Implementation of HolguinHistogramManager
 */

// Own includes
#include "HolguinHistogramManager.h"

/*
 * Default constructor
 */
HolguinHistogramManager::HolguinHistogramManager() :
  fKValueVector(0),
  fnKValues(0),
  fKValuesFound(false)
{  
  // Energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
        for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){ 
          for(int iKValue = 0; iKValue < kMaxKValues; iKValue++){   
            fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue] = NULL;
            fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue] = NULL;
          } // k-value loop
        } // Centrality loop
      } // Energy weight loop
    } // Track pT loop
  } // Jet pT loop
}

/*
 * Constructor with input file
 */
HolguinHistogramManager::HolguinHistogramManager(TString inputDirectory) :
  HolguinHistogramManager()
{  
  // Load all the graphs from the given directory
  LoadGraphs(inputDirectory); 
}

/*
 * Copy constructor
 */
HolguinHistogramManager::HolguinHistogramManager(const HolguinHistogramManager& in) :
  fKValueVector(in.fKValueVector),
  fnKValues(in.fnKValues),
  fKValuesFound(in.fKValuesFound)
{
  // Copy constructor
    
  // Energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
        for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){  
          for(int iKValue = 0; iKValue < kMaxKValues; iKValue++){         
            fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue] = in.fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue];
            fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue] = in.fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue];
          } // k-value loop
        } // Centrality loop
      } // Energy weight loop
    } // Track pT loop
  } // Jet pT loop
}

/*
 *  Load all the graphs from input directory. We assume that files for each configuration are present in the given directory.
 *
 *   Arguments: TString inputDirectory = Path to directory where .dat files are located
 */
void HolguinHistogramManager::LoadGraphs(TString inputDirectory){

  // If there is a "/" sign at the end of the inputDirectory string, remove it
  if(inputDirectory.EndsWith("/")){
    inputDirectory.Remove(inputDirectory.Capacity()-2);
  }
   
  // First, read all the calculations for the PbPb distirbution itself
  TString currentFile;
  int iKValue;
  std::vector<std::pair<double, TGraph*>> encodedGraphs;
  for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){
     for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
        for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){

          // Read the distribution file from the input directory
          currentFile = Form("%s/ForFitting_%s_b3p2GLV.dat", inputDirectory.Data(), fJetPtName[iJetPt]);

          // Read the graphs from this file
          encodedGraphs = GetGraphsFromDatFile(currentFile);

          // Find the grahps from the encoded data
          iKValue = 0;

          for(auto thisGraph : encodedGraphs){
            
            // The graph is provided as the second component of the encoded message
            fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue] = thisGraph.second;

            // If we have not yet determined k-values for this directory, read the k-values from file to vector
            if(!fKValuesFound){
              fKValueVector.push_back(thisGraph.first);
            }

            // Increment the k-value counter
            iKValue++;
          }

          // We have now successfuly determined the k-values from directory, raise a flag to celebrate!
          if(!fKValuesFound){
            fKValuesFound = true;
            fnKValues = iKValue;
          }
          
          // Next, read the ratio file from the input directory
          currentFile = Form("%s/ForFitting_%s_b3p2GLVratioAA_pp.dat", inputDirectory.Data(), fJetPtName[iJetPt]);

          // Read the graphs from this file
          encodedGraphs = GetGraphsFromDatFile(currentFile);

          // Find the grahps from the encoded data
          iKValue = 0;
          for(auto thisGraph : encodedGraphs){
            
            // The graph is provided as the second component of the encoded message
            fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue] = thisGraph.second;

            // Increment the k-value counter
            iKValue++;
          }

        } // Energy weight loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
}

/*
 * Get a vector of graphs that extract the data points from a .dat file
 * The file is expected to have the following format:
 * First line of the file defined the meaning of each columns
 * The first column has the x-axis points
 * The other columns are the y-axis points for different curves
 * The other rows give different points for the curves
 *
 * TString fileName = Name of the .dat file from which the data is read
 *
 * return: Structure that gives a vector of graphs, and a k-value associated to each graph
 */
std::vector<std::pair<double, TGraph*>> HolguinHistogramManager::GetGraphsFromDatFile(TString fileName){
  
  // Set up the file names file for reading
  std::ifstream file_stream(fileName);
  std::string line;
  std::vector<TString> lineVector;
  lineVector.clear();

  // Open the file names file for reading
  if( file_stream.is_open() ) {
    
    // Loop over the lines in the file
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      TString lineString(line);
      
      // Put all non-empty lines to file names vector
      if( lineString.CompareTo("", TString::kExact) != 0 ) {
        lineVector.push_back(lineString);
      } // Empty line if
      
    } // Loop over lines in the file
    
  // If cannot read the file, give error and end program
  } else {
    std::cout << "Error, could not open " << fileName.Data() << " for reading" << std::endl;
    return std::vector<std::pair<double, TGraph*>>();
  }

  // Find the number of predictions and points within prediction from the input
  TObjArray* parameterArray = lineVector.at(0).Tokenize(" ");  // Tokenize the string from every ' ' character
  const int numberOfPredictions = parameterArray->GetEntries()-1;
  const int numberOfPoints = lineVector.size()-1;

  // Create vectors to hold the information about the predictions
  std::vector<double> deltaRPoints;
  std::vector<std::vector<double>> eecPredictions;
  std::vector<double> valueOfK;

  // Resize the vectors to have the desired dimensions
  deltaRPoints.resize(numberOfPoints);
  eecPredictions.resize(numberOfPredictions);
  valueOfK.resize(numberOfPredictions);
  for(auto& prediction : eecPredictions){
    prediction.resize(numberOfPoints);
  }

  // Loop over the input and collect the information into arrays
  int numberOfParameters;
  TObjString* currentParameterObject;
  TString currentParameter;
  for(int iLine = 0; iLine < lineVector.size(); iLine++){
    parameterArray = lineVector.at(iLine).Tokenize(" ");
    numberOfParameters = parameterArray->GetEntries();
    for(int iParameter = 0; iParameter < numberOfParameters; iParameter++){   // Loop over all the files in the array
      currentParameterObject = (TObjString*)parameterArray->At(iParameter);
      currentParameter = currentParameterObject->String();

      // The first line is special, as it gives k-values for each following line
      if(iLine == 0){

        // Skip the first parameter in the first line. This is just DeltaR
        if(iParameter > 0){
          valueOfK[iParameter-1] = atof(currentParameter);
        }

      // The following lines give the data that we want to gather to the graphs
      } else {
       
        // The first column is the deltaR value, other columns show points for EECs
        if(iParameter == 0){
          deltaRPoints[iLine-1] = atof(currentParameter);
        } else {
          eecPredictions[iParameter-1][iLine-1] = atof(currentParameter);
        }
      }
    } // Parameter loop loop
  } // Line loop

  // Make graphs out of the information arrays and put them into a vector
  std::vector<std::pair<double, TGraph*>> graphedPrediction;
  for(int iPrediction = 0; iPrediction < numberOfPredictions; iPrediction++){
    graphedPrediction.push_back(std::make_pair(valueOfK[iPrediction], new TGraph(numberOfPoints, deltaRPoints.data(), eecPredictions[iPrediction].data())));
  }

  // Return the vector of predictions
  return graphedPrediction;

}

// ******************************************* //
// Getters for energy-energy correlator graphs //
// ******************************************* //

// Getter for energy-energy correlators from PbPb collisions using bin indices
TGraph* HolguinHistogramManager::GetEnergyEnergyCorrelatorPbPb(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iKValue) const{

  // If any of the bins is negative, indicating the histogram does not exist, return null
  if(iCentrality < 0 || iJetPt < 0 || iTrackPt < 0 || iEnergyWeight < 0 || iKValue < 0) return NULL;

  return fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue];
}

// Getter for energy-energy correlators from PbPb collisions using bin borders
TGraph* HolguinHistogramManager::GetEnergyEnergyCorrelatorPbPb(std::pair<int,int> centralityBin, std::pair<double,double> jetPtBin, double trackPtCut, double energyWeight, double valueOfK) const{
  int iCentrality = FindCentralityBinIndex(centralityBin);
  int iJetPt = FindJetPtBinIndex(jetPtBin);
  int iTrackPt = FindTrackPtBinIndex(trackPtCut);
  int iEnergyWeight = FindEnergyWeightIndex(energyWeight);
  int iKValue = FindKValueIndex(valueOfK);
  return GetEnergyEnergyCorrelatorPbPb(iCentrality, iJetPt, iTrackPt, iEnergyWeight, iKValue);
}

// Getter for PbPb to pp energy-energy correlator ratios using bin indices
TGraph* HolguinHistogramManager::GetEnergyEnergyCorrelatorPbPbToPpRatio(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iKValue) const{

  // If any of the bins is negative, indicating the histogram does not exist, return null
  if(iCentrality < 0 || iJetPt < 0 || iTrackPt < 0 || iEnergyWeight < 0 || iKValue < 0) return NULL;

  return fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iKValue];

}

// Getter for PbPb to pp energy-energy correlator ratios using bin borders
TGraph* HolguinHistogramManager::GetEnergyEnergyCorrelatorPbPbToPpRatio(std::pair<int,int> centralityBin, std::pair<double,double> jetPtBin, double trackPtCut, double energyWeight, double valueOfK) const{
  int iCentrality = FindCentralityBinIndex(centralityBin);
  int iJetPt = FindJetPtBinIndex(jetPtBin);
  int iTrackPt = FindTrackPtBinIndex(trackPtCut);
  int iEnergyWeight = FindEnergyWeightIndex(energyWeight);
  int iKValue = FindKValueIndex(valueOfK);
  return GetEnergyEnergyCorrelatorPbPbToPpRatio(iCentrality, iJetPt, iTrackPt, iEnergyWeight, iKValue);
}

// Find bin indices from bin borders

// Get a centrality bin index from a given centrality bin borders
int HolguinHistogramManager::FindCentralityBinIndex(std::pair<int,int> centralityBin) const{
  for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){
    if(centralityBin.first == fCentralityBinBorders[iCentrality].first){
      if(centralityBin.second == fCentralityBinBorders[iCentrality].second){
        return iCentrality;
      }
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get a jet pT bin index from a given jet pT bin borders
int HolguinHistogramManager::FindJetPtBinIndex(std::pair<double,double> jetPtBin) const{
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    if(TMath::Abs(jetPtBin.first - fJetPtBinBorders[iJetPt].first) < 0.0001){
      if(TMath::Abs(jetPtBin.second - fJetPtBinBorders[iJetPt].second) < 0.0001){
        return iJetPt;
      }
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get a track pT bin index from a given track pT cut
int HolguinHistogramManager::FindTrackPtBinIndex(double trackPtBin) const{
  for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
    if(TMath::Abs(trackPtBin-fTrackPtCuts[iTrackPt]) < 0.0001){
      return iTrackPt;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get an energy weight bin index from a given energy weight value
int HolguinHistogramManager::FindEnergyWeightIndex(double energyWeight) const{
  for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
    if(TMath::Abs(energyWeight-fEnergyWeights[iEnergyWeight]) < 0.0001){
      return iEnergyWeight;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get an k-value index from a given k-value
int HolguinHistogramManager::FindKValueIndex(double valueOfK) const{
  for(int iKValue = 0; iKValue < fnKValues; iKValue++){
    if(TMath::Abs(valueOfK - fKValueVector.at(iKValue)) < 0.0001){
      return iKValue;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get the number of k-values included in the input folder
int HolguinHistogramManager::GetNKValues() const{
  return fnKValues;
}
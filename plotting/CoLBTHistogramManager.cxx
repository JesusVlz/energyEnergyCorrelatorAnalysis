/*
 * Implementation of CoLBTHistogramManager
 */

// Own includes
#include "CoLBTHistogramManager.h"

/*
 * Default constructor
 */
CoLBTHistogramManager::CoLBTHistogramManager()
{  
  // Energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
        for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){ 
          for(int iQValue = 0; iQValue < kQValues; iQValue++){   
            fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iQValue] = NULL;
          } // k-value loop
        } // Centrality loop
      } // Energy weight loop
    } // Track pT loop
  } // Jet pT loop
}

/*
 * Constructor with input file
 */
CoLBTHistogramManager::CoLBTHistogramManager(TString inputDirectory) :
  CoLBTHistogramManager()
{  
  // Load all the graphs from the given directory
  LoadGraphs(inputDirectory); 
}

/*
 * Copy constructor
 */
CoLBTHistogramManager::CoLBTHistogramManager(const CoLBTHistogramManager& in)
{
  // Copy constructor
    
  // Energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
        for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){  
          for(int iQValue = 0; iQValue < kQValues; iQValue++){         
            fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iQValue] = in.fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iQValue];
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
void CoLBTHistogramManager::LoadGraphs(TString inputDirectory){

  // If there is a "/" sign at the end of the inputDirectory string, remove it
  if(inputDirectory.EndsWith("/")){
    inputDirectory.Remove(inputDirectory.Capacity()-2);
  }

  // Helper variables for creating graphs
  std::vector<std::vector<double>> parameterValue;
  std::vector<double> xValue;
  std::vector<double> yValue;
  std::vector<double> xError;
  std::vector<double> yError;
   
  // First, read all the calculations for the PbPb distirbution itself
  TString currentFile;
  std::vector<std::pair<double, TGraph*>> encodedGraphs;
  for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){
     for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
        for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
          for(int iQValue = 0; iQValue < kQValues; iQValue++){

            // Read the desired file from the input directory
            currentFile = Form("%s/eec%s_%s_%s.dat", inputDirectory.Data(), fJetPtName[iJetPt], fTrackPtName[iTrackPt], fQValueName[iQValue]);

            // Read all the lines from the current input file
            std::ifstream file_stream(currentFile);
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
    
            // If cannot read the file, give error and warn about resulting NULL graph
            } else {
              std::cout << "Error! Could not open " << currentFile.Data() << " for reading" << std::endl;
              std::cout << "TGraph corresponding to this file will be NULL!" << std::endl;
              continue;
            }

            // After we have all the information in a vector, we can dissect it into arrays
            // The first line is a comment, so it needs to be subtracted from the number of points
            const int nPoints = lineVector.size()-1;

            // Initialize the vectors to have nPoint zeros
            xValue.clear();
            yValue.clear();
            xError.clear();
            yError.clear();
            for(int iPoint = 0; iPoint < nPoints; iPoint++){
              xValue.push_back(0);
              yValue.push_back(0);
              xError.push_back(0);
              yError.push_back(0);
            }

            // The parameters in the input file are in this order:
            // 0: DeltaR value for the points
            // 1: EEC value for points
            // 2: Error for DeltaR
            // 3: Error for EEC
            parameterValue.clear();
            for(int i = 0; i < 4; i++){
              parameterValue.push_back(xValue);
            }

            // Loop over the lines in the input files and collect the numbers to arrays
            int numberOfParameters;
            TObjArray* parameterArray;
            TObjString* currentParameterObject;
            TString currentParameter;

            // Skip the first line in the file. It only contains comments
            for(int iLine = 1; iLine < lineVector.size(); iLine++){

              // The different parameters are separated by spaces, so tokenize strings to get individual parameters
              parameterArray = lineVector.at(iLine).Tokenize(" ");
              numberOfParameters = parameterArray->GetEntries();
              if(numberOfParameters != 4){
                std::cout << "Error! File " << currentFile.Data() << " does not have the correct format!" << std::endl;
                std::cout << "Cannot read the parameters." << std::endl;
                break;
              }

              // Loop over the parameters and collect them to the parameter array
              for(int iParameter = 0; iParameter < numberOfParameters; iParameter++){   
                currentParameterObject = (TObjString*)parameterArray->At(iParameter);
                currentParameter = currentParameterObject->String();
                parameterValue[iParameter][iLine-1] = atof(currentParameter);
              } // Parameter loop
            } // Point loop

            // The deltaR points in the x-axis are common for all graphs
            for(int iPoint = 0; iPoint < nPoints; iPoint++){
              xValue[iPoint] = parameterValue[0][iPoint];
              yValue[iPoint] = parameterValue[1][iPoint];
              xError[iPoint] = parameterValue[2][iPoint];
              yError[iPoint] = parameterValue[3][iPoint];
            }

            // Create the energy-energy correlator graph for PbPb for this bin
            fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iQValue] = new TGraphErrors(nPoints, xValue.data(), yValue.data(), xError.data(), yError.data());

          } // q-value loop
        } // Energy weight loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop  
}


// ******************************************* //
// Getters for energy-energy correlator graphs //
// ******************************************* //

// Getter for PbPb to pp energy-energy correlator ratios using bin indices
TGraphErrors* CoLBTHistogramManager::GetEnergyEnergyCorrelatorPbPbToPpRatio(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iQValue) const{

  // If any of the bins is negative, indicating the histogram does not exist, return null
  if(iCentrality < 0 || iJetPt < 0 || iTrackPt < 0 || iEnergyWeight < 0 || iQValue < 0) return NULL;

  return fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iQValue];

}

// Getter for PbPb to pp energy-energy correlator ratios using bin borders
TGraphErrors* CoLBTHistogramManager::GetEnergyEnergyCorrelatorPbPbToPpRatio(std::pair<int,int> centralityBin, std::pair<double,double> jetPtBin, double trackPtCut, double energyWeight, double qValue) const{
  int iCentrality = FindCentralityBinIndex(centralityBin);
  int iJetPt = FindJetPtBinIndex(jetPtBin);
  int iTrackPt = FindTrackPtBinIndex(trackPtCut);
  int iEnergyWeight = FindEnergyWeightIndex(energyWeight);
  int iQValue = FindQValueIndex(qValue);
  return GetEnergyEnergyCorrelatorPbPbToPpRatio(iCentrality, iJetPt, iTrackPt, iEnergyWeight, iQValue);
}

// Find bin indices from bin borders

// Get a centrality bin index from a given centrality bin borders
int CoLBTHistogramManager::FindCentralityBinIndex(std::pair<int,int> centralityBin) const{
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
int CoLBTHistogramManager::FindJetPtBinIndex(std::pair<double,double> jetPtBin) const{
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
int CoLBTHistogramManager::FindTrackPtBinIndex(double trackPtBin) const{
  for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
    if(TMath::Abs(trackPtBin-fTrackPtCuts[iTrackPt]) < 0.0001){
      return iTrackPt;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get an energy weight bin index from a given energy weight value
int CoLBTHistogramManager::FindEnergyWeightIndex(double energyWeight) const{
  for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
    if(TMath::Abs(energyWeight-fEnergyWeights[iEnergyWeight]) < 0.0001){
      return iEnergyWeight;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get an q-value index from a given q-value
int CoLBTHistogramManager::FindQValueIndex(double qValue) const{
  for(int iQValue = 0; iQValue < kQValues; iQValue++){
    if(TMath::Abs(qValue - fQValues[iQValue]) < 0.0001){
      return iQValue;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

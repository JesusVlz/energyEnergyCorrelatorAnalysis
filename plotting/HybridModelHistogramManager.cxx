/*
 * Implementation of HybridModelHistogramManager
 */

// Own includes
#include "HybridModelHistogramManager.h"

/*
 * Default constructor
 */
HybridModelHistogramManager::HybridModelHistogramManager()
{  
  // Energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
        for(int iWake = 0; iWake < kWakeConfigurations; iWake++){
          fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight][iWake] = NULL;
          for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){ 
            fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake] = NULL;
            fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake] = NULL;
          } // Centrality loop
        } // Wake type loop
      } // Energy weight loop
    } // Track pT loop
  } // Jet pT loop
}

/*
 * Constructor with input file
 */
HybridModelHistogramManager::HybridModelHistogramManager(TString inputDirectory) :
  HybridModelHistogramManager()
{  
  // Load all the graphs from the given directory
  LoadGraphs(inputDirectory); 
}

/*
 * Copy constructor
 */
HybridModelHistogramManager::HybridModelHistogramManager(const HybridModelHistogramManager& in)
{
  // Copy constructor
    
  // Energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
        for(int iWake = 0; iWake < kWakeConfigurations; iWake++){
          fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight][iWake] = in.fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight][iWake];
          for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){           
            fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake] = in.fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake];
            fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake] = in.fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake];
          } // Centrality loop
        } // Wake type loop
      } // Energy weight loop
    } // Track pT loop
  } // Jet pT loop
}

/*
 *  Load all the graphs from input directory. We assume that files for each configuration are present in the given directory.
 *
 *   Arguments: TString inputDirectory = Path to directory where .dat files are located
 */
void HybridModelHistogramManager::LoadGraphs(TString inputDirectory){

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

  // Loop over all the files that are present in the input directory and create graphs from the information inside the files
  TString currentFile;
  for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){
     for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
        for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
          for(int iWake = 0; iWake < kWakeConfigurations; iWake++){

            // Read the desired file from the input directory
            currentFile = Form("%s/HYBRID_Hadrons_%s_5020_inc_pthat50_%s_%s_%s_%s_%s_EEC.dat", inputDirectory.Data(), fSubtractionName[iWake], fCentralityName[iCentrality], fWakeName[iWake], fJetPtName[iJetPt], fTrackPtName[iTrackPt], fEnergyWeightName[iEnergyWeight]);

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
            const int nPoints = lineVector.size();

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
            // 1: Central value of the pp graph
            // 2: Error for the pp graph
            // 3: Upper error band value for PbPb graph
            // 4: Lower error band value for PbPb graph
            // 5: Upper error band value for PbPb/pp ratio
            // 6: Lower error band value for PbPb/pp ratio
            parameterValue.clear();
            for(int i = 0; i < 7; i++){
              parameterValue.push_back(xValue);
            }

            // Loop over the lines in the input files and collect the numbers to arrays
            int numberOfParameters;
            TObjArray* parameterArray;
            TObjString* currentParameterObject;
            TString currentParameter;
            for(int iPoint = 0; iPoint < lineVector.size(); iPoint++){

              // The different parameters are separated by tab character, so tokenize strings to get individual parameters
              parameterArray = lineVector.at(iPoint).Tokenize("\t");
              numberOfParameters = parameterArray->GetEntries();
              if(numberOfParameters != 7){
                std::cout << "Error! File " << currentFile.Data() << " does not have the correct format!" << std::endl;
                std::cout << "Cannot read the parameters." << std::endl;
                break;
              }

              // Loop over the parameters and collect them to the parameter array
              for(int iParameter = 0; iParameter < numberOfParameters; iParameter++){   
                currentParameterObject = (TObjString*)parameterArray->At(iParameter);
                currentParameter = currentParameterObject->String();
                parameterValue[iParameter][iPoint] = atof(currentParameter);
              } // Parameter loop
            } // Point loop

            // The deltaR points in the x-axis are common for all graphs
            for(int iPoint = 0; iPoint < nPoints; iPoint++){
              xValue[iPoint] = parameterValue[0][iPoint];
              xError[iPoint] = 0;
            }

            // Only fill pp histograms once, the same information is contained in several data files
            if(iCentrality == 0){

              // Collect the y-values and errors for pp graph
              for(int iPoint = 0; iPoint < nPoints; iPoint++){
                yValue[iPoint] = parameterValue[1][iPoint];
                yError[iPoint] = parameterValue[2][iPoint];
              }

              // Create the energy-energy correlator graph for pp for this bin
              fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight][iWake] = new TGraphErrors(nPoints, xValue.data(), yValue.data(), xError.data(), yError.data());

            }

            // Collect the y-values and errors for PbPb graph
            for(int iPoint = 0; iPoint < nPoints; iPoint++){
              yValue[iPoint] = (parameterValue[3][iPoint] + parameterValue[4][iPoint]) / 2;
              yError[iPoint] = (parameterValue[3][iPoint] - parameterValue[4][iPoint]) / 2;
            }

            // Create the energy-energy correlator graph for PbPb for this bin
            fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake] = new TGraphErrors(nPoints, xValue.data(), yValue.data(), xError.data(), yError.data());

            // Collect the y-values and errors for PbPb to pp ratio graph
            for(int iPoint = 0; iPoint < nPoints; iPoint++){
              yValue[iPoint] = (parameterValue[5][iPoint] + parameterValue[6][iPoint]) / 2;
              yError[iPoint] = (parameterValue[5][iPoint] - parameterValue[6][iPoint]) / 2;
            }

            // Create the energy-energy correlator graph for PbPb for this bin
            fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake] = new TGraphErrors(nPoints, xValue.data(), yValue.data(), xError.data(), yError.data());


          } // Wake type loop
        } // Energy weight loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
}

// Getters for energy-energy correlator graphs

// Getter for energy-energy correlators from pp collisions using bin indices
TGraphErrors* HybridModelHistogramManager::GetEnergyEnergyCorrelatorPp(const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iWake) const{

  // If any of the bins is negative, indicating the histogram does not exist, return null
  if(iJetPt < 0 || iTrackPt < 0 || iEnergyWeight < 0 || iWake < 0) return NULL;

  return fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight][iWake];
}

// Getter for energy-energy correlators from pp collisions using bin borders
TGraphErrors* HybridModelHistogramManager::GetEnergyEnergyCorrelatorPp(std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight, int iWake) const{
  int iJetPt = FindJetPtBinIndex(jetPtBin);
  int iTrackPt = FindTrackPtBinIndex(trackPtCut);
  int iEnergyWeight = FindEnergyWeightIndex(energyWeight);
  return GetEnergyEnergyCorrelatorPp(iJetPt, iTrackPt, iEnergyWeight, iWake);
}

// Getter for energy-energy correlators from PbPb collisions using bin indices
TGraphErrors* HybridModelHistogramManager::GetEnergyEnergyCorrelatorPbPb(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iWake) const{

  // If any of the bins is negative, indicating the histogram does not exist, return null
  if(iCentrality < 0 || iJetPt < 0 || iTrackPt < 0 || iEnergyWeight < 0 || iWake < 0) return NULL;

  return fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake];
}

// Getter for energy-energy correlators from PbPb collisions using bin borders
TGraphErrors* HybridModelHistogramManager::GetEnergyEnergyCorrelatorPbPb(std::pair<int,int> centralityBin, std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight, int iWake) const{
  int iCentrality = FindCentralityBinIndex(centralityBin);
  int iJetPt = FindJetPtBinIndex(jetPtBin);
  int iTrackPt = FindTrackPtBinIndex(trackPtCut);
  int iEnergyWeight = FindEnergyWeightIndex(energyWeight);
  return GetEnergyEnergyCorrelatorPbPb(iCentrality, iJetPt, iTrackPt, iEnergyWeight, iWake);
}

// Getter for PbPb to pp energy-energy correlator ratios using bin indices
TGraphErrors* HybridModelHistogramManager::GetEnergyEnergyCorrelatorPbPbToPpRatio(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iWake) const{

  // If any of the bins is negative, indicating the histogram does not exist, return null
  if(iCentrality < 0 || iJetPt < 0 || iTrackPt < 0 || iEnergyWeight < 0 || iWake < 0) return NULL;

  return fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight][iWake];

}

// Getter for PbPb to pp energy-energy correlator ratios using bin borders
TGraphErrors* HybridModelHistogramManager::GetEnergyEnergyCorrelatorPbPbToPpRatio(std::pair<int,int> centralityBin, std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight, int iWake) const{
  int iCentrality = FindCentralityBinIndex(centralityBin);
  int iJetPt = FindJetPtBinIndex(jetPtBin);
  int iTrackPt = FindTrackPtBinIndex(trackPtCut);
  int iEnergyWeight = FindEnergyWeightIndex(energyWeight);
  return GetEnergyEnergyCorrelatorPbPbToPpRatio(iCentrality, iJetPt, iTrackPt, iEnergyWeight, iWake);
}

// Find bin indices from bin borders

// Get a centrality bin index from a given centrality bin borders
int HybridModelHistogramManager::FindCentralityBinIndex(std::pair<int,int> centralityBin) const{
  for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){
    if(centralityBin.first == fCentralityBinBorders[iCentrality].first){
      if(centralityBin.second == fCentralityBinBorders[iCentrality].second){
        return iCentrality;
      }

      // The last bin 50-70 differs from analysis. Here accept only matching lower bin border
      if(centralityBin.first == 50){
        return iCentrality;
      }
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get a jet pT bin index from a given jet pT bin borders
int HybridModelHistogramManager::FindJetPtBinIndex(std::pair<int,int> jetPtBin) const{
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    if(jetPtBin.first == fJetPtBinBorders[iJetPt].first){
      if(jetPtBin.second == fJetPtBinBorders[iJetPt].second){
        return iJetPt;
      }
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get a track pT bin index from a given track pT cut
int HybridModelHistogramManager::FindTrackPtBinIndex(double trackPtBin) const{
  for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
    if(TMath::Abs(trackPtBin-fTrackPtCuts[iTrackPt]) < 0.0001){
      return iTrackPt;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get an energy weight bin index from a given energy weight value
int HybridModelHistogramManager::FindEnergyWeightIndex(double energyWeight) const{
  for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
    if(TMath::Abs(energyWeight-fEnergyWeights[iEnergyWeight]) < 0.0001){
      return iEnergyWeight;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Getter for a nice legend name for the wake configuration
const char* HybridModelHistogramManager::GetWakeName(const int iWake) const{
  if(iWake < 0 || iWake >= kWakeConfigurations) return "NonSense";
  return fWakeLegendName[iWake];
}

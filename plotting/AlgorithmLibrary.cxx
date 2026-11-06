/*
 * Library of useful algorithms, easily accessible from a class
 */

// Own includes
#include "AlgorithmLibrary.h"

/*
 *  Constructor
 */
AlgorithmLibrary::AlgorithmLibrary()
{
  // Constructor
}

/*
 *  Copy constructor
 */
AlgorithmLibrary::AlgorithmLibrary(const AlgorithmLibrary& in)
{
  // Copy constructor
}

/*
 * Assingment operator
 */
AlgorithmLibrary& AlgorithmLibrary::operator=(const AlgorithmLibrary& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  return *this;
}

/*
 *  Destructor
 */
AlgorithmLibrary::~AlgorithmLibrary()
{
  // Destructor
}

/*
 * Find the minimum and maximum values from a histogram. They must be more extreme than the current values
 *
 *  Arguments: TH1D* histogram = Histogram from which the minimum and maximum values are searched
 *             std::pair<double,double> currentMinMax = The found values need to be more extreme than these to be accepted
 */
std::pair<double,double> AlgorithmLibrary::FindHistogramMinMax(TH1D* histogram, std::pair<double,double> currentMinMax){
  
  // Use the whole histogram as a search range
  std::pair<int,int> searchRange = std::make_pair(1, histogram->GetNbinsX());
  return FindHistogramMinMax(histogram, currentMinMax, searchRange);
  
}

/*
 * Find the minimum and maximum values from a histogram. They must be more extreme than the current values
 *
 *  Arguments: TH1D* histogram = Histogram from which the minimum and maximum values are searched
 *             std::pair<double,double> currentMinMax = The found values need to be more extreme than these to be accepted
 *             std::pair<double,double> searchRange = bin value range from which the minimum and maximum values are searched
 */
std::pair<double,double> AlgorithmLibrary::FindHistogramMinMax(TH1D* histogram, std::pair<double,double> currentMinMax, std::pair<double,double> searchRange){
  
  // Find the bin range from which the minimum and maximum values are searched
  double epsilon = 0.00001;
  int firstBin = histogram->GetXaxis()->FindBin(searchRange.first+epsilon);
  int lastBin = histogram->GetXaxis()->FindBin(searchRange.second-epsilon);

  std::pair<int,int> searchBins = std::make_pair(firstBin, lastBin);

  // Return the minimum and maximum values from the histogram
  return FindHistogramMinMax(histogram, currentMinMax, searchBins);
  
}

/*
 * Find the minimum and maximum values from a histogram. They must be more extreme than the current values
 *
 *  Arguments: TH1D* histogram = Histogram from which the minimum and maximum values are searched
 *             std::pair<double,double> currentMinMax = The found values need to be more extreme than these to be accepted
 *             std::pair<int,int> searchRange = Bin index range from which the minimum and maximum values are searched
 */
std::pair<double,double> AlgorithmLibrary::FindHistogramMinMax(TH1D* histogram, std::pair<double,double> currentMinMax, std::pair<int,int> searchRange){
  
  // As initial guess, take the given minimum and maximum values
  std::pair<double,double> newMinMax = std::make_pair(currentMinMax.first, currentMinMax.second);
  
  // Loop through all the bins in the histogram and update the minimum and maximum values
  double currentValue, currentError;
  for(int iBin = searchRange.first; iBin <= searchRange.second; iBin++){
    currentValue = histogram->GetBinContent(iBin);
    currentError = histogram->GetBinError(iBin);
    if((currentValue-currentError) < newMinMax.first) newMinMax.first = currentValue-currentError;
    if((currentValue+currentError) > newMinMax.second) newMinMax.second = currentValue+currentError;
  }
  
  // Return the minimum and maximum values from the histogram
  return newMinMax;
  
}

/*
 * Rebin one dimensional histogram with asymmetric bin edges
 *
 * Arguments:
 *  TH1D* histogramInNeedOfRebinning = Histogram to be rebinned
 *  const int nBins = Number of bins for the rebinned histogram
 *  const double* binEdges = Bin edges for the rebinned histogram
 */
TH1D* AlgorithmLibrary::RebinAsymmetric(TH1D* histogramInNeedOfRebinning, const int nBins, const double* binEdges){
  
  // First, check that the new bin boundaries are also bin boundaries in the original histogram
  bool binsGood = CheckBinBoundaries(nBins,binEdges,histogramInNeedOfRebinning->GetXaxis());
  if(!binsGood){
    std::cout << "Cannot rebin histogram " << histogramInNeedOfRebinning->GetName() << " because given bin borders do not match with the bin borders of the original histogram!" << std::endl;
    return histogramInNeedOfRebinning;
  }
  
  // Clone the original histogram
  TH1D *clonedHistogram = (TH1D*) histogramInNeedOfRebinning->Clone("_rebinned");
  
  // Set the new binning for histogram (destroys content in each bin)
  clonedHistogram->SetBins(nBins,binEdges);
  
  // Make sure that each bin is set to zero
  for(int iBin = 1; iBin <= clonedHistogram->GetNbinsX(); iBin++){
    clonedHistogram->SetBinContent(iBin,0);
    clonedHistogram->SetBinError(iBin,0);
  }
  
  // Add the contents back to the histogram that was rebinned
  double binContent, binError, binCenter, oldContent, oldError;
  int newBin;
  double binWidth;
  for(int iBin = 1; iBin <= histogramInNeedOfRebinning->GetNbinsX(); iBin++){
    
    // Read the contents from the non-rebinned histogram
    binWidth = histogramInNeedOfRebinning->GetBinWidth(iBin);
    binContent = histogramInNeedOfRebinning->GetBinContent(iBin)*binWidth;  // Remove previous bin width normalization
    binError = histogramInNeedOfRebinning->GetBinError(iBin)*binWidth;      // Remove previous bin width normalization
    binCenter = histogramInNeedOfRebinning->GetBinCenter(iBin);
    
    // Add the contents to the rebinned histgram
    newBin = clonedHistogram->FindBin(binCenter);
    oldContent = clonedHistogram->GetBinContent(newBin);
    oldError = clonedHistogram->GetBinError(newBin);
    clonedHistogram->SetBinContent(newBin,binContent+oldContent);
    clonedHistogram->SetBinError(newBin,TMath::Sqrt(binError*binError+oldError*oldError));
  }
  
  // Normalize the bin contents to bin width
  for(int iBin = 1; iBin <= clonedHistogram->GetNbinsX(); iBin++){
    binWidth = clonedHistogram->GetBinWidth(iBin);
    binContent = clonedHistogram->GetBinContent(iBin);
    binError = clonedHistogram->GetBinError(iBin);
    clonedHistogram->SetBinContent(iBin,binContent/binWidth);
    clonedHistogram->SetBinError(iBin,binError/binWidth);
  }
  
  // Return the rebinned histogram
  return clonedHistogram;
}

/*
 * Rebin a two-dimensional histogram
 *
 *  Arguments:
 *   TH2D* histogramInNeedOfRebinning = The two dimensional histogram that is going to be rebinned
 *   const int nBinsX = Number of bins after rebinning in the x-axis
 *   const double* binBordersX = Bin borders for rebinning for the x-axis
 *   const int nBinsY = Number of bins after rebinning in the y-axis
 *   const double* binBordersY = Bin borders for rebinning for the y-axis
 *   const bool undoBinArea = Undo bin area normalization from the input histogram
 *   const bool normalizeBinArea = Normalize the output histogram contents with bin area
 *
 *   return: Rebinned histogram
 */
TH2D* AlgorithmLibrary::RebinHistogram(TH2D* histogramInNeedOfRebinning, const int nBinsX, const double* binBordersX, const int nBinsY, const double* binBordersY, const bool undoBinArea, const bool normalizeBinArea){

  // First, check that the new bin boundaries are also bin boundaries in the original histogram
  bool binsGood = CheckBinBoundaries(nBinsX, binBordersX, histogramInNeedOfRebinning->GetXaxis());
  if(!binsGood){
    std::cout << "Cannot rebin histogram " << histogramInNeedOfRebinning->GetName() << " because of a bin edge problem in x-axis!" << std::endl;
    return histogramInNeedOfRebinning;
  }
  
  binsGood = CheckBinBoundaries(nBinsY, binBordersY, histogramInNeedOfRebinning->GetYaxis());
  if(!binsGood){
    std::cout << "Cannot rebin histogram " << histogramInNeedOfRebinning->GetName() << " because of a bin edge problem in y-axis!" << std::endl;
    return histogramInNeedOfRebinning;
  }
  
  // Root does not offer a method to directly rebin a two-dimensional histogram, so I have implemented my own
  // Helper variables for rebinning
  double currentBinContent;
  double currentBinError;
  double nonRebinnedContent;
  double nonRebinnedError;
  int newHistogramIndex;
  double xValue;
  double yValue;
  double binWidthX;
  double binWidthY;
  double binArea = 1;
  
  // Create the histogram with new binning
  TString newName = Form("%sRebinned",histogramInNeedOfRebinning->GetName());
  TH2D *rebinnedHistogram = new TH2D(newName,newName,nBinsX,binBordersX,nBinsY,binBordersY);
  
  // Loop over all the bins in the old histogram and insert the content to the new histogram
  for(int iBinX = 1; iBinX <= histogramInNeedOfRebinning->GetNbinsX(); iBinX++){
    xValue = histogramInNeedOfRebinning->GetXaxis()->GetBinCenter(iBinX);
    binWidthX = histogramInNeedOfRebinning->GetXaxis()->GetBinWidth(iBinX);
    
    for(int iBinY = 1; iBinY <= histogramInNeedOfRebinning->GetNbinsY(); iBinY++){
      yValue = histogramInNeedOfRebinning->GetYaxis()->GetBinCenter(iBinY);
      binWidthY = histogramInNeedOfRebinning->GetYaxis()->GetBinWidth(iBinY);
      
      // If requested, undo bin area normalization from the input histogram
      if(undoBinArea) binArea = binWidthX*binWidthY;
      
      // Find the global bin index from the new histogram correcponding to the bin in the old histogram
      newHistogramIndex = rebinnedHistogram->FindBin(xValue,yValue);
      
      // Add the bin content from the old histogram to the new histogram, adding errors in quadrature
      nonRebinnedContent = histogramInNeedOfRebinning->GetBinContent(iBinX,iBinY) * binArea;
      nonRebinnedError = histogramInNeedOfRebinning->GetBinError(iBinX,iBinY) * binArea;
      currentBinContent = rebinnedHistogram->GetBinContent(newHistogramIndex);
      currentBinError = rebinnedHistogram->GetBinError(newHistogramIndex);
      rebinnedHistogram->SetBinContent(newHistogramIndex,currentBinContent+nonRebinnedContent);
      rebinnedHistogram->SetBinError(newHistogramIndex,TMath::Sqrt(TMath::Power(currentBinError,2)+TMath::Power(nonRebinnedError,2)));
    } // DeltaEta loop
  } // DeltaPhi loop
  
  // After rebinning, do bin area normalization for each bin if requested
  if(normalizeBinArea){
    for(int iBinX = 1; iBinX <= rebinnedHistogram->GetNbinsX(); iBinX++){
      binWidthX = rebinnedHistogram->GetXaxis()->GetBinWidth(iBinX);
      for(int iBinY = 1; iBinY <= rebinnedHistogram->GetNbinsY(); iBinY++){
        binWidthY = rebinnedHistogram->GetYaxis()->GetBinWidth(iBinY);
        currentBinContent = rebinnedHistogram->GetBinContent(iBinX,iBinY);
        currentBinError = rebinnedHistogram->GetBinError(iBinX,iBinY);
        rebinnedHistogram->SetBinContent(iBinX,iBinY,currentBinContent/(binWidthY*binWidthX));
        rebinnedHistogram->SetBinError(iBinX,iBinY,currentBinError/(binWidthY*binWidthX));
      } // y-axis loop
    } // x-axis loop
  } // Normalizing by bin area
  
  return rebinnedHistogram;
  
}

/*
 * Checker that new bin boundaries are compatible with the old ones
 *
 *  Arguments:
 *   const int nCheckedBins = Number of new bins to be checked
 *   const double* checkedBins = Bin boundaries to be checked
 *   const TAxis* originalAxis = Original axis against with the new bins are checked
 *
 *   return: True, if all the new bin boundaries can be found from the original axis. False if not.
 */
bool AlgorithmLibrary::CheckBinBoundaries(const int nCheckedBins, const double* checkedBins, const TAxis* originalAxis){
  
  // Flag, if the current bin is a bin boundary in the histogram to be rebinned
  bool binOK = false;
  
  // First, check that the bin boundaries for the rebinned histogram match with the old bin boundaries
  for(int iCheckedBin = 0; iCheckedBin < nCheckedBins + 1; iCheckedBin++){
    binOK = false;
    for(int iOldBin = 1; iOldBin <= originalAxis->GetNbins()+1; iOldBin++){
      
      // We the bin edge is close enough to one original bin, accept the bin
      if(TMath::Abs(originalAxis->GetBinLowEdge(iOldBin)-checkedBins[iCheckedBin]) < 1e-4){
        binOK = true;
        break;
      }
      
    } // Loop over bins in the original axis
    if(!binOK){ // If the bin is not in original histogram, print error message and return false
      std::cout << "The bin boundary " << checkedBins[iCheckedBin] << " is not a bin boundary in the original histogram!" << std::endl;
      return false;
    }
  } // Loop over bins to be checked
  
  // If all is good, return true
  return true;
}

// Normalize all the columns of a 2-D histogram to a given value
void AlgorithmLibrary::NormalizeMatrix(TH2D* histogramInNeedOfNormalization, const double value, const int direction){
  
  if(direction == 1){
    NormalizeColumns(histogramInNeedOfNormalization, value);
  } else {
    NormalizeRows(histogramInNeedOfNormalization, value);
  }
  
}

// Normalize all the columns of a 2-D histogram to a given value
void AlgorithmLibrary::NormalizeColumns(TH2D* histogramInNeedOfNormalization, const double value){
  
  // Helper variables
  double binContent;
  double binError;
  double binSum;
  double normalizationFactor;
  
  // Loop over all the y-bins for a given x-bin and normalize the content
  for(int iBinX = 1; iBinX <= histogramInNeedOfNormalization->GetNbinsX(); iBinX++){
    
    binSum = 0;
    
    // First calculate the sum of contents in the y-bins
    for(int iBinY = 1; iBinY <= histogramInNeedOfNormalization->GetNbinsY(); iBinY++){
      binSum += histogramInNeedOfNormalization->GetBinContent(iBinX, iBinY);
    }
    
    // If there is no content, no need to normalize
    if(binSum == 0) continue;
    
    // Get the normalization factor
    normalizationFactor = value / binSum;
    
    // Normalize the bin content in each bin
    for(int iBinY = 1; iBinY <= histogramInNeedOfNormalization->GetNbinsY(); iBinY++){
      binContent = histogramInNeedOfNormalization->GetBinContent(iBinX, iBinY);
      binError = histogramInNeedOfNormalization->GetBinError(iBinX, iBinY);
      binContent *= normalizationFactor;
      binError *= normalizationFactor;
      histogramInNeedOfNormalization->SetBinContent(iBinX, iBinY, binContent);
      histogramInNeedOfNormalization->SetBinError(iBinX, iBinY, binError);
    }
  }
  
}

// Normalize all the rows of a 2-D histogram to a given value
void AlgorithmLibrary::NormalizeRows(TH2D* histogramInNeedOfNormalization, const double value){
  
  // Helper variables
  double binContent;
  double binError;
  double binSum;
  double normalizationFactor;
  
  // Loop over all the x-bins for a given y-bin and normalize the content
  for(int iBinY = 1; iBinY <= histogramInNeedOfNormalization->GetNbinsY(); iBinY++){
    
    binSum = 0;
    
    // First calculate the sum of contents in the y-bins
    for(int iBinX = 1; iBinX <= histogramInNeedOfNormalization->GetNbinsX(); iBinX++){
      binSum += histogramInNeedOfNormalization->GetBinContent(iBinX, iBinY);
    }
    
    // If there is no content, no need to normalize
    if(binSum == 0) continue;
    
    // Get the normalization factor
    normalizationFactor = value / binSum;
    
    // Normalize the bin content in each bin
    for(int iBinX = 1; iBinX <= histogramInNeedOfNormalization->GetNbinsX(); iBinX++){
      binContent = histogramInNeedOfNormalization->GetBinContent(iBinX, iBinY);
      binError = histogramInNeedOfNormalization->GetBinError(iBinX, iBinY);
      binContent *= normalizationFactor;
      binError *= normalizationFactor;
      histogramInNeedOfNormalization->SetBinContent(iBinX, iBinY, binContent);
      histogramInNeedOfNormalization->SetBinError(iBinX, iBinY, binError);
    }
  }
  
}

/*
 * Method for rotating a two dimensional histogram. Assumes constant bin size.
 *
 *  Arguments:
 *   TH2D* originalHistogram = Histogram to be rotated
 *
 *   return: Histogram where x and y axes have been switched
 */
TH2D* AlgorithmLibrary::RotateHistogram(TH2D* originalHistogram){
  
  // Find the binning information
  double nBinsX = originalHistogram->GetNbinsX();
  double nBinsY = originalHistogram->GetNbinsY();
  double lowBinX = originalHistogram->GetXaxis()->GetBinLowEdge(1);
  double highBinX = originalHistogram->GetXaxis()->GetBinUpEdge(nBinsX);
  double lowBinY = originalHistogram->GetYaxis()->GetBinLowEdge(1);
  double highBinY = originalHistogram->GetYaxis()->GetBinUpEdge(nBinsY);
  
  // Create a new histogram with inverted axes
  TString newName = Form("%sRotated", originalHistogram->GetName());
  TString newTitle = Form("%sRotated", originalHistogram->GetTitle());
  TH2D *rotatedHistogram = new TH2D(newName, newTitle, nBinsY, lowBinY, highBinY, nBinsX, lowBinX, highBinX);
  
  // Fill the new histogram
  double binValue, binError;
  for(int iX = 1; iX <= nBinsX; iX++){
    for(int iY = 1; iY <= nBinsY; iY++){
      binValue = originalHistogram->GetBinContent(iX, iY);
      binError = originalHistogram->GetBinError(iX, iY);
      rotatedHistogram->SetBinContent(iY, iX, binValue);
      rotatedHistogram->SetBinError(iY, iX, binError);
    }
  }
  
  // Return the new histogram
  return rotatedHistogram;
  
}

/* Remove uncertainties from a histogram by making them negligibly small.
 *
 *  TH1* histogramInNeedForAccuracy = Histogram from which the uncertainties are removed
 *  const double epsilon = Small uncertainty with which the original uncertianty is replaced with
 */
void AlgorithmLibrary::RemoveUncertainties(TH1* histogramInNeedForAccuracy, const double epsilon){

  double negligibleUncertainty;
  double currentValue;
  for(int iBin = 1; iBin <= histogramInNeedForAccuracy->GetNbinsX(); iBin++){
    currentValue = histogramInNeedForAccuracy->GetBinContent(iBin);
    negligibleUncertainty = currentValue * epsilon;
    histogramInNeedForAccuracy->SetBinError(iBin, negligibleUncertainty);
  }

}

/*
 * Transform histogram to another one that shows the relative uncertainty of the original histogram in each bin
 *
 *  TH1D* transformedHistogram = Histogram in need of transformation
 *  const bool centerAtOne = True: Put the central value at one and relative systematic uncertainty around it.
 *                           False: Show relative systematic uncertainty as the histogram value. Disable errors.
 */
void AlgorithmLibrary::TransformToRelativeUncertainty(TH1D* transformedHistogram, const bool centerAtOne){
  
  double binContent;
  double binError;
  double relativeError;

  for(int iBin = 1; iBin <= transformedHistogram->GetNbinsX(); iBin++){
    binContent = transformedHistogram->GetBinContent(iBin);
    binError = transformedHistogram->GetBinError(iBin);
    relativeError = binError / binContent;

    if(centerAtOne){
      transformedHistogram->SetBinContent(iBin, 1);
      transformedHistogram->SetBinError(iBin, relativeError);
    } else {
      transformedHistogram->SetBinContent(iBin, relativeError);
      transformedHistogram->SetBinError(iBin, 0);
    }
  }

}

/*
 * Transform histogram describing relative uncertainties to one describing absolute uncertainties
 *
 *  TH1D* transformedHistogram = Histogram in need of transformation
 *  TH1D* absoluteScaleHistogram = Histogram giving the absolute scale for transformation
 *  const bool centerAtOne = True: The central value  in the relative uncertainty histogram is at one and relative systematic uncertainty around it.
 *                           False: The relative uncertainty is desribed in the histogram value in the relative uncertainty histogram
 */
void AlgorithmLibrary::TransformToAbsoluteUncertainty(TH1D* transformedHistogram, TH1D* absoluteScaleHistogram, const bool centerAtOne){

  double binContent;
  double binError;
  double relativeError;

  for(int iBin = 1; iBin <= transformedHistogram->GetNbinsX(); iBin++){

    if(centerAtOne){
      relativeError = transformedHistogram->GetBinError(iBin);
    } else {
      relativeError = transformedHistogram->GetBinContent(iBin);
    }

    binContent = absoluteScaleHistogram->GetBinContent(iBin);
    binError = binContent*relativeError;

    transformedHistogram->SetBinContent(iBin, binContent);
    transformedHistogram->SetBinError(iBin, binError);
  }

}

/*
 * Suppress single bin fluctuations in the fluctuating histogram
 *
 *  TH1D* fluctuatingHistogram = Histogram suffering from fluctuations
 *  const double lowRange = Lowest value in x-axis taken into account in fluctuation suppression
 *  const double highRange = Highest value on x-axis taken into account in fluctuation suppression
 *  const double threshold = Definition of a fluctuation that is too big for a single bin
 *  const double suppressionLevel = Definition how much fluctuating bins are suppressed
 */
void AlgorithmLibrary::SuppressSingleBinFluctuations(TH1D* fluctuatingHistogram, const double lowRange, const double highRange, const double threshold, const double suppressionLevel){

  double epsilon = 0.00001;
  int firstBin = fluctuatingHistogram->GetXaxis()->FindBin(lowRange+epsilon);
  int lastBin = fluctuatingHistogram->GetXaxis()->FindBin(highRange-epsilon);

  // Fluctuations are not suppressed in the first and last bin of the range, but they are in all the bins in between
  double previousValue, currentValue, nextValue, biggestNeighbor;
  for(int iBin = firstBin+1; iBin < lastBin; iBin++){
    previousValue = fluctuatingHistogram->GetBinContent(iBin-1);
    currentValue = fluctuatingHistogram->GetBinContent(iBin);
    nextValue = fluctuatingHistogram->GetBinContent(iBin+1);

    if(currentValue > threshold*previousValue && currentValue > threshold*nextValue){
      biggestNeighbor = previousValue;
      if(nextValue > previousValue) biggestNeighbor = nextValue;
      fluctuatingHistogram->SetBinContent(iBin, biggestNeighbor+biggestNeighbor*suppressionLevel);
    }
  }

}

/*
 * Square the contents of a histogram
 *
 *  TH1* transformedHistogram = Histogram that will get its contents exponentiated
 *  const double exponent = Amount of exponentiation applied to the histogram bins
 */
void AlgorithmLibrary::ExponentiateHistogram(TH1* transformedHistogram, const double exponent){

  // Helper variables for the transformation
  double binContent;
  double binError;

  // Loop over all the bins in the histogram
  for(int iBin = 1; iBin <= transformedHistogram->GetNbinsX(); iBin++){
    
    // Find current bin content and error in each bin
    binContent = transformedHistogram->GetBinContent(iBin);
    binError = transformedHistogram->GetBinError(iBin);

    // Square the bin content and properly scale the bin error
    transformedHistogram->SetBinContent(iBin, TMath::Power(binContent, exponent));
    transformedHistogram->SetBinError(iBin, exponent * TMath::Power(binContent, exponent-1) * binError);

  } // Bin loop

}

/*
 * Square the contents of a histogram
 *
 *  TH1* transformedHistogram = Histogram that will get its contents squared
 */
void AlgorithmLibrary::SquareHistogram(TH1* transformedHistogram){

  // Use the exponentiation function to do the dirty work
  ExponentiateHistogram(transformedHistogram, 2);

}

/*
 * Transform histogram to cumulant
 *
 *  TH1D* originalHistogram = Histogram from which cumulant distribution is calculated
 *  const int lowestBin = Bin index of the lowest bin included in the cumulant calculation
 *
 *  return: Cumulant distribution corresponding to the original histogram
 */
TH1D* AlgorithmLibrary::GetCumulant(TH1* originalHistogram, const int lowestBin){
  
  // First, clone the original histogram to get the correct binning
  TH1D* cumulantHistogram = (TH1D*) originalHistogram->Clone(Form("%sCumulant", originalHistogram->GetName()));

  // For the cumulant distribution, set the bins before the lowest bin border to zero
  for(int iBin = 1; iBin < lowestBin; iBin++){
    cumulantHistogram->SetBinContent(iBin, 0);
    cumulantHistogram->SetBinError(iBin, 0);
  }

  // Helper variables for integration
  double binIntegral;
  double binIntegralError;

  // For the other bins, calculate the cumulant by integrating from the lowest bin to the current bin
  for(int iBin = lowestBin; iBin <= originalHistogram->GetNbinsX(); iBin++){
    binIntegral = originalHistogram->IntegralAndError(lowestBin, iBin, binIntegralError, "width");
    cumulantHistogram->SetBinContent(iBin, binIntegral);
    cumulantHistogram->SetBinError(iBin, binIntegralError);
  }

  // Return the cumulant histogram
  return cumulantHistogram;

}

/*
 * Change the normalization of the histogram from bin width based to bin area based
 *
 *  TH1* transformedHistogram = Histogram that will get its contents squared
 */
void AlgorithmLibrary::ChangeBinWidthToBinAreaNormalization(TH1* normalizedHistogram){

  // Helper variables
  double binContent;
  double binError;
  double binWidth;
  double binArea;
  double lowBinBorder;
  double highBinBorder;

  // Loop over all the bins in the histogram
  for(int iBin = 1; iBin <= normalizedHistogram->GetNbinsX(); iBin++){
    
    // Find current bin content, error, and width information
    binContent = normalizedHistogram->GetBinContent(iBin);
    binError = normalizedHistogram->GetBinError(iBin);
    lowBinBorder = normalizedHistogram->GetXaxis()->GetBinLowEdge(iBin);
    highBinBorder = normalizedHistogram->GetXaxis()->GetBinUpEdge(iBin);
    binWidth = highBinBorder - lowBinBorder;
    binArea = TMath::Pi() * (highBinBorder*highBinBorder - lowBinBorder*lowBinBorder);

    // Square the bin content and properly scale the bin error
    normalizedHistogram->SetBinContent(iBin, binContent * binWidth / binArea);
    normalizedHistogram->SetBinError(iBin, binError * binWidth / binArea);

  } // Bin loop

}

/*
 * Get a TString with today's date
 */
TString AlgorithmLibrary::GetToday(){

  // Get the current time_point
  std::chrono::system_clock::time_point now = std::chrono::system_clock::now();

  // Convert the time_point to time_t
  std::time_t now_t = std::chrono::system_clock::to_time_t(now);

  // Convert the time_t to tm in local timezone
  std::tm* now_tm = std::localtime(&now_t);

  // Create an output string stream
  std::ostringstream oss;

  // Write the date in year-month-day format to the stream
  oss << std::put_time(now_tm, "%Y-%m-%d");

  // Get the string from the stream
  std::string date_str = oss.str();

  // Return the string as TString
  return date_str;

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
 * return: Structure that gives a vector of graphs, and a label associated to each graph
 */
std::vector<std::pair<TString, TGraph*>> AlgorithmLibrary::GetGraphsFromDatFile(TString fileName){
  
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
    return std::vector<std::pair<TString,TGraph*>>();
  }

  // Find the number of predictions and points within prediction from the input
  TObjArray* parameterArray = lineVector.at(0).Tokenize(" ");  // Tokenize the string from every ' ' character
  const int numberOfPredictions = parameterArray->GetEntries()-1;
  const int numberOfPoints = lineVector.size()-1;

  // Create vectors to hold the information about the predictions
  std::vector<double> deltaRPoints;
  std::vector<std::vector<double>> eecPredictions;
  std::vector<TString> predictionLabel;

  // Resize the vectors to have the desired dimensions
  deltaRPoints.resize(numberOfPoints);
  eecPredictions.resize(numberOfPredictions);
  predictionLabel.resize(numberOfPredictions);
  for(auto& prediction : eecPredictions){
    prediction.resize(numberOfPoints);
  }

  // Loop over the input and collect the information into arrays
  int numberOfParameters;
  TObjString* currentParameterObject;
  TString currentParameter;
  for(std::vector<TString>::size_type iLine = 0; iLine < lineVector.size(); iLine++){
    parameterArray = lineVector.at(iLine).Tokenize(" ");
    numberOfParameters = parameterArray->GetEntries();
    for(int iParameter = 0; iParameter < numberOfParameters; iParameter++){   // Loop over all the files in the array
      currentParameterObject = (TObjString*)parameterArray->At(iParameter);
      currentParameter = currentParameterObject->String();

      // The first line is special, as it gives labels for each following line
      if(iLine == 0){

        // Ship the first parameter in the first line. This would give the x-axis label
        if(iParameter > 0){
          predictionLabel[iParameter-1] = currentParameter;
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
  std::vector<std::pair<TString, TGraph*>> graphedPrediction;
  for(int iPrediction = 0; iPrediction < numberOfPredictions; iPrediction++){
    graphedPrediction.push_back(std::make_pair(predictionLabel[iPrediction], new TGraph(numberOfPoints, deltaRPoints.data(), eecPredictions[iPrediction].data())));
  }

  // Return the vector of predictions
  return graphedPrediction;

}

/*
 * Create a histogram contents of a graph. Linear interpolation between graph points are used
 */
TH1D* AlgorithmLibrary::Histogrammify(TGraph* gSource, TH1D* hBinning){

  // Read the binning for the histogram from the example histogram
  TH1D* histogrammifiedGraph = (TH1D*) hBinning->Clone(Form("histogrammified%s", gSource->GetName()));
  histogrammifiedGraph->Reset();

  // For each of the histogram bins, evaluate the graph content using linear interpolation
  double binCenter;
  for(int iBin = 1; iBin <= hBinning->GetNbinsX(); iBin++){
    binCenter = histogrammifiedGraph->GetBinCenter(iBin);
    histogrammifiedGraph->SetBinContent(iBin, gSource->Eval(binCenter));
  }

  // Return the histogrammified graph
  return histogrammifiedGraph;

}

/*
 * Create a histogram contents of a TGraphErrors. Linear interpolation between graph points are used.
 * Errors in the y-axis are taken into account in the created histogram. Errors in the x-axis are ignored.
 */
TH1D* AlgorithmLibrary::HistogrammifyWithErrors(TGraphErrors* gSource, TH1D* hBinning){

  // Create graph with upper and lower uncertainty bands such that they can be nicely interpolated
  TGraph* upperErrorGraph = (TGraph*) gSource->Clone("upperErrorGraph");
  TGraph* lowerErrorGraph = (TGraph*) gSource->Clone("lowerErrorGraph");

  // Shift the points in the upper and lower error graphs up and down by the error of the point
  double yValue, yError;
  for(int iPoint = 0; iPoint < gSource->GetN(); iPoint++){
    yValue = gSource->GetPointY(iPoint);
    yError = gSource->GetErrorY(iPoint);
    upperErrorGraph->SetPointY(iPoint, yValue+yError);
    lowerErrorGraph->SetPointY(iPoint, yValue-yError);
  }

  // Histogrammify all graphs
  TH1D* histogrammifiedGraph = Histogrammify(gSource, hBinning);
  TH1D* histogrammifiedUpperErrorGraph = Histogrammify(upperErrorGraph, hBinning);
  TH1D* histogrammifiedLowerErrorGraph = Histogrammify(lowerErrorGraph, hBinning);

  // For each of the histogram bins, set the bin error to be the largest between upper and lower errors
  // They should be the same by construction, but do this to be safe
  double binContent, upperbinContent, lowerBinContent;
  double binError;
  for(int iBin = 1; iBin <= hBinning->GetNbinsX(); iBin++){
    binContent = histogrammifiedGraph->GetBinContent(iBin);
    upperbinContent = histogrammifiedUpperErrorGraph->GetBinContent(iBin);
    lowerBinContent = histogrammifiedLowerErrorGraph->GetBinContent(iBin);
    binError = upperbinContent - binContent;
    if(binContent - lowerBinContent > binError) binError = binContent - lowerBinContent;
    histogrammifiedGraph->SetBinError(iBin, binError);
  }

  // Return the histogrammified graph
  return histogrammifiedGraph;

}

/*
 * Get a string out of a number
 *
 *  double number = Number from which we want to dig a string out
 *  int accuracy = Number of numbers after the decimal point to show
 *  char decimalPoint = Character used to describe the decimal point in the string
 */
TString AlgorithmLibrary::StringifyNumber(double number, int accuracy, char decimalPoint){

  // Make a string for the number
  TString numberString;
  double epsilon = 1e-10;

  // Check if the number is an integer
  if(TMath::Abs(number - static_cast<int>(number)) < epsilon){
    // If the number is integer already, just return integer string 
    numberString = Form("%.0f", number);
  } else {
    // If the number is not integer, return the value with the defined accuracy
    TString accuracyString = Form("%%.%df", accuracy);
    numberString = Form(accuracyString, number);
    if(decimalPoint != '.'){
      numberString.ReplaceAll('.',decimalPoint);
    }
  }

  // Return the stringified number
  return numberString;

}

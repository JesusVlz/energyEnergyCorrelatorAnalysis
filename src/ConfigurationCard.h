// Class for reading necessary configuration for an analysis from a card file
// Original code written by Jyväskylä ALICE group
// Modified for CMS by Jussi Viinikainen

#ifndef CONFIGURATIONCARD_H
#define CONFIGURATIONCARD_H

// C++ includes
#include <iostream>
#include <fstream>
#include <vector>

// Root includes
#include <TString.h>
#include <TVector.h>
#include <TObjString.h>
#include <TFile.h>
#include <TF1.h>
#include <THashList.h>
#include <TNamed.h>

using namespace std;
#ifndef MAXDIMBUFFER
#define MAXDIMBUFFER
const int kMaxDimBuffer = 300; // Maximum length of a line read to a buffer
#endif

class ConfigurationCard {

  //====   M e m b e r    F u n c t i o n s   ========

public:

  ConfigurationCard(); // constructor
  ConfigurationCard(const char *filename); // constructor
  ConfigurationCard& operator=(const ConfigurationCard& obj);

  virtual ~ConfigurationCard();

  void AddToKeyTable( TString key, int index );

  bool HasKey(TString keyword) const; // Check optional settings without requiring a value
  float  Get(TString keyword, int VectorComponent=0) const; //get TVector component
  TString  GetStr(TString keyword ) const; //get TVector component
  TVector* GetVector( TString keyword ) ;
  int GetN(TString keyword) const;       //get TVector dimension
  int GetBin(TString keyword, double value) const;  // Find the bin for value from keyword vector
  int GetNBin(TString keyword) const;   // Get number of bins related to keyword
  void PrintOut();
  void WriteCard(TDirectory *file) const;
  void ReadInputLine( const char* buffer );
  void SetGitHash(const char* hash);

protected:
  
  void    ReadInputCard();
  int     GetNwithIndex(int i) const{ return fValuesVector[i].GetNrows(); }
  unsigned int GetTVectorIndex(TString keyword, int tol=0) const;

  //====   D a t a    M e m b e r s  ========

  char fCardName[255];                       // File name for card
  int  fnEntry;                              // Number of lines in configuration file
  std::vector< TString > fKeyWordVector;     // Array of key words
  std::vector< TVector > fValuesVector;      // Array of float number config parameter vectors
  std::vector< TString > fValueString;       // Storage of raw input string for each item
  TObjString fGitHash;                       // String for git hash
  THashList fKeyTable;                       // key map with hash algorithm

};

#endif























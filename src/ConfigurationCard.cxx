// Implementation of the card to read the configuration from a file

// Header include
#include "ConfigurationCard.h"

ConfigurationCard::ConfigurationCard() :
fnEntry(0),
fKeyWordVector(0),
fValuesVector(0),
fValueString(0),
fGitHash("NotSet"),
fKeyTable(0)
{   
  //constructor
}

ConfigurationCard::ConfigurationCard(const char *filename) :
fnEntry(0),
fKeyWordVector(0),
fValuesVector(0),
fValueString(0),
fGitHash("NotSet"),
fKeyTable(0)
{  
  //constructor
  
  strcpy(fCardName, filename);  //needed in PrintOut()
  
  if( strlen( filename ) > 0 ){
    
    //----  r e a d   t h e   c a r d ----
    ReadInputCard(); //read config file fill TVectors
  }
}

/*
 * Equal sign operator
 */
ConfigurationCard& ConfigurationCard::operator=(const ConfigurationCard& obj){
  // equal sign operator
  
  // TODO: Implement this
  
  return *this;
}

ConfigurationCard::~ConfigurationCard(){
  // destructor
}


bool ConfigurationCard::HasKey(TString keyword) const{
  return fKeyTable.FindObject(keyword.Data()) != nullptr;
}

unsigned int ConfigurationCard::GetTVectorIndex(TString keyword, int tol) const{
  
  // Returns fIndex of a TVector according to its position in std::hash_map
  
  UInt_t i = 0;
  Int_t ind = -1;
  
  //cout<<"ConfigurationCard_SEARCH_MODE_HASHLIST"<<endl;
  TNamed * ko = (TNamed*)fKeyTable.FindObject( keyword.Data() );
  //cout<<ko<<endl;
  if(ko) ind = ko->GetUniqueID();
  if( ind == -1 ){
    for( UInt_t ii=0;ii<fKeyWordVector.size();ii++ ){
      if( fKeyWordVector[ii] == keyword ) return i;
    }
    if( tol == 0 ){
      cout << "ERROR: \""<<keyword.Data()<<"\" must be defined  "<< endl;
      exit(1);
    }else if ( tol == 1 ){
      cout << "Warning: \""<<keyword.Data()<<"\" is not exist. return default value  "<< endl;
      return -1;
    }else{
      return -1;
    }
  }
  
  return ind;
}

/*
 * Get the size of the vector corresponding to keyword
 */
int ConfigurationCard::GetN(TString keyword) const{
  //returns size of TVector
  unsigned int findex = GetTVectorIndex(keyword);
  return (int) fValuesVector[findex].GetNrows();
}

/*
 * Get the number of bins corresponding to keyword
 */
int ConfigurationCard::GetNBin(TString keyword) const{
  return GetN(keyword)-1;
}

/*
 * Get the vector corresponding to a keyword
 */
TVector *  ConfigurationCard::GetVector(TString keyword ){
  int findex = GetTVectorIndex(keyword);
  return &fValuesVector[findex];
}

/*
 * Get the vector component corresponding to a keyword and index in vector
 */
float ConfigurationCard::Get(TString keyword, int VectorComponent) const{
  //returns VectorComponent Component of  fValuesVector TVector for given keyword
  int findex = GetTVectorIndex(keyword);
  if(0<=VectorComponent && VectorComponent<GetNwithIndex(findex)){
    return fValuesVector[findex](VectorComponent+1);
  }else{
    cout<<"ERROR: fValuesVector findex out of range "<<keyword.Data()<<endl;
    cout << "   Max findex: " << GetN(keyword) -  1<< " Asked: " <<  VectorComponent << endl;
    exit(1);
  }
}

/*
 * Get a string corresponding to keyword
 */
TString ConfigurationCard::GetStr(TString keyword) const{
  int findex = GetTVectorIndex(keyword, 1);
  if( findex < 0  ) return TString("");
  return fValueString[findex];
}

/*
 * Get a bin index for a value in keyword
 */
int ConfigurationCard::GetBin(TString keyword, double value) const{
  
  // First find the dimension of vector and check that it is reasonable
  int dimension = GetN(keyword);
  if(dimension < 2) return -1; // We need to have at least two values in a vector so that it represents a bin
  
  // If the given value is below the lowest bin range, return -1
  if(value < Get(keyword)) return -1;
  
  // Check to which bin the value falls
  for(int i = 1; i < dimension; i++){
    if(value < Get(keyword,i)) return i-1;
  }
  
  // If value is exactly on the high limit, put it in the last bin
  if(value == Get(keyword,dimension-1)) return dimension-1;
  
  // If we are here, the value is larger than the highest bin limit. Return -1 in this case.
  return -1;
}

/*
 * Adding key to keytable
 */
void ConfigurationCard::AddToKeyTable(TString key, int index){
  TNamed * ko = new TNamed(key.Data(),"");
  ko->SetUniqueID(index);
  fKeyTable.Add(ko); // TODO check duplicate
}

/*
 * Read the input card
 */
void ConfigurationCard::ReadInputCard(){
  // read card
  
  char buffer[kMaxDimBuffer];
  ifstream incard;
  
  cout << "Reading configuration card from file: " << fCardName << endl;
  incard.open(fCardName,ios::in);
  
  if(!incard){
    cout<<"ERROR: Config file <"<<fCardName<<"> not found!"<<endl;
    exit(1);
  }
  
  while(!incard.eof()){ // Loop over the input configuration line by line
    
    incard.getline(buffer,kMaxDimBuffer); //read a line
    
    if(fnEntry > 1000){//is the file reasonably long?
      cout<<"Maximum number of 1000 lines reached in ConfigurationCard.C"<<endl;
      exit(1);
    }
    
    ReadInputLine( buffer );
    
    fnEntry++;
  }//while eof

}

/*
 *  Read a line from the input card
 */
void ConfigurationCard::ReadInputLine( const char *buffer ){
  // parse a line
  
  TString tstr(buffer); // Convert the line in the buffer to TString
  
  if( tstr.BeginsWith("#") ) return; //Skip comments
  tstr.ReplaceAll("\t"," "); // Get rid of tabulators
  
  // Remove comments in line
  Ssiz_t startOFcomment = tstr.First('#');
  if(startOFcomment>0){
    tstr.Remove(startOFcomment,tstr.Length() - startOFcomment);
  }
  
  // Remove white spaces from the beginning
  if(tstr.BeginsWith(" ")){
    Ssiz_t startOFkeyword = 0;
    while(1){
      TString s = tstr[startOFkeyword];
      if(s.CompareTo(" ")) break;
      startOFkeyword++;
    }
    tstr.Replace(0,startOFkeyword,"",0);
  }
  
  // Separate inputs
  TObjArray *lineContents = tstr.Tokenize(" ");
  
  if(lineContents->GetEntriesFast() < 1) return; // Skip empty lines
  
  //----- Read a keyword -----
  TString entryname = ((TObjString*)(lineContents->At(0)))->String(); // Read a keyword
  
  // Give a warning if no parameters set for a keyword
  if(lineContents->GetEntriesFast() == 1){
    cout<<"WARNING: No parameters given for keyword "<<entryname.Data()<<" on configuration card!"<<endl;
  } else {
    
    //----- Read parameters -----
    vector< float > items; // Auxiliary vector
    
    for(int i=1; i<lineContents->GetEntriesFast(); i++){ // Loop over the numbers
      TString token = ((TObjString*)(lineContents->At(i)))->String(); // Read number as a string
      
      if(token.IsFloat()){
        items.push_back(token.Atof()); // If string is float number, store it to vector
      }else{
        items.push_back(0); // TODO: Not sure if this is the wanted behavior...
        // cout<<"ERROR: char "<<token.Data()<<" among numbers"<<endl;
        // exit(1);
      }
    } // End of the line contents loop
    
    
    // Fill TVectors and Map
    fKeyWordVector.push_back( entryname.Data() ); // Put the new keyword at the end of the array
    
    
    fValuesVector.push_back( TVector( 1, items.size(), &items[0]) ); // Store TVector to array
    fValueString.push_back( ((TObjString*)(lineContents->At(1)))->String() );
    
    AddToKeyTable( entryname, fValuesVector.size()-1 );
    
  }//else
  
  lineContents->~TObjArray(); // Remove array from heap
}

/*
 *  Set a value for the git hash
 */
void ConfigurationCard::SetGitHash( const char *hash ){
  fGitHash = hash;
}

/*
 * Print the card to console
 */
void ConfigurationCard::PrintOut(){
  // echo
  cout<<endl<<"======== "<<fCardName<<" ========="<<endl;
  cout << "GitHash: " << fGitHash.String().Data() << endl;
  for(unsigned int i=0; i < fValuesVector.size(); i++){
    cout << Form("%15s",fKeyWordVector[i].Data()); //print keyword
    cout << " (dim = " << fValuesVector[i].GetNrows() << ") "; //print size of TVector
    for(int j = 1; j <= fValuesVector[i].GetNrows(); j++){
      cout << fValuesVector[i][j] << " "; //TVector components
    }
    cout << endl;
  }
}


void ConfigurationCard::WriteCard(TDirectory *file) const{
  // write
  cout<<endl<<"====== Writing into file ========="<<endl;
  
  if(!file->GetDirectory("JCard")) {
    file->mkdir("JCard");//directory to store input parameters
  }
  file->cd("JCard");
  fGitHash.Write("GitHash");
  for(unsigned int i=0;i<fValuesVector.size();i++){
    fValuesVector[i].Write(fKeyWordVector[i]);
    if(fKeyWordVector[i] == "MixingFileList") TObjString(fValueString[i]).Write("MixingFileListPath");
    if(fKeyWordVector[i] == "JetVetoMapFile" || fKeyWordVector[i] == "JetVetoMapHistogram" ||
       fKeyWordVector[i] == "JetVetoPtConvention")
      TObjString(fValueString[i]).Write(fKeyWordVector[i]+"String");
  }
}

// SingleJetCorrector
// v3.0
// Author: Yi Chen
// 
// This class applies JEC for any given level using TF1 as the workhorse
// Supposedly runs faster than v1.0
// v3.0: one can add list of text files to apply them one by one

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <limits>

#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

class SingleJetCorrector
{
private:
   enum Type { TypeNone, TypeJetPT, TypeJetEta, TypeJetPhi, TypeJetArea, TypeRho };
   bool Initialized;
   //bool IsFunction; // Edit: Removed this because compiler complains that it is not used
   double JetPT = std::numeric_limits<double>::quiet_NaN();
   double JetEta = std::numeric_limits<double>::quiet_NaN();
   double JetPhi = std::numeric_limits<double>::quiet_NaN();
   double JetArea = std::numeric_limits<double>::quiet_NaN();
   double Rho = std::numeric_limits<double>::quiet_NaN();
   std::vector<std::string> Formulas;
   std::vector<std::vector<double>> Parameters;
   std::vector<std::vector<Type>> BinTypes;
   std::vector<std::vector<double>> BinRanges;
   std::vector<std::vector<Type>> Dependencies;
   std::vector<std::vector<double>> DependencyRanges;
   std::vector<TF1 *> Functions;
public:
   SingleJetCorrector()                  { Initialized = false; }
   SingleJetCorrector(std::string File)  { Initialized = false; Initialize(File); }
   ~SingleJetCorrector()                 { for(auto P : Functions) if(P != nullptr) delete P; }
   void SetJetPT(double value)     { JetPT = value; }
   void SetJetEta(double value)    { JetEta = value; }
   void SetJetPhi(double value)    { JetPhi = value; }
   void SetJetArea(double value)   { JetArea = value; }
   void SetRho(double value)       { Rho = value; }
   void Initialize(std::string FileName);
   std::vector<std::string> BreakIntoParts(std::string Line);
   bool CheckDefinition(std::string Line);
   std::string StripBracket(std::string Line);
   SingleJetCorrector::Type ToType(std::string Line);
   double GetCorrection();
   double GetCorrectedPT();
   double GetValue(Type T);
private:
   std::string Hack4(std::string Formula, char V, int N);
};

class JetCorrector
{
private:
   std::vector<SingleJetCorrector> JEC;
   double JetPT = std::numeric_limits<double>::quiet_NaN();
   double JetEta = std::numeric_limits<double>::quiet_NaN();
   double JetPhi = std::numeric_limits<double>::quiet_NaN();
   double JetArea = std::numeric_limits<double>::quiet_NaN();
   double Rho = std::numeric_limits<double>::quiet_NaN();
public:
   JetCorrector()                               {}
   JetCorrector(std::string File)               { Initialize(File); }
   JetCorrector(std::vector<std::string> Files) { Initialize(Files); }
   void Initialize(std::string File)            { std::vector<std::string> X; X.push_back(File); Initialize(X); }
   void Initialize(std::vector<std::string> Files);
   void SetJetPT(double value)     { JetPT = value; }
   void SetJetEta(double value)    { JetEta = value; }
   void SetJetPhi(double value)    { JetPhi = value; }
   void SetJetArea(double value)   { JetArea = value; }
   void SetRho(double value)       { Rho = value; }
   double GetCorrection();
   double GetCorrectedPT();
   std::vector<double> GetCorrectedPTStages(); // raw followed by each actual sequential level
};

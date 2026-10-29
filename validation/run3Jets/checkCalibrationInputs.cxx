// Synthetic software guard, built with -fno-access-control only for this fixture.
#include "../../src/EECAnalyzer.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
void Require(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
void Set(ConfigurationCard& c,const char* key,float value){auto* v=c.GetVector(key);(*v)[1]=value;}
int main(int argc,char** argv){
 try{
  if(argc!=2)throw std::runtime_error("Pass a scratch directory");
  const std::filesystem::path out=argv[1];std::filesystem::create_directories(out);
  {
   ConfigurationCard input("cardEECOOInclusiveR04.input");Set(input,"BackgroundMethods",0);
   EECAnalyzer a({},&input,1,-1);a.fJetUncertainty=new JetUncertainty();
   for(int system: {1,3,9,10,11,12}){
    a.fIsOOData=system==9||system==10;a.fIsPp2024Data=system==11||system==12;
    a.ConfigureJetUncertainty(37,53,.123,-2.3);
    const bool run3=a.fIsOOData||a.fIsPp2024Data;
    Require(a.fJetUncertainty->JetPT==(run3?53:37),"JES pT input changed legacy or failed future Run-3 convention");
    Require(a.fJetUncertainty->JetEta==.123 && a.fJetUncertainty->JetPhi==-2.3,"JES changed caller's lookup coordinates");
   }
  }
  const auto one=out/"one.txt",two=out/"two.txt";
  {std::ofstream f(one);f<<"{1 JetEta 1 JetPt [0] Correction L2Relative}\n-5 5 3 1 1000 2\n";}
  {std::ofstream f(two);f<<"{1 JetEta 1 JetPt x/10 Correction L2Residual}\n-5 5 2 1 1000\n";}
  {
   JetCorrector sequential(std::vector<std::string>{one.string(),two.string()});
   sequential.SetJetPT(10);sequential.SetJetEta(0);sequential.SetJetPhi(0);sequential.SetRho(0);sequential.SetJetArea(0);
   Require(std::abs(sequential.GetCorrectedPT()-40)<1e-10,"Correction levels were not evaluated sequentially");
  }
  for(const char* name: {"cardEECOOInclusiveR04.input","cardEECOOMCInclusiveR04.input","cardEECpp5v36TeVR04.input","cardEECpp5v36TeVMCR04.input"}){
   ConfigurationCard card(name);Set(card,"BackgroundMethods",0);
   auto a=std::make_unique<EECAnalyzer>(std::vector<TString>{},&card,1,-1);
   a->RunAnalysis(); // No input files: inspect actual configured level sequence only.
   Require(a->fJetCorrector->JEC.size()==(card.Get("DataType")==11?2:1),"Residual was enabled for MC/OO or omitted from pp data");
   Require(card.Get("JetAxis")==1 && a->fJetAxis==1,"Legacy WTA convention was changed");
   Require(a->fJetUncertaintyMode==0,"Unsupported JES/JER was enabled");
  }
  std::cout<<"PASS: legacy WTA preserved; sequential raw-pT levels; data-only residuals; disabled JES/JER; future Run-3 corrected-pT JES input with legacy input unchanged\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}

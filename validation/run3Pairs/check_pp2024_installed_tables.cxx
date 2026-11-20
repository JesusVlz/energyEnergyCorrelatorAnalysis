// Native lookup software test, not physics closure.
#include "../../src/TrackPairEfficiencyCorrector.h"
#include <TFile.h>
#include <TH1.h>
#include <TVectorT.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv){
 if(argc!=2){std::cerr<<"Usage: check_pp2024_installed_tables TABLE.root\n";return 2;}
 try{
  TFile f(argv[1],"READ");if(f.IsZombie())throw std::runtime_error("Unreadable ROOT table");
  TrackPairEfficiencyCorrector r(argv[1],false,false);r.SetDisableCorrection(false);
  auto* pt=dynamic_cast<TVectorT<float>*>(f.Get("JCard/TrackPairPtBinEdges"));
  auto* jp=dynamic_cast<TVectorT<float>*>(f.Get("JCard/JetPtBinEdgesEEC"));
  if(!pt||!jp||!r.IsRun3Table())throw std::runtime_error("Missing Run-3 axes");
  auto mid=[](TVectorT<float>* v,int i){return .5*((*v)[v->GetLwb()+i]+(*v)[v->GetLwb()+i+1]);};
  long checks=0,fallbacks=0;
  for(int t=0;t<pt->GetNrows()-1;++t)for(int a=0;a<=t;++a){
   auto* inc=dynamic_cast<TH1*>(f.Get(Form("trackPairEfficiencyCorrection_C0T%dA%d",t,a)));
   if(!inc)throw std::runtime_error("Missing inclusive histogram");
   for(int j=-1;j<jp->GetNrows()-1;++j){
    auto* h=j<0?inc:dynamic_cast<TH1*>(f.Get(Form("trackPairEfficiencyCorrectionCloseToJet_C0T%dA%dJ%d",t,a,j)));
    if(!h)throw std::runtime_error("Missing differential histogram");
    for(int b=1;b<=h->GetNbinsX();++b){
     double ratio=h->GetBinContent(b);
     if(!std::isfinite(ratio)||ratio<=0){ratio=inc->GetBinContent(b);++fallbacks;}
     const double dr=h->GetXaxis()->GetBinCenter(b),jet=j<0?-1:mid(jp,j);
     if(!std::isfinite(ratio)||ratio<=0){
      bool rejected=false;
      try{r.GetTrackPairEfficiencyCorrection(dr,-.5,mid(pt,t),mid(pt,a),jet);}
      catch(const std::out_of_range&){rejected=true;}
      if(!rejected)throw std::runtime_error("Unsupported lookup accepted");
     }else{
      const auto got=r.GetTrackPairEfficiencyCorrection(dr,-.5,mid(pt,t),mid(pt,a),jet);
      const double expected=1.0/ratio;
      if(!std::isfinite(got.first)||std::abs(got.first-expected)>1e-12*std::max(1.0,std::abs(expected)))
       throw std::runtime_error("Native inverse differs from stored histogram");
     }
     ++checks;
    }
   }
  }
  for(int p=0;p<3;++p){
   bool rejected=false;
   try{r.GetTrackPairEfficiencyCorrection(p==0?2.:.1,-.5,4.,p==1?1.:3.,p==2?1.:50.);}
   catch(const std::out_of_range&){rejected=true;}
   if(!rejected)throw std::runtime_error("Out-of-domain lookup accepted");
   ++checks;
  }
  std::cout<<"PASS native lookup checks="<<checks<<" fallback_probes="<<fallbacks<<'\n';
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
 return 0;
}

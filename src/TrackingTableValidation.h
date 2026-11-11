#ifndef EEC_TRACKING_TABLE_VALIDATION_H
#define EEC_TRACKING_TABLE_VALIDATION_H
#include <cmath>
#include <stdexcept>
#include <string>
#include <initializer_list>
#include <TAxis.h>
#include <TH1.h>

namespace EECRun3Tracking {
inline bool InsideAxis(const TAxis* axis, double x) {
  return axis && std::isfinite(x) && x >= axis->GetXmin() && x < axis->GetXmax();
}
inline void SameAxis(const TAxis* a, const TAxis* b, const std::string& context) {
  if(!a || !b || a->GetNbins() != b->GetNbins())
    throw std::runtime_error(context + ": tracking table axis dimensions differ");
  for(int i=1;i<=a->GetNbins()+1;++i) {
    const double x=a->GetBinLowEdge(i), y=b->GetBinLowEdge(i);
    if(!std::isfinite(x) || !std::isfinite(y) || std::abs(x-y)>1e-9 ||
       (i>1 && x<=a->GetBinLowEdge(i-1)))
      throw std::runtime_error(context + ": tracking table axis edges differ or are invalid");
  }
}
inline void ValidateTables(std::initializer_list<const TH1*> tables, int dimension, const std::string& context) {
  const TH1* first = *tables.begin();
  if(!first || first->GetDimension()!=dimension)
    throw std::runtime_error(context + ": invalid tracking histogram dimension");
  for(const TH1* h:tables) {
    if(!h || h->GetDimension()!=dimension) throw std::runtime_error(context + ": missing or incompatible tracking histogram");
    SameAxis(first->GetXaxis(),h->GetXaxis(),context);
    SameAxis(first->GetYaxis(),h->GetYaxis(),context);
    if(dimension==3) SameAxis(first->GetZaxis(),h->GetZaxis(),context);
  }
  // Verify the coordinate ordering used by the wrappers can cover the stated
  // acceptance. No interpolation, extrapolation or statistical smoothing.
  if(first->GetXaxis()->GetXmin()>-2.4 || first->GetXaxis()->GetXmax()<2.4 ||
     !InsideAxis(first->GetYaxis(),2.0) || !InsideAxis(first->GetYaxis(),6.1))
    throw std::runtime_error(context + ": expected X=detector eta and Y=track pT coverage");
  if(dimension==3 && (!InsideAxis(first->GetZaxis(),0.) || !InsideAxis(first->GetZaxis(),199.)))
    throw std::runtime_error(context + ": expected Z=hiBin coverage 0 through 199");
}
inline bool ValidContents(double efficiency,double fake,double secondary,double multiple) {
  // Efficiencies slightly above one are not clipped; invalid rates are rejected.
  return std::isfinite(efficiency) && efficiency>0 &&
         std::isfinite(fake) && fake>=0 && fake<1 &&
         std::isfinite(secondary) && secondary>=0 && secondary<1 &&
         std::isfinite(multiple) && multiple>=0;
}
} // namespace EECRun3Tracking
#endif

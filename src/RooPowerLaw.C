#include "RooPowerLaw.h"
/// Advertise analytical integral over fX
Int_t RooPowerLaw::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analyticVars, const char * rangeName) const
{
  if(matchArgs(allVars,analyticVars,fX) && fX.min(rangeName) > 0.0){
    return 1;
  }
  else{
    return 0;
  }
}

/// Evaluate analytical integral over fX.
Double_t RooPowerLaw::analyticalIntegral(Int_t code, const char * rangeName) const {
  if(code == 1){

    double integral(0.0);
    double minX(std::max(0.0, fX.min(rangeName)));
    double maxX(std::max(0.0, fX.max(rangeName)));

    if(minX != maxX){
  //std::cout << "MinX = " <<  minX << ", MaxX = " << maxX << std::endl;
  if(fIndex == -1.0){
    integral = fNormalisation*(std::log(maxX) - std::log(minX));
  }
  else{
    double newIndex(fIndex + 1.0);
    integral = (std::pow(maxX, newIndex) - std::pow(minX, newIndex))/newIndex;
  }
    }
    return integral;
  }
  assert(0);
  return 0;
}

void RooPowerLaw::generateEvent(Int_t code){
  // should never be called unless code == 1
  assert(code == 1);

  double index = fIndex;
  double uni = RooRandom::randomGenerator()->Uniform();

  double xgen = std::pow((std::pow(fX.max(), index + 1.0) - std::pow(fX.min(), index + 1))*uni + std::pow(fX.min(), index + 1.0), 1.0/(index+1.0));

  fX = xgen;
}

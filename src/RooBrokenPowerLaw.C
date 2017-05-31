#include "RooBrokenPowerLaw.h"
/// Advertise analytical integral over fX
Int_t RooBrokenPowerLaw::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analyticVars, const char * rangeName) const
{
  if(matchArgs(allVars,analyticVars,fX) && fX.min(rangeName) > 0.0){
    return 1;
  }
  else{
    return 0;
  }
}

Double_t RooBrokenPowerLaw::CalcPLIntegral(double min, double max, double index) const {

  double integral(0.0);
  if(index == -1){
    integral = std::log(max) - std::log(min);
  }
  else{
    double newIndex(index + 1.0);
    integral = (std::pow(max, newIndex) - std::pow(min, newIndex))/newIndex;
  }

  return integral;
}

/// Evaluate analytical integral over fX.
Double_t RooBrokenPowerLaw::analyticalIntegral(Int_t code, const char * rangeName) const {
  if(code == 1){

    double integral(0.0);
    double minX(std::max(0.0, fX.min(rangeName)));
    double maxX(std::max(0.0, fX.max(rangeName)));

    if(minX != maxX){
  //std::cout << "MinX = " <<  minX << ", MaxX = " << maxX << std::endl;
  if(maxX < fXBreak){ // simple power law with first index
    integral = CalcPLIntegral(minX, maxX, fIndex1);
  }
  else if(minX > fXBreak){ // scaled power law with second index
    integral = std::pow(fXBreak, fIndex1-fIndex2)*CalcPLIntegral(minX, maxX, fIndex2);
  }
  else{ // broken power law
    integral = CalcPLIntegral(minX, fXBreak, fIndex1) + std::pow(fXBreak, fIndex1-fIndex2)*CalcPLIntegral(fXBreak, maxX, fIndex2);
  }
    }
    return integral;
  }
  assert(0);
  return 0;
}

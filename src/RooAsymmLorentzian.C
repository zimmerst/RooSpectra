#include "RooAsymmLorentzian.h"

/// Advertise analytical integral over fX.
Int_t RooAsymmLorentzian::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analyticVars, const char * ) const {
  if(matchArgs(allVars,analyticVars,fX)){
    return 1;
  }
  else{
    return 0;
  }
}

  /// Evaluate analytical integral over fX.
Double_t RooAsymmLorentzian::analyticalIntegral(Int_t code, const char * rangeName) const {
  if(code == 1){

    double integral(0.0);

    if(fX.max(rangeName) < fMean){ // we are only integrating left of the peak
  integral = M_PI*fGammaL*(ROOT::Math::cauchy_cdf(fX.max(rangeName), fGammaL, fMean) - ROOT::Math::cauchy_cdf(fX.min(rangeName), fGammaL, fMean));
    }
    else if(fX.min(rangeName) > fMean){ // we are only integrating above the peak
  integral = M_PI*fGammaR*(ROOT::Math::cauchy_cdf(fX.max(rangeName), fGammaR, fMean) - ROOT::Math::cauchy_cdf(fX.min(rangeName), fGammaR, fMean));
    }
    else{ // the integration region spans the peak
  integral = M_PI*fGammaL*(ROOT::Math::cauchy_cdf(0.0, fGammaL) - ROOT::Math::cauchy_cdf(fX.min(rangeName), fGammaL, fMean)); // left hand part
  integral += M_PI*fGammaR*(ROOT::Math::cauchy_cdf(fX.max(rangeName), fGammaR, fMean) - ROOT::Math::cauchy_cdf(0.0, fGammaR)); // right hand part
    }
    return integral;
  }

  assert(0);
  return 0;
}

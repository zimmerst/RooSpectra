#ifndef ROOASYMMLORENTZIAN
#define ROOASYMMLORENTZIAN

#include <RooAbsPdf.h>
#include <RooRealProxy.h>
#include <RooAbsReal.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/PdfFuncMathCore.h>

#include <cassert>
#include <iostream>

class RooAsymmLorentzian : public RooAbsPdf {

protected :

  RooRealProxy fX;
  RooRealProxy fMean;
  RooRealProxy fGammaL;
  RooRealProxy fGammaR;

  Double_t evaluate() const {

    double leftOrRight = fX-fMean;

    if(leftOrRight < 0.0){ // left
      return M_PI*fGammaL*ROOT::Math::cauchy_pdf(fX, fGammaL, fMean);
    }
    else{ // right
      return M_PI*fGammaR*ROOT::Math::cauchy_pdf(fX, fGammaR, fMean);
    }
  }

public :

  /// Default constructor
  RooAsymmLorentzian():
    RooAbsPdf("RooAsymmLorentzian", "Asymmetric Lorentzian Pdf")
  {}

  /// Constructor
  RooAsymmLorentzian(const  char * name, const char * title, RooAbsReal& x, RooAbsReal& mean, RooAbsReal& gammaL, RooAbsReal& gammaR):
    RooAbsPdf(name, title),
    fX("fX", "Dependent", this, x),
    fMean("fMean", "Mean", this, mean),
    fGammaL("fGammaL", "Gamma left", this, gammaL),
    fGammaR("fGammaR", "Gamma right", this, gammaR)
  {}

  /// Copy constructor
  RooAsymmLorentzian(const RooAsymmLorentzian& other, const char* name = 0):
    RooAbsPdf(other, name),
    fX("fX", this, other.fX),
    fMean("fMean", this, other.fMean),
    fGammaL("fGammaL", this, other.fGammaL),
    fGammaR("fGammaR", this, other.fGammaR)
  {}

  /// Destructor
  virtual ~RooAsymmLorentzian(){}

  /// Advertise analytical integral over fX.
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analyticVars, const char * ) const;

  /// Evaluate analytical integral over fX.
  Double_t analyticalIntegral(Int_t code, const char * rangeName) const;


  Double_t peakArea(){
    return 0.5*M_PI*(fGammaR + fGammaL);
  }

  /// Clone method.
  virtual TObject * clone(const char * newname = "") const { 
    return new RooAsymmLorentzian(*this, newname); 
  }

  ClassDef(RooAsymmLorentzian, 1);

};

ClassImp(RooAsymmLorentzian)

#endif /* ROOASYMMLORENTZIAN */

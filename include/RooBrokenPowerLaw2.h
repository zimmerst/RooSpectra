#ifndef RooBrokenPowerLaw2
#define RooBrokenPowerLaw2

#include <RooAbsPdf.h>
#include <RooRealProxy.h>
#include <RooAbsReal.h>
#include <RooRandom.h>

#include <cassert>
#include <cmath>

class RooBrokenPowerLaw2 : public RooAbsPdf {

protected :

  RooRealProxy fX;
  RooRealProxy fIndex1;
  RooRealProxy fIndex2;
  RooRealProxy fXBreak;
  RooRealProxy falpha;
  RooRealProxy fNormalisation
  Double_t evaluate() const {
    double value(0.0);
    if (fX > 0.0 && falpha != 0){
      value = fNormalization * std::pow(fX/fXBreak,-fIndex1) * std::pow((1+std::pow(E/fXBreak, 1/falpha)), -(fIndex2-fIndex1)*falpha);
    }
    return value;
  }

public :

  /// Default constructor.
  RooBrokenPowerLaw2():
    RooAbsPdf("RooBrokenPowerLaw2", "Broken Power Law Pdf, HESS parametrization")
  {}

  /// Constructor.
  RooBrokenPowerLaw2(const  char * name, const char * title, RooAbsReal& x, RooAbsReal& norm, RooAbsReal& index1,
                     RooAbsReal& index2, RooAbsReal& xBreak, RooAbsReal& alpha):
    RooAbsPdf(name, title),
    fX("fX", "Dependent", this, x),
    fNormalization("fNormalization", "Normalization", this, norm),
    fIndex1("fIndex1", "First Index", this, index1),
    fIndex2("fIndex2", "Second Index", this, index2),
    fXBreak("fXBreak", "Break Ordinate", this, xBreak),
    falpha("falpha", "sharpness indicator", this, alpha)
  {}

  /// Copy constructor.
  RooBrokenPowerLaw2(const RooBrokenPowerLaw2& other, const char* name = 0):
    RooAbsPdf(other, name),
    fX("fX", this, other.fX),
    fIndex1("fNormalization", this, other.fNormalization),
    fIndex1("fIndex1", this, other.fIndex1),
    fIndex2("fIndex2", this, other.fIndex2),
    fXBreak("fXBreak", this, other.fXBreak),
    falpha("falpha", this, other.falpha)
  {}

  /// Destructor.
  virtual ~RooBrokenPowerLaw2(){}

  /// Advertise analytical integral over fX
  //Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analyticVars, const char * rangeName) const ;

  Double_t CalcPLIntegral(double min, double max, double index) const ;

  /// Evaluate analytical integral over fX.
  //Double_t analyticalIntegral(Int_t code, const char * rangeName) const;
#ifdef EXTRAS

  // Advertise event generation capability for fX.
  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const {
    if (matchArgs(directVars,generateVars,fX)) return 1 ;  
    return 0 ;
  }

  /** Implement event genration for fX. See http://mathworld.wolfram.com/RandomNumber.html for 
   * method details.
   */
/*  void generateEvent(Int_t code){
    // should never be called unless code == 1
    assert(code == 1);

    double index = fIndex;
    double uni = RooRandom::randomGenerator()->Uniform();

    double xgen = std::pow((std::pow(fX.max(), index + 1.0) - std::pow(fX.min(), index + 1))*uni + std::pow(fX.min(), index + 1.0), 1.0/(index+1.0));

    fX = xgen;
  }
*/
#endif
    
  /// Clone method.
  virtual TObject * clone(const char * newname = "") const { 
    return new RooBrokenPowerLaw2(*this, newname); 
  }

  ClassDef(RooBrokenPowerLaw2, 1);

};

#endif /* RooBrokenPowerLaw2 */

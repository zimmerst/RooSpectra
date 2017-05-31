#ifndef ROOPOWERLAW
#define ROOPOWERLAW

#include <RooAbsPdf.h>
#include <RooRealProxy.h>
#include <RooAbsReal.h>
#include <RooRandom.h>

#include <cassert>
#include <cmath>

class RooPowerLaw : public RooAbsPdf {

protected :

  RooRealProxy fX;
  RooRealProxy fIndex;
  RooRealProxy fNormalisation;

  Double_t evaluate() const {
    return fX > 0.0 ? fNormalisation*std::pow(fX, fIndex) : 0.0;
  }

public :

  /// Default constructor.
  RooPowerLaw():
    RooAbsPdf("RooPowerLaw", "Power Law Pdf")
  {}

  /// Constructor.
  RooPowerLaw(const  char * name, const char * title, RooAbsReal& x, RooAbsReal& index, RooAbsReal& normalisation):
    RooAbsPdf(name, title),
    fX("fX", "Dependent", this, x),
    fIndex("fIndex", "Index", this, index),
    fNormalisation("fNormalisation", "Normalisation", this, normalisation)
  {}

  /// Copy constructor.
  RooPowerLaw(const RooPowerLaw& other, const char* name = 0):
    RooAbsPdf(other, name),
    fX("fX", this, other.fX),
    fIndex("fIndex", this, other.fIndex),
    fNormalisation("fNormalisation", this, other.fNormalisation)
  {}

  /// Destructor.
  virtual ~RooPowerLaw(){}

  /// Advertise analytical integral over fX
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analyticVars, const char * rangeName) const;

  /// Evaluate analytical integral over fX.
  Double_t analyticalIntegral(Int_t code, const char * rangeName) const;
  // Advertise event genration capability for fX.
  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const {
    if (matchArgs(directVars,generateVars,fX)) return 1 ;  
    return 0 ;
  }

  /** Implement event genration for fX. See http://mathworld.wolfram.com/RandomNumber.html for 
   * method details.
   */
  void generateEvent(Int_t code);
    
  /// Clone method.
  virtual TObject * clone(const char * newname = "") const { 
    return new RooPowerLaw(*this, newname); 
  }

  ClassDef(RooPowerLaw, 1);

};

#endif /* ROOPOWERLAW */

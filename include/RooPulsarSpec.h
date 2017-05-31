#ifndef ROOPULSARSPEC
#define ROOPULSARSPEC

#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooRealProxy.h>

#include <cmath>

/** \class RooPulsarSpec
 *
 * \brief Class which represents a generic parameterisation of Pulsar-Like photon 
 * spectrum.
 */
class RooPulsarSpec: public RooAbsPdf {

public :
  /** \brief Default constructor 
   */
  RooPulsarSpec(){}

  /** \brief Constructor.
   *
   * \param name - A unique textual label for this isnatnce of the PDF.
   * \param title - A title (used in plots) for the PDF.
   * \param energy - The emitted photon energy in TeV.
   * \param gamma - Power Law component index.
   * \param beta - Exponential component index.
   * \param cutoff - Exponential component cutoff energy.
   */
 RooPulsarSpec(const char * name, const char * title, RooAbsReal& energy, RooAbsReal& gamma, RooAbsReal& beta, RooAbsReal& cutoff):
  RooAbsPdf(name, title),
    fEnergy("fEnergy", "Photon energy", this, energy),
    fGamma("fGamma", "Power Law component index", this, gamma),
    fBeta("fBeta", "Exponential component index", this, beta),
    fCutoff("fCutoff", "Exponential component cutoff energy", this, cutoff)
  {}

  /** \brief Copy constructor.
   */
  RooPulsarSpec(const RooPulsarSpec& other, const char * name):
  RooAbsPdf(other, name),
    fEnergy("fEnergy", this, other.fEnergy),
    fGamma("fGamma", this, other.fGamma),
    fBeta("fBeta", this, other.fBeta),
    fCutoff("fCutoff", this, other.fCutoff)
  {}

  /** \brief Destructor.
   */
  virtual ~RooPulsarSpec(){}

  /** \brief Clone function must be defined so that RooPulsarSpec is not pure virtual!
   */
  virtual TObject* clone(const char* newname) const { return new RooPulsarSpec(*this,newname); }

protected : 

  /// The photon energy in TeV.
  RooRealProxy fEnergy;
  /// The Power Law component index.
  RooRealProxy fGamma;
  /// The Exponential component index.
  RooRealProxy fBeta;
  /// The Exponential component cutoff energy in TeV.
  RooRealProxy fCutoff;

  /** \brief Evaluate the value of the PDF.
   *
   * \return The differerential number of emitted photons per energy interval (\f$ dN/dE \f$)
   */
  Double_t evaluate() const {

    double gamma = fGamma;
    double beta = fBeta;
    double cutoff = fCutoff;
    double energy = fEnergy;

    if(energy < 0){
      return 1e-3;
    }

    double value = std::pow(energy, -gamma)*std::exp(-std::pow(energy/cutoff, beta));

    return value;

  }

private :

  ClassDef(RooPulsarSpec,1) // Pulsar Spectrum PDF

};

ClassImp(RooPulsarSpec);

#endif /* ROOPULSARSPEC */


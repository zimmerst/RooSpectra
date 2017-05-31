#ifndef ROODMSPEC
#define ROODMSPEC

#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooRealProxy.h>

#include <cmath>

/** \class RooDMSpec
 *
 * \brief Class which represents a generic parameterisation of a Dark-Matter-Like photon 
 * spectrum as described by Bergstrom.
 */
class RooDMSpec: public RooAbsPdf {

public :
  /** \brief Default constructor 
   */
  RooDMSpec(){}

  /** \brief Constructor.
   *
   * \param name - A unique textual label for this isnatnce of the PDF.
   * \param title - A title (used in plots) for the PDF.
   * \param energy - The emitted photon energy in TeV.
   * \param mass - The DM particle mass in TeV.
   */
  RooDMSpec(const char * name, const char * title, RooAbsReal& energy, RooAbsReal& mass):
    RooAbsPdf(name, title),
    fEnergy("fEnergy", "Photon energy", this, energy),
    fMass("fMass", "Particle mass", this, mass)
  {}

  /** \brief Copy constructor.
   */
  RooDMSpec(const RooDMSpec& other, const char * name):
    RooAbsPdf(other, name),
    fEnergy("fEnergy", this, other.fEnergy),
    fMass("fMass", this, other.fMass)
  {}

  /** \brief Destructor.
   */
  virtual ~RooDMSpec(){}

  /** \brief Clone function must be defined so that RooDMSpec is not pure virtual!
   */
  virtual TObject* clone(const char* newname) const { return new RooDMSpec(*this,newname); }

protected : 

  /// The photon energy in TeV.
  RooRealProxy fEnergy;
  /// The DM particle mass in TeV.
  RooRealProxy fMass;

  /** \brief Evaluate the value of the PDF.
   *
   * \return The differerential number of emitted photons per energy interval (\f$ dN/dE \f$)
   */
  Double_t evaluate() const {

    double mass = fMass;
    double energy = fEnergy;
    double energyDatum = 2.0e-4;

    double value = 0.73*std::exp(-7.8*energy/mass);
    value /= (mass*(std::pow(energy/mass, 1.5) + energyDatum));

    return value;

  }

private :

  ClassDef(RooDMSpec,1) // Dark Matter Spectrum PDF

};

ClassImp(RooDMSpec);

#endif /* ROODMSPEC */


#ifndef ROODMSPEC2
#define ROODMSPEC2

#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooRealProxy.h>

#include <cmath>

/** \class RooDMSpec2
 *
 * \brief Class which represents a generic parameterisation of a Dark-Matter-Like photon 
 * spectrum as described by Fornengo, Pieri and Scopel (2004) arXiv:hep-ph/0407342v4.
 */
class RooDMSpec2: public RooAbsPdf {
    
    public :
    /** \brief Default constructor 
     */
    RooDMSpec2(){}
    
    /** \brief Constructor.
     *
     * \param name - A unique textual label for this instance of the PDF.
     * \param title - A title (used in plots) for the PDF.
     * \param energy - The emitted photon energy in TeV.
     * \param mass - The DM particle mass in TeV.
     */
    RooDMSpec2(const char * name, const char * title, RooAbsReal& energy, RooAbsReal& mass):
    RooAbsPdf(name, title),
    fEnergy("fEnergy", "Photon energy", this, energy),
    fMass("fMass", "Particle mass", this, mass)
    {}
    
    /** \brief Copy constructor.
     */
    RooDMSpec2(const RooDMSpec2& other, const char * name):
    RooAbsPdf(other, name),
    fEnergy("fEnergy", this, other.fEnergy),
    fMass("fMass", this, other.fMass)
    {}
    
    /** \brief Destructor.
     */
    virtual ~RooDMSpec2(){}
    
    /** \brief Clone function must be defined so that RooDMSpec2 is not pure virtual!
     */
    virtual TObject* clone(const char* newname) const { return new RooDMSpec2(*this,newname); }
    
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
        
        double x = energy/mass;
        
        // parameters are as defined for a b-quark final state in Fornengo, Pieri and Scopel (2004) arXiv:hep-ph/0407342v4
        // "Neutralino annihilation into gamma-rays in the Milky Way and in external galaxies"
        double eta = 1.0;
        double a = -1.5;
        double b = 0.48;
        double c = -16.87;
        double d = 21.09;
        double e = -22.49;
                
        return eta*std::pow(x, a)*std::exp(b + c*x + d*x*x + e*x*x*x);;
        
    }
    
    private :
    
    ClassDef(RooDMSpec2,1) // Dark Matter Spectrum PDF
    
};

ClassImp(RooDMSpec2);

#endif /* ROODMSPEC2 */


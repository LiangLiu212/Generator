#ifndef G4INCLParticleTable_hh
#define G4INCLParticleTable_hh 1

#include <string>
#include <vector>
#include <cassert>

#include "G4INCLParticleType.hh"
#include "G4INCLParticleSpecies.hh"
#include "G4INCLLogger.hh"
#include "G4INCLConfig.hh"
#include "G4INCLHFB.hh"

#ifdef INCLXX_IN_GEANT4_MODE
#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
#endif
#include "G4INCLGlobals.hh"
#include "G4INCLNaturalIsotopicDistributions.hh"

namespace G4INCL {

  namespace ParticleTable {

    const int maxClusterMass = 12;
    const int maxClusterCharge = 8;

    const int clusterTableZSize = maxClusterCharge+1;
    const int clusterTableASize = maxClusterMass+1;
    const int clusterTableSSize = 4;

    const double effectiveNucleonMass = 938.2796;
    const double effectiveNucleonMass2 = 8.8036860777616e5;
    const double effectiveDeltaMass = 1232.0;
    const double effectiveDeltaWidth = 130.0;
    const double effectivePionMass = 138.0;
    const double effectiveLambdaMass = 1115.683;
    const double effectiveSigmaMass = 1197.45; // max value
    const double effectiveXiMass = 1321.71; // max value
    const double effectiveKaonMass = 497.614; // max value
    const double effectiveAntiKaonMass = 497.614; // max value
    const double effectiveEtaMass = 547.862;
    const double effectiveOmegaMass = 782.65;
    const double effectiveEtaPrimeMass = 957.78;
    const double effectivePhotonMass = 0.0;
    extern G4ThreadLocal double minDeltaMass;
    extern G4ThreadLocal double minDeltaMass2;
    extern G4ThreadLocal double minDeltaMassRndm;

    /// \brief Initialize the particle table
    void initialize(Config const * const theConfig = 0);

    /// \brief Get the isospin of a particle
    int getIsospin(const ParticleType t);

    /// \brief Get the native INCL name of the particle
    std::string getName(const ParticleType t);

    /// \brief Get the short INCL name of the particle
    std::string getShortName(const ParticleType t);

    /// \brief Get the native INCL name of the particle
    std::string getName(const ParticleSpecies &s);

    /// \brief Get the short INCL name of the particle
    std::string getShortName(const ParticleSpecies &s);

    /// \brief Get the native INCL name of the ion
    std::string getName(const int A, const int Z);

    /// \brief Get the native INCL name of the ion
    std::string getName(const int A, const int Z, const int S);

    /// \brief Get the short INCL name of the ion
    std::string getShortName(const int A, const int Z);

    /// \brief Get INCL nuclear mass (in MeV/c^2)
    double getINCLMass(const int A, const int Z, const int S);

    /// \brief Get INCL particle mass (in MeV/c^2)
    double getINCLMass(const ParticleType t);

#ifndef INCLXX_IN_GEANT4_MODE
    /// \brief Do we have this particle mass?
    double hasMassTable(const unsigned int A, const unsigned int Z);

    /** \brief Weizsaecker mass formula
     *
     * Return the nuclear mass, as calculated from Weizsaecker's mass formula.
     * Adapted from the Geant4 source.
     *
     * \param A the mass number
     * \param Z the charge number
     * \return the nuclear mass [MeV/c^2]
     */
    double getWeizsaeckerMass(const int A, const int Z);
#endif

    ///\brief Get particle mass (in MeV/c^2)
    double getRealMass(const G4INCL::ParticleType t);
    ///\brief Get nuclear mass (in MeV/c^2)
    double getRealMass(const int A, const int Z, const int S = 0);

    /**\brief Get Q-value (in MeV/c^2)
     *
     * Uses the getTableMass function to compute the Q-value for the
     * following reaction:
     * \f[ (A_1,Z_1) + (A_2, Z_2) --> (A_1+A_2,Z_1+Z_2) \f]
     */
    double getTableQValue(const int A1, const int Z1, const int S1, const int A2, const int Z2, const int S2);

    /**\brief Get Q-value (in MeV/c^2)
     *
     * Uses the getTableMass function to compute the Q-value for the
     * following reaction:
     * \f[ (A_1,Z_1) + (A_2, Z_2) --> (A_3,Z_3) + (A1+A2-A3,Z1+Z2-Z3) \f]
     */
    double getTableQValue(const int A1, const int Z1, const int S1, const int A2, const int Z2, const int S2, const int A3, const int Z3, const int S3);

    double getTableSpeciesMass(const ParticleSpecies &p);

    /// \brief Get mass number from particle type
    int getMassNumber(const ParticleType t);

    /// \brief Get charge number from particle type
    int getChargeNumber(const ParticleType t);
    
    /// \brief Get strangeness number from particle type
    int getStrangenessNumber(const ParticleType t);

    double getNuclearRadius(const ParticleType t, const int A, const int Z);
    double getLargestNuclearRadius(const int A, const int Z);
    double getRadiusParameter(const ParticleType t, const int A, const int Z);
    double getMaximumNuclearRadius(const ParticleType t, const int A, const int Z);
    double getSurfaceDiffuseness(const ParticleType t, const int A, const int Z);

    /// \brief Return the RMS of the momentum distribution (light clusters)
    double getMomentumRMS(const int A, const int Z);

    /// \brief Return INCL's default separation energy
    double getSeparationEnergyINCL(const ParticleType t, const int /*A*/, const int /*Z*/);

    /// \brief Return the real separation energy
    double getSeparationEnergyReal(const ParticleType t, const int A, const int Z);

    /// \brief Return the real separation energy only for light nuclei
    double getSeparationEnergyRealForLight(const ParticleType t, const int A, const int Z);

    /// \brief Getter for protonSeparationEnergy
    double getProtonSeparationEnergy();

    /// \brief Getter for neutronSeparationEnergy
    double getNeutronSeparationEnergy();

    /// \brief Getter for antiprotonSeparationEnergy
    double getantiProtonSeparationEnergy();

    /// \brief Getter for antineutronSeparationEnergy
    double getantiNeutronSeparationEnergy();

     /// \brief Getter for antilambdaSeparationEnergy
    double getantiLambdaSeparationEnergy();

     /// \brief Getter for lambdaSeparationEnergy
    double getLambdaSeparationEnergy();

    /// \brief Setter for protonSeparationEnergy
    void setProtonSeparationEnergy(const double s);

    /// \brief Setter for protonSeparationEnergy
    void setNeutronSeparationEnergy(const double s);

     /// \brief Setter for lambdaSeparationEnergy
    void setLambdaSeparationEnergy(const double s);

    /// \brief Setter for antilambdaSeparationEnergy
    void setantiLambdaSeparationEnergy(const double s);

    /// \brief Get the name of the element from the atomic number
    std::string getElementName(const int Z);

    /// \brief Get the name of an unnamed element from the IUPAC convention
    std::string getIUPACElementName(const int Z);

    /// \brief Get the name of the element from the atomic number
    int parseElement(std::string pS);

    /** \brief Parse a IUPAC element name
     *
     * Note: this function is UGLY. Look at it at your own peril.
     *
     * \param pS a normalised string (lowercase)
     * \return the charge number of the nuclide, or zero on fail
     */
    int parseIUPACElement(std::string const &pS);

    IsotopicDistribution const &getNaturalIsotopicDistribution(const int Z);

    int drawRandomNaturalIsotope(const int Z);

    // Typedefs and pointers for transparent handling of mass functions
    //typedef double (*NuclearMassFn)(const int, const int);
    typedef double (*NuclearMassFn)(const int, const int, const int);
    typedef double (*ParticleMassFn)(const ParticleType);
    /// \brief Static pointer to the mass function for nuclei
    extern G4ThreadLocal NuclearMassFn getTableMass;
    /// \brief Static pointer to the mass function for particles
    extern G4ThreadLocal ParticleMassFn getTableParticleMass;

    // Typedefs and pointers for transparent handling of separation energies
    typedef double (*SeparationEnergyFn)(const ParticleType, const int, const int);
    /// \brief Static pointer to the separation-energy function
    extern G4ThreadLocal SeparationEnergyFn getSeparationEnergy;

    // Typedefs and pointers for transparent handling of Fermi momentum
    typedef double (*FermiMomentumFn)(const int, const int);
    extern G4ThreadLocal FermiMomentumFn getFermiMomentum;

    /// \brief Return the constant value of the Fermi momentum
    double getFermiMomentumConstant(const int /*A*/, const int /*Z*/);

    /** \brief Return the constant value of the Fermi momentum - special for light
     *
     * This function should always return PhysicalConstants::Pf for heavy
     * nuclei, and values from the momentumRMS table for light nuclei.
     *
     * \param A mass number
     * \param Z charge number
     */
    double getFermiMomentumConstantLight(const int A, const int Z);

    /** \brief Return the value Fermi momentum from a fit
     *
     * This function returns a fitted Fermi momentum, based on data from Moniz
     * et al., Phys. Rev. Lett. 26 (1971) 445. The fitted functional form is
     * \f[
     * p_F(A)=\alpha-\beta\cdot e^{(-A\cdot\gamma)}
     * \f]
     * with \f$\alpha=259.416\f$ MeV/\f$c\f$, \f$\beta=152.824\f$ MeV/\f$c\f$
     * and \f$\gamma=9.5157\cdot10^{-2}\f$.
     *
     * \param A mass number
     */
    double getFermiMomentumMassDependent(const int A, const int /*Z*/);

    /** \brief Get the value of the r-p correlation coefficient
     *
     * \param t the type of the particle (Proton or Neutron)
     * \return the value of the r-p correlation coefficient
     */
    double getRPCorrelationCoefficient(const ParticleType t);
    
    /// \brief Get the configuration of src-pair correlations
    bool getsrcPairConfig();
    
    /// \brief Get the distance between src nucleons
    float getsrcPairDistance();

    /// \brief Get the thickness of the neutron skin
    double getNeutronSkin();

    /// \brief Get the size of the neutron halo
    double getNeutronHalo();

    /// \brief Get the type of pion
    ParticleType getPionType(const int isosp);

    /// \brief Get the type of nucleon
    ParticleType getNucleonType(const int isosp);

    /// \brief Get the type of delta
    ParticleType getDeltaType(const int isosp);

    /// \brief Get the type of sigma
    ParticleType getSigmaType(const int isosp);
    
    /// \brief Get the type of kaon
    ParticleType getKaonType(const int isosp);

    /// \brief Get the type of antikaon
    ParticleType getAntiKaonType(const int isosp);
       
    /// \brief Get the type of xi
    ParticleType getXiType(const int isosp);
    
    /// \brief Get the type of antinucleon
    ParticleType getAntiNucleonType(const int isosp);

    /// \brief Get the type of antidelta
    ParticleType getAntiXiType(const int isosp);

    /// \brief Get the type of antisigma
    ParticleType getAntiSigmaType(const int isosp); 

    /// \brief Get particle width (in s)
    double getWidth(const ParticleType t);
  }
}

#endif


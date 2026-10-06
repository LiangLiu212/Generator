#ifndef G4INCLNuclearDensity_hh
#define G4INCLNuclearDensity_hh 1

#include <vector>
#include <map>
#include <cassert>
#include "G4INCLThreeVector.hh"
#include "G4INCLIFunction1D.hh"
#include "G4INCLParticle.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLRandom.hh"
#include "G4INCLINuclearPotential.hh"
#include "G4INCLInterpolationTable.hh"

namespace G4INCL {

  class NuclearDensity {
  public:
    NuclearDensity(const int A, const int Z, const int S, InterpolationTable const * const rpCorrelationTableProton, InterpolationTable const * const rpCorrelationTableNeutron, InterpolationTable const * const rpCorrelationTableLambda);
    ~NuclearDensity();

    /// \brief Copy constructor
    NuclearDensity(const NuclearDensity &rhs);

    /// \brief Assignment operator
    NuclearDensity &operator=(const NuclearDensity &rhs);

    /// \brief Helper method for the assignment operator
    void swap(NuclearDensity &rhs);

    /** \brief Get the maximum allowed radius for a given momentum.
     *  \param t type of the particle
     *  \param p absolute value of the particle momentum, divided by the
     *           relevant Fermi momentum.
     *  \return maximum allowed radius.
     */
    double getMaxRFromP(const ParticleType t, const double p) const;

    double getMinPFromR(const ParticleType t, const double r) const;

    double getMaximumRadius() const { return theMaximumRadius; };

    /** \brief The radius used for calculating the transmission coefficient.
     *
     * \return the radius
     */
    double getTransmissionRadius(Particle const * const p) const {
      const ParticleType t = p->getType();
      assert(t!= antiLambda && t!=antiNeutron && t!=Neutron && t!=PiZero && t!=DeltaZero && t!=Eta && t!=Omega && t!=EtaPrime && t!=Photon && t!= Lambda && t!=SigmaZero && t!=KZero && t!=KZeroBar && t!=KShort && t!=KLong); // no neutral particles here
      if(t==Composite) {
        return transmissionRadius[t] +
          ParticleTable::getNuclearRadius(t, p->getA(), p->getZ());
      } else if(t==antiComposite){
        return transmissionRadius[t] + ParticleTable::getNuclearRadius(t,-(p->getA()),-(p->getZ()));
      } else
        return transmissionRadius[t];
    };

    /** \brief The radius used for calculating the transmission coefficient.
     *
     * \return the radius
     */
    double getTransmissionRadius(ParticleType type) const {
      assert(type!=Composite);
      assert(type!=antiComposite);
      return transmissionRadius[type];
    };

    /// \brief Get the mass number.
    int getA() const { return theA; }

    /// \brief Get the charge number.
    int getZ() const { return theZ; }

    /// \brief Get the strange number.
    int getS() const { return theS; }

    double getProtonNuclearRadius() const { return theProtonNuclearRadius; }
    void setProtonNuclearRadius(const double r) { theProtonNuclearRadius = r; }

  private:

    /** \brief Initialize the transmission radius. */
    void initializeTransmissionRadii();

    int theA, theZ, theS;
    double theMaximumRadius;
    /// \brief Represents INCL4.5's R0 variable
    double theProtonNuclearRadius;

    /* \brief map of transmission radii per particle type */
    double transmissionRadius[UnknownParticle];

    InterpolationTable const *rFromP[UnknownParticle];
    InterpolationTable const *pFromR[UnknownParticle];
  };

}

#endif

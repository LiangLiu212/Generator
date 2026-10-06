#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"

#ifndef TransmissionChannel_hh
#define TransmissionChannel_hh 1
namespace G4INCL {
  class FinalState;

  class TransmissionChannel : public IChannel {
  public:
    TransmissionChannel(Nucleus *n, Particle *p);
    TransmissionChannel(Nucleus *n, Particle *p, const double TOut);
    TransmissionChannel(Nucleus *n, Particle *p, const double kOut, const double cosR);
    virtual ~TransmissionChannel();

    void fillFinalState(FinalState *fs);

  private:
    /** \brief Modify particle that leaves the nucleus.
     *
     * Modify the particle momentum and/or position when the particle leaves
     * the nucleus.
     */
    void particleLeaves();

    /** \brief Kinetic energy of the transmitted particle
     *
     * Calculate the kinetic energy of the particle outside the nucleus, if the
     * value has not been provided as a pre-calculated argument to the
     * constructor.
     */
    double initializeKineticEnergyOutside();

    Nucleus * const theNucleus;
    Particle * const theParticle;

    /// \brief True if refraction should be applied
    const bool refraction;

    /// \brief Momentum of the particle outside the nucleus
    const double pOutMag;

    /// \brief Kinetic energy of the particle outside the nucleus
    const double kineticEnergyOutside;

    /// \brief Cosine of the refraction angle
    const double cosRefractionAngle;

    INCL_DECLARE_ALLOCATION_POOL(TransmissionChannel)
  };
}
#endif // TransmissionChannel_hh

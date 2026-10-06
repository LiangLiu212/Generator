#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"

#ifndef REFLECTIONCHANNEL_HH_
#define REFLECTIONCHANNEL_HH_
namespace G4INCL {
  class FinalState;

  class ReflectionChannel : public IChannel {
  public:
    ReflectionChannel(Nucleus *n, Particle *p);
    virtual ~ReflectionChannel();

    void fillFinalState(FinalState *fs);

  private:
    /** \brief Sine^2 of the smallest acceptable reflection angle / 4
     *
     * Particles impinging almost tangentially on the surface generate a large
     * number of reflections. If the impinging angle is larger than a certain
     * limit, we move the particle a bit towards the inside of the nucleus.
     * The limit angle is arbitrarily chosen to correspond to 100 reflections
     * per orbit. The position scaling factor is given by the
     * positionScalingFactor member.
     */
    static const double sinMinReflectionAngleSquaredOverFour;
    /// \brief Scaling factor for excessively tangential reflection
    static const double positionScalingFactor;
    Nucleus *theNucleus;
    Particle *theParticle;

    INCL_DECLARE_ALLOCATION_POOL(ReflectionChannel)
  };
}
#endif //REFLECTIONCHANNEL_HH_

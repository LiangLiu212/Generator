#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLAllocationPool.hh"

#ifndef G4INCLParticleEntry_hh
#define G4INCLParticleEntry_hh 1

namespace G4INCL {
  class FinalState;

  class ParticleEntryChannel : public IChannel {
  public:
    ParticleEntryChannel(Nucleus *n, Particle *p);
    virtual ~ParticleEntryChannel();

    void fillFinalState(FinalState *fs);

  private:
    /** \brief Modify particle that enters the nucleus.
     *
     * Modify the particle momentum and/or position when the particle enters
     * the nucleus.
     *
     * \return false if the particle enters below 0 total energy.
     */
    bool particleEnters(const double theQValueCorrection);

    Nucleus *theNucleus;
    Particle *theParticle;

    INCL_DECLARE_ALLOCATION_POOL(ParticleEntryChannel)
  };
}

#endif

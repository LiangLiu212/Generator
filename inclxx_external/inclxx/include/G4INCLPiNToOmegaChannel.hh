#ifndef G4INCLPiNToOmegaChannel_hh
#define G4INCLPiNToOmegaChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class PiNToOmegaChannel : public IChannel {
    public:
      PiNToOmegaChannel(Particle *, Particle *);
      virtual ~PiNToOmegaChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(PiNToOmegaChannel);
  };
}

#endif

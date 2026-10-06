#ifndef G4INCLOmegaNToPiNChannel_hh
#define G4INCLOmegaNToPiNChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class OmegaNToPiNChannel : public IChannel {
    public:
      OmegaNToPiNChannel(Particle *, Particle *);
      virtual ~OmegaNToPiNChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(OmegaNToPiNChannel);
  };
}

#endif

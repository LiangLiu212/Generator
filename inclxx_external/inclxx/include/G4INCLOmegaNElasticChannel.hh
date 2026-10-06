#ifndef G4INCLOmegaNElasticChannel_hh
#define G4INCLOmegaNElasticChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class OmegaNElasticChannel : public IChannel {
    public:
      OmegaNElasticChannel(Particle *, Particle *);
      virtual ~OmegaNElasticChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(OmegaNElasticChannel);
  };
}

#endif

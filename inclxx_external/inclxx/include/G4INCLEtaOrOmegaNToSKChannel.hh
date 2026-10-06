#ifndef G4INCLEtaOrOmegaNToSKChannel_hh
#define G4INCLEtaOrOmegaNToSKChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class EtaOrOmegaNToSKChannel : public IChannel {
    public:
      EtaOrOmegaNToSKChannel(Particle *, Particle *);
      virtual ~EtaOrOmegaNToSKChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(EtaOrOmegaNToSKChannel);
  };
}

#endif

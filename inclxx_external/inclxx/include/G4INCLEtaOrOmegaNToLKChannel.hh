#ifndef G4INCLEtaOrOmegaNToLKChannel_hh
#define G4INCLEtaOrOmegaNToLKChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class EtaOrOmegaNToLKChannel : public IChannel {
    public:
      EtaOrOmegaNToLKChannel(Particle *, Particle *);
      virtual ~EtaOrOmegaNToLKChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(EtaOrOmegaNToLKChannel);
  };
}

#endif

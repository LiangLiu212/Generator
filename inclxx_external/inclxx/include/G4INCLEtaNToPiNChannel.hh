#ifndef G4INCLEtaNToPiNChannel_hh
#define G4INCLEtaNToPiNChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class EtaNToPiNChannel : public IChannel {
    public:
      EtaNToPiNChannel(Particle *, Particle *);
      virtual ~EtaNToPiNChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(EtaNToPiNChannel);
  };
}

#endif

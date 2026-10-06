#ifndef G4INCLEtaNElasticChannel_hh
#define G4INCLEtaNElasticChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class EtaNElasticChannel : public IChannel {
    public:
      EtaNElasticChannel(Particle *, Particle *);
      virtual ~EtaNElasticChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(EtaNElasticChannel);
  };
}

#endif

#ifndef G4INCLPiNElasticChannel_hh
#define G4INCLPiNElasticChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class PiNElasticChannel : public IChannel {
    public:
      PiNElasticChannel(Particle *, Particle *);
      virtual ~PiNElasticChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(PiNElasticChannel)
  };
}

#endif

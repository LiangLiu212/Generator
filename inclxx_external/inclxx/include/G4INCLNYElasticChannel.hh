#ifndef G4INCLNYElasticChannel_hh
#define G4INCLNYElasticChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NYElasticChannel : public IChannel {
    public:
      NYElasticChannel(Particle *, Particle *);
      virtual ~NYElasticChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NYElasticChannel);
  };
}

#endif

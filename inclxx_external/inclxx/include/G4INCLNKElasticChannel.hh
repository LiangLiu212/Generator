#ifndef G4INCLNKElasticChannel_hh
#define G4INCLNKElasticChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NKElasticChannel : public IChannel {
    public:
      NKElasticChannel(Particle *, Particle *);
      virtual ~NKElasticChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NKElasticChannel);
  };
}

#endif

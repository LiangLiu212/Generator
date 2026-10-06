#ifndef G4INCLNNbarElasticChannel_hh
#define G4INCLNNbarElasticChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNbarElasticChannel : public IChannel {
    public:
      NNbarElasticChannel(Particle *, Particle *);
      virtual ~NNbarElasticChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NNbarElasticChannel);
  };
}

#endif

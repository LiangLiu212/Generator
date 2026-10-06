#ifndef G4INCLNNbarToLLbarChannel_hh
#define G4INCLNNbarToLLbarChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNbarToLLbarChannel : public IChannel {
    public:
      NNbarToLLbarChannel(Particle *, Particle *);
      virtual ~NNbarToLLbarChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NNbarToLLbarChannel);
  };
}

#endif

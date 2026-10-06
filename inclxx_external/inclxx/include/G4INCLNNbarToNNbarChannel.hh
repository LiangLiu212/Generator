#ifndef G4INCLNNbarToNNbarChannel_hh
#define G4INCLNNbarToNNbarChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNbarToNNbarChannel : public IChannel {
    public:
      NNbarToNNbarChannel(Particle *, Particle *);
      virtual ~NNbarToNNbarChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NNbarToNNbarChannel);
  };
}

#endif

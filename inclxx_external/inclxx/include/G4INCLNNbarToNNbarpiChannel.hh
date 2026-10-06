#ifndef G4INCLNNbarToNNbarpiChannel_hh
#define G4INCLNNbarToNNbarpiChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNbarToNNbarpiChannel : public IChannel {
    public:
      NNbarToNNbarpiChannel(Particle *, Particle *);
      virtual ~NNbarToNNbarpiChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NNbarToNNbarpiChannel);
  };
}

#endif

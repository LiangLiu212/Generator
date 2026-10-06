#ifndef G4INCLNNbarToNNbar3piChannel_hh
#define G4INCLNNbarToNNbar3piChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNbarToNNbar3piChannel : public IChannel {
    public:
      NNbarToNNbar3piChannel(Particle *, Particle *);
      virtual ~NNbarToNNbar3piChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NNbarToNNbar3piChannel);
  };
}

#endif

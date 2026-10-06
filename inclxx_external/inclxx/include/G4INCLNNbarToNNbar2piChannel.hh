#ifndef G4INCLNNbarToNNbar2piChannel_hh
#define G4INCLNNbarToNNbar2piChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNbarToNNbar2piChannel : public IChannel {
    public:
      NNbarToNNbar2piChannel(Particle *, Particle *);
      virtual ~NNbarToNNbar2piChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NNbarToNNbar2piChannel);
  };
}

#endif

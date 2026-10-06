#ifndef G4INCLNNbarCEXChannel_hh
#define G4INCLNNbarCEXChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNbarCEXChannel : public IChannel {
    public:
      NNbarCEXChannel(Particle *, Particle *);
      virtual ~NNbarCEXChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NNbarCEXChannel);
  };
}

#endif

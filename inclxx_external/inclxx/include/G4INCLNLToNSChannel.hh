#ifndef G4INCLNLToNSChannel_hh
#define G4INCLNLToNSChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NLToNSChannel : public IChannel {
    public:
      NLToNSChannel(Particle *, Particle *);
      virtual ~NLToNSChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NLToNSChannel);
  };
}

#endif

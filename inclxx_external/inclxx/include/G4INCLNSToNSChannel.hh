#ifndef G4INCLNSToNSChannel_hh
#define G4INCLNSToNSChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NSToNSChannel : public IChannel {
    public:
      NSToNSChannel(Particle *, Particle *);
      virtual ~NSToNSChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NSToNSChannel);
  };
}

#endif

#ifndef G4INCLNSToNLChannel_hh
#define G4INCLNSToNLChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NSToNLChannel : public IChannel {
    public:
      NSToNLChannel(Particle *, Particle *);
      virtual ~NSToNLChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NSToNLChannel);
  };
}

#endif

#ifndef G4INCLNDeltaToNSKChannel_hh
#define G4INCLNDeltaToNSKChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NDeltaToNSKChannel : public IChannel {
    public:
      NDeltaToNSKChannel(Particle *, Particle *);
      virtual ~NDeltaToNSKChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NDeltaToNSKChannel);
  };
}

#endif

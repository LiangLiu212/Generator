#ifndef G4INCLNDeltaToNNKKbChannel_hh
#define G4INCLNDeltaToNNKKbChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NDeltaToNNKKbChannel : public IChannel {
    public:
      NDeltaToNNKKbChannel(Particle *, Particle *);
      virtual ~NDeltaToNNKKbChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NDeltaToNNKKbChannel);
  };
}

#endif

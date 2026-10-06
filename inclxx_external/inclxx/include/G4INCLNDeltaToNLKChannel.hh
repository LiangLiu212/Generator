#ifndef G4INCLNDeltaToNLKChannel_hh
#define G4INCLNDeltaToNLKChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NDeltaToNLKChannel : public IChannel {
    public:
      NDeltaToNLKChannel(Particle *, Particle *);
      virtual ~NDeltaToNLKChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NDeltaToNLKChannel);
  };
}

#endif

#ifndef G4INCLNNToNNKKbChannel_hh
#define G4INCLNNToNNKKbChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNToNNKKbChannel : public IChannel {
    public:
      NNToNNKKbChannel(Particle *, Particle *);
      virtual ~NNToNNKKbChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NNToNNKKbChannel);
  };
}

#endif

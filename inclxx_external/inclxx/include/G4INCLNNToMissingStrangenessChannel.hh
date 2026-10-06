#ifndef G4INCLNNToMissingStrangenessChannel_hh
#define G4INCLNNToMissingStrangenessChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNToMissingStrangenessChannel : public IChannel {
    public:
      NNToMissingStrangenessChannel(Particle *, Particle *);
      virtual ~NNToMissingStrangenessChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NNToMissingStrangenessChannel);
  };
}

#endif

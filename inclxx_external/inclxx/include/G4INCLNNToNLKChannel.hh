#ifndef G4INCLNNToNLKChannel_hh
#define G4INCLNNToNLKChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNToNLKChannel : public IChannel {
    public:
      NNToNLKChannel(Particle *, Particle *);
      virtual ~NNToNLKChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NNToNLKChannel);
  };
}

#endif

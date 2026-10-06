#ifndef G4INCLNNToNLK2piChannel_hh
#define G4INCLNNToNLK2piChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNToNLK2piChannel : public IChannel {
    public:
      NNToNLK2piChannel(Particle *, Particle *);
      virtual ~NNToNLK2piChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NNToNLK2piChannel);
  };
}

#endif

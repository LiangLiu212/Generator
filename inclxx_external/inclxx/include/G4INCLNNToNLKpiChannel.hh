#ifndef G4INCLNNToNLKpiChannel_hh
#define G4INCLNNToNLKpiChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NNToNLKpiChannel : public IChannel {
    public:
      NNToNLKpiChannel(Particle *, Particle *);
      virtual ~NNToNLKpiChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NNToNLKpiChannel);
  };
}

#endif

#ifndef G4INCLNKToNKpiChannel_hh
#define G4INCLNKToNKpiChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NKToNKpiChannel : public IChannel {
    public:
      NKToNKpiChannel(Particle *, Particle *);
      virtual ~NKToNKpiChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NKToNKpiChannel);
  };
}

#endif

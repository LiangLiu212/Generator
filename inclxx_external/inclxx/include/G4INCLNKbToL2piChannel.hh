#ifndef G4INCLNKbToL2piChannel_hh
#define G4INCLNKbToL2piChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NKbToL2piChannel : public IChannel {
    public:
      NKbToL2piChannel(Particle *, Particle *);
      virtual ~NKbToL2piChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NKbToL2piChannel);
  };
}

#endif

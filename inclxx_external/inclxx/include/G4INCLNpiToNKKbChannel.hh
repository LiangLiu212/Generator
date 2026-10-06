#ifndef G4INCLNpiToNKKbChannel_hh
#define G4INCLNpiToNKKbChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NpiToNKKbChannel : public IChannel {
    public:
      NpiToNKKbChannel(Particle *, Particle *);
      virtual ~NpiToNKKbChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NpiToNKKbChannel);
  };
}

#endif

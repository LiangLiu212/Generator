#ifndef G4INCLNKbToNKbpiChannel_hh
#define G4INCLNKbToNKbpiChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NKbToNKbpiChannel : public IChannel {
    public:
      NKbToNKbpiChannel(Particle *, Particle *);
      virtual ~NKbToNKbpiChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NKbToNKbpiChannel);
  };
}

#endif

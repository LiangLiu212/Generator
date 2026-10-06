#ifndef G4INCLNKbToNKbChannel_hh
#define G4INCLNKbToNKbChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NKbToNKbChannel : public IChannel {
    public:
      NKbToNKbChannel(Particle *, Particle *);
      virtual ~NKbToNKbChannel();

      void fillFinalState(FinalState *fs);

      ThreeVector KaonMomentum(Particle const * const kaon, Particle const * const nucleon);
      
    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NKbToNKbChannel);
  };
}

#endif

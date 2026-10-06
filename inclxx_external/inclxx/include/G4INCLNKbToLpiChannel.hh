#ifndef G4INCLNKbToLpiChannel_hh
#define G4INCLNKbToLpiChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NKbToLpiChannel : public IChannel {
    public:
      NKbToLpiChannel(Particle *, Particle *);
      virtual ~NKbToLpiChannel();

      void fillFinalState(FinalState *fs);

      ThreeVector KaonMomentum(Particle const * const kaon, Particle const * const nucleon);
      
    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NKbToLpiChannel);
  };
}

#endif

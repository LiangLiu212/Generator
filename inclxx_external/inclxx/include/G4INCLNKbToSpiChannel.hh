#ifndef G4INCLNKbToSpiChannel_hh
#define G4INCLNKbToSpiChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NKbToSpiChannel : public IChannel {
    public:
      NKbToSpiChannel(Particle *, Particle *);
      virtual ~NKbToSpiChannel();

      void fillFinalState(FinalState *fs);

      ThreeVector KaonMomentum(Particle const * const kaon, Particle const * const nucleon, const int WhichChannel);
      
    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NKbToSpiChannel);
  };
}

#endif

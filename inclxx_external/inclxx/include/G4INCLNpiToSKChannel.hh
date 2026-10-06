#ifndef G4INCLNpiToSKChannel_hh
#define G4INCLNpiToSKChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NpiToSKChannel : public IChannel {
    public:
      NpiToSKChannel(Particle *, Particle *);
      virtual ~NpiToSKChannel();

      void fillFinalState(FinalState *fs);
      
      ThreeVector KaonMomentum(Particle const * const pion, Particle const * const nucleon, int WhichChannel);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NpiToSKChannel);
  };
}

#endif

#ifndef G4INCLNpiToLKChannel_hh
#define G4INCLNpiToLKChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NpiToLKChannel : public IChannel {
    public:
      NpiToLKChannel(Particle *, Particle *);
      virtual ~NpiToLKChannel();

      void fillFinalState(FinalState *fs);
      
      ThreeVector KaonMomentum(Particle const * const pion, Particle const * const nucleon);

    private:
      Particle *particle1, *particle2;
      
      INCL_DECLARE_ALLOCATION_POOL(NpiToLKChannel);
  };
}

#endif

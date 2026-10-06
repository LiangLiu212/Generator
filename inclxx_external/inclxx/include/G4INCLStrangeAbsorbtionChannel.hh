#ifndef G4INCLStrangeAbsorbtionChannel_hh
#define G4INCLStrangeAbsorbtionChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class StrangeAbsorbtionChannel : public IChannel {
    public:
      StrangeAbsorbtionChannel(Particle *, Particle *);
      virtual ~StrangeAbsorbtionChannel();

      void fillFinalState(FinalState *fs);

    private:
      void sampleAngles(double*, double*, double*);
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(StrangeAbsorbtionChannel);
  };
}

#endif

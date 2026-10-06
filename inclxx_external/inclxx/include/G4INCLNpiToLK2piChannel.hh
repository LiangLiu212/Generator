#ifndef G4INCLNpiToLK2piChannel_hh
#define G4INCLNpiToLK2piChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NpiToLK2piChannel : public IChannel {
    public:
      NpiToLK2piChannel(Particle *, Particle *);
      virtual ~NpiToLK2piChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NpiToLK2piChannel);
  };
}

#endif

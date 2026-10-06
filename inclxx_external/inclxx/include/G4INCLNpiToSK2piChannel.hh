#ifndef G4INCLNpiToSK2piChannel_hh
#define G4INCLNpiToSK2piChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NpiToSK2piChannel : public IChannel {
    public:
      NpiToSK2piChannel(Particle *, Particle *);
      virtual ~NpiToSK2piChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NpiToSK2piChannel);
  };
}

#endif

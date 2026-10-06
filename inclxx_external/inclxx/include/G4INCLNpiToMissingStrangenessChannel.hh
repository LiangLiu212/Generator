#ifndef G4INCLNpiToMissingStrangenessChannel_hh
#define G4INCLNpiToMissingStrangenessChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NpiToMissingStrangenessChannel : public IChannel {
    public:
      NpiToMissingStrangenessChannel(Particle *, Particle *);
      virtual ~NpiToMissingStrangenessChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NpiToMissingStrangenessChannel);
  };
}

#endif

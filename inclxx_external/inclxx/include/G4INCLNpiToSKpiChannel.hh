#ifndef G4INCLNpiToSKpiChannel_hh
#define G4INCLNpiToSKpiChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NpiToSKpiChannel : public IChannel {
    public:
      NpiToSKpiChannel(Particle *, Particle *);
      virtual ~NpiToSKpiChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NpiToSKpiChannel);
  };
}

#endif

#ifndef G4INCLNpiToLKpiChannel_hh
#define G4INCLNpiToLKpiChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NpiToLKpiChannel : public IChannel {
    public:
      NpiToLKpiChannel(Particle *, Particle *);
      virtual ~NpiToLKpiChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NpiToLKpiChannel);
  };
}

#endif

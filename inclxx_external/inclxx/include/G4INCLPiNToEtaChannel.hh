#ifndef G4INCLPiNToEtaChannel_hh
#define G4INCLPiNToEtaChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class PiNToEtaChannel : public IChannel {
    public:
      PiNToEtaChannel(Particle *, Particle *);
      virtual ~PiNToEtaChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(PiNToEtaChannel);
  };
}

#endif

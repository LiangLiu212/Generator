#ifndef G4INCLPiNToDeltaChannel_hh
#define G4INCLPiNToDeltaChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class PiNToDeltaChannel : public IChannel {
    public:
      PiNToDeltaChannel(Particle *, Particle *);
      virtual ~PiNToDeltaChannel();

      void fillFinalState(FinalState *fs);

    private:
      Particle *particle1, *particle2;

      INCL_DECLARE_ALLOCATION_POOL(PiNToDeltaChannel)
  };
}

#endif

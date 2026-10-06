#ifndef G4INCLOmegaNToPiPiNChannel_hh
#define G4INCLOmegaNToPiPiNChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class OmegaNToPiPiNChannel : public IChannel {
    public:
      OmegaNToPiPiNChannel(Particle *, Particle *);
      virtual ~OmegaNToPiPiNChannel();

      void fillFinalState(FinalState *fs);

    private:
      int ind2; // can be changed
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(OmegaNToPiPiNChannel);
  };
}

#endif

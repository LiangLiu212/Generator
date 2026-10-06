#ifndef G4INCLEtaNToPiPiNChannel_hh
#define G4INCLEtaNToPiPiNChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class EtaNToPiPiNChannel : public IChannel {
    public:
      EtaNToPiPiNChannel(Particle *, Particle *);
      virtual ~EtaNToPiPiNChannel();

      void fillFinalState(FinalState *fs);

    private:
      int ind2; // can be changed
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(EtaNToPiPiNChannel);
  };
}

#endif

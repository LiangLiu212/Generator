#ifndef G4INCLPiNToMultiPionsChannel_hh
#define G4INCLPiNToMultiPionsChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class PiNToMultiPionsChannel : public IChannel {
    public:
      PiNToMultiPionsChannel(const int, Particle *, Particle *);
      virtual ~PiNToMultiPionsChannel();

      void fillFinalState(FinalState *fs);

    private:
      int npion;
      int ind2; // like isosp, can be changed in isospinRepartition
      int isosp[4];
      Particle *particle1, *particle2;

      static const double angularSlope;

      void isospinRepartition(int ipi);

      INCL_DECLARE_ALLOCATION_POOL(PiNToMultiPionsChannel)
  };
}

#endif

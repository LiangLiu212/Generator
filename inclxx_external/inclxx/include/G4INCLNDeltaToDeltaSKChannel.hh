#ifndef G4INCLNDeltaToDeltaSKChannel_hh
#define G4INCLNDeltaToDeltaSKChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NDeltaToDeltaSKChannel : public IChannel {
    public:
      NDeltaToDeltaSKChannel(Particle *, Particle *);
      virtual ~NDeltaToDeltaSKChannel();

      void fillFinalState(FinalState *fs);

    private:
      double sampleDeltaMass(double ecm);
      
      Particle *particle1, *particle2;

      static const double angularSlope;

      INCL_DECLARE_ALLOCATION_POOL(NDeltaToDeltaSKChannel);
  };
}

#endif

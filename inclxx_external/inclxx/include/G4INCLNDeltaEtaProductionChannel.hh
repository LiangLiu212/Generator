#ifndef G4INCLNDeltaEtaProductionChannel_hh
#define G4INCLNDeltaEtaProductionChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NDeltaEtaProductionChannel : public IChannel {
  public:
    NDeltaEtaProductionChannel(Particle *, Particle *);
    virtual ~NDeltaEtaProductionChannel();

    void fillFinalState(FinalState *fs);

  private:
    double sampleDeltaMass(double ecm);

    Particle *particle1, *particle2;
   
    static const double angularSlope;

    static const int maxTries;
    INCL_DECLARE_ALLOCATION_POOL(NDeltaEtaProductionChannel);
  };
}

#endif

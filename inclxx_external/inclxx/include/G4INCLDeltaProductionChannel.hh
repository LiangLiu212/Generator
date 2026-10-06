#ifndef G4INCLDeltaProductionChannel_hh
#define G4INCLDeltaProductionChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"
#include "G4INCLSrcChannel.hh"

namespace G4INCL {
  class DeltaProductionChannel : public IChannel {
  public:
    DeltaProductionChannel(Particle *, Particle *, Nucleus *n = nullptr);
    virtual ~DeltaProductionChannel();

    void fillFinalState(FinalState *fs);

  private:
    double sampleDeltaMass(double ecm);

    Particle *particle1, *particle2;
    Nucleus *thenucleus;
    SrcChannel *srcChannel;

    static const int maxTries;
    INCL_DECLARE_ALLOCATION_POOL(DeltaProductionChannel)
  };
}

#endif

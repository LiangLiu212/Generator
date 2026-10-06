#ifndef G4INCLDeltaDecayChannel_hh
#define G4INCLDeltaDecayChannel_hh 1

#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class DeltaDecayChannel : public IChannel {
  public:
    DeltaDecayChannel(Particle *, ThreeVector const &);
    virtual ~DeltaDecayChannel();

    static double computeDecayTime(Particle *p);
    void fillFinalState(FinalState *fs);

  private:
    void sampleAngles(double*, double*, double*);

    Particle *theParticle;
    ThreeVector const incidentDirection;

    INCL_DECLARE_ALLOCATION_POOL(DeltaDecayChannel)
  };
}

#endif

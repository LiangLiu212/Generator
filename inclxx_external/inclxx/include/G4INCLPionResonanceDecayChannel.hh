#ifndef G4INCLPionResonanceDecayChannel_hh
#define G4INCLPionResonanceDecayChannel_hh 1

#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class PionResonanceDecayChannel : public IChannel {
  public:
    PionResonanceDecayChannel(Particle *, ThreeVector const &);
    virtual ~PionResonanceDecayChannel();

    static double computeDecayTime(Particle *p);
    void fillFinalState(FinalState *fs);

  private:
    void sampleAngles(double*, double*, double*);
	  
	static const double angularSlope;

    Particle *theParticle;
    ThreeVector const incidentDirection;

    INCL_DECLARE_ALLOCATION_POOL(PionResonanceDecayChannel);
  };
}

#endif

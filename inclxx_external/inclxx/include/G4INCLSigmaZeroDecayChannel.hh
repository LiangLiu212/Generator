#ifndef G4INCLSigmaZeroDecayChannel_hh
#define G4INCLSigmaZeroDecayChannel_hh 1

#include "G4INCLIChannel.hh"
#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class SigmaZeroDecayChannel : public IChannel {
  public:
    SigmaZeroDecayChannel(Particle *, ThreeVector const &);
    virtual ~SigmaZeroDecayChannel();
    
	static double computeDecayTime(Particle *p);
    void fillFinalState(FinalState *fs);

  private:
    void sampleAngles(double*, double*, double*);
    
    Particle *theParticle;
    ThreeVector const incidentDirection;
    
    INCL_DECLARE_ALLOCATION_POOL(SigmaZeroDecayChannel);
  };
}

#endif

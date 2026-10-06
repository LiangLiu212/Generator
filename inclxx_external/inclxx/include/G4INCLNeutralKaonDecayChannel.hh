#ifndef G4INCLNeutralKaonDecayChannel_hh
#define G4INCLNeutralKaonDecayChannel_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

namespace G4INCL {
  class NeutralKaonDecayChannel : public IChannel {
  public:
    NeutralKaonDecayChannel(Particle *);
    virtual ~NeutralKaonDecayChannel();
    
    void fillFinalState(FinalState *fs);

  private:
    
    Particle *theParticle;
    
    INCL_DECLARE_ALLOCATION_POOL(NeutralKaonDecayChannel);
  };
}

#endif

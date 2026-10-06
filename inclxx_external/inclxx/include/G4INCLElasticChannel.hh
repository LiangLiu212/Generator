
#include "G4INCLAllocationPool.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLParticle.hh"
#include "G4INCLSrcChannel.hh"

#ifndef G4INCLElasticChannel_HH
#define G4INCLElasticChannel_HH 1

namespace G4INCL {
class ElasticChannel : public IChannel {

public:
  ElasticChannel(Particle *p1, Particle *p2, Nucleus *n = nullptr);
  virtual ~ElasticChannel();

  void fillFinalState(FinalState *fs);

private:
  Particle *particle1, *particle2;
  Nucleus *thenucleus;
  SrcChannel *srcChannel;

  INCL_DECLARE_ALLOCATION_POOL(ElasticChannel)
};

} // namespace G4INCL

#endif

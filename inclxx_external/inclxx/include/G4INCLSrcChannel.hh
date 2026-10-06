/*
 * \file G4INCLSrcChannel.hh
 *
 * \date Feb 24, 2022
 * \author Jose Luis Rodriguez-Sanchez
 */

#include "G4INCLAllocationPool.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLParticle.hh"

#include "G4INCLEventInfo.hh"

#ifndef G4INCLSrcChannel_HH
#define G4INCLSrcChannel_HH 1

namespace G4INCL {
class SrcChannel : public IChannel {

public:
  SrcChannel(Particle *p1, Particle *p2, Nucleus *n);
  virtual ~SrcChannel();

  void fillFinalState(FinalState *fs);
  void fillFinalState(FinalState *fs, ParticleType , ParticleType);

private:
  Particle *particle1, *particle2;
  ParticleType ftype1, ftype2;
  Particle *srcpartner;
  Nucleus *thenucleus;
  double fDistSrc;
  
  EventInfo theEventInfo;

  /**
   *  Compute the current number of src pairs.
   */
  Particle *findpairpartner(Particle *pt);

  INCL_DECLARE_ALLOCATION_POOL(SrcChannel)
};

} // namespace G4INCL

#endif /* G4INCLSrcChannel_HH */

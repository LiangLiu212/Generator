/** \file G4INCLRecombinationChannel.hh
 * \brief Delta-nucleon recombination channel.
 *
 *  \date 25 March 2011
 * \author Davide Mancusi
 */

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIChannel.hh"
#include "G4INCLFinalState.hh"
#include "G4INCLAllocationPool.hh"

#ifndef G4INCLRECOMBINATIONCHANNEL_HH_
#define G4INCLRECOMBINATIONCHANNEL_HH_

namespace G4INCL {
  class RecombinationChannel : public IChannel {

  public:
    RecombinationChannel(Particle *p1, Particle *p2);
    virtual ~RecombinationChannel();

    void fillFinalState(FinalState *fs);

  private:
    Particle *theNucleon, *theDelta;

    INCL_DECLARE_ALLOCATION_POOL(RecombinationChannel)
  };

}

#endif /* G4INCLRECOMBINATIONCHANNEL_HH_ */

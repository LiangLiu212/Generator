/*
 * G4INCLIChannel.hh
 *
 *  \date Jun 5, 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLICHANNEL_H_
#define G4INCLICHANNEL_H_

#include "G4INCLAllocationPool.hh"

namespace G4INCL {

  class FinalState;

  /**
   * Channel generates a final state of an avatar.
   */
  class IChannel {
  public:
    IChannel() {}
    virtual ~IChannel() {}

    FinalState *getFinalState();
    virtual void fillFinalState(FinalState *fs) = 0;

    INCL_DECLARE_ALLOCATION_POOL(IChannel)
  };

}

#endif /* G4INCLICHANNEL_H_ */

/*
 * G4INCLIRandomGenerator.hh
 *
 *  \date 7 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLIRANDOMGENERATOR_HH_
#define G4INCLIRANDOMGENERATOR_HH_

#include "G4INCLRandomSeedVector.hh"

namespace G4INCL {

  class IRandomGenerator {
  public:
    IRandomGenerator() {}
    virtual ~IRandomGenerator() {}

    virtual Random::SeedVector getSeeds() = 0;
    virtual void setSeeds(const Random::SeedVector &) = 0;
    virtual double flat() = 0;
  };

}

#endif /* G4INCLIRANDOMGENERATOR_HH_ */

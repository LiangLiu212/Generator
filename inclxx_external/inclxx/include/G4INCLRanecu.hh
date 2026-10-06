/*
 * G4INCLRanecu.hh
 *
 *  \date 7 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLRANECU_HH_
#define G4INCLRANECU_HH_

#include "G4INCLIRandomGenerator.hh"
#include <cassert>

namespace G4INCL {

  class Ranecu: public G4INCL::IRandomGenerator {
  public:
    Ranecu();
    Ranecu(const Random::SeedVector &sv);
    virtual ~Ranecu();

    Random::SeedVector getSeeds() {
      Random::SeedVector sv;
      sv.push_back(iseed1);
      sv.push_back(iseed2);
      return sv;
    }

    void setSeeds(const Random::SeedVector &sv) {
      assert(sv.size()>=2);
      iseed1 = sv.at(0);
      iseed2 = sv.at(1);
    }

    double flat();

  private:
    int iseed1, iseed2;
  };

}

#endif /* G4INCLRANECU_HH_ */

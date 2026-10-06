/*
 * G4INCLRanecu.hh
 *
 *  \date 7 June 2009
 * \author Pekka Kaitaniemi
 */

#ifndef G4INCLCONSTANT_RANDOM_HH
#define G4INCLCONSTANT_RANDOM_HH 1

#include "G4INCLIRandomGenerator.hh"

namespace G4INCL {

  class ConstantRandom : public G4INCL::IRandomGenerator {
  public:
    ConstantRandom() {};
    virtual ~ConstantRandom() {};

    SeedVector getSeeds() const {
      return SeedVector;
    }

    void setSeeds(const SeedVector &) {};
    double flat() { return 0.4; }
  };

}

#endif

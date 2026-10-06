#ifndef G4INCLGeant4Random_hh
#define G4INCLGeant4Random_hh 1

#ifdef INCLXX_IN_GEANT4_MODE

#include "globals.hh"
#include "Randomize.hh"

#include "G4INCLLogger.hh"
#include "G4INCLIRandomGenerator.hh"

namespace G4INCL {

  class Geant4RandomGenerator : public G4INCL::IRandomGenerator {
  public:
    Geant4RandomGenerator() {};
    Geant4RandomGenerator(const Random::SeedVector &) {};
    virtual ~Geant4RandomGenerator() {};

    Random::SeedVector getSeeds() {
      INCL_WARN("getSeeds not supported.");
      Random::SeedVector sv;
      return sv;
    }

    void setSeeds(const Random::SeedVector &) {
      INCL_WARN("setSeeds not supported.");
    }

    double flat() {
      return G4UniformRand();
    }
  };

}

#endif // INCLXX_IN_GEANT4_MODE

#endif

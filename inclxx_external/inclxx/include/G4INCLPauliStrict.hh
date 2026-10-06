#ifndef G4INCLPauliStrict_hh
#define G4INCLPauliStrict_hh 1

#include "G4INCLIPauli.hh"

namespace G4INCL {
  class PauliStrict : public IPauli {
  public:
    PauliStrict();
    ~PauliStrict();

    bool isBlocked(ParticleList const &, Nucleus const * const);
  };
}

#endif

#ifndef G4INCLPauliGlobal_hh
#define G4INCLPauliGlobal_hh 1

#include "G4INCLIPauli.hh"

namespace G4INCL {
  class PauliGlobal : public IPauli {
  public:
    PauliGlobal();
    ~PauliGlobal();

    bool isBlocked(ParticleList const &, Nucleus const * const);
  };
}

#endif

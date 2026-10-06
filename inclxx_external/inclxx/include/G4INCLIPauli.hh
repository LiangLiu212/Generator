#ifndef G4INCLIPauli_hh
#define G4INCLIPauli_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {
  class IPauli {
  public:
    IPauli() {};
    virtual ~IPauli() {};

    virtual bool isBlocked(ParticleList const &, Nucleus const * const) = 0;

  protected:
  };
}

#endif

#ifndef G4INCLPauliStrictStandard_hh
#define G4INCLPauliStrictStandard_hh 1

#include "G4INCLIPauli.hh"
#include "G4INCLNucleus.hh"

namespace G4INCL {

  class PauliStrictStandard : public IPauli {
  public:
    PauliStrictStandard();
    ~PauliStrictStandard();

    /// \brief Dummy copy constructor to silence Coverity warning
    PauliStrictStandard(const PauliStrictStandard &rhs);

    /// \brief Dummy assignment operator to silence Coverity warning
    PauliStrictStandard &operator=(const PauliStrictStandard &rhs);

    bool isBlocked(ParticleList const &, Nucleus const * const);

  private:
    IPauli *theStrictBlocker, *theStandardBlocker;
  };
}

#endif

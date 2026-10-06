#ifndef G4INCLPauliStandard_hh
#define G4INCLPauliStandard_hh 1

#include "G4INCLIPauli.hh"

namespace G4INCL {
  class PauliStandard : public IPauli {
  public:
    PauliStandard();
    ~PauliStandard();

    bool isBlocked(ParticleList const &, Nucleus const * const);
    double getBlockingProbability(Particle const * const, Nucleus const * const) const;

  private:
    /** \brief Phase-space blocking cell
     *
     * This variable is the third root of volume of the blocking phase-space
     * cell. Its value is
     * \f[ {\left(2.38\cdot\frac{9\pi}{2}\right)}^{1/6}\sqrt\hbar \f]
     * which yields a phase-space test volume of \f$2.38h\f$.
     */
    const double cellSize;
  };
}

#endif

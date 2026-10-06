#ifndef G4INCLCDPP_hh
#define G4INCLCDPP_hh 1

#include "G4INCLParticle.hh"
#include "G4INCLNucleus.hh"
#include "G4INCLIPauli.hh"
#include "G4INCLINuclearPotential.hh"

namespace G4INCL {
  class CDPP : public IPauli {
  public:
    CDPP();
    ~CDPP();

    bool isBlocked(ParticleList const &, Nucleus const * const);

    inline void processOneParticle(Particle const * const p) {
      if(p->isBaryon()) {
        const double Tf = thePotential->getFermiEnergy(p);
        const double T = p->getKineticEnergy();

        if(T > Tf) {
          const double sep = thePotential->getSeparationEnergy(p);
          Sk += sep;
        } else {
          TbelowTf += T - p->getPotentialEnergy();
        }
      } else if(p->isMeson() || p->isResonance()) {
        const double sep = thePotential->getSeparationEnergy(p);
        Sk += sep;
      }
    }

  protected:
    double Sk;
    double TbelowTf;
    NuclearPotential::INuclearPotential const * thePotential;

  };
}

#endif

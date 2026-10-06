/** \file G4INCLNuclearPotentialConstant.hh
 * \brief Isospin- and energy-independent nuclear potential.
 *
 * Provides a constant nuclear potential (V0).
 *
 * \date 17 January 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLNUCLEARPOTENTIALCONSTANT_HH
#define G4INCLNUCLEARPOTENTIALCONSTANT_HH 1

#include "G4INCLINuclearPotential.hh"

namespace G4INCL {

  namespace NuclearPotential {

    class NuclearPotentialConstant : public INuclearPotential {

      public:
        NuclearPotentialConstant(const int A, const int Z, const bool pionPotential);
        virtual ~NuclearPotentialConstant();

        double getNucleonPotential() const { return vNucleon; }
        double getDeltaPotential() const { return vDelta; }

        virtual double computePotentialEnergy(const Particle * const p) const;

      private:
        double vNucleon, vDelta, vSigma, vLambda;

        void initialize();

    };

  }
}

#endif /* G4INCLNUCLEARPOTENTIALCONSTANT_HH */

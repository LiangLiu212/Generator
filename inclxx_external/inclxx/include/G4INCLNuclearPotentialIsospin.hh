/** \file G4INCLNuclearPotentialIsospin.hh
 * \brief Isospin-dependent nuclear potential.
 *
 * Provides an isospin-dependent nuclear potential.
 *
 * \date 28 February 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLNUCLEARPOTENTIALISOSPIN_HH
#define G4INCLNUCLEARPOTENTIALISOSPIN_HH 1

#include "G4INCLINuclearPotential.hh"

namespace G4INCL {

  namespace NuclearPotential {

    class NuclearPotentialIsospin : public INuclearPotential {

      public:
        NuclearPotentialIsospin(const int A, const int Z, const bool pionPotential);
        virtual ~NuclearPotentialIsospin();

        virtual double computePotentialEnergy(const Particle * const p) const;

      private:
        double vProton, vNeutron;
        double vDeltaPlusPlus, vDeltaPlus, vDeltaZero, vDeltaMinus;
        double vSigmaPlus, vSigmaZero, vSigmaMinus, vLambda;
        double vantiProton, vantiNeutron;

        void initialize();

    };

  }
}

#endif /* G4INCLNUCLEARPOTENTIALISOSPIN_HH */

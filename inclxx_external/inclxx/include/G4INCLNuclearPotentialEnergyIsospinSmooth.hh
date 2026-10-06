/** \file G4INCLNuclearPotentialEnergyIsospinSmooth.hh
 * \brief Isospin- and energy-dependent nuclear potential.
 *
 * Provides an isospin- and energy-dependent nuclear potential. The cusp at 200
 * MeV is replaced by a smooth exponential.
 *
 * \date 16 February 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLNUCLEARPOTENTIALENERGYISOSPINSMOOTH_HH
#define G4INCLNUCLEARPOTENTIALENERGYISOSPINSMOOTH_HH 1

#include "G4INCLNuclearPotentialIsospin.hh"

namespace G4INCL {

  namespace NuclearPotential {

    class NuclearPotentialEnergyIsospinSmooth : public NuclearPotentialIsospin {

      public:
        NuclearPotentialEnergyIsospinSmooth(const int A, const int Z, const bool pionPotential);
        virtual ~NuclearPotentialEnergyIsospinSmooth();

        virtual double computePotentialEnergy(const Particle * const p) const;

      private:
        /// Slope of the V(T) curve
        static const double alpha;

        /// Distance from the cusp where the exponential kicks in
        static const double deltaE;

    };

  }
}

#endif /* G4INCLNUCLEARPOTENTIALENERGYISOSPINSMOOTH_HH */

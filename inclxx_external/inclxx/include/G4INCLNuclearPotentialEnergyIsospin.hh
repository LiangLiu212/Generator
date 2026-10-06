/** \file G4INCLNuclearPotentialEnergyIsospin.hh
 * \brief Isospin- and energy-dependent nuclear potential.
 *
 * Provides an isospin- and energy-dependent nuclear potential.
 *
 * \date 21 March 2011
 * \author Davide Mancusi
 */

#ifndef G4INCLNUCLEARPOTENTIALENERGYISOSPIN_HH
#define G4INCLNUCLEARPOTENTIALENERGYISOSPIN_HH 1

#include "G4INCLNuclearPotentialIsospin.hh"

namespace G4INCL {

  namespace NuclearPotential {

    class NuclearPotentialEnergyIsospin : public NuclearPotentialIsospin {

      public:
        NuclearPotentialEnergyIsospin(const int A, const int Z, const bool pionPotential);
        virtual ~NuclearPotentialEnergyIsospin();

        virtual double computePotentialEnergy(const Particle * const p) const;

      private:
        /// Slope of the V(T) curve
      static const double alpha;

    };

  }
}

#endif /* G4INCLNUCLEARPOTENTIALENERGYISOSPIN_HH */

/** \file G4INCLIPhaseSpaceGenerator.hh
 * \brief Abstract interface for the phase-space generators
 *
 * \date 15th November 2013
 * \author Davide Mancusi
 */

#ifndef G4INCLIPHASESPACEGENERATOR_HH
#define G4INCLIPHASESPACEGENERATOR_HH

#include "G4INCLParticle.hh"
#include "G4INCLThreeVector.hh"

namespace G4INCL {
  /// \brief Abstract interface for the phase-space generators
  class IPhaseSpaceGenerator {
    public:

      IPhaseSpaceGenerator() {}
      virtual ~IPhaseSpaceGenerator() {}

      /// \brief Generate an event in the CM frame
      virtual void generate(const double sqrtS, ParticleList &particles) = 0;

  };
}

#endif

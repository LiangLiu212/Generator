#ifndef G4INCLPHASESPACEKOPYLOV_HH
#define G4INCLPHASESPACEKOPYLOV_HH

#include "G4INCLThreeVector.hh"
#include "G4INCLParticle.hh"
#include "G4INCLIPhaseSpaceGenerator.hh"
#include <vector>

namespace G4INCL {
  /// \brief Generate momenta using the Kopylov method
  class PhaseSpaceKopylov : public IPhaseSpaceGenerator {

    public:
      /** \brief Generate momenta according to a uniform, non-Lorentz-invariant phase-space model
       *
       * This function will assign momenta to the particles in the list that is
       * passed as an argument. The event is generated in the CM frame.
       *
       * \param sqrtS total centre-of-mass energy of the system
       * \param particles list of particles
       */
      void generate(const double sqrtS, ParticleList &particles);

    private:
      /// \brief Internal function used by the Kopylov algorithm
      double betaKopylov(int K) const;

      std::vector<double> masses;
      std::vector<double> sumMasses;
      ThreeVector PFragCM, PRestCM;
      ThreeVector boostV;
  };
}

#endif // G4INCLPHASESPACEKOPYLOV_HH
